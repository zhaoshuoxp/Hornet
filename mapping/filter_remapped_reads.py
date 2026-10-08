#!/usr/bin/env python3
"""filter_remapped_reads

Keep the reads (pairs) from to.remap.bam for which every allelic version
remapped exactly like the original read.

Default (WASP-like) criteria, for each remapped record of a read (pair):

* not a secondary or supplementary alignment (any such record discards the
  read: the swapped read maps to more than one place);
* mapped to the same chromosome, start position and strand as the original,
  with the same CIGAR (hard clips ignored);
* paired-end: flagged as a proper pair, and its mate is at the original
  mate position;
* the number of such records equals the number of versions written by
  find_intersecting_snps (times two for pairs).

The remapped BAM can be in any order.

--legacy runs the original Hornet algorithm instead (start position only,
CIGAR restricted to M/N/S operations; needs a remapped BAM sorted by read
name), reproducing its output exactly.

"""
from __future__ import print_function
import sys
import gzip
import re
import time
from array import array

import pysam

_HARD_CLIP = re.compile(r'[0-9]+H')


def _cigar_key(cigarstring):
    if cigarstring is None:
        return 0
    if 'H' in cigarstring:
        cigarstring = _HARD_CLIP.sub('', cigarstring)
    return hash(cigarstring)


def _now():
    return time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))


def run_strict(to_remap_path, remap_path, keep_path, is_paired_end, threads=1):
    """WASP-like filter; see the module docstring."""
    sys.stderr.write("{} Starting run\n".format(_now()))
    threads = max(1, threads)
    n = 2 if is_paired_end else 1

    # Pass 1: the original reads, in remap-ID order (n records per ID)
    to_remap = pysam.AlignmentFile(to_remap_path, "rb", threads=threads)
    orig_tid = array('i', [-1])
    orig_pos = array('q', [0] * n)
    orig_rev = bytearray(n)
    orig_cig = array('q', [0] * n)
    k = 0
    for rec in to_remap:
        if k % n == 0:
            orig_tid.append(rec.reference_id)
        orig_pos.append(rec.reference_start)
        orig_rev.append(1 if rec.is_reverse else 0)
        orig_cig.append(_cigar_key(rec.cigarstring))
        k += 1
    if k % n:
        sys.exit("ERROR: %s has an odd number of reads" % to_remap_path)
    n_ids = k // n
    names = list(to_remap.references)
    to_remap.close()

    # Pass 2: remapped reads, any order
    remap = pysam.AlignmentFile(remap_path, "rb", threads=threads)
    name_tid = {c: i for i, c in enumerate(names)}
    tidmap = [name_tid.get(c, -2) for c in remap.references]
    expected = array('i', [0] * (n_ids + 1))
    correct = array('i', [0] * (n_ids + 1))
    moved = bytearray(n_ids + 1)
    multi = bytearray(n_ids + 1)
    for rec in remap:
        fields = rec.query_name.split(':')
        last = fields[-1]
        if '/' in last:
            last = last.split('/', 1)[0]
        try:
            rid = int(fields[0])
            num = int(last)
        except ValueError:
            sys.exit("ERROR: unexpected read name %s in %s (not from "
                     "find_intersecting_snps?)" % (rec.query_name, remap_path))
        if rid < 1 or rid > n_ids:
            sys.exit("ERROR: read %s has no counterpart in %s"
                     % (rec.query_name, to_remap_path))
        expected[rid] = num * n
        flag = rec.flag
        if flag & 0x900:
            multi[rid] = 1
            continue
        if flag & 4:
            continue
        m = 1 if (n == 2 and flag & 0x80) else 0
        j = rid * n + m
        tid = tidmap[rec.reference_id]
        ok = (tid == orig_tid[rid] and rec.reference_start == orig_pos[j]
              and (1 if flag & 16 else 0) == orig_rev[j]
              and _cigar_key(rec.cigarstring) == orig_cig[j])
        if ok and n == 2:
            ok = (flag & 2 and not flag & 8 and rec.next_reference_id >= 0
                  and tidmap[rec.next_reference_id] == tid
                  and rec.next_reference_start == orig_pos[rid * n + 1 - m])
        if ok:
            correct[rid] += 1
        else:
            moved[rid] = 1
    remap.close()

    # Pass 3: write the kept reads
    to_remap = pysam.AlignmentFile(to_remap_path, "rb", threads=threads)
    keep = pysam.AlignmentFile(keep_path, "wb", template=to_remap,
                               threads=threads)
    n_keep = n_moved = n_multi = n_missing = 0
    k = 0
    keep_write = keep.write
    for rec in to_remap:
        rid = k // n + 1
        k += 1
        if multi[rid]:
            if k % n == 0:
                n_multi += 1
        elif moved[rid]:
            if k % n == 0:
                n_moved += 1
        elif expected[rid] == 0 or correct[rid] != expected[rid]:
            if k % n == 0:
                n_missing += 1
        else:
            keep_write(rec)
            if k % n == 0:
                n_keep += 1
    keep.close()
    to_remap.close()
    sys.stderr.write("{} Finished\n\n".format(_now()))
    sys.stderr.write(
        ("RUN STATISTICS:\n\tTotal remapped read (pair)s: {}\n\t"
         "Read (pair)s remapped to the correct position: {} {:.2%}\n"
         "\tDiscarded, remapped to another position, strand or CIGAR: {}\n"
         "\tDiscarded, secondary/supplementary alignments: {}\n"
         "\tDiscarded, unmapped or missing remapped reads: {}\n")
        .format(n_ids, n_keep, n_keep / float(n_ids) if n_ids else 0.0,
                n_moved, n_multi, n_missing))


def _remapped_ok(remap_bam, is_paired_end):
    """First pass over the name-sorted remapped reads.

    Returns (correct_maps, nums) exactly as built by the original code:
    correct_maps -- remap IDs of every record that mapped back to the original
                    position without indels (one entry per record)
    nums         -- expected number of such records for each remap ID
    """
    refnames = list(remap_bam.references)
    # Index of the position field to compare with, for forward/reverse reads
    rev_field = 3 if is_paired_end else 2
    mult = 2 if is_paired_end else 1

    correct_maps = []
    nums = []
    append_correct = correct_maps.append
    append_num = nums.append
    prev_read_num = 0
    # Keep track of remapped read index
    counter = 1

    it = iter(remap_bam)
    # Get a list of reads that remapped correctly
    remap_read = next(it)
    while True:
        fields = remap_read.query_name.strip().split(":")
        read_num = int(fields[0])

        # For each original read, keep track of how many alternate reads were mapped
        if read_num != prev_read_num:

            # If none of the alternate versions of the read got mapped, set num to 1
            if read_num != counter:
                skipped = read_num - counter
                nums += ([1]*skipped)
                counter = read_num

            append_num(int(fields[-1]) * mult)

            prev_read_num = read_num
            counter += 1

        tid = remap_read.reference_id
        if tid != -1:
            pos = int(fields[rev_field] if remap_read.is_reverse else fields[2])
            if remap_read.reference_start == pos and refnames[tid] == fields[1]:
                # Throw out the remapped read if it remapped with a deletion
                ok = True
                for cig in remap_read.cigartuples or ():
                    op = cig[0]
                    if op != 0 and op != 3 and op != 4:
                        ok = False
                        break
                if ok:
                    append_correct(read_num)

        try:
            remap_read = next(it)
        except StopIteration:
            break
        except Exception:
            break
    return correct_maps, nums


def run(to_remap_bam, remap_bam, keep_bam, is_paired_end, threads=1,
        legacy=False):
    """Core function."""

    sys.stderr.write(
        "{} Starting run\n".format(
            time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))
        )
    )

    threads = max(1, threads)
    to_remap_bam = pysam.AlignmentFile(to_remap_bam, "rb", threads=threads)
    remap_bam = pysam.AlignmentFile(remap_bam, "rb", threads=threads)
    keep_bam = pysam.AlignmentFile(keep_bam, "wb", template=to_remap_bam,
                                   threads=threads)

    # List of correctly mapped reads, represented by read ID
    # (e.g. 42 in @42:chr1:14536:3)
    # List of how many alternate versions of each read there are
    # (e.g. 3 in @42:chr1:14536:3)
    correct_maps, nums = _remapped_ok(remap_bam, is_paired_end)

    correct_maps.sort()

    # Original aligned reads
    orig_read = next(to_remap_bam)

    # Number of different reads generated from the original read
    orig_num = int(nums[0])

    # Line number of the remap_bam file (if single end data) or read pair
    # number if paired end data.
    line_num = 1

    # Index for walking through correct_maps.
    map_indx = 0

    # Number of correctly mapped reads for the current read (pair).
    correct = 0

    # Total number of correctly mapped read (pairs).
    total_correct = 0

    keep_write = keep_bam.write
    n_correct_maps = len(correct_maps)
    last_correct = correct_maps[-1] if correct_maps else None

    end_of_file = False
    while (not end_of_file and
           (map_indx < n_correct_maps) and
           (line_num <= last_correct)):
        if line_num != correct_maps[map_indx]:

            # If we saw the correct number of remaps for the last read, keep it.
            if correct == orig_num:
                total_correct += 1
                keep_write(orig_read)

                # If the data is paired end, write out the paired read.
                if is_paired_end:
                    try:
                        orig_read = next(to_remap_bam)
                    except:
                        raise ValueError(
                            "File ended unexpectedly (no pair found)."
                        )
                    keep_write(orig_read)
            elif is_paired_end or legacy:
                # Skip the mate of the discarded read.  (The original code
                # also did this for single-end data, silently dropping the
                # read that follows every discarded read: --legacy.)
                try:
                    second_read = next(to_remap_bam)
                except:
                    end_of_file=True
                    break
            try:
                orig_read = next(to_remap_bam)
                orig_num = nums[line_num]
            except StopIteration:
                end_of_file = True
            line_num += 1
            correct = 0
        else:
            correct += 1
            map_indx += 1

    if correct == orig_num:
        total_correct += 1
        keep_write(orig_read)

        # If the data is paired end, write out the paired read.
        if is_paired_end:
            try:
                orig_read = next(to_remap_bam)
            except:
                sys.stderr.write("File ended unexpectedly (no pair found).")
                keep_bam.close()
                exit()
            keep_write(orig_read)
    keep_bam.close()
    sys.stderr.write(
        "{} Finished\n\n".format(
            time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))
        )
    )
    sys.stderr.write(
        ("RUN STATISTICS:\n\tTotal remapped read (pair)s: {}\n\t"
         "Read (pair)s remapped to the correct position: {} {:.2%}\n")
        .format(
            line_num, total_correct, total_correct/float(line_num)
        )
    )

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', dest='is_paired_end',
                        default=False, help=('Indicates that reads are '
                                             'paired-end (default is single).'))
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help='Threads for BAM (de)compression (default: 4)')
    parser.add_argument("--legacy", action='store_true', default=False,
                        help=('Use the original Hornet algorithm and reproduce '
                              'its output exactly, including its single-end '
                              'bug (the read after each discarded read was '
                              'skipped). Needs a remapped BAM sorted by name.'))
    h = ('to.remap.bam file from find_intersecting_snps.py.')
    parser.add_argument("to_remap_bam", help=h)
    parser.add_argument("remap_bam", help='BAM from remapping the fastqs (any order; sorted by name with --legacy)')
    parser.add_argument("keep_bam", help=(
        'File to write correctly remapped reads to.'
    ))

    options = parser.parse_args()

    if options.legacy:
        run(options.to_remap_bam, options.remap_bam, options.keep_bam,
            options.is_paired_end, options.threads, legacy=True)
    else:
        run_strict(options.to_remap_bam, options.remap_bam, options.keep_bam,
                   options.is_paired_end, options.threads)

if __name__ == '__main__':
    main()
