#!/usr/bin/env python3
"""filter_remapped_reads

This is a rewrite of the official WASP code, which no longer requires
a .num.gz file. This code is compatable with the rewritten version of
the first step (find_intersecting_snps), which no longer generates the
.num.gz file. The bam file of remapped reads must be sorted by read name.

High-performance version: identical output to the original Hornet code, with
a single parse of each read name, cached reference names, a C-level CIGAR
test and multi-threaded BAM (de)compression (``--threads``).

"""
from __future__ import print_function
import sys
import gzip
import time

import pysam


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


def run(to_remap_bam, remap_bam, keep_bam, is_paired_end, threads=1):
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
            else:
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
    h = ('to.remap.bam file from find_intersecting_snps.py.')
    parser.add_argument("to_remap_bam", help=h)
    parser.add_argument("remap_bam", help='Remapped bam file.')
    parser.add_argument("keep_bam", help=(
        'File to write correctly remapped reads to.'
    ))

    options = parser.parse_args()

    run(options.to_remap_bam, options.remap_bam, options.keep_bam,
        options.is_paired_end, options.threads)

if __name__ == '__main__':
    main()
