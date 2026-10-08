#!/usr/bin/env python3
""" find_intersecting_snps

Find the reads that overlap heterozygous SNPs and write, for each of them,
the versions of the read with the other allele(s) for remapping (step 2 of
the Hornet/WASP pipeline; see README.rst).

Outputs, next to the input BAM (<base> = input name without extension):
  <base>.keep.bam       reads that overlap no SNP (kept without remapping)
  <base>.to.remap.bam   original reads that must be remapped
  <base>.remap.fq.gz    allele-swapped reads to remap (.fq1.gz/.fq2.gz for -p)
  <base>.dropped.bam    reads overlapping indels / too many SNPs / with
                        inconsistent mates

Implementation notes: SNPs are kept in sorted arrays (~12 bytes per SNP)
and found by binary search over the aligned blocks of the CIGAR; BAM and
fastq (de)compression is multi-threaded (-t); -j processes the chromosomes of
a sorted, indexed BAM in parallel with an order-preserving merge (output
identical to a serial run). --legacy reproduces the original Hornet code
exactly. The C program c/hornet implements the same tool, with identical
output, several times faster.
"""
from __future__ import print_function
import sys
import argparse
import gzip
import os
import shutil
import subprocess
import tempfile
import time
import zlib
import struct
from array import array
from bisect import bisect_left
from collections import defaultdict, Counter
from glob import glob
from os import path
from pysam import AlignmentFile as Samfile

try:
    from functools import reduce
    from operator import mul
except ImportError as exc:
    # We better hope we're in Python 2.
    print(exc)

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

# Maximum number of allelic combinations (including the original) per read
# or read pair; set by --max-seqs (the original code used 32).
MAX_SEQS_PER_READ = 64
DRAW_PROGRESS = False
# Reproduce the original Hornet output exactly, including its bugs
# (--legacy). See README.rst, "Fixes".
LEGACY = False


def _now():
    return time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))


def product(iterable):
    "Returns the product of all items in the iterable"
    return reduce(mul, iterable, 1)


###############################################################################
# SNP loading
###############################################################################

class SNPTable(object):
    """Compact storage of the SNPs.

    For each chromosome:  ``pos[chrom]``  -- sorted array of 0-based positions
                          ``code[chrom]`` -- array of indices into ``alleles``
    ``alleles[code]`` is ``ref + alt`` exactly as the original snp_dict values,
    and ``is_indel[code]`` reproduces the original get_indels() test.

    When a position occurs more than once the last occurrence wins, matching
    the dict assignment semantics of the original code.
    """

    def __init__(self):
        self.pos = {}
        self.code = {}
        self.alleles = []
        self.is_indel = []
        self._allele_index = {}

    def allele_code(self, alleles):
        code = self._allele_index.get(alleles)
        if code is None:
            code = len(self.alleles)
            self._allele_index[alleles] = code
            self.alleles.append(alleles)
            if LEGACY:
                # Original: ('-' in alleles) or (max(len(i) for i in alleles)
                # > 1), where iterating a string yields single characters, so
                # only '-' marked an indel.
                self.is_indel.append(('-' in alleles) or
                                     (max(len(i) for i in alleles) > 1))
            else:
                # Like WASP: anything but a single-base ref and alt (e.g.
                # 'A AT' or 'C G,T') is treated as an indel/unsupported site.
                self.is_indel.append(('-' in alleles) or len(alleles) != 2)
        return code

    def __len__(self):
        return len(self.pos)

    def __contains__(self, chrom):
        return chrom in self.pos

    def chroms(self):
        return list(self.pos)

    def set_chrom(self, chrom, positions, codes):
        """positions/codes: numpy arrays in file order (may be unsorted and
        contain duplicates)."""
        if len(positions) == 0:
            return
        order = np.argsort(positions, kind='stable')
        positions = positions[order]
        codes = codes[order]
        last = np.ones(len(positions), dtype=bool)
        last[:-1] = positions[1:] != positions[:-1]
        positions = positions[last]
        codes = codes[last]
        if positions[0] >= -2**31 and positions[-1] < 2**31:
            parr = array('i')
            parr.frombytes(positions.astype(np.int32).tobytes())
        else:
            parr = array('q')
            parr.frombytes(positions.astype(np.int64).tobytes())
        carr = array('I')
        carr.frombytes(codes.astype(np.uint32).tobytes())
        self.pos[chrom] = parr
        self.code[chrom] = carr

    def to_dict(self):
        "The original dict-of-dicts representation (for tests / debugging)."
        snp_dict = defaultdict(dict)
        for chrom in self.pos:
            alleles = self.alleles
            snp_dict[chrom] = {p: alleles[c] for p, c in
                               zip(self.pos[chrom], self.code[chrom])}
        return snp_dict


def _parse_python(fname, has_chrom, chrom_only=None):
    """Pure python parser with exactly the original semantics.

    Returns {chrom: (positions list, allele string list)} for consolidated
    files, or (positions list, allele string list) for per-chrom files.
    """
    if has_chrom:
        out = defaultdict(lambda: ([], []))
        for line in gzip.open(fname, 'rt', encoding='ascii'):
            chrom, pos, ref, alt = line.split()
            if chrom_only is not None and chrom != chrom_only:
                continue
            ent = out[chrom]
            ent[0].append(int(pos) - 1)
            ent[1].append("".join([ref, alt]))
        return dict(out)
    poss, alls = [], []
    for line in gzip.open(fname, 'rt', encoding='ascii'):
        pos, ref, alt = line.split()
        poss.append(int(pos) - 1)
        alls.append("".join([ref, alt]))
    return poss, alls


def _parse_pandas(fname, has_chrom, chrom_only=None):
    """Fast parser. Returns the same structure as _parse_python but with
    numpy arrays: positions (int64) and (unique allele strings, codes)."""
    import csv
    import pandas as pd
    names = ['chrom', 'pos', 'ref', 'alt'] if has_chrom else ['pos', 'ref', 'alt']
    dtype = {'pos': np.int64, 'ref': 'category', 'alt': 'category'}
    if has_chrom:
        dtype['chrom'] = 'category'
    df = pd.read_csv(fname, sep=r'\s+', header=None, names=names,
                     dtype=dtype, compression='gzip', na_filter=False,
                     keep_default_na=False, quoting=csv.QUOTE_NONE,
                     engine='c', comment=None, skip_blank_lines=False,
                     encoding='ascii')
    for col in ('ref', 'alt') + (('chrom',) if has_chrom else ()):
        cats = df[col].cat.categories
        if len(cats) and (cats == '').any():
            raise ValueError("malformed line in %s" % fname)

    def pack(sub):
        rc = sub['ref'].cat.codes.to_numpy().astype(np.int64)
        ac = sub['alt'].cat.codes.to_numpy().astype(np.int64)
        rcats = sub['ref'].cat.categories
        acats = sub['alt'].cat.categories
        pair = rc * len(acats) + ac
        uniq, inv = np.unique(pair, return_inverse=True)
        strs = ["".join([rcats[u // len(acats)], acats[u % len(acats)]])
                for u in uniq]
        return (sub['pos'].to_numpy() - 1, strs, inv.astype(np.int64))

    if not has_chrom:
        return pack(df)
    out = {}
    ccodes = df['chrom'].cat.codes.to_numpy()
    for ci, chrom in enumerate(df['chrom'].cat.categories):
        if chrom_only is not None and chrom != chrom_only:
            continue
        sub = df[ccodes == ci]
        if len(sub):
            out[chrom] = pack(sub)
    return out


def _load_file(args):
    fname, has_chrom, chrom_only = args
    try:
        import pandas  # noqa: F401
        return _parse_pandas(fname, has_chrom, chrom_only)
    except ImportError:
        pass
    except Exception:  # fall back to reference parser to get its error
        pass
    res = _parse_python(fname, has_chrom, chrom_only)
    if has_chrom:
        return {c: _pack_python(*v) for c, v in res.items()}
    return _pack_python(*res)


def _pack_python(poss, alls):
    index = {}
    codes = np.empty(len(alls), dtype=np.int64)
    for i, a in enumerate(alls):
        c = index.get(a)
        if c is None:
            c = index[a] = len(index)
        codes[i] = c
    strs = [None] * len(index)
    for a, c in index.items():
        strs[c] = a
    return np.asarray(poss, dtype=np.int64), strs, codes


def _map_pool(func, jobs, processes):
    if processes > 1 and len(jobs) > 1:
        import multiprocessing as mp
        try:
            ctx = mp.get_context('fork')
        except ValueError:  # pragma: no cover
            ctx = mp.get_context()
        with ctx.Pool(min(processes, len(jobs))) as pool:
            return pool.map(func, jobs, chunksize=1)
    return [func(j) for j in jobs]


def get_snps(snpdir, chrom_only=None, processes=1):
    """Get SNPs from a single file or directory of files

    Returns a SNPTable (see above), equivalent to the original dictionary of
    dictionaries:
    snp_dict = {
                'chrom1' : {
                            pos1 : [ref, alt],
                            pos2 : [ref, alt],
                           },
                'chrom2' : {...},
                ...
                }
    where positions are 0-based
    """
    table = SNPTable()

    def add(chrom, parts):
        # parts: list of (positions, allele_strings, codes) in file order
        poss, codes = [], []
        for p, strs, c in parts:
            remap = np.array([table.allele_code(s) for s in strs],
                             dtype=np.int64)
            poss.append(p)
            codes.append(remap[c] if len(c) else c)
        table.set_chrom(chrom, np.concatenate(poss), np.concatenate(codes))

    if path.exists(path.join(snpdir, 'all.txt.gz')):
        print(_now(), ".... Loading snps from consolidated file.")
        res = _load_file((path.join(snpdir, 'all.txt.gz'), True, chrom_only))
        for chrom, part in res.items():
            add(chrom, [part])
        return table

    jobs, chroms = [], []
    for fname in glob(path.join(snpdir, '*.txt.gz')):
        chrom = path.basename(fname).split('.snps.txt.gz')[0]
        if chrom_only is not None and chrom != chrom_only:
            continue
        print(_now(), ".... Loading snps from", fname)
        jobs.append((fname, False, None))
        chroms.append(chrom)
    results = _map_pool(_load_file, jobs, processes)
    by_chrom = {}
    for chrom, part in zip(chroms, results):
        by_chrom.setdefault(chrom, []).append(part)
    for chrom, parts in by_chrom.items():
        add(chrom, parts)
    if len(table) == 0:
        print("ERROR: Could not find appropriate SNP files")
        sys.exit(1)
    return table


def get_snps_vcf(vcf, sample=None, chrom_only=None):
    """Heterozygous sites of one sample from a VCF/BCF file.

    SNVs (single-base alleles) are stored with the sample's two alleles.
    Every reference base covered by a heterozygous indel or other non-SNV
    allele is marked as an indel site, so that overlapping reads are dropped.
    Homozygous and missing genotypes are ignored.
    """
    import pysam
    verbosity = pysam.set_verbosity(0)   # no 'index not found' message
    vf = pysam.VariantFile(vcf)
    pysam.set_verbosity(verbosity)
    samples = list(vf.header.samples)
    if sample is None:
        if len(samples) != 1:
            sys.exit("ERROR: %s has %d samples; choose one with --sample"
                     % (vcf, len(samples)))
        sample = samples[0]
    elif sample not in samples:
        sys.exit("ERROR: sample %s not found in %s" % (sample, vcf))
    vf.subset_samples([sample])
    print(_now(), ".... Loading heterozygous sites of sample", sample, "from", vcf)
    table = SNPTable()
    snv = defaultdict(lambda: ([], []))
    indel = defaultdict(list)
    if chrom_only is not None:
        try:
            it = vf.fetch(chrom_only)
        except (ValueError, KeyError):
            it = (r for r in vf if r.chrom == chrom_only)
    else:
        it = vf
    VALID = set('ACGT')
    for rec in it:
        gt = rec.samples[0]['GT']
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None \
                or gt[0] == gt[1]:
            continue
        al = rec.alleles
        i, j = sorted(gt)
        a, b = al[i], al[j]
        if len(a) == 1 and len(b) == 1 and a in VALID and b in VALID:
            ent = snv[rec.chrom]
            ent[0].append(rec.start)
            ent[1].append(a + b)
        else:
            indel[rec.chrom].extend(range(rec.start, max(rec.stop, rec.start + 1)))
    indel_code = table.allele_code('-')
    for chrom in set(snv) | set(indel):
        pos, alls = snv.get(chrom, ([], []))
        codes = [table.allele_code(x) for x in alls]
        ipos = indel.get(chrom, [])
        # indel marks come last so that they win over SNVs at the same base
        table.set_chrom(chrom,
                        np.asarray(list(pos) + list(ipos), dtype=np.int64),
                        np.asarray(codes + [indel_code] * len(ipos),
                                   dtype=np.int64))
    if len(table) == 0:
        print("ERROR: no heterozygous sites found for sample", sample)
        sys.exit(1)
    return table


###############################################################################
# Read processing
###############################################################################

RC_TABLE = {
    ord('A'):ord('T'),
    ord('T'):ord('A'),
    ord('C'):ord('G'),
    ord('G'):ord('C'),
}

def reverse_complement(seq):
    "Reverse complements the string input"
    return seq.translate(RC_TABLE)[::-1]


_NO_SITES = ()


def find_sites(read, tab):
    """Variant sites hit by the aligned (M/=/X) bases of ``read``.

    Equivalent to scanning ``read.get_aligned_pairs(matches_only=True)`` and
    testing every reference position against the SNP table, but uses binary
    search over the CIGAR blocks.

    Returns (sites, n_aligned) where sites is a list of tuples
    (aligned_index, read_pos, ref_pos, allele_code) in alignment order.
    """
    posarr, codes = tab
    n = len(posarr)
    start = read.reference_start
    if start < 0:
        # get_aligned_pairs uses an unsigned position; nothing ever matches.
        return _NO_SITES, 0
    end = read.reference_end
    i = bisect_left(posarr, start)
    if i >= n:
        return _NO_SITES, 0
    if end is not None and posarr[i] >= end:
        return _NO_SITES, 0
    cigar = read.cigartuples
    if not cigar:
        return _NO_SITES, 0
    out = []
    pos = start
    qpos = 0
    aidx = 0
    for op, l in cigar:
        if op == 0 or op == 7 or op == 8:
            e = pos + l
            if i < n and posarr[i] < pos:
                i = bisect_left(posarr, pos, i)
            while i < n and posarr[i] < e:
                p = posarr[i]
                d = p - pos
                out.append((aidx + d, qpos + d, p, codes[i]))
                i += 1
            pos = e
            qpos += l
            aidx += l
        elif op == 1 or op == 4 or op == 6:
            qpos += l
        elif op == 2 or op == 3:
            pos += l
    return out, aidx


def get_read_seqs_sites(read, sites, n_aligned, table, dispositions):
    """ For each read, get all possible SNP substitutions

    for N biallelic snps in the read, will return 2^N reads

    Same semantics as the original get_read_seqs(), but only visits the
    aligned positions that are variant sites.  ``sites`` must be non-empty.
    """
    num_snps = 0
    seq0 = read.seq
    seqs = [seq0]
    alleles_of = table.alleles
    is_indel = table.is_indel
    prev = -2
    for aidx, read_pos, ref_pos, code in sites:
        if len(seqs) > MAX_SEQS_PER_READ:
            # The original checks the size at *every* aligned position, after
            # the indel test.  Was the check triggered at an earlier,
            # non-variant position?
            if aidx != prev + 1:
                dispositions['toss_manysnps'] += 1
                return []
            if is_indel[code]:
                dispositions['toss_indel'] += 1
                return []
            dispositions['toss_manysnps'] += 1
            return []
        if is_indel[code]:
            dispositions['toss_indel'] += 1
            return []
        prev = aidx

        dispositions['total_snps'] += 1
        alleles = alleles_of[code]
        read_base = seq0[read_pos]
        if read_base in alleles:
            dispositions['ref_match'] += 1
            num_snps += 1
            for new_allele in alleles:
                if new_allele == read_base:
                    continue
                for seq in list(seqs):
                    # Note that we make a copy up-front to avoid modifying
                    # the list we're iterating over
                    new_seq = seq[:read_pos] + new_allele + seq[read_pos+1:]
                    seqs.append(new_seq)
        else:
            dispositions['no_match'] += 1
    if len(seqs) > MAX_SEQS_PER_READ and prev + 1 < n_aligned:
        dispositions['toss_manysnps'] += 1
        return []
    if len(seqs) == 1:
        dispositions['no_snps'] += 1
    else:
        dispositions['has_snps'] += 1
    return seqs


def get_read_seqs_checked(read, sites, table, dispositions, phased=False):
    """Single-end allelic versions of a read (default, non-legacy mode).

    Reads overlapping an indel are dropped, then reads with more than
    MAX_SEQS_PER_READ allelic combinations. Unphased: every combination of
    the alleles (bases matching neither allele are left alone). Phased (-P):
    one version with every SNP switched to the other allele.
    """
    alleles_of = table.alleles
    is_indel = table.is_indel
    for site in sites:
        if is_indel[site[3]]:
            dispositions['toss_indel'] += 1
            return []
    if product(len(alleles_of[site[3]]) for site in sites) > MAX_SEQS_PER_READ:
        dispositions['toss_manysnps'] += 1
        return []
    seq0 = read.seq
    if phased:
        new = list(seq0)
        for _, read_pos, _, code in sites:
            alleles = alleles_of[code]
            dispositions['total_snps'] += 1
            k = alleles.find(seq0[read_pos])
            if k < 0:
                dispositions['no_match'] += 1
                dispositions['non_refalt_base'] += 1
                return []
            dispositions['ref_match'] += 1
            new[read_pos] = alleles[1 - k]
        dispositions['has_snps'] += 1
        return [seq0, ''.join(new)]
    seqs = [seq0]
    for _, read_pos, _, code in sites:
        dispositions['total_snps'] += 1
        alleles = alleles_of[code]
        read_base = seq0[read_pos]
        if read_base in alleles:
            dispositions['ref_match'] += 1
            for new_allele in alleles:
                if new_allele == read_base:
                    continue
                for seq in list(seqs):
                    seqs.append(seq[:read_pos] + new_allele + seq[read_pos+1:])
        else:
            dispositions['no_match'] += 1
    if len(seqs) == 1:
        dispositions['no_snps'] += 1
    else:
        dispositions['has_snps'] += 1
    return seqs


def _dual_read_seqs_core(snps, read_posns, seq1, seq2, dispositions, phased):
    """Body of the original get_dual_read_seqs() after the SNP lookup.

    Kept verbatim (including its variable shadowing in the unphased loop) so
    that the generated sequences and statistics are identical.
    """
    if product(len(i) for i in snps.values()) > MAX_SEQS_PER_READ:
        dispositions['toss_manysnps'] += 1
        return [[], []]

    if len(snps) == 0:
        dispositions['no_snps'] += 1
        return [[seq1], [seq2]]
    num_alleles = len(next(iter(snps.values())))
    if num_alleles > 2:
        # This happens if the SNP dict has multiple rows with the same position
        # so we just toss the read.
        dispositions['toss_manysnps'] += 1
        return [[], []]
        #  raise NotImplementedError("We can't yet do multiple phased genomes")

    if phased:
        reads1 = [[seq1], []]
        reads2 = [[seq2], []]
        last_pos1 = 0
        last_pos2 = 0
        for ref_pos in sorted(snps):
            pos1, pos2 = read_posns[ref_pos]
            alleles = snps[ref_pos]
            if not LEGACY:
                # The original never counted SNPs in phased mode (and then
                # divided by zero when printing the statistics).
                dispositions['total_snps'] += 1
                if ((pos1 is not None and seq1[pos1] not in alleles) or
                        (pos2 is not None and seq2[pos2] not in alleles)):
                    dispositions['no_match'] += 1
                    dispositions['non_refalt_base'] += 1
                    return [[], []]
                dispositions['ref_match'] += 1
            if pos1 is not None:
                is_alt = alleles.find(seq1[pos1])
                if is_alt == -1:
                    dispositions['non_refalt_base']
                    return [[], []]
                allele = alleles[1-is_alt]
                reads1[1].append(seq1[last_pos1:pos1])
                reads1[1].append(allele)
                last_pos1 = pos1 + 1

            if pos2 is not None:
                is_alt = alleles.find(seq2[pos2])
                if is_alt == -1:
                    dispositions['non_refalt_base']
                    return [[], []]
                allele = alleles[1-is_alt]
                reads2[1].append(seq2[last_pos2:pos2])
                reads2[1].append(allele)
                last_pos2 = pos2 + 1
            if pos1 is not None and pos2 is not None and seq1[pos1] != seq2[pos2]:
                dispositions['toss_anomalous_phase'] += 1
                return [[], []]

        reads1[1].append(seq1[last_pos1:])
        reads2[1].append(seq2[last_pos2:])
        if len(snps) == 0:
            dispositions['no_snps'] += 1
        else:
            dispositions['has_snps'] += 1
        return [[''.join(r) for r in reads1], [''.join(r) for r in reads2]]


    seqs1, seqs2 = [seq1], [seq2]


    for ref_pos in snps:
        match = False
        dispositions['total_snps'] += 1
        alleles = snps[ref_pos]
        pos1, pos2 = read_posns[ref_pos]
        new_seqs1 = []
        new_seqs2 = []
        if pos1 is None:
            if seq2[pos2] in alleles:
                dispositions['ref_match'] += 1
            else:
                dispositions['no_match'] += 1

            for allele in alleles:
                if allele == seq2[pos2]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(seq1)
                    new_seqs2.append(''.join([seq2[:pos2], allele, seq2[pos2+1:]]))

        elif pos2 is None:
            if seq1[pos1] in alleles:
                dispositions['ref_match'] += 1
            else:
                dispositions['no_match'] += 1

            for allele in alleles:
                if allele == seq1[pos1]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(''.join([seq1[:pos1], allele, seq1[pos1+1:]]))
                    new_seqs2.append(seq2)
        else:
            if seq1[pos1] != seq2[pos2]:
                dispositions['toss_anomalous_phase'] += 1
                return [[], []]
            if seq2[pos2] in alleles:
                dispositions['ref_match'] += 1
            else:
                dispositions['no_match'] += 1
            for allele in alleles:
                if allele == seq2[pos2]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(''.join([seq1[:pos1], allele, seq1[pos1+1:]]))
                    new_seqs2.append(''.join([seq2[:pos2], allele, seq2[pos2+1:]]))
        seqs1.extend(new_seqs1)
        seqs2.extend(new_seqs2)

    if len(seqs1) == 1:
        dispositions['no_snps'] += 1
    else:
        dispositions['has_snps'] += 1
    return seqs1, seqs2


def get_dual_read_seqs_sites(read1, read2, sites1, sites2, table,
                             dispositions, phased=False):
    """ For each pair of reads, get all concordant SNP substitutions

    Note that if the reads overlap, the matching positions in read1 and read2
    will get the same subsitution as each other.

    ``sites1``/``sites2`` come from find_sites(); at least one is non-empty.
    """
    is_indel = table.is_indel
    for site in sites1:
        if is_indel[site[3]]:
            dispositions['toss_indel'] += 1
            return [[], []]
    for site in sites2:
        if is_indel[site[3]]:
            dispositions['toss_indel'] += 1
            return [[], []]
    alleles_of = table.alleles
    snps = {}
    read_posns = defaultdict(lambda: [None, None])
    for _, read_pos1, ref_pos, code in sites1:
        snps[ref_pos] = alleles_of[code]
        read_posns[ref_pos][0] = read_pos1
    for _, read_pos2, ref_pos, code in sites2:
        snps[ref_pos] = alleles_of[code]
        read_posns[ref_pos][1] = read_pos2
    return _dual_read_seqs_core(snps, read_posns, read1.seq, read2.seq,
                                dispositions, phased)


###############################################################################
# Output
###############################################################################

class GzipTextWriter(object):
    """Buffered gzip writer; uses pigz when available."""

    def __init__(self, fname, threads=1, level=6, bufsize=1 << 20):
        self._buf = []
        self._size = 0
        self._bufsize = bufsize
        self._proc = None
        pigz = shutil.which('pigz')
        if pigz:
            self._fh = open(fname, 'wb')
            self._proc = subprocess.Popen(
                [pigz, '-c', '-%d' % level, '-p', str(max(1, threads))],
                stdin=subprocess.PIPE, stdout=self._fh)
            self._out = self._proc.stdin
        else:
            self._fh = None
            self._out = gzip.open(fname, 'wb', compresslevel=level)

    def write(self, text):
        self._buf.append(text)
        self._size += len(text)
        if self._size >= self._bufsize:
            self.flush()

    def flush(self):
        if self._buf:
            self._out.write(''.join(self._buf).encode('ascii'))
            self._buf = []
            self._size = 0

    def close(self):
        self.flush()
        self._out.close()
        if self._proc is not None:
            if self._proc.wait() != 0:
                raise IOError("pigz failed")
            self._fh.close()


def write_read_seqs(both_read_seqs, keep, remap_bam, fastqs, dropped=None, remap_num=0):
    """Write the given reads out to the appropriate file


    If there are no SNPs in the read, write it to the BAM file `keep`

    If there are too many SNPs, and `dropped` is provided, write the original
    read out to `dropped`

    Otherwise, output all possible substitutions to the fastqs for remapping,
    as well as a bam file containing the original read.
    """
    reads, seqs = zip(*both_read_seqs)
    assert len(reads) == len(fastqs)
    num_seqs = len(both_read_seqs[0][1])
    if num_seqs == 0:
        if dropped is not None:
            for read in reads:
                dropped.write(read)
            return 0
        else:
            return 0
    elif num_seqs == 1:
        for read, seqs in both_read_seqs:
            keep.write(read)
    else:
        assert len(reads) > 0
        for read in reads:
            remap_bam.write(read)
        left_pos = min(r.pos for r in reads)
        right_pos = max(r.pos for r in reads)
        loc_line = '{}:{}:{}:{}:{}'.format(
            remap_num,
            read.reference_name,
            left_pos,
            right_pos,
            len(seqs[0])-1,
        )
        rev = [read.is_reverse for read in reads]
        if LEGACY:
            quals = [read.qual for read in reads]
        else:
            # The sequence of a reverse-strand read is reverse complemented
            # back to the sequenced orientation; its qualities must be
            # reversed too.
            quals = [(read.qual[::-1] if (is_rev and read.qual is not None)
                      else read.qual) for read, is_rev in zip(reads, rev)]

        out = []
        first = True
        # Some python fanciness to deal with single or paired end reads (or
        # n-ended reads, if such technology ever happens.
        for read_seqs in zip(*seqs):
            if first:
                first = False
                continue
            for seq, is_rev, qual, fastq in zip(read_seqs, rev, quals, fastqs):
                fastq.write(
                    "@{loc_line}\n{seq}\n+{loc_line}\n{qual}\n"
                    .format(
                        loc_line=loc_line,
                        seq=reverse_complement(seq) if is_rev else seq,
                        qual=qual)
                    )
        return 1
    return 0


###############################################################################
# Main loop
###############################################################################

# Reads dropped before anything else (default mode): unmapped, secondary,
# QC-fail, supplementary; for pairs also reads whose mate is unmapped.
FLAG_DROP_SINGLE = 0x4 | 0x100 | 0x200 | 0x800
FLAG_DROP_PAIRED = FLAG_DROP_SINGLE | 0x8


class PlainTextWriter(object):
    """Uncompressed buffered writer with the GzipTextWriter interface."""

    def __init__(self, fname, bufsize=1 << 20):
        self._out = open(fname, 'wb')
        self._buf = []
        self._size = 0
        self._bufsize = bufsize

    write = GzipTextWriter.write
    flush = GzipTextWriter.flush

    def close(self):
        self.flush()
        self._out.close()


class _CrossChromosome(Exception):
    pass


def process_reads(insam, table, out_prefix, is_paired=True, phased=False,
                  keep_only=None, threads=1, region=None, worker=False):
    """Classify reads and write the outputs with the given prefix.

    Returns a dict with the read_results Counter and bookkeeping used by the
    parallel mode.

    With ``region`` only the reads of one reference (or '*' for the reads
    without coordinates) are processed. With ``worker`` the fastqs are written
    uncompressed, and the information needed to prove that the per-chromosome
    split cannot change the pairing of reads is collected.
    """
    chroms = insam.references
    bam_threads = max(1, threads)
    keep = Samfile(out_prefix + 'keep.bam', 'wb', template=insam,
                   threads=bam_threads)
    remap_bam = Samfile(out_prefix + 'to.remap.bam', 'wb', template=insam,
                        threads=bam_threads)
    dropped_bam = Samfile(out_prefix + 'dropped.bam', 'wb', template=insam,
                          threads=bam_threads)
    names = ['remap.fq1.gz', 'remap.fq2.gz'] if is_paired else ['remap.fq.gz']
    if worker:
        fastqs = [PlainTextWriter(out_prefix + n) for n in names]
    else:
        fastqs = [GzipTextWriter(out_prefix + n, threads) for n in names]

    # Per-tid SNP tables
    empty = (array('i'), array('I'))
    tabs = [((table.pos[c], table.code[c]) if c in table else empty)
            for c in chroms]
    nonempty = [len(t[0]) > 0 for t in tabs]

    unpaired_reads = [{}, {}]
    stored = array('q')
    stored_append = stored.append
    read_results = Counter()
    remap_num = 1
    n_total = n_skipped = n_nosnps = n_not_proper = n_unmapped = n_flagged = 0
    keep_write = keep.write
    dropped_write = dropped_bam.write

    if region is None:
        iterator = iter(insam)
        if DRAW_PROGRESS:
            try:
                from tqdm import tqdm
                iterator = tqdm(iterator, total=insam.mapped, mininterval=5,
                                unit_scale=True, unit='Rd')
            except Exception:
                pass
    else:
        iterator = insam.fetch(region)

    try:
        for read in iterator:
            tid = read.reference_id
            if keep_only is not None and chroms[tid] != keep_only:
                n_skipped += 1
                continue
            n_total += 1
            if not is_paired:
                if LEGACY:
                    if phased:
                        raise NotImplementedError("Not yet implemented. It should work with paired end reads, though")
                elif read.flag & FLAG_DROP_SINGLE:
                    n_flagged += 1
                    continue
                if tid >= 0 and nonempty[tid]:
                    sites, n_al = find_sites(read, tabs[tid])
                else:
                    sites = _NO_SITES
                if not sites:
                    n_nosnps += 1
                    keep_write(read)
                    continue
                if LEGACY:
                    read_seqs = get_read_seqs_sites(read, sites, n_al, table,
                                                    read_results)
                else:
                    read_seqs = get_read_seqs_checked(read, sites, table,
                                                      read_results, phased)
                remap_num += write_read_seqs([(read, read_seqs)], keep, remap_bam,
                                             fastqs, dropped_bam, remap_num)
                continue

            flag = read.flag
            if not LEGACY:
                if flag & FLAG_DROP_PAIRED:
                    n_flagged += 1
                    continue
                if not flag & 2 or read.next_reference_id != tid:
                    n_not_proper += 1
                    continue
            if flag & 2:  # proper pair
                slot_self = 1 if flag & 128 else 0 # 0 if is_read1, 1 if is_read2
                slot_other = 1 if flag & 64 else 0
                qname = read.query_name
                if worker and read.next_reference_id != tid:
                    raise _CrossChromosome()
                other = unpaired_reads[slot_other].pop(qname, None)
                if other is None:
                    unpaired_reads[slot_self][qname] = read
                    if worker:
                        stored_append(hash(qname))
                    continue
                both_reads = [None, None]
                both_reads[slot_self] = read
                both_reads[slot_other] = other
                read1, read2 = both_reads
                if read1.flag & 4 or read2.flag & 4:
                    n_unmapped += 1
                    dropped_write(read1)
                    dropped_write(read2)
                    continue
                tid1 = read1.reference_id
                if tid1 >= 0 and nonempty[tid1]:
                    tab = tabs[tid1]
                    sites1, _ = find_sites(read1, tab)
                    sites2, _ = find_sites(read2, tab)
                else:
                    sites1 = sites2 = _NO_SITES
                if not sites1 and not sites2:
                    n_nosnps += 1
                    keep_write(read1)
                    keep_write(read2)
                    continue
                both_seqs = get_dual_read_seqs_sites(read1, read2, sites1, sites2,
                                                     table, read_results, phased)
                remap_num += write_read_seqs([(read1, both_seqs[0]),
                                              (read2, both_seqs[1])],
                                             keep, remap_bam, fastqs, dropped_bam,
                                             remap_num)
            else:
                n_not_proper += 1
                # Most tools assume reads are paired and do not check IDs. Drop it out.
                continue

    finally:
        keep.close()
        remap_bam.close()
        dropped_bam.close()
        for fq in fastqs:
            fq.close()

    read_results['total'] += n_total
    read_results['skipped'] += n_skipped
    read_results['no_snps'] += n_nosnps
    read_results['not_proper_pair'] += n_not_proper
    read_results['unmapped read'] += n_unmapped
    read_results['flagged'] += n_flagged

    out = {
        'results': read_results,
        'n_unpaired': len(unpaired_reads[0]) + len(unpaired_reads[1]),
        'n_remap': remap_num - 1,
        'n_seen': n_total + n_skipped,
    }
    if worker:
        left = array('q', [hash(q) for d in unpaired_reads for q in d])
        out['leftover'] = left
        out['stored'] = stored
    return out


def _pct(num, den):
    try:
        return "%.2f" % ((num / den)*100)
    except ZeroDivisionError:
        return "nan"


def print_stats(read_results, n_unpaired, is_paired, keep_only, phased=False):
    print()
    print(_now(), ".... Finished!")
    print()
    print("RUN STATISTICS:")

    if is_paired:
        total_pairs = read_results['total']//2
        print("  Total input reads:", total_pairs, "pairs.")
        print("  Unpaired reads:", n_unpaired, "(" + \
            _pct(n_unpaired, total_pairs) + "%)")
    else:
        total_pairs = read_results['total']
        print("  Total input reads:", total_pairs)
    if keep_only is not None:
        print("  Reads on skipped chromosomes", read_results['skipped'])
    if not LEGACY:
        print("  Reads dropped [unmapped/secondary/supplementary/QC-fail"
              + ("/mate unmapped" if is_paired else "") + "]:",
              read_results['flagged'], "(" + _pct(read_results['flagged'],
                                                  total_pairs) + "%)")

    print("  Reads with no SNPs:", read_results['no_snps'], "(" + _pct(read_results\
        ['no_snps'], total_pairs) + "%)")
    print("  Reads overlapping SNPs:", read_results['has_snps'], "(" + _pct(read_results\
        ['has_snps'], total_pairs) + "%)")

    total_snps = read_results['total_snps']
    print("  Total SNPs covered:", total_snps)
    print("\tReference SNP matches:", read_results['ref_match'], "(" + _pct(read_results\
            ['ref_match'], total_snps) + "%)")
    print("\tNon-reference SNP matches:", read_results['no_match'], "(" + _pct(read_results\
            ['no_match'], total_snps) + "%)")

    print("  Reads dropped [INDEL]:", read_results['toss_indel'], "(" + _pct(read_results\
        ['toss_indel'], total_pairs) + "%)")
    print("  Reads dropped [too many SNPs]:", read_results['toss_manysnps'], "(" + \
        _pct(read_results['toss_manysnps'], total_pairs) + "%)")

    if is_paired:
        print("  Reads dropped [anomalous pair]:", read_results['toss_anomalous_phase'], "(" + \
            _pct(read_results['toss_anomalous_phase'], total_pairs) + "%)")
    if phased and not LEGACY:
        print("  Reads dropped [non ref/alt base]:", read_results['non_refalt_base'], "(" + \
            _pct(read_results['non_refalt_base'], total_pairs) + "%)")
    if is_paired:
        print("  Reads dropped [not proper pair]:", read_results['not_proper_pair'], "(" + \
            _pct(read_results['not_proper_pair'], total_pairs) + "%)")

    print()
    sys.stdout.flush()


def _basename(insam):
    fname = insam.filename
    if isinstance(fname, bytes):
        fname = fname.decode()
    return fname.rsplit('.', 1)[0]


###############################################################################
# Per-chromosome parallel mode
###############################################################################

_BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')


def _bgzf_blocks(fh):
    """Yield raw BGZF blocks (gzip header included)."""
    while True:
        head = fh.read(12)
        if not head:
            return
        if len(head) < 12 or head[:4] != b'\x1f\x8b\x08\x04':
            raise IOError("not a BGZF file")
        xlen = struct.unpack('<H', head[10:12])[0]
        extra = fh.read(xlen)
        bsize = None
        j = 0
        while j + 4 <= len(extra):
            slen = struct.unpack('<H', extra[j+2:j+4])[0]
            if extra[j] == 66 and extra[j+1] == 67 and slen == 2:
                bsize = struct.unpack('<H', extra[j+4:j+6])[0]
            j += 4 + slen
        if bsize is None:
            raise IOError("not a BGZF file")
        body = fh.read(bsize + 1 - 12 - xlen)  # cdata + crc32 + isize
        yield head + extra + body


def _bgzf_block_data(block):
    xlen = struct.unpack('<H', block[10:12])[0]
    return zlib.decompress(block[12 + xlen:-8], -15)


def _bam_header_len(data):
    "Length of the binary BAM header in data, or None if incomplete."
    if len(data) < 8:
        return None
    o = 8 + struct.unpack('<i', data[4:8])[0]
    if len(data) < o + 4:
        return None
    n_ref = struct.unpack('<i', data[o:o+4])[0]
    o += 4
    for _ in range(n_ref):
        if len(data) < o + 4:
            return None
        o += 4 + struct.unpack('<i', data[o:o+4])[0] + 4
    return o if len(data) >= o else None


def bam_cat(out_name, in_names):
    """Concatenate BAM files that share a header (like ``samtools cat``)."""
    with open(out_name, 'wb') as out:
        for k, name in enumerate(in_names):
            with open(name, 'rb') as fh:
                blocks = _bgzf_blocks(fh)
                if k > 0:
                    # Skip the header; htslib flushes the BGZF stream after
                    # writing it, so it occupies whole blocks.
                    data = b''
                    need = None
                    while need is None:
                        data += _bgzf_block_data(next(blocks))
                        need = _bam_header_len(data)
                    if len(data) != need:
                        raise IOError("BAM header does not end on a block "
                                      "boundary: " + name)
                for block in blocks:
                    if block == _BGZF_EOF:
                        continue
                    out.write(block)
        out.write(_BGZF_EOF)


_WORKER_STATE = {}


def _worker(job):
    """Process one reference in a forked worker."""
    region, prefix = job
    st = _WORKER_STATE
    insam = Samfile(st['infile'])
    try:
        res = process_reads(insam, st['table'], prefix, st['is_paired'],
                            st['phased'], st['keep_only'], 1, region=region,
                            worker=True)
    except _CrossChromosome:
        return None
    stored = res.pop('stored')
    with open(prefix + 'stored.bin', 'wb') as fh:
        stored.tofile(fh)
    return res


def _renumber_job(job):
    """Compress a worker fastq into a gzip member, shifting the remap IDs."""
    src, dst, offset, level = job
    with open(src, 'rb') as fh, open(dst, 'wb') as raw:
        comp = zlib.compressobj(level, zlib.DEFLATED, 31)
        out = []
        size = 0
        for i, line in enumerate(fh):
            if offset and i % 2 == 0:
                # '@N:...' (header) or '+N:...' (separator)
                head, rest = line[1:].split(b':', 1)
                line = b'%c%d:%s' % (line[0], int(head) + offset, rest)
            out.append(line)
            size += len(line)
            if size > (1 << 22):
                raw.write(comp.compress(b''.join(out)))
                out = []
                size = 0
        raw.write(comp.compress(b''.join(out)))
        raw.write(comp.flush())
    os.unlink(src)


def assign_reads_parallel(infile, table, is_paired, phased, keep_only,
                          processes, threads, basename):
    """Per-chromosome parallel processing of a sorted & indexed BAM.

    Returns (Counter, n_unpaired), or None if the input is not suitable or the
    split could change the result (the caller then runs serially).
    """
    try:
        insam = Samfile(infile)
        if not insam.has_index():
            return None
        if insam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            return None
        expected = {}
        for s in insam.get_index_statistics():
            if s.total > 0:
                expected[insam.get_tid(s.contig)] = s.total
        n_nocoor = insam.nocoordinate
    except Exception:
        return None
    refs = insam.references
    regions = [(refs[t], expected[t]) for t in sorted(expected)]
    if n_nocoor:
        regions.append(('*', n_nocoor))
    if len(regions) < 2:
        return None

    tmpdir = tempfile.mkdtemp(prefix='.hornet_tmp_',
                              dir=path.dirname(path.abspath(basename)))
    try:
        jobs = [(reg, path.join(tmpdir, '%d.' % k))
                for k, (reg, _) in enumerate(regions)]
        _WORKER_STATE.update(infile=infile, table=table, is_paired=is_paired,
                             phased=phased, keep_only=keep_only)
        import multiprocessing as mp
        ctx = mp.get_context('fork')
        with ctx.Pool(min(processes, len(jobs))) as pool:
            # Largest references first for load balancing.
            order = sorted(range(len(jobs)), key=lambda k: -regions[k][1])
            pending = {k: pool.apply_async(_worker, (jobs[k],)) for k in order}
            results = [pending[k].get() for k in range(len(jobs))]
            _WORKER_STATE.clear()

            if any(r is None for r in results):
                sys.stderr.write("find_intersecting_snps: proper pair with a "
                                 "mate on another chromosome; running "
                                 "serially\n")
                return None
            for (reg, n_exp), res in zip(regions, results):
                if res['n_seen'] != n_exp:
                    sys.stderr.write("find_intersecting_snps: index read "
                                     "counts do not match; running "
                                     "serially\n")
                    return None
            if is_paired:
                # A serial scan pairs reads by name across chromosomes. Make
                # sure no read left unpaired on one reference shares its name
                # with a read stored on a later reference.
                seen_left = np.empty(0, dtype=np.int64)
                for k, res in enumerate(results):
                    if len(seen_left):
                        stored = np.fromfile(jobs[k][1] + 'stored.bin',
                                             dtype=np.int64)
                        if np.isin(stored, seen_left).any():
                            sys.stderr.write(
                                "find_intersecting_snps: read names shared "
                                "between chromosomes; running serially\n")
                            return None
                    left = np.frombuffer(res['leftover'], dtype=np.int64)
                    seen_left = np.union1d(seen_left, left)

            # Stitch the outputs together in file order.
            fq_names = (['remap.fq1.gz', 'remap.fq2.gz'] if is_paired
                        else ['remap.fq.gz'])
            fq_jobs = []
            for name in fq_names:
                offset = 0
                for k, res in enumerate(results):
                    src = jobs[k][1] + name
                    fq_jobs.append((src, src + '.gz', offset, 6))
                    offset += res['n_remap']
            pool.map(_renumber_job, fq_jobs, chunksize=1)

        for name in ['keep.bam', 'to.remap.bam', 'dropped.bam']:
            bam_cat('.'.join([basename, name]),
                    [jobs[k][1] + name for k in range(len(jobs))])
        for name in fq_names:
            with open('.'.join([basename, name]), 'wb') as out:
                for k in range(len(jobs)):
                    with open(jobs[k][1] + name + '.gz', 'rb') as fh:
                        shutil.copyfileobj(fh, out, 1 << 22)

        total = Counter()
        for res in results:
            total.update(res['results'])
        n_unpaired = sum(res['n_unpaired'] for res in results)
        return total, n_unpaired
    finally:
        _WORKER_STATE.clear()
        shutil.rmtree(tmpdir, ignore_errors=True)


def assign_reads(insam, snp_dict, indel_dict=None, is_paired=True, phased=False,
                 keep_only=None, threads=1, processes=1):
    """ Loop through all the reads in insam and output them to the appropriate file


    """
    basename = _basename(insam)
    result = None
    if processes > 1 and not (phased and not is_paired):
        fname = insam.filename
        if isinstance(fname, bytes):
            fname = fname.decode()
        result = assign_reads_parallel(fname, snp_dict, is_paired, phased,
                                       keep_only, processes, threads,
                                       basename)
    if result is None:
        res = process_reads(insam, snp_dict, basename + '.', is_paired,
                            phased, keep_only, threads)
        result = res['results'], res['n_unpaired']
    print_stats(result[0], result[1], is_paired, keep_only, phased)



def get_indels(snp_dict):
    """Kept for API compatibility; indel status lives in the SNPTable."""
    return None



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paired_end",
                        action='store_true',
                        dest='is_paired_end', default=False,
                        help=('Indicates that reads are '
                              'paired-end (default is single).'))
    parser.add_argument('-C', '--limit-to-chrom',
                        default=None,
                        help="Limit loading of SNPs to the specified chromosome")

    parser.add_argument("-s", "--sorted",
                        action='store_true', dest='is_sorted', default=False,
                        help=('Indicates that the input bam file'
                              ' is coordinate sorted (default is False)'))
    parser.add_argument('-P', "--phased-snp-data",
                       dest='is_phased', action='store_true', default=False,
                       help=('SNP data is phased---useful if there is a '
                             + 'defined "reference" and "alternate", such as '
                             + 'in  a hybrid. If enabled, will only generate '
                             + '2 reads per input read, rather than 2^N_snps.')
                      )
    parser.add_argument('--progressbar', action='store_true', default=False,
                        help='Show a progress bar if possible')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help=('Threads for BAM/fastq (de)compression '
                              '(default: 4)'))
    parser.add_argument('--vcf', default=None,
                        help=('VCF/BCF file with genotypes; the heterozygous '
                              'sites of --sample are used instead of snp_dir.'))
    parser.add_argument('--sample', default=None,
                        help='Sample in --vcf (default: the only sample)')
    parser.add_argument('--max-seqs', type=int, default=64,
                        help=('Drop reads (pairs) with more than this many '
                              'allelic combinations, i.e. more than '
                              'log2(N) SNPs (default: 64)'))
    parser.add_argument('--legacy', action='store_true', default=False,
                        help=('Reproduce the output of the original Hornet '
                              'code exactly, including its bugs (phased-mode '
                              'statistics, qualities of reverse-strand reads '
                              'in the remap fastqs), with no read flag '
                              'filtering and --max-seqs 32.'))
    parser.add_argument('-j', '--processes', type=int, default=1,
                        help=('Process chromosomes in parallel with this many '
                              'worker processes. Requires a coordinate sorted '
                              'and indexed BAM; otherwise falls back to a '
                              'serial run. Output is identical either way. '
                              'Also used to load SNP files in parallel. '
                              '(default: 1)'))


    parser.add_argument("infile", help="BAM/SAM/CRAM file from the first mapping")
    snp_dir_help = ('Directory containing the SNPs segregating within the '
                    'sample in question (which need to be checked for '
                    'mappability issues).  This directory should contain '
                    'sorted files of SNPs separated by chromosome and named: '
                    '<chrname>.snps.txt.gz (where chrname is the name of the '
                    'chromosomes in the bam file). These files should contain '
                    '3 columns: position RefAllele AltAllele. Not needed '
                    'with --vcf.')

    parser.add_argument("snp_dir", nargs='?', default=None, help=snp_dir_help)

    options = parser.parse_args()
    DRAW_PROGRESS = False
    if options.progressbar:
        DRAW_PROGRESS = True
    LEGACY = options.legacy
    if (options.vcf is None) == (options.snp_dir is None):
        parser.error("give either snp_dir or --vcf")
    if LEGACY:
        if options.vcf or options.max_seqs != 64:
            parser.error("--legacy cannot be combined with --vcf or --max-seqs")
        MAX_SEQS_PER_READ = 32
    else:
        MAX_SEQS_PER_READ = options.max_seqs

    INSAM = Samfile(options.infile, threads=max(1, options.threads))
    if options.vcf:
        SNP_DICT = get_snps_vcf(options.vcf, options.sample,
                                options.limit_to_chrom)
    else:
        SNP_DICT = get_snps(options.snp_dir, options.limit_to_chrom,
                            options.processes)

    print(_now(), ".... Done with SNPs.")
    sys.stdout.flush()

    assign_reads(INSAM, SNP_DICT, None, options.is_paired_end,
                 options.is_phased, options.limit_to_chrom,
                 options.threads, options.processes)
