#!/usr/bin/env python3
"""Slow, simple reference implementation of find_intersecting_snps (default
mode), used to check the optimised Python and C versions.

usage: ref_find.py [-p] [-P] [-C chrom] [--max-seqs N]
                   (--vcf FILE [--sample NAME] | --snp-dir DIR) in.bam out_prefix

Writes out_prefix.{keep,to.remap,dropped}.sam and out_prefix.remap.fq
(single-end) or out_prefix.remap.fq1/.fq2 (paired-end), uncompressed.
"""
import argparse
import glob
import gzip
import os
import pysam

RC = str.maketrans('ACGT', 'TGCA')


def load_text(d, chrom_only):
    snps = {}
    allf = os.path.join(d, 'all.txt.gz')
    if os.path.exists(allf):
        for line in gzip.open(allf, 'rt'):
            c, p, r, a = line.split()
            if chrom_only is None or c == chrom_only:
                snps.setdefault(c, {})[int(p) - 1] = r + a
        return snps
    for f in glob.glob(os.path.join(d, '*.txt.gz')):
        c = os.path.basename(f).split('.snps.txt.gz')[0]
        if chrom_only is not None and c != chrom_only:
            continue
        for line in gzip.open(f, 'rt'):
            p, r, a = line.split()
            snps.setdefault(c, {})[int(p) - 1] = r + a
    return snps


def load_vcf(f, sample, chrom_only):
    verbosity = pysam.set_verbosity(0)   # no 'index not found' message
    vf = pysam.VariantFile(f)
    pysam.set_verbosity(verbosity)
    if sample is None:
        sample = list(vf.header.samples)[0]
    snv, indel = {}, {}
    for rec in vf:
        if chrom_only is not None and rec.chrom != chrom_only:
            continue
        gt = rec.samples[sample]['GT']
        if gt is None or len(gt) != 2 or None in gt or gt[0] == gt[1]:
            continue
        a, b = (rec.alleles[i] for i in sorted(gt))
        if len(a) == 1 and len(b) == 1 and a in 'ACGT' and b in 'ACGT':
            snv.setdefault(rec.chrom, {})[rec.start] = a + b
        else:
            for p in range(rec.start, max(rec.stop, rec.start + 1)):
                indel.setdefault(rec.chrom, set()).add(p)
    for c, ps in indel.items():
        for p in ps:
            snv.setdefault(c, {})[p] = '-'
    return snv


def is_indel(alleles):
    return '-' in alleles or len(alleles) != 2


def product(xs):
    r = 1
    for x in xs:
        r *= x
    return r


def sites_of(read, table):
    return [(q, p, table[p]) for q, p in read.get_aligned_pairs(matches_only=True)
            if p in table]


def single(read, sites, phased, max_seqs):
    if any(is_indel(a) for _, _, a in sites):
        return []
    if product(len(a) for _, _, a in sites) > max_seqs:
        return []
    seq = read.query_sequence
    if phased:
        new = list(seq)
        for q, _, a in sites:
            k = a.find(seq[q])
            if k < 0:
                return []
            new[q] = a[1 - k]
        return [seq, ''.join(new)]
    seqs = [seq]
    for q, _, a in sites:
        b = seq[q]
        if b in a:
            for na in a:
                if na != b:
                    seqs = seqs + [s[:q] + na + s[q + 1:] for s in seqs]
    return seqs


def dual(r1, r2, s1, s2, phased, max_seqs):
    if any(is_indel(a) for _, _, a in s1 + s2):
        return [[], []]
    snps, q1s, q2s = {}, {}, {}
    for q, p, a in s1:
        snps[p] = a
        q1s[p] = q
    for q, p, a in s2:
        snps[p] = a
        q2s[p] = q
    if product(len(a) for a in snps.values()) > max_seqs:
        return [[], []]
    seq1, seq2 = r1.query_sequence, r2.query_sequence
    if phased:
        n1, n2 = list(seq1), list(seq2)
        for p in sorted(snps):
            a = snps[p]
            q1, q2 = q1s.get(p), q2s.get(p)
            for q, seq, new in ((q1, seq1, n1), (q2, seq2, n2)):
                if q is not None:
                    k = a.find(seq[q])
                    if k < 0:
                        return [[], []]
                    new[q] = a[1 - k]
            if q1 is not None and q2 is not None and seq1[q1] != seq2[q2]:
                return [[], []]
        return [[seq1, ''.join(n1)], [seq2, ''.join(n2)]]
    S1, S2 = [seq1], [seq2]
    for p, a in snps.items():
        q1, q2 = q1s.get(p), q2s.get(p)
        if q1 is not None and q2 is not None and seq1[q1] != seq2[q2]:
            return [[], []]
        base = seq2[q2] if q2 is not None else seq1[q1]
        new1, new2 = [], []
        for na in a:
            if na == base:
                continue
            for x, y in zip(S1, S2):
                new1.append(x if q1 is None else x[:q1] + na + x[q1 + 1:])
                new2.append(y if q2 is None else y[:q2] + na + y[q2 + 1:])
        S1 += new1
        S2 += new2
    return [S1, S2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-p', action='store_true')
    ap.add_argument('-P', action='store_true')
    ap.add_argument('-C', default=None)
    ap.add_argument('--max-seqs', type=int, default=64)
    ap.add_argument('--vcf')
    ap.add_argument('--sample')
    ap.add_argument('--snp-dir')
    ap.add_argument('bam')
    ap.add_argument('out')
    o = ap.parse_args()
    snps = load_vcf(o.vcf, o.sample, o.C) if o.vcf else load_text(o.snp_dir, o.C)
    bam = pysam.AlignmentFile(o.bam)
    refs = bam.references
    keep = pysam.AlignmentFile(o.out + '.keep.sam', 'w', template=bam)
    remap = pysam.AlignmentFile(o.out + '.to.remap.sam', 'w', template=bam)
    dropped = pysam.AlignmentFile(o.out + '.dropped.sam', 'w', template=bam)
    if o.p:
        fqs = [open(o.out + '.remap.fq1', 'w'), open(o.out + '.remap.fq2', 'w')]
    else:
        fqs = [open(o.out + '.remap.fq', 'w')]
    state = {'num': 1}

    def write(reads, seqs):
        n = len(seqs[0])
        if n == 0:
            for r in reads:
                dropped.write(r)
        elif n == 1:
            for r in reads:
                keep.write(r)
        else:
            for r in reads:
                remap.write(r)
            loc = '%d:%s:%d:%d:%d' % (state['num'], reads[-1].reference_name,
                                      min(r.reference_start for r in reads),
                                      max(r.reference_start for r in reads), n - 1)
            for i in range(1, n):
                for r, s, fq in zip(reads, seqs, fqs):
                    seq, qual = s[i], r.qual
                    if r.is_reverse:
                        seq = seq.translate(RC)[::-1]
                        if qual is not None:
                            qual = qual[::-1]
                    fq.write('@%s\n%s\n+%s\n%s\n' % (loc, seq, loc, qual))
            state['num'] += 1

    unpaired = [{}, {}]
    for r in bam:
        if o.C is not None and refs[r.reference_id] != o.C:
            continue
        f = r.flag
        if not o.p:
            if f & (4 | 256 | 512 | 2048):
                continue
            s = sites_of(r, snps.get(r.reference_name, {}))
            if not s:
                keep.write(r)
                continue
            write([r], [single(r, s, o.P, o.max_seqs)])
            continue
        if f & (4 | 8 | 256 | 512 | 2048):
            continue
        if not f & 2 or r.next_reference_id != r.reference_id:
            continue
        me, other = (1 if f & 128 else 0), (1 if f & 64 else 0)
        if r.query_name not in unpaired[other]:
            unpaired[me][r.query_name] = r
            continue
        both = [None, None]
        both[me] = r
        both[other] = unpaired[other].pop(r.query_name)
        r1, r2 = both
        table = snps.get(r1.reference_name, {})
        s1, s2 = sites_of(r1, table), sites_of(r2, table)
        if not s1 and not s2:
            keep.write(r1)
            keep.write(r2)
            continue
        write([r1, r2], dual(r1, r2, s1, s2, o.P, o.max_seqs))
    for f in (keep, remap, dropped):
        f.close()
    for f in fqs:
        f.close()


if __name__ == '__main__':
    main()
