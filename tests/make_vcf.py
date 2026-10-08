#!/usr/bin/env python3
"""Make a VCF from a directory of Hornet SNP files, plus the equivalent SNP
directory that Hornet should derive from it.

usage: make_vcf.py snp_dir bam_for_header out.vcf.gz expected_snp_dir seed

SNV rows get random genotypes (het, phased het, hom-ref, hom-alt, missing,
haploid) and a few become multi-allelic 1/2 records; rows with '-' become
heterozygous deletions spanning two bases. Positions that occur twice are
skipped. A second sample with other genotypes is added.
"""
import sys
import os
import glob
import gzip
import random
import pysam

snp_dir, bam, out, exp_dir, seed = sys.argv[1:6]
rng = random.Random(int(seed))
hdr = pysam.AlignmentFile(bam).header
os.makedirs(exp_dir, exist_ok=True)
rows = {}
for f in glob.glob(os.path.join(snp_dir, '*.snps.txt.gz')):
    c = os.path.basename(f).split('.snps.txt.gz')[0]
    seen = {}
    for line in gzip.open(f, 'rt'):
        p, r, a = line.split()
        seen.setdefault(int(p), []).append((r, a))
    rows[c] = {p: v[0] for p, v in seen.items() if len(v) == 1}

lines = []
for c in [s['SN'] for s in hdr.to_dict()['SQ']]:
    if c not in rows:
        continue
    exp = []
    taken = set()
    for p in sorted(rows[c]):
        r, a = rows[c][p]
        if p in taken:
            continue
        other = rng.choice(['0/0', '0/1', '1/1'])
        if r == '-' or a == '-':
            if p + 1 in rows[c] or p + 1 in taken:
                continue
            taken.update((p, p + 1))
            lines.append((c, p, 'AC', 'A', '0/1', other))
            exp += ['%d\tA\t-\n' % p, '%d\tA\t-\n' % (p + 1)]
            continue
        taken.add(p)
        if len(r) != 1 or len(a) != 1 or r not in 'ACGT' or a not in 'ACGT' or r == a:
            continue
        x = rng.random()
        if x < 0.6:
            gt = '0/1'
        elif x < 0.7:
            gt = '1|0'
        elif x < 0.75:
            gt = '0/0'
        elif x < 0.8:
            gt = '1/1'
        elif x < 0.85:
            gt = './.'
        elif x < 0.88:
            gt = '1'
        elif x < 0.93:
            b = rng.choice([z for z in 'ACGT' if z not in (r, a)])
            lines.append((c, p, r, a + ',' + b, '1/2', other))
            exp.append('%d\t%s\t%s\n' % (p, a, b))
            continue
        else:
            gt = '0|1'
        lines.append((c, p, r, a, gt, other))
        if gt in ('0/1', '1|0', '0|1'):
            exp.append('%d\t%s\t%s\n' % (p, r, a))
    with gzip.open(os.path.join(exp_dir, c + '.snps.txt.gz'), 'wt') as fh:
        fh.write(''.join(exp))

tmp = out + '.tmp.vcf'
with open(tmp, 'w') as fh:
    fh.write('##fileformat=VCFv4.2\n')
    for s in hdr.to_dict()['SQ']:
        fh.write('##contig=<ID=%s,length=%d>\n' % (s['SN'], s['LN']))
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tother\tsample1\n')
    for c, p, r, a, gt, other in lines:
        fh.write('%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\t%s\n' % (c, p, r, a, other, gt))
pysam.tabix_compress(tmp, out, force=True)
os.unlink(tmp)
