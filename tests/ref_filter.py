#!/usr/bin/env python3
"""Straightforward reference implementation of filter_remapped_reads.

usage: ref_filter.py [-p] to.remap.bam remap.bam out.sam

A read (pair) from to.remap.bam is kept if every one of its remapped
versions mapped back to the original position: the number of remapped
records that map correctly must equal the number expected from the read name
(times two for pairs). "Mapping correctly" uses Hornet's criteria: same
chromosome, same start (left position for forward reads, right position for
reverse reads of pairs) and a CIGAR made only of M, N and S operations.
"""
import sys
from collections import Counter
import pysam

args = sys.argv[1:]
paired = args[0] == '-p'
if paired:
    args = args[1:]
to_remap, remap, out = args

correct = Counter()
expected = {}
f = pysam.AlignmentFile(remap)
for r in f:
    fields = r.query_name.split(':')
    rid, chrom, num = int(fields[0]), fields[1], int(fields[-1])
    expected[rid] = num * (2 if paired else 1)
    pos = int(fields[3] if (paired and r.is_reverse) else fields[2])
    if (r.reference_id != -1 and r.reference_start == pos and
            r.reference_name == chrom and
            all(op in (0, 3, 4) for op, _ in (r.cigartuples or ()))):
        correct[rid] += 1

t = pysam.AlignmentFile(to_remap)
o = pysam.AlignmentFile(out, 'w', template=t)
recs = list(t)
n = 2 if paired else 1
for k in range(0, len(recs), n):
    rid = k // n + 1
    if rid in expected and correct[rid] == expected[rid]:
        for r in recs[k:k + n]:
            o.write(r)
o.close()
