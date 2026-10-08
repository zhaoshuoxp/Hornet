#!/usr/bin/env python3
"""Slow, simple reference implementation of filter_remapped_reads (default,
WASP-like mode).  usage: ref_filter_strict.py [-p] to.remap.bam remap.bam out.sam
"""
import sys
import pysam

args = sys.argv[1:]
paired = args[0] == '-p'
if paired:
    args = args[1:]
to_remap, remap, out = args
n = 2 if paired else 1

t = pysam.AlignmentFile(to_remap)
orig = list(t)
units = {i // n + 1: orig[i:i + n] for i in range(0, len(orig), n)}


def cigar(r):
    return [c for c in (r.cigartuples or []) if c[0] != 5]


expected, correct, bad = {}, {}, set()
for r in pysam.AlignmentFile(remap):
    fields = r.query_name.split(':')
    rid = int(fields[0])
    expected[rid] = int(fields[-1].split('/')[0]) * n
    if r.is_secondary or r.is_supplementary:
        bad.add(rid)
        continue
    if r.is_unmapped:
        continue
    mates = units[rid]
    o = mates[1] if (paired and r.is_read2) else mates[0]
    ok = (r.reference_name == o.reference_name and r.reference_start == o.reference_start
          and r.is_reverse == o.is_reverse and cigar(r) == cigar(o))
    if paired:
        other = mates[0] if r.is_read2 else mates[1]
        ok = ok and r.is_proper_pair and not r.mate_is_unmapped and \
            r.next_reference_name == o.reference_name and \
            r.next_reference_start == other.reference_start
    if ok:
        correct[rid] = correct.get(rid, 0) + 1
    else:
        bad.add(rid)

o = pysam.AlignmentFile(out, 'w', template=t)
for rid, mates in units.items():
    if rid not in bad and rid in expected and correct.get(rid, 0) == expected[rid]:
        for r in mates:
            o.write(r)
o.close()
