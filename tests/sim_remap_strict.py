#!/usr/bin/env python3
"""Fake an aligner's output for the remap fastqs, with realistic records.

usage: sim_remap_strict.py to.remap.bam out.bam seed fq1 [fq2]

Most versions map back exactly like the original read (same position,
strand, CIGAR, mate). Some are perturbed: dropped, unmapped, shifted, given a
different CIGAR or strand, a wrong mate position, not flagged proper, or get
an extra secondary/supplementary record. The output is shuffled (not sorted
by name).
"""
import sys
import gzip
import random
import pysam

tmpl, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
fqs = sys.argv[4:]
paired = len(fqs) == 2
n = 2 if paired else 1
rng = random.Random(seed)
t = pysam.AlignmentFile(tmpl)
orig = list(t)
units = {i // n + 1: orig[i:i + n] for i in range(0, len(orig), n)}
recs = []
handles = [gzip.open(f, 'rt') for f in fqs]
while True:
    entries = []
    for h in handles:
        lines = [h.readline() for _ in range(4)]
        if not lines[0]:
            break
        entries.append(lines)
    if not entries:
        break
    if rng.random() < 0.03:
        continue
    name = entries[0][0][1:].strip()
    rid = int(name.split(':')[0])
    x = rng.random()
    for m, lines in enumerate(entries):
        o = units[rid][m]
        a = pysam.AlignedSegment(t.header)
        a.query_name = name
        seq = lines[1].strip()
        qual = lines[3].strip()
        if o.is_reverse:   # aligner reports the reverse strand sequence
            seq = seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
            if qual != 'None':
                qual = qual[::-1]
        a.query_sequence = seq
        if qual != 'None':
            a.query_qualities = pysam.qualitystring_to_array(qual)
        a.flag = o.flag & (1 | 2 | 16 | 32 | 64 | 128)
        a.reference_id = o.reference_id
        a.reference_start = o.reference_start
        a.cigartuples = o.cigartuples
        a.mapping_quality = 60
        if paired:
            other = units[rid][1 - m]
            a.next_reference_id = other.reference_id
            a.next_reference_start = other.reference_start
        if x < 0.02:
            a.flag = (a.flag | 4) & ~2
            a.cigartuples = None
        elif x < 0.04:
            a.reference_start += rng.randint(1, 30)
        elif x < 0.05:
            # same start, different CIGAR: soft-clip the first two bases
            # (or, if the read does not start with a match, move it)
            ops = list(o.cigartuples)
            if ops[0][0] == 0 and ops[0][1] > 3:
                a.cigartuples = [(4, 2), (0, ops[0][1] - 2)] + ops[1:]
            else:
                a.reference_start += 1
        elif x < 0.06:
            a.flag ^= 16
        elif x < 0.07 and paired:
            a.next_reference_start += 5
        elif x < 0.08 and paired:
            a.flag &= ~2
        recs.append(a)
        if 0.08 <= x < 0.085 and m == 0:
            s = pysam.AlignedSegment(t.header)
            s.query_name = name
            s.query_sequence = seq
            s.flag = a.flag | (256 if rng.random() < 0.5 else 2048)
            s.reference_id = o.reference_id
            s.reference_start = o.reference_start + 1000
            s.cigartuples = [(0, len(seq))]
            recs.append(s)
rng.shuffle(recs)
with pysam.AlignmentFile(out, 'wb', template=t) as fh:
    for a in recs:
        fh.write(a)
