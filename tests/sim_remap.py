"""Fake a remapping of Hornet remap fastqs. usage: sim_remap.py template.bam out.bam seed fq1 [fq2]"""
import sys, gzip, random
import pysam
tmpl, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
fqs = sys.argv[4:]
rng = random.Random(seed)
t = pysam.AlignmentFile(tmpl)
refs = list(t.references)
o = pysam.AlignmentFile(out + '.unsorted.bam', 'wb', template=t)
handles = [gzip.open(f, 'rt') for f in fqs]
paired = len(fqs) == 2
while True:
    recs = []
    for h in handles:
        lines = [h.readline() for _ in range(4)]
        if not lines[0]:
            break
        recs.append(lines)
    if not recs:
        break
    if rng.random() < 0.03:
        continue  # aligner dropped this alternate
    name = recs[0][0][1:].strip()
    f = name.split(':')
    chrom, left, right = f[1], int(f[2]), int(f[3])
    for m, lines in enumerate(recs):
        seq = lines[1].strip()
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = seq
        rev = (m == 1) if paired else (rng.random() < 0.5)
        x = rng.random()
        if x < 0.05:
            a.flag = 4 | (1 if paired else 0) | ((64, 128)[m] if paired else 0)
            a.reference_id = -1; a.reference_start = -1
        else:
            a.flag = (16 if rev else 0) | ((1 | 2 | (64, 128)[m]) if paired else 0)
            a.reference_id = refs.index(chrom) if chrom in refs else 0
            p = (right if rev else left) if paired else left
            if x < 0.12:
                p += rng.randint(1, 50)
            elif x < 0.14:
                a.reference_id = (a.reference_id + 1) % len(refs)
            a.reference_start = p
            if x > 0.97:
                h = len(seq) // 2
                a.cigartuples = [(0, h), (2, 2), (0, len(seq) - h)]
            elif x > 0.95:
                a.cigartuples = [(4, 3), (0, len(seq) - 3)]
            else:
                a.cigartuples = [(0, len(seq))]
            a.mapping_quality = 60
        o.write(a)
o.close()
pysam.sort('-n', '-o', out, out + '.unsorted.bam')
