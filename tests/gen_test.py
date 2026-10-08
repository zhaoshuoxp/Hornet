"""Generate a synthetic dataset that exercises the corner cases of Hornet.

usage: gen_test.py outdir seed n_pairs
"""
import sys, os, random, gzip
import pysam

out, seed, n_pairs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
collide = float(sys.argv[4]) if len(sys.argv) > 4 else 0.002
rng = random.Random(seed)
os.makedirs(out, exist_ok=True)
os.makedirs(os.path.join(out, 'snps'), exist_ok=True)

chroms = [('chr1', 400000, 0.01), ('chr2', 200000, 0.004), ('chr3', 60000, 0.12),
          ('chrM', 16000, 0.0), ('chrUn', 5000, 0.01)]
genome = {c: ''.join(rng.choice('ACGT') for _ in range(L)) for c, L, _ in chroms}
snps = {}
for c, L, dens in chroms:
    if dens == 0:
        continue
    lines = []
    d = {}
    for p in range(L):
        if rng.random() < dens:
            ref = genome[c][p]
            r = rng.random()
            if r < 0.03:
                alt = '-'
            elif r < 0.05:
                ref, alt = '-', rng.choice('ACGT')
            elif r < 0.06:
                alt = rng.choice('ACGT') + rng.choice('ACGT')   # multi-base, no '-'
            elif r < 0.07:
                alt = ref                                        # ref == alt
            elif r < 0.075:
                ref = rng.choice('ACGT')                         # wrong ref
                alt = rng.choice('ACGT')
            else:
                alt = rng.choice([b for b in 'ACGT' if b != ref])
            lines.append((p + 1, ref, alt))
            d[p] = ref + alt
            if rng.random() < 0.01:   # duplicate row, later wins
                alt2 = rng.choice('ACGT')
                lines.append((p + 1, ref, alt2))
                d[p] = ref + alt2
    if c == 'chr2':
        rng.shuffle(lines)  # unsorted file
    snps[c] = d
    with gzip.open(os.path.join(out, 'snps', c + '.snps.txt.gz'), 'wt') as fh:
        for p, r, a in lines:
            fh.write('%d\t%s %s\n' % (p, r, a))

hdr = {'HD': {'VN': '1.6', 'SO': 'unsorted'},
       'SQ': [{'SN': c, 'LN': L} for c, L, _ in chroms]}


def make_seq(c, pos, cigar):
    g = genome[c]
    s = []
    p = pos
    for op, l in cigar:
        if op in (0, 7, 8):
            for k in range(l):
                b = g[p + k] if p + k < len(g) else 'A'
                al = snps.get(c, {}).get(p + k)
                if al and rng.random() < 0.5:
                    b = rng.choice(al)
                    if b == '-':
                        b = 'N'
                if rng.random() < 0.005:
                    b = rng.choice('ACGTN')
                s.append(b)
            p += l
        elif op in (1, 4):
            s.append(''.join(rng.choice('ACGT') for _ in range(l)))
        elif op in (2, 3):
            p += l
    return ''.join(s), p


def rand_cigar(length):
    r = rng.random()
    if r < 0.6:
        return [(0, length)]
    ops = []
    remaining = length
    if rng.random() < 0.3:
        s = rng.randint(1, 10); ops.append((4, s)); remaining -= s
    if rng.random() < 0.05:
        ops.insert(0, (5, rng.randint(1, 20)))
    while remaining > 0:
        m = min(remaining, rng.randint(5, 60))
        op = rng.choice([0, 0, 0, 7, 8]) if rng.random() < 0.2 else 0
        ops.append((op, m)); remaining -= m
        if remaining <= 0:
            break
        g = rng.random()
        if g < 0.3:
            ops.append((3, rng.randint(50, 3000)))
        elif g < 0.5:
            ops.append((2, rng.randint(1, 5)))
        elif g < 0.7:
            i = min(remaining, rng.randint(1, 4)); ops.append((1, i)); remaining -= i
    return ops


def qual(n):
    return ''.join(chr(33 + rng.randint(2, 40)) for _ in range(n))


reads = []
mapped_chroms = [c for c in ('chr1', 'chr2', 'chr3', 'chrM')]
weights = [5, 3, 3, 1]
for i in range(n_pairs):
    c = rng.choices(mapped_chroms, weights)[0]
    L = dict((x[0], x[1]) for x in chroms)[c]
    l1 = rng.choice([50, 75, 100, 150])
    l2 = l1 if rng.random() < 0.9 else rng.choice([36, 50, 101])
    p1 = rng.randint(0, L - 5000)
    p2 = p1 + rng.randint(-50, 400)
    p2 = max(0, min(L - 5000, p2))
    name = 'read%d' % i
    if rng.random() < collide and i > 10:
        name = 'read%d' % rng.randint(0, i - 1)  # name collision
    c1 = rand_cigar(l1); c2 = rand_cigar(l2)
    recs = []
    for mate, (p, cig) in enumerate([(p1, c1), (p2, c2)]):
        a = pysam.AlignedSegment()
        a.query_name = name
        seq, end = make_seq(c, p, cig)
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(qual(len(seq)))
        a.reference_id = [x[0] for x in chroms].index(c)
        a.reference_start = p
        a.cigartuples = cig
        a.mapping_quality = 255
        recs.append(a)
    r1, r2 = recs
    rev1 = rng.random() < 0.5
    f1 = 1 | 64 | (16 if rev1 else 32)
    f2 = 1 | 128 | (32 if rev1 else 16)
    if rng.random() < 0.93:
        f1 |= 2; f2 |= 2
    r1.flag, r2.flag = f1, f2
    for a, b in ((r1, r2), (r2, r1)):
        a.next_reference_id = b.reference_id
        a.next_reference_start = b.reference_start
    x = rng.random()
    if x < 0.01:   # mate unmapped but proper flag kept (bad aligner)
        r2.flag |= 4; r1.flag |= 8
    elif x < 0.015:  # one mate missing (MAPQ filtered)
        r2 = None
    elif x < 0.02:   # secondary alignment of r1
        s = pysam.AlignedSegment()
        s.query_name = name; s.query_sequence = r1.query_sequence
        s.query_qualities = r1.query_qualities
        s.reference_id = r1.reference_id
        s.reference_start = min(r1.reference_start + 1000, L - 5000)
        s.cigartuples = [(0, len(r1.query_sequence))]
        s.flag = r1.flag | 256
        s.next_reference_id = r1.reference_id; s.next_reference_start = r2.reference_start
        reads.append(s)
    elif x < 0.022:  # mate on other chromosome, not proper
        r2.flag &= ~2; r1.flag &= ~2
    elif x < 0.024:  # no qualities
        r1.query_qualities = None
    reads.append(r1)
    if r2 is not None:
        reads.append(r2)

# some fully unmapped pairs
for i in range(n_pairs // 100):
    for mate in (64, 128):
        a = pysam.AlignedSegment()
        a.query_name = 'unm%d' % i
        a.query_sequence = ''.join(rng.choice('ACGT') for _ in range(50))
        a.query_qualities = pysam.qualitystring_to_array(qual(50))
        a.flag = 1 | 4 | 8 | mate
        a.reference_id = -1; a.reference_start = -1
        a.next_reference_id = -1; a.next_reference_start = -1
        reads.append(a)

# group by name (like aligner output), shuffle pair order a bit
rng.shuffle(reads)
reads.sort(key=lambda a: a.query_name)
with pysam.AlignmentFile(os.path.join(out, 'pe.bam'), 'wb', header=hdr) as fh:
    for a in reads:
        fh.write(a)
with pysam.AlignmentFile(os.path.join(out, 'se.bam'), 'wb', header=hdr) as fh:
    for a in reads:
        if a.flag & 128:
            continue
        if a.flag & 1:
            a.flag = a.flag & (4 | 16 | 256)
        fh.write(a)
for n in ('pe', 'se'):
    pysam.sort('-o', os.path.join(out, n + '.sorted.bam'), os.path.join(out, n + '.bam'))
    pysam.index(os.path.join(out, n + '.sorted.bam'))
print('records', len(reads))
