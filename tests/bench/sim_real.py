"""Simulate a heterozygous individual on hg38 and paired reads from it.

usage: sim_real.py genome.fa outdir n_pairs snp_rate seed
"""
import sys, os, gzip
import numpy as np
import pysam

fa, out, n_pairs, rate, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed)
os.makedirs(os.path.join(out, 'snps'), exist_ok=True)
chroms = ['chr%d' % i for i in range(1, 23)] + ['chrX']
F = pysam.FastaFile(fa)
B = np.frombuffer(b'ACGT', dtype=np.uint8)
haps = {}
for c in chroms:
    g = np.frombuffer(F.fetch(c).upper().encode(), dtype=np.uint8).copy()
    L = len(g)
    n = rng.binomial(L, rate)
    pos = np.unique(rng.integers(0, L, n))
    pos = pos[np.isin(g[pos], B)]
    ref = g[pos]
    alt = B[(np.searchsorted(B, ref) + rng.integers(1, 4, len(pos))) % 4]
    kind = rng.random(len(pos))
    h1 = g.copy()
    het = kind >= 0.05
    h1[pos[het]] = alt[het]
    haps[c] = (g, h1)
    with gzip.open(os.path.join(out, 'snps', c + '.snps.txt.gz'), 'wt', compresslevel=3) as fh:
        lines = []
        for p, r, a, k in zip(pos.tolist(), ref.tobytes().decode(), alt.tobytes().decode(), kind.tolist()):
            if k < 0.03:
                lines.append('%d\t%s\t-\n' % (p + 1, r))
            elif k < 0.05:
                lines.append('%d\t-\t%s\n' % (p + 1, a))
            else:
                lines.append('%d\t%s\t%s\n' % (p + 1, r, a))
        fh.write(''.join(lines))
    print(c, L, len(pos), flush=True)

lens = np.array([len(haps[c][0]) for c in chroms], dtype=float)
which = rng.choice(len(chroms), n_pairs, p=lens / lens.sum())
RC = bytes.maketrans(b'ACGTN', b'TGCAN')
RL = 100
f1 = gzip.open(os.path.join(out, 'r1.fq.gz'), 'wt', compresslevel=1)
f2 = gzip.open(os.path.join(out, 'r2.fq.gz'), 'wt', compresslevel=1)
buf1, buf2 = [], []
qual = 'I' * RL
k = 0
for ci in range(len(chroms)):
    c = chroms[ci]
    m = int((which == ci).sum())
    g0, g1 = haps[c]
    L = len(g0)
    flen = np.clip(rng.normal(300, 50, m).astype(int), RL, 1000)
    start = rng.integers(0, L - 1001, m)
    hap = rng.integers(0, 2, m)
    flip = rng.integers(0, 2, m)
    for s, fl, h, fp in zip(start.tolist(), flen.tolist(), hap.tolist(), flip.tolist()):
        g = g1 if h else g0
        frag = g[s:s + fl].tobytes()
        if frag.count(b'N') > 10:
            continue
        a = bytearray(frag[:RL])
        b = bytearray(frag[-RL:].translate(RC)[::-1])
        for r in (a, b):
            for e in rng.integers(0, RL, 1 if rng.random() < 0.4 else 0):
                r[e] = b'ACGT'[rng.integers(0, 4)]
        if fp:
            a, b = b, a
        name = 'sim%d' % k
        k += 1
        buf1.append('@%s/1\n%s\n+\n%s\n' % (name, a.decode(), qual))
        buf2.append('@%s/2\n%s\n+\n%s\n' % (name, b.decode(), qual))
        if len(buf1) > 100000:
            f1.write(''.join(buf1)); f2.write(''.join(buf2)); buf1, buf2 = [], []
f1.write(''.join(buf1)); f2.write(''.join(buf2))
f1.close(); f2.close()
print('pairs', k)
