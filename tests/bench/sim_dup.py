"""ATAC/ChIP-like duplicate-rich simulation on hg38.

usage: sim_dup.py genome.fa outdir n_regions region_len n_molecules seed

Fragments come from a few hundred 'peaks' at high depth, so many independent
molecules share start positions; each molecule is PCR-amplified (geometric
number of copies) and every copy gets its own sequencing errors and random
base qualities. Read names: m<molecule>_c<copy>. Het SNVs (1/150 bp) on one
haplotype; molecules are drawn from either haplotype with p=0.5.
"""
import sys, os, gzip
import numpy as np
import pysam

fa, out, n_reg, reg_len, n_mol, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
rng = np.random.default_rng(seed)
os.makedirs(os.path.join(out, 'snps_snv'), exist_ok=True)
F = pysam.FastaFile(fa)
chroms = ['chr%d' % i for i in range(1, 23)]
B = np.frombuffer(b'ACGT', dtype=np.uint8)
regions = []
snps = {}
while len(regions) < n_reg:
    c = chroms[rng.integers(0, len(chroms))]
    L = F.get_reference_length(c)
    s = int(rng.integers(1000000, L - 1000000))
    seq = F.fetch(c, s, s + reg_len).upper()
    if seq.count('N') or seq.count('A') + seq.count('T') > 0.7 * reg_len:
        continue
    g0 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    g1 = g0.copy()
    pos = np.unique(rng.integers(0, reg_len, reg_len // 150))
    alt = B[(np.searchsorted(B, g0[pos]) + rng.integers(1, 4, len(pos))) % 4]
    g1[pos] = alt
    regions.append((c, s, g0, g1))
    for p, r, a in zip(pos.tolist(), g0[pos].tobytes().decode(), alt.tobytes().decode()):
        snps.setdefault(c, []).append((s + p + 1, r, a))
for c, lst in snps.items():
    lst.sort()
    with gzip.open(os.path.join(out, 'snps_snv', c + '.snps.txt.gz'), 'wt') as fh:
        fh.write(''.join('%d\t%s\t%s\n' % x for x in lst))

RL = 75
RC = bytes.maketrans(b'ACGT', b'TGCA')
f1 = gzip.open(os.path.join(out, 'r1.fq.gz'), 'wt', compresslevel=1)
f2 = gzip.open(os.path.join(out, 'r2.fq.gz'), 'wt', compresslevel=1)
reg_idx = rng.integers(0, n_reg, n_mol)
n_reads = 0
buf1, buf2 = [], []
for m in range(n_mol):
    c, s, g0, g1 = regions[reg_idx[m]]
    fl = int(np.clip(rng.normal(200, 40), RL + 10, reg_len - 1))
    st = int(rng.integers(0, reg_len - fl))
    g = g1 if rng.random() < 0.5 else g0
    frag = g[st:st + fl].tobytes()
    a0, b0 = frag[:RL], frag[-RL:].translate(RC)[::-1]
    flip = rng.random() < 0.5
    copies = int(rng.geometric(0.6))     # mean 1.67 copies -> ~40% duplicates
    for k in range(copies):
        reads = []
        for r in (a0, b0):
            r = bytearray(r)
            q = rng.integers(15, 41, RL)
            for e in np.nonzero(rng.random(RL) < 0.004)[0]:
                r[e] = b'ACGT'[rng.integers(0, 4)]
            reads.append((r.decode(), ''.join(chr(33 + x) for x in q)))
        if flip:
            reads.reverse()
        name = 'm%d_c%d' % (m, k)
        buf1.append('@%s/1\n%s\n+\n%s\n' % (name, reads[0][0], reads[0][1]))
        buf2.append('@%s/2\n%s\n+\n%s\n' % (name, reads[1][0], reads[1][1]))
        n_reads += 1
    if len(buf1) > 50000:
        f1.write(''.join(buf1)); f2.write(''.join(buf2)); buf1, buf2 = [], []
f1.write(''.join(buf1)); f2.write(''.join(buf2)); f1.close(); f2.close()
print('molecules', n_mol, 'read pairs', n_reads, 'snvs', sum(len(v) for v in snps.values()))
