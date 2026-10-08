"""Count ref/alt bases at SNV sites in a BAM.  usage: count_alleles.py in.bam snp_dir [label]

Primary, mapped reads only. Prints reads overlapping SNVs, ref/alt/other base
counts and the reference fraction ref/(ref+alt) (truth in the simulation: 0.5).
"""
import sys, io, contextlib
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'mapping'))
import find_intersecting_snps as F
import pysam

bam, snpdir = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else bam
with contextlib.redirect_stdout(io.StringIO()):
    T = F.get_snps(snpdir, None, 8)
f = pysam.AlignmentFile(bam, threads=4)
tabs = [((T.pos[c], T.code[c]) if c in T else None) for c in f.references]
ref = alt = other = nreads = nreads_snp = 0
for r in f:
    if r.flag & (4 | 256 | 2048):
        continue
    nreads += 1
    tab = tabs[r.reference_id]
    if tab is None:
        continue
    sites, _ = F.find_sites(r, tab)
    if not sites:
        continue
    nreads_snp += 1
    seq = r.query_sequence
    for _, q, _, code in sites:
        a = T.alleles[code]
        b = seq[q]
        if b == a[0]:
            ref += 1
        elif b == a[1]:
            alt += 1
        else:
            other += 1
print("%s\treads=%d\treads_with_snv=%d\tref=%d\talt=%d\tother=%d\tref_fraction=%.5f"
      % (label, nreads, nreads_snp, ref, alt, other, ref / max(1, ref + alt)))
