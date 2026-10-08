"""Duplicate-removal metrics on sim_dup.py data. usage: dedup_metrics.py out.bam input.bam label

pairs      : read pairs kept (counted by primary read1)
molecules  : distinct original molecules among them
resid_dups : pairs - molecules (PCR duplicates that survived)
lost       : molecules present in the input but absent from the output
"""
import sys
import pysam


def mols(bam):
    pairs = 0
    seen = set()
    for r in pysam.AlignmentFile(bam):
        if r.flag & (4 | 256 | 2048) or (r.flag & 1 and not r.flag & 64):
            continue
        pairs += 1
        seen.add(r.query_name.split('_c')[0])
    return pairs, seen


out, inp, label = sys.argv[1:4]
p_in, m_in = mols(inp)
p, m = mols(out)
print("%s\tpairs=%d\tmolecules=%d\tresid_dups=%d\tlost=%d\t(input: pairs=%d molecules=%d)"
      % (label, p, len(m), p - len(m), len(m_in - m), p_in, len(m_in)))
