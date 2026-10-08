#!/usr/bin/env python3
"""Check the reverse-strand quality fix in the remap fastqs.

usage: check_fastq_quals.py to.remap.bam legacy.fq[1].gz fixed.fq[1].gz [legacy.fq2.gz fixed.fq2.gz]

The fixed fastqs must equal the legacy (original) ones, except that the
qualities of reads that are on the reverse strand in to.remap.bam are
reversed. Exits non-zero on any difference.
"""
import sys
import gzip
import pysam

bam = sys.argv[1]
pairs = [(sys.argv[i], sys.argv[i + 1]) for i in range(2, len(sys.argv), 2)]
n_mates = len(pairs)

# to.remap.bam holds the original reads in remap-ID order (1 or 2 per ID)
is_rev = {}
recs = list(pysam.AlignmentFile(bam))
for k in range(0, len(recs), n_mates):
    is_rev[k // n_mates + 1] = [r.is_reverse for r in recs[k:k + n_mates]]

n_checked = n_reversed = 0
for mate, (legacy, fixed) in enumerate(pairs):
    with gzip.open(legacy, 'rt') as a, gzip.open(fixed, 'rt') as b:
        while True:
            ra = [a.readline() for _ in range(4)]
            rb = [b.readline() for _ in range(4)]
            if not ra[0] and not rb[0]:
                break
            if ra[:3] != rb[:3]:
                sys.exit("header/sequence differs: %r vs %r" % (ra[0], rb[0]))
            rid = int(ra[0][1:].split(':', 1)[0])
            qa, qb = ra[3].rstrip('\n'), rb[3].rstrip('\n')
            if is_rev[rid][mate] and qa != 'None':   # 'None': no qualities
                expect = qa[::-1]
                n_reversed += 1
            else:
                expect = qa
            if qb != expect:
                sys.exit("quality mismatch for %s" % ra[0].strip())
            n_checked += 1
if n_reversed == 0:
    sys.exit("no reverse-strand reads checked")
print("checked %d fastq records, %d reverse-strand" % (n_checked, n_reversed))
