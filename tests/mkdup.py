import sys, random, pysam
src, out = sys.argv[1], sys.argv[2]; rng = random.Random(5)
f = pysam.AlignmentFile(src); recs = list(f)
byname = {}
for r in recs: byname.setdefault(r.query_name, []).append(r)
o = pysam.AlignmentFile(out + '.u.bam', 'wb', template=f)
for name, rs in byname.items():
    for k in range(rng.choice([1,1,1,2,3,6])):
        for r in rs:
            a = pysam.AlignedSegment.from_dict(r.to_dict(), f.header); a.query_name = '%s_d%d' % (name, k); o.write(a)
o.close(); pysam.sort('-o', out, out + '.u.bam')
