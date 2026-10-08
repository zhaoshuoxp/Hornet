#!/usr/bin/env python3
import gzip
from itertools import islice
from sys import argv

if __name__ == "__main__":
    fastq = gzip.open(argv[1], 'rt')
    base, _ = argv[1].rsplit('.remap', 1)
    num = gzip.open(base+'.to.remap.num.gz', 'wt')
    last = ''
    out = []
    # Header lines are every 4th line of the fastq.
    for line in islice(fastq, 0, None, 4):
        curr, counts = line.strip().rsplit(':', 1)
        if curr != last:
            counts = int(counts)
            out.append('{}\n'.format(2*counts))
            if len(out) >= 100000:
                num.write(''.join(out))
                out = []
        last = curr
    num.write(''.join(out))
    num.close()
