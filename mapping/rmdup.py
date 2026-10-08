#!/usr/bin/env python3
"""Unbiased removal of duplicate single-end reads.

High-performance version: same algorithm and same sequence of random choices
as the original (so ``--seed`` reproduces the original script seeded the same
way), with multi-threaded BAM (de)compression.
"""
import random
from random import choice
import pysam
import os
import sys
import argparse


def rmdup(infile, outfile):
    write = outfile.write
    chr = ""
    pos = 0
    linelistplus = []
    linelistminus = []
    for line in infile:
        rname = line.reference_id
        lpos = line.reference_start
        if rname != chr or lpos != pos:
            if linelistplus:
                write(choice(linelistplus))
            if linelistminus:
                write(choice(linelistminus))
            chr = rname
            pos = lpos
            linelistplus = []
            linelistminus = []
        flag = line.flag
        if flag == 0:
            linelistplus.append(line)
        elif flag == 16:
            linelistminus.append(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help="input BAM or SAM file (must be sorted!)")
    parser.add_argument("output_bam", help="output BAM or SAM file")
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Threads for BAM (de)compression (default: 4)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Seed for the random choice of duplicates')

    options = parser.parse_args()
    if options.seed is not None:
        random.seed(options.seed)
    threads = max(1, options.threads)

    if options.input_bam.endswith(".sam") or options.input_bam.endswith("sam.gz"):
        infile = pysam.AlignmentFile(options.input_bam, "r", threads=threads)
    else:
        # assume binary BAM file
        infile = pysam.AlignmentFile(options.input_bam, "rb", threads=threads)

    if options.output_bam.endswith(".sam"):
        # output in text SAM format
        outfile = pysam.AlignmentFile(options.output_bam, "w", template=infile)
    elif options.output_bam.endswith(".bam"):
        # output in binary compressed BAM format
        outfile = pysam.AlignmentFile(options.output_bam, "wb", template=infile,
                                      threads=threads)
    else:
        raise ValueError("name of output file must end with .bam or .sam")

    rmdup(infile, outfile)
    infile.close()
    outfile.close()


if __name__ == '__main__':
    main()
