#!/usr/bin/env python3
"""Unbiased removal of duplicate paired-end reads.

High-performance version: same algorithm and same sequence of random choices
as the original (so ``--seed`` reproduces the original script seeded the same
way), with mate lookup by hash instead of a linear scan of the mate's
position, and multi-threaded BAM (de)compression.

The original script called ``infile.next()``, which does not exist for pysam
under Python 3; this version uses the iterator protocol instead.
"""
import random
from random import choice
import pysam
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help="input BAM or SAM file (must be sorted!)")
    parser.add_argument("output_bam", help="output BAM or SAM file")
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Threads for BAM (de)compression (default: 4)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Seed for the random choice of duplicates')

    max_window=10000
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

    readf=Read_filter(infile,outfile,max_window)
    infile.close()
    outfile.close()

class Read_filter:
    def __init__(self,infile,outfile,max_window):
        self.read_table=[[] for x in range(max_window)]
        # occupied[i] == 1 iff read_table[i] is non-empty; lets us find the
        # next non-empty slot with bytearray.find instead of a python loop.
        self.occupied=bytearray(max_window)
        self.num_reads=0
        self.cur_pos=0
        self.chr=""
        self._iter=iter(infile)
        self.cur_read=next(self._iter)
        self.infile=infile
        self.outfile=outfile
        self.max_window=max_window
        self.finished=False
        while not self.finished:
            self.fill_table()
            self.empty_slot()
        self.empty_table()

    def fill_table(self):
        if self.cur_read.reference_id != self.chr:
            self.empty_table()
            self.chr=self.cur_read.reference_id

        if self.num_reads==0:
            self.cur_pos=self.cur_read.reference_start
        chr = self.chr
        limit = self.cur_pos + self.max_window
        W = self.max_window
        table = self.read_table
        occupied = self.occupied
        nxt = self._iter.__next__
        cur_read = self.cur_read
        n = 0
        while True:
            rpos = cur_read.reference_start
            if cur_read.reference_id != chr or rpos >= limit:
                break
            s = rpos % W
            table[s].append(cur_read)
            occupied[s] = 1
            n += 1
            try:
                cur_read = nxt()
            except:
                self.finished=True
                break
        self.num_reads += n
        self.cur_read = cur_read

    def empty_table(self):
        while self.num_reads>0:
            self.empty_slot()

    def empty_slot(self):
        W = self.max_window
        table = self.read_table
        cur_pos = self.cur_pos
        slot = table[cur_pos % W]
        ends=dict()

        for read in slot:
            mate_pos=read.next_reference_start
            if mate_pos in ends:
                ends[mate_pos].append(read)
            else:
                ends[mate_pos]=[read]

        self.num_reads-=len(slot)
        table[cur_pos % W] = []
        self.occupied[cur_pos % W] = 0

        write = self.outfile.write
        # qname -> first read with that name, for each mate slot looked at.
        # The table is not modified while the slot is emptied.
        by_name = {}
        for end_key in ends:
            read_list=ends[end_key]
            if end_key < cur_pos:
                continue

            mate_slot = end_key % W
            names = by_name.get(mate_slot)
            if names is None:
                names = {}
                for mate_read in table[mate_slot]:
                    names.setdefault(mate_read.query_name, mate_read)
                by_name[mate_slot] = names
            while len(read_list)>0:
                keep_indx=choice(range(len(read_list)))
                keep_read=read_list.pop(keep_indx)
                mate_read = names.get(keep_read.query_name)
                if mate_read is not None:
                    write(keep_read)
                    write(mate_read)
        if self.num_reads>0:
            # while len(table[self.cur_pos % W])==0: self.cur_pos+=1
            start = self.cur_pos % W
            i = self.occupied.find(1, start)
            if i < 0:
                i = self.occupied.find(1, 0, start) + W
            self.cur_pos += i - start



if __name__ == '__main__':
    main()
