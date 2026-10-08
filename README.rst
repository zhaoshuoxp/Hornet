######
Hornet
######

Unbiased allele-specific read mapping, based on the
`WASP <https://github.com/bmvdgeijn/WASP>`_ remap-and-filter strategy
(`van de Geijn et al., Nature Methods 2015
<https://www.nature.com/articles/nmeth.3582>`_).

Reads that carry the non-reference allele of a SNP map less well than reads
carrying the reference allele, which inflates reference allele counts. Hornet
finds the reads that overlap heterozygous SNPs, swaps their alleles, has you
remap the swapped reads with the same aligner, and keeps only the reads whose
every allelic version maps exactly like the original. The remaining reads
can be counted without reference mapping bias.

This fork (`zhaoshuoxp/Hornet <https://github.com/zhaoshuoxp/Hornet>`_ of
`TheFraserLab/Hornet <https://github.com/TheFraserLab/Hornet>`_) provides

* a fast C program, ``c/hornet``, and equivalent Python scripts in
  ``mapping/`` (identical output; the C version is about 6-10x faster);
* WASP-like rigour: read flag filters, a remap filter that also checks
  CIGAR, strand and mate positions, a configurable limit on allelic
  combinations;
* VCF/BCF input (heterozygous sites of one sample);
* fixes for three bugs of the original code, and ``--legacy`` to reproduce
  the original output exactly.

.. contents::
   :local:
   :depth: 1

Installation
############

C program (recommended), needs `htslib <https://github.com/samtools/htslib>`_
>= 1.10:

.. code:: shell

   git clone https://github.com/zhaoshuoxp/Hornet.git
   cd Hornet/c
   make HTSLIB_PREFIX=$CONDA_PREFIX     # htslib from conda: conda install -c bioconda htslib
   # or: make HTSLIB_DIR=/path/to/htslib-1.x   (compiled source tree, static)
   # or: make                                  (htslib found by pkg-config)
   cp hornet ~/bin/                      # anywhere on your PATH

An htslib built with libdeflate (as on bioconda) roughly halves the time
spent compressing output.

Python scripts: Python 3 with ``pysam`` and ``numpy``; ``pandas`` (faster SNP
loading) and ``pigz`` (faster fastq compression) are used if available. Run
them directly: ``python mapping/find_intersecting_snps.py ...``.

The two implementations take the same options and produce the same files;
below, ``hornet find_intersecting_snps`` can be replaced by
``python mapping/find_intersecting_snps.py`` and ``hornet
filter_remapped_reads`` by ``python mapping/filter_remapped_reads.py``.

Quick start
###########

Paired-end example (bowtie2, e.g. ATAC-seq or ChIP-seq):

.. code:: shell

   # 1. map, keep primary alignments with MAPQ >= 10
   bowtie2 -p 16 -X 2000 -x genome -1 R1.fq.gz -2 R2.fq.gz \
     | samtools view -b -q 10 -F 0x904 -o sample.bam -

   # 2. reads overlapping heterozygous SNPs of the sample -> fastqs to remap
   hornet find_intersecting_snps -p -t 16 --vcf genotypes.vcf.gz --sample S1 sample.bam
   #   writes sample.keep.bam, sample.to.remap.bam, sample.dropped.bam,
   #          sample.remap.fq1.gz, sample.remap.fq2.gz

   # 3. remap with exactly the same aligner, options and filters as in step 1
   bowtie2 -p 16 -X 2000 -x genome -1 sample.remap.fq1.gz -2 sample.remap.fq2.gz \
     | samtools view -b -q 10 -F 0x904 -o sample.remap.bam -

   # 4. keep the reads whose every allelic version remapped like the original
   #    (the remapped BAM does not need to be sorted)
   hornet filter_remapped_reads -p -t 16 sample.to.remap.bam sample.remap.bam sample.remap.keep.bam

   # 5. merge, sort, index
   samtools merge -f -@ 16 -o - sample.keep.bam sample.remap.keep.bam \
     | samtools sort -@ 16 -o sample.wasp.bam -
   samtools index sample.wasp.bam

   # 6. remove duplicates without allelic bias (see "Duplicate removal")
   samtools addreplacerg -r 'ID:S1' -r 'SM:S1' -o sample.wasp.rg.bam sample.wasp.bam
   java -jar picard.jar MarkDuplicates -I sample.wasp.rg.bam -O sample.wasp.dedup.bam \
     -M sample.dup_metrics.txt --REMOVE_DUPLICATES true --DUPLICATE_SCORING_STRATEGY RANDOM

   # 7. count alleles at heterozygous sites (e.g. in peaks), for example with
   #    wasp2-count count-variants, GATK ASEReadCounter or bcftools mpileup

Single-end data: drop ``-p`` (``find_intersecting_snps`` writes
``sample.remap.fq.gz``) and remap a single fastq.

Instead of a VCF you can give a directory of SNP files (see `SNP input`_):
``hornet find_intersecting_snps -p sample.bam snp_dir/``.

Notes by data type
==================

ATAC-seq, ChIP-seq, DNA-seq
  Works with any DNA aligner (bowtie2, BWA). Use the same options for both
  mappings (``-X 2000`` for ATAC-seq), run Hornet before Tn5 shifting or any
  other read modification, remove duplicates as described below. Single-end
  ChIP-seq is supported.

RNA-seq
  Remap with STAR using the same options (including ``--sjdbFileChrStartEnd``
  with the ``SJ.out.tab`` of the first pass if you used 2-pass mapping).
  Spliced alignments are handled; the remapped CIGAR must reproduce the
  splice junctions. STAR's built-in ``--waspOutputMode`` is an alternative
  for RNA-seq.

Hybrids and phased genotypes
  With ``-P`` each read is remapped once, with every SNP switched to the
  other allele (instead of all 2^N combinations). This is the natural test
  when the two genomes of the sample are known, e.g. F1 hybrids of inbred
  strains, and is valid whenever a read comes from one haplotype.

SNP input
#########

Either a VCF/BCF file:

``--vcf FILE [--sample NAME]``
  Heterozygous sites of the sample (``0/1``, ``1|0``, ``1/2``, ...) are used;
  homozygous, missing and haploid genotypes are ignored, so a multi-sample
  population VCF can be used directly. SNVs keep the sample's two alleles
  (for ``1/2`` the two alternative alleles). Heterozygous indels and other
  non-SNV alleles are not swapped: every reference base they cover is marked,
  and reads overlapping such a base are dropped (as in WASP). ``--sample``
  can be omitted if the file has one sample. The VCF does not need an index.

or a directory with one file per chromosome named ``<chrom>.snps.txt.gz``
(chromosome names as in the BAM), with three whitespace-separated columns::

  position(1-based)  ref_allele  alt_allele

or a single ``all.txt.gz`` in that directory with a leading chromosome
column. List only the sample's heterozygous sites. Use ``-`` for an indel
allele; any allele that is not a single base (e.g. ``A AT``) is also treated
as an indel. When a position is listed more than once the last line wins.

Commands
########

find_intersecting_snps
======================

.. code:: text

   hornet find_intersecting_snps [options] in.bam snp_dir
   hornet find_intersecting_snps [options] --vcf FILE [--sample NAME] in.bam

   -p, --paired_end       paired-end reads
   -P, --phased-snp-data  one swapped read per read (see "Hybrids and phased genotypes")
   -C, --limit-to-chrom   only load SNPs of, and process reads on, this chromosome
   --vcf FILE, --sample NAME
                          use the heterozygous sites of a sample in a VCF/BCF
   --max-seqs N           drop reads (pairs) with more than N allelic combinations
                          (default 64, i.e. at most 6 SNPs per read or pair)
   -t, --threads N        threads for BAM/fastq (de)compression (default 4)
   -j, --processes N      Python only: process chromosomes in parallel (needs a
                          sorted, indexed BAM; same output as a serial run)
   --legacy               reproduce the original Hornet output exactly

The BAM does not need to be sorted. Outputs, named after the input
(``<base>`` = name without extension):

``<base>.keep.bam``
  reads (pairs) that overlap no SNP; kept without remapping.
``<base>.to.remap.bam``
  original reads (pairs) that overlap SNPs.
``<base>.remap.fq.gz`` or ``<base>.remap.fq1.gz``/``.fq2.gz``
  allele-swapped versions to remap, in the sequenced orientation. Read names
  are ``ID:chrom:left:right:N`` (ID numbers the read/pair in
  ``to.remap.bam``, N is the number of versions).
``<base>.dropped.bam``
  reads overlapping indels, with too many SNPs, or with mates that disagree
  at a shared SNP; these should be excluded from allele counting.

Reads that are unmapped, secondary, supplementary or QC-fail are discarded
first; for paired-end data also reads with an unmapped mate, a mate on
another chromosome, or no proper-pair flag. Duplicate-flagged reads are kept
(with their flag). Statistics are printed at the end.

filter_remapped_reads
=====================

.. code:: text

   hornet filter_remapped_reads [-p] [-t N] [--legacy] to.remap.bam remap.bam keep.bam

``remap.bam`` is the output of remapping the fastqs; it can be in any order.
A read (pair) is kept only if

* none of its remapped versions has a secondary or supplementary alignment
  (the swapped read maps to more than one place);
* every remapped record is on the original chromosome, at the original
  start, on the original strand and with the original CIGAR (hard clips
  ignored);
* for pairs, every record is a proper pair with its mate at the original
  mate position;
* all versions were remapped (the number of records equals the number of
  versions, times two for pairs).

Duplicate removal
#################

Duplicates must not be removed by keeping the read with the best alignment
score or mapping quality, since those are lower for reads carrying the
non-reference allele. Use Picard with random selection (it keeps the read
with the highest hash of its name, which is independent of the allele):

.. code:: shell

   java -jar picard.jar MarkDuplicates -I in.bam -O dedup.bam -M metrics.txt \
     --REMOVE_DUPLICATES true --DUPLICATE_SCORING_STRATEGY RANDOM

Picard needs a read group (``samtools addreplacerg``). ``samtools markdup
-r`` (after ``samtools collate``/``fixmate -m``/``sort``) is faster and keeps
the read with the highest base-quality sum, which is also not allele
dependent. Duplicates can be removed after step 5 (as in WASP) or before
step 2 (less to remap). For data with UMIs use a UMI-aware tool.

The original ``rmdup.py``/``rmdup_pe.py`` have been removed: ``rmdup_pe.py``
removed almost no duplicates (it wrote every pair whose mate it found in its
window), and ``rmdup.py`` skipped the last position and all reads with flags
other than 0 or 16. On an ATAC-like simulation with PCR duplicates (900,357
pairs from 541,132 molecules, ``tests/bench/sim_dup.py``), Picard (default
and ``RANDOM``), ``samtools markdup``, ``samtools rmdup`` and WASP's
``rmdup_pe.py`` found the same duplicates (6-7 residual duplicate pairs) and
none changed the reference-allele fraction beyond sampling error, while
``rmdup_pe.py`` kept 900,355 pairs.

Benchmarks
##########

Simulated hg38 data (``tests/bench``): a heterozygous individual, reads
drawn from either haplotype with p = 0.5, so the true reference-allele
fraction is 0.5; bowtie2, MAPQ >= 10, 256-core server.

Bias removal (reads at SNVs kept; fraction of reference alleles among them):

================================ ================== ================== ==================
data set / method                reference fraction reads at SNVs kept find + filter time
================================ ================== ================== ==================
**A**: 3.5 M pairs, 1 SNV/500 bp
raw mapping                      0.50137            1,227,389
Hornet                           0.50040            1,223,141 (99.7%)  7 s
WASP 0.3 (original)              0.50040            1,223,149 (99.7%)  524 s
**B**: 1.8 M pairs, 1 SNV/100 bp
raw mapping                      0.50388            2,164,071
Hornet                           0.50014            2,135,130 (98.7%)  11 s
WASP 0.3 (original)              0.50022            2,134,213 (98.6%)  852 s
original Hornet                  0.50011            2,104,709 (97.3%)
================================ ================== ================== ==================

Both Hornet and WASP bring the reference fraction to 0.5 within sampling
error (about +-0.0003).

Speed on data set A (``-p``):

====================== ======================= ============================= ============================
step                   original Hornet         Python (``mapping/``)         C (``hornet``)
====================== ======================= ============================= ============================
find_intersecting_snps 367 s, 716 MB           62 s (-t 4); 14 s (-j 24)     10.6 s (-t 4); 5.9 s (-t 16)
filter_remapped_reads  28.5 s (+ name sort)    15.7 s (-t 4)                 2.4 s (-t 4)
====================== ======================= ============================= ============================

SNPs take ~12 bytes each in memory (the original used ~150 bytes per SNP,
about 10 GB for 70 M SNPs).

Comparison with WASP and WASP2
##############################

* *WASP 0.3* removes the bias equally well (table above) and applies the same
  kind of read and remap filters, but needs a sorted BAM and is ~80x slower.
  Its HDF5 mode can restrict the swapped versions to the haplotypes observed
  in a panel.
* *WASP2* (1.4.1) adds indel-aware remapping, allele counting, imbalance
  tests and single-cell support, but its ``make-reads`` supports only
  paired-end reads with phased genotypes (``Single-End not Implemented``,
  ``Unphased not Implemented``). In our test it also kept almost no
  SNP-overlapping reads: ``make-reads`` writes only the haplotype that
  differs from the read but names it as one of two (``total_seqs = 2`` is
  fixed in ``rust/src/bam_remapper.rs``), so ``filter-remapped`` waits for a
  second version that never comes and discards the read. The same code is on
  the current master branch.

A practical combination is Hornet for the mapping-bias filter and WASP2 (or
another tool) for counting and testing allelic imbalance.

Differences from the original Hornet
####################################

Bugs fixed (``--legacy`` restores the original behaviour):

1. Phased mode (``-P``) never counted covered SNPs and crashed with
   ``ZeroDivisionError`` when printing the statistics.
2. The fastq qualities of reverse-strand reads were not reversed along with
   the reverse-complemented sequence.
3. In single-end mode ``filter_remapped_reads`` skipped the read following
   every discarded read.

New default behaviour (``--legacy`` turns all of it off):

* read flag filters in ``find_intersecting_snps`` (see above; the original
  only checked the proper-pair flag, and nothing for single-end data);
* ``--max-seqs`` (default 64; the original used a fixed 32 and, for
  single-end reads, checked it only at the next aligned base);
* alleles that are not single bases are treated as indels (the original
  only recognised ``-``);
* ``-P`` also works for single-end data;
* the WASP-like ``filter_remapped_reads`` (the original compared only the
  start of each read, rejected any remapped read with an insertion or
  deletion, and needed a remapped BAM sorted by name);
* ``--vcf`` input; the duplicate-removal scripts are gone.

With ``--legacy`` (and the original command lines) both implementations
reproduce the original output exactly: the same BAM records in the same
order, the same fastq content and statistics (the remapped BAM must then be
sorted by name, as before).

Testing
#######

``tests/run_equivalence.sh`` (needs samtools and git) checks, on synthetic
data with many corner cases:

* ``--legacy`` against the original scripts (taken from upstream commit
  ``e642ac0``), Python and C;
* the default mode of both implementations against the simple reference
  implementations ``tests/ref_find.py`` and ``tests/ref_filter_strict.py``,
  including VCF input, ``-P``, ``-C``, ``--max-seqs`` and the parallel mode.

``tests/bench`` contains the simulations used for the benchmarks.
