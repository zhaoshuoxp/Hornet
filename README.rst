######
Hornet
######

A pipeline for unbiased allele-specific read mapping based on the
`WASP pipeline <https://github.com/bmvdgeijn/WASP>`_ by the
`Pritchard lab <http://pritchardlab.stanford.edu/home.html>`_.

Introduction
############

Hornet is a simple set of tools for unbiased allele-specific read mapping
that is originally described in
`van de Geijn B, McVicker G, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci" <https://www.nature.com/articles/nmeth.3582>`_

High-performance version
########################

This fork contains two faster implementations of every step. Both produce
**exactly the same results as the original code**: the same BAM records in
the same order, the same fastq content (byte for byte after decompression) and
the same run statistics. This is checked by ``tests/run_equivalence.sh``,
which runs the original scripts (taken from upstream commit ``e642ac0``) next
to the new ones on synthetic data covering the corner cases (indels,
duplicated and multi-allelic SNP rows, ``N``/``D``/``I``/``S``/``H``/``=``/``X``
CIGARs, >32 SNP combinations, unmapped mates, secondary alignments, improper
pairs, phased mode, ``-C``, consolidated ``all.txt.gz``, ...).

``c/hornet`` (recommended)
  A single C program built on htslib with one sub-command per script::

    hornet find_intersecting_snps [-p] [-P] [-C chrom] [-t threads] in.bam SNP_dir
    hornet filter_remapped_reads  [-p] [-t threads] to.remap.bam remap.sort.bam out.bam
    hornet rmdup    [-t threads] [--seed N] in.bam out.bam
    hornet rmdup_pe [-t threads] [--seed N] in.bam out.bam

  Build it with::

    cd c
    make HTSLIB_PREFIX=$CONDA_PREFIX     # htslib from conda/bioconda, or
    make HTSLIB_DIR=/path/to/htslib-1.x  # a compiled htslib source tree, or
    make                                 # htslib found by pkg-config

  An htslib built with libdeflate (e.g. ``conda install -c bioconda htslib``)
  roughly halves the time spent compressing output.

``mapping/*.py``
  Drop-in replacements for the original Python scripts (same command lines,
  plus ``-t/--threads``; ``find_intersecting_snps.py`` also has
  ``-j/--processes``, which processes the chromosomes of a sorted and indexed
  BAM in parallel). Requires pysam and numpy; pandas and pigz are used if
  available.

``rmdup.py``/``rmdup_pe.py`` pick duplicates at random. Both new versions make
the same sequence of calls to Python's Mersenne Twister as the original, so
``--seed N`` gives the same reads as the original script run after
``random.seed(N)`` (the C version re-implements CPython's generator).

Benchmark
---------

Simulated hg38 data: 3.55 M read pairs (bowtie2, MAPQ >= 10), 6.0 M SNPs
(5% indels). Wall-clock time, 256-core Linux server; the C program is linked
to htslib 1.21 with libdeflate.

======================== ================== =================================== ===============================
step                     original           Python (new)                        C
======================== ================== =================================== ===============================
find_intersecting_snps   366.8 s, 716 MB    137.0 s (-t 1), 55.2 s (-t 4),      41.7 s (-t 1), 10.4 s (-t 4),
-p                                          13.7 s (-j 24, sorted input),       4.6 s (-t 16), 4.2 s (-t 32),
                                            ~270 MB                             ~190 MB
filter_remapped_reads -p 28.5 s, 154 MB     23.1 s (-t 1), 10.1 s (-t 4)        8.1 s (-t 1), 2.3 s (-t 4), 36 MB
rmdup_pe                 359.3 s [1]_       23.1 s                              5.5 s
rmdup                    28.1 s             6.7 s                               3.0 s
======================== ================== =================================== ===============================

.. [1] The original ``rmdup_pe.py`` calls ``AlignmentFile.next()``, which does
   not exist in Python 3; it was run through a small compatibility shim
   (``tests/shim_run.py``).

SNP memory use is ~12 bytes per SNP in both new versions, against ~150 bytes
per SNP in the original nested dictionaries (about 10 GB for 70 M SNPs).

Behaviour kept for compatibility
--------------------------------

Because the output has to be identical, a few quirks of the original code are
reproduced on purpose:

* Phased mode (``-P``) never counts covered SNPs, so the original crashes with
  ``ZeroDivisionError`` while printing the statistics (after all output files
  are written). The new versions print ``nan`` for that percentage and finish.
* The fastq qualities of reverse-strand reads are written in reference
  orientation although the sequence is reverse-complemented.
* In single-end mode ``filter_remapped_reads`` also skips the read that
  follows a read which failed to remap.
* ``rmdup.py`` keeps only reads with flag 0 or 16 and never writes the last
  position of the file; ``rmdup_pe.py`` writes every read whose mate is found
  in its window (so it removes few duplicates) and drops pairs whose mates
  start at the same position.

Overview
########

Reads are mapped normally using a mapper chosen by the user (must output
BAM or SAM format).  Then mapped reads that overlap single nucleotide
polymorphisms (SNPs) are identified. For each read that overlaps a SNP, its
genotype is swapped with that of the other allele and the read is re-mapped.
Re-mapped reads that fail to map to exactly the same location in the genome are
discarded.

Step 1:
-------

Map the fastq files using your favorite mapper/options and filter for quality using a cutoff of your choice

Example:
~~~~~~~~

.. code:: shell
 
	STAR \
            --genomeDir genome_dir/STAR \
            --outFileNamePrefix analysis_dir/STAR1/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis_dir/STAR1/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn input.R1.fastq.gz input.R2.fastq.gz \
	    --outFilterMultimapNmax 1

Step 2:
-------

Use `find_intersecting_snps` to identify reads that may have mapping biases. 
It is recommended to remove duplicate reads prior to this step to save time 
provessing reads that will ultimately be discarded.

Usage::

	find_intersecting_snps.py [-p][-C][-P][-s] <input.bam> <SNP_directory>
	   -p indicates that reads are paired-end (default is single)
	   -C limits the list of SNPs to those found on a specific chromosome
	   -P indicates that the list of SNPs are phased; if so, only 2 reads will generated 
	      per input read that overlaps a SNP (as in a hybrid)
	   -s indicates that input data are coordinate sorted
	   <input.bam> is the bamfile from the initial mapping process
	   <SNP_directory> is the directory containing the SNPs segregating within the
	      sample in question (which need to be checked for mappability issues).  This directory
	      should contain coordinate sorted files of SNPs separated by chromosome and named:
	         chr<#>.snps.txt.gz
	      These files should contain 3 columns: position RefAllele AltAllele


Output::

	input.keep.bam - bamfile with reads that did not intersect SNPs or indels and therefore can
	   be kept without remapping
	input.to.remap.bam - bamfile with original reads that overlapped SNPs that need to be remapped
	input.remap.fq.gz - fastq file containing the reads with the new variants to remap. If the
	    paired-end option is used two files ending with .fq1.gz and .fq2.gz will be output.
	    
	Run statistics: run statistics will be reported once this step is complete, and include:
	  - Total input reads
	  - Reads with no SNPs
	  - Reads overlapping SNPs
	  - Total SNPs covered (refers to the number of instances where a read was found to overlap a SNP)
	  - Reference SNP matches (refers to the % of 'Total SNPs covered' that match EITHER allele in the SNP file; 
	  	should be as close to 100% as possible if SNPs have been called correctly)
	  - Non-reference SNP matches (should be a low %, caused by mutations or sequencing errors)
	  - Reads dropped [INDEL]
	  - Reads dropped [too many SNPs]
	  - Reads dropped [multivalent SNPs] (indicates a single read had SNPs matching both the ref and alt alleles)
	
	Note: Reads that overlap indels are currently excluded and will not be present in any of the 'remap' files
	or the input.keep.bam file. For this reason the total number of reads will not add up to the
	number of reads provided in the input.sort.bam file.


Example (after removing duplicates):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   python find_intersecting_snps.py \
        -p -P -s\
        rmdup.bam SNP_directory

Step 3
------
Map the input.remap.fq.gz file(s) using the same mapping arguments used in Step 1. 
The arguments should be exactly the same as those in Step 1 EXCEPT for arguments that
directly modify the reads that are used by the aligner. For example the read trimming
arguments to bowtie (-3 and -5 arguments) should be used in Step 1 ONLY because
they modify the reads that are output by bowtie. If you used 2-pass mapping in STAR 
where splice junctions were included, you can use the same SJ.out.tab file for remapping 
in this step.

Example:
~~~~~~~~

.. code:: shell

	STAR \
            --genomeDir genome_dir/STAR \
            --outFileNamePrefix analysis_dir/remap/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis_dir/remap/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn rmdup.remap.fq1.gz rmdup.remap.fq2.gz \
	    --outFilterMultimapNmax 1
	    
	    mv analysis_dir/remap/remapAligned.out.bam rmdup.remap.bam

Step 4
------
Use filter_remapped_reads.py to retrieve reads that remapped correctly.
The remapped bam file MUST be sorted by read name for this step. The read names
contain the original mapping information, which is needed in this step to 
determine whether the reads mapped to the same location. 

Usage::

	filter_remapped_reads.py [-p] <to.remap.bam> <remapped_reads.bam> <output.bam>
	   -p option indicates that the reads are paired-end
	   <to.remap.bam> output from find_intersecting_snps.py which contains
	      the original aligned reads that were remapped
	   <remapped_reads.bam> output from the second mapping step (Step 3)
	   <output.bam> file where reads that are kept after remapping are stored
	   

Example:
~~~~~~~~

.. code:: shell

 samtools sort -n rmdup.remap.bam -o rmdup.remap.sort.bam

 python filter_remapped_reads.py \
        -p \
        rmdup.to.remap.bam rmdup.remap.sort.bam \
        rmdup.remap.kept.bam

At the end of the pipeline, rmdup.remap.keep.bam and rmdup.remap.keep.bam
can be sorted and merged for a complete set of mappability filtered aligned reads. 
The merged file should then be indexed:

.. code:: shell

  samtools sort rmdup.remap.kept.bam -o rmdup.remap.kept.sort.bam
  samtools sort rmdup.keep.bam -o rmdup.keep.sort.bam
  samtools merge rmdup.remap.kept.merged.bam rmdup.keep.sort.bam rmdup.remap.kept.sort.bam
  samtools index rmdup.remap.kept.merged.bam



Dependencies
############

The Python scripts need Python 3,
`pysam <https://github.com/pysam-developers/pysam>`_ and
`numpy <http://www.numpy.org>`_. `pandas <https://pandas.pydata.org>`_
(faster SNP loading) and ``pigz`` (faster fastq compression) are used when
available.

The C program needs `htslib <https://github.com/samtools/htslib>`_ (>= 1.10)
and zlib.

Installation
############

.. code:: shell

   git clone https://github.com/zhaoshuoxp/Hornet.git
   cd Hornet/c && make HTSLIB_PREFIX=$CONDA_PREFIX   # C version
   pip install pysam numpy pandas                    # Python version

The Python scripts in ``mapping/`` can then be run directly; copy
``c/hornet`` to a directory on your ``PATH``. To check both against the
original code: ``tests/run_equivalence.sh`` (needs samtools).
