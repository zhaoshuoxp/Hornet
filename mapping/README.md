# Hornet Python scripts

Python implementation of the Hornet steps; same options and output as the C
program `../c/hornet` (which is faster). See `../README.rst` for the full
pipeline, input formats, options and duplicate removal.

* `find_intersecting_snps.py [-p] [-P] [--vcf FILE --sample NAME | snp_dir] in.bam`
  -- reads overlapping heterozygous SNPs, allele-swapped fastqs to remap.
  `-j N` processes the chromosomes of a sorted, indexed BAM in parallel.
* `filter_remapped_reads.py [-p] to.remap.bam remap.bam keep.bam`
  -- keep reads whose every allelic version remapped like the original.
* `make_num_from_fq.py` -- writes the `.to.remap.num.gz` file of the old WASP
  workflow from a remap fastq (not needed by `filter_remapped_reads.py`).

`find_intersecting_snps_2.py` and `filter_remapped_reads_2.py` are links kept
for backward compatibility. `--legacy` reproduces the original Hornet output
exactly.

Requires Python 3 with pysam and numpy (pandas and pigz are used if
available). Tests: `../tests/run_equivalence.sh`.
