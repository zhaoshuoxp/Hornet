# Benchmark data used in README.rst

* `sim_real.py genome.fa outdir n_pairs snp_rate seed` -- heterozygous individual on
  hg38 and paired 100 bp reads (data sets A: `4000000 0.002 42`, B: `2000000 0.01 7`).
* `sim_dup.py genome.fa outdir 400 3000 600000 11` -- ATAC-like, duplicate-rich data
  (read names `m<molecule>_c<copy>`).
* `count_alleles.py in.bam snp_dir` -- reference/alternative bases at SNVs.
* `dedup_metrics.py out.bam in.bam label` -- kept pairs, distinct molecules and residual
  duplicates after duplicate removal.

Reads were mapped with `bowtie2 -x hg38 -1 r1.fq.gz -2 r2.fq.gz | samtools view -b -q 10`;
remapping used the same command.
