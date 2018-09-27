# CUT&RUN Pipeline
## Author: Spencer Nystrom, Chris Uyehara

### Instructions:
Edit `basename_cols` in Snakefile to create & rename fastq files into sensible format


Requires config file as tsv in form:
additional columns can be added to describe contents of file
```
sample	rep	fastq_r1	fastq_r2
sample1	rep1	SLN02_1-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R1_001.fastq.gz	SLN02_1-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R2_001.fastq.gz
sample1	rep2	SLN02_1-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R1_001.fastq.gz	SLN02_1-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R2_001.fastq.gz
sample2	rep1	SLN02_2-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R1_001.fastq.gz	SLN02_2-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R2_001.fastq.gz
sample2	rep2	SLN02_2-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R1_001.fastq.gz	SLN02_2-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R2_001.fastq.gz
sample3	rep1	SLN02_3-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R1_001.fastq.gz	SLN02_3-Geno-Time-Assay-Condition-Rep1_TAATA_L001_R2_001.fastq.gz
sample3	rep2	SLN02_3-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R1_001.fastq.gz	SLN02_3-Geno-Time-Assay-Condition-Rep2_TAATA_L001_R2_001.fastq.gz
```

### TODO:
 - add step that merges techincal replicate fastq files from config
 - find better solution to copying fastq files & merging them
 - test on real data
