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

1. Run `preProcessSampleConfig.py ss.tsv sample rep` to copy & rename raw fastq files to current directory.
2. Edit Snakefile `basename_cols = ['sample', 'rep']` (use the same column names used in the call to `preProcessSampleConfig.py`
3. Edit Snakefile `file_info_path` to point to your sample sheet (in this example `ss.tsv`)
 - Configure additional parameters if needed
4. Edit slurmConfig.json to configure default parameters if necessary (i.e. account to match lab group)
5. Run `sh slurmSubmission.sh`
6. you can monitor job progress in another terminal using `sacct -S now`


