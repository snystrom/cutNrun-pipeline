mkdir Fastq/

touch Fastq/WildType-Rep1_L{1,2}_R{1,2}.fastq.gz

echo sample$'\t'rep$'\t'fastq_r1$'\t'fastq_r2 > sampleInfo.tsv
echo wt$'\t'rep1$'\t'Fastq/WildType-Rep1_L1_R1.fastq.gz$'\t'Fastq/WildType-Rep1_L1_R2.fastq.gz >> sampleInfo.tsv
echo wt$'\t'rep1$'\t'Fastq/WildType-Rep1_L2_R1.fastq.gz$'\t'Fastq/WildType-Rep1_L2_R2.fastq.gz >> sampleInfo.tsv

snakemake -n

rm sampleInfo.tsv
rm -r Fastq/
