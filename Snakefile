import re, os, sys, shutil
import pandas as pd
import preProcessSampleConfig as pre

##############################
# Module Versions:

bowtie2Ver = str('bowtie2/2.2.8') 
samtoolsVer = str('samtools/1.3.1') 
bedtoolsVer = str('bedtools/2.25.0') 
picardVer = str('2.2.4') 
picardPath = str('/nas02/apps/picard-' + picardVer + '/picard-tools-' + picardVer + '/picard.jar') 
deeptoolsVer = str('deeptools/2.4.1')
macsVer = str('macs/2016-02-15')
ucscVer = str('ucsctools/320')
rVer = str('r/3.3.1')

python3Ver = str('python/3.5.1')

bbmapVer = str('bbmap/37.50')

##############################

file_info_path = "ss.tsv"
basename_columns = ['sample', 'rep']
pool_basename_columns = ['sample']

# http://metagenomic-methods-for-microbial-ecottttttt.readthedocs.io/en/latest/day-1/
REFGENOMEPATH = '/proj/mckaylab/genomeFiles/dm3/RefGenome/dm3'
SPIKEGENOMEPATH = '/proj/seq/data/sacCer3_UCSC/Sequence/Bowtie2Index/genome'

REFGENOME = 'dm3'
SPIKEGENOME = 'sacCer3'

chromSize_Path = '/proj/mckaylab/genomeFiles/dm3/dm3.chrom.sizes'

speciesList  = [REFGENOME, SPIKEGENOME]
indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}
fragTypes    = ['allFrags', '20to120', '150to700']
normTypeList = ['', '_spikeNorm', '_rpgcNorm']

sampleSheet, pool_sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-")

# rename fastq read1 and read2 files to basename_R1.fq.gz basename_R2.fq.gz
def move_fastq(read1, read2, baseNames):
	# vectorized rename
	# requires .fastq.gz type
	# TODO: make generic & use input type as output type
	fastq_dir = 'Fastq/'

	if os.path.exists(fastq_dir):
		#os.rmdir(fastq_dir)
		shutil.rmtree(fastq_dir)

	os.mkdir(fastq_dir)

	for r1, r2, baseName in zip(read1, read2, baseNames):
		#os.symlink(r1, fastq_dir + baseName + "_R1.fastq.gz")
		#os.symlink(r2, fastq_dir + baseName + "_R2.fastq.gz")
		shutil.copyfile(r1, fastq_dir + baseName + "_R1.fastq.gz")
		shutil.copyfile(r2, fastq_dir + baseName + "_R2.fastq.gz")

## TODO:
# Make better solution to this:
# -- add column to sampleSheet with old fastq names so they can be un-renamed,
# -- change fastq_r1 and fastq_r2 to NEW filename in sampleSheet
move_fastq(sampleSheet.fastq_r1, sampleSheet.fastq_r2, sampleSheet.baseName)

rule all:
	input:
		expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1','2']),
		expand("Sam/{sample}_{species}_trim.sam", sample = sampleSheet.baseName, species = speciesList),
		expand("Bam/{sample}_{species}_trim_q5_dupsRemoved_noYUHet.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam", "bam.bai"}),
		expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_noYUHet_{fragType}{normType}.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes, normType = normTypeList, ftype = {"bw", "bg"}),
		expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_noYUHet_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes)
		#expand("Bam/{sample}.bam", sample = sampleSheet.baseName.unique()),
		#expand("Bed/{sample}_{size}.bed", sample = sampleSheet.baseName, size = {"allFrags", "20to120", "150to170"}),
		#expand("Peaks/{sample}_{size}_peaks.narrowPeak", sample = sampleSheet.baseName, size = {"allFrags", "20to120", "150to170"}),
		#expand("BigWigz/{sample}_{size}.bw", sample = sampleSheet.baseName, size = {"allFrags", "20to120", "150to170"})


rule trim_adapter:
	input:
		r1 = "Fastq/{sample}_R1.fastq.gz",
		r2 = "Fastq/{sample}_R2.fastq.gz"
	output:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "Fastq/{sample}_R2_trim.fastq.gz",
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	params:
		module  = bbmapVer
	shell:
		"""
		module purge && module load {params.module}
		bbduk.sh in1={input.readOne} in2={input.readTwo} out1={output.readOne} out2={output.readTwo} stats={output.adapterStats} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 > {output.trimStats}
		"""

rule align:
	input:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "Fastq/{sample}_R2_trim.fastq.gz"
	output:
		sam = "Sam/{sample}_{species}_trim.sam",
		logInfo = "logs/{sample}_{species}_bowtie2.txt"
	threads: 8
	params:
		refgenome = lambda wildcards: indexDict[wildcards.species],
		module = bowtie2Ver
	shell:
		"""
		module purge && module load {params.module}
		bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {params.refgenome} -1 {input.readOne} -2 {input.readTwo} -S {output.sam} 2> {output.logInfo}
		"""

rule convertToBam:
	input:
		'Sam/{sample}_{species}_trim.sam'
	output:
		'Bam/{sample}_{species}_trim.bam'
	params:
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools view -b {input} > {output}
		"""

rule qFilter:
	input:
		'Bam/{sample}_{species}_trim.bam'
	output:
		temp('Bam/{sample}_{species}_trim_q5.bam')
	params:
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools view -@ 4 -bq 5 {input} > {output}
		"""

rule markDups:
	input:
		'Bam/{sample}_{species}_trim_q5.bam'
	output:
		sorted = temp('Bam/{sample}_{species}_trim_q5_sorted.bam'),
		markedDups = temp('Bam/{sample}_{species}_trim_q5_dupsMarked.bam'),
		PCRdups = "PCRdups/{sample}_{species}_trim_PCR_duplicates"
	params:
		module = picardVer,
		picardPath = picardPath
	shell:
		"""
		module purge && module load {params.module}
		java -Xmx8g -jar {params.picardPath} SortSam INPUT= {input} OUTPUT= {output.sorted} SORT_ORDER=coordinate
		java -Xmx8g -jar {params.picardPath} MarkDuplicates INPUT= {output.sorted} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true
		"""

rule removeDups:
	input:
		'Bam/{sample}_{species}_trim_q5_dupsMarked.bam'
	output:
		bam = 'Bam/{sample}_{species}_trim_q5_dupsRemoved.bam',
		index = 'Bam/{sample}_{species}_trim_q5_dupsRemoved.bam.bai'
	params:
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools view -@ 4 -bF 0x400 {input} > {output.bam} 
		samtools index {output.bam}
		"""

rule sortBam:
	input:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved.bam'
	output:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam'
	params:
		module = samtoolsVer
	threads: 4
	shell:
		"""
		module purge && module load {params.module}
		samtools sort -@ {threads} -o {output} {input}
		"""

rule indexSortedBam:
	input:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam'
	output:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam.bai'
	params:
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools index {input}	
		"""
