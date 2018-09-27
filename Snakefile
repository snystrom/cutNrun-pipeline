import re, os, sys
import pandas as pd

##############################
# Module Versions:

bowtie2Ver = str('bowtie2/2.2.8') samtoolsVer = str('samtools/1.3.1') bedtoolsVer = str('bedtools/2.25.0') picardVer = str('2.2.4') picardPath = str('/nas02/apps/picard-' + picardVer + '/picard-tools-' + picardVer + '/picard.jar') 
deeptoolsVer = str('deeptools/2.4.1')
macsVer = str('macs/2016-02-15')
ucscVer = str('ucsctools/320')
rVer = str('r/3.3.1')

python3Ver = str('python/3.5.1')

##############################
# add install script that does pybedtools and pysam


#import py_sam_2_spikenormbg


# http://metagenomic-methods-for-microbial-ecologists.readthedocs.io/en/latest/day-1/
REFGENOMEPATH = '/proj/mckaylab/genomeFiles/dm3/RefGenome/dm3'
SPIKEGENOMEPATH = '/proj/seq/data/sacCer3_UCSC/Sequence/Bowtie2Index/genome'

REFGENOME = 'dm3'
SPIKEGENOME = 'sacCer3'

chromSize_Path = '/proj/mckaylab/genomeFiles/dm3/dm3.chrom.sizes'

# Reconfigure for json or figure out what yaml structure is supposed to be
#configfile: 'sampleConfig.yaml'

# TODO:
# set filetype to be any of {fastq.gz, fastq, fq, fq.gz}
fqs = glob.glob("*.fastq.gz")
fastq_regex = re.compile("(.+\d+_\d+)-(?P<rep>(?P<id>.+)-(Rep\d+))_[A,T,G,C]+_L\d+_(R\d)_001\.fastq\.gz")

def extract_group(strings, groupID, regex):
	return([regex.match(s).group(groupID) for s in strings])

sampleList = extract_group(fqs, "rep", fastq_regex)
groupList  = extract_group(fqs, "id", fastq_regex)

baseName = fqs.split(".")[0]

sampleSheet = pd.read_table("ss.tsv", delimiter = "\t")

ss['sample_output'] = ss.[['sample', 'rep']].apply(lambda x : "{}-{}".format(*x), axis = 1)

rule all:
	expand("Sam/{sample}.sam", sample = baseName),
	expand("Bam/{sample}.bam", sample = baseName),
	expand("Bed/{sample}_{size}.bed", sample = baseName, size = {"allFrags", "20to120", "150to170"}),
	expand("Peaks/{sample}_{size}_peaks.narrowPeak", sample = baseName, size = {"allFrags", "20to120", "150to170"}),
	expand("BigWigz/{sample}_{size}.bw", sample = baseName, size = {"allFrags", "20to120", "150to170"})

rule align_ref:
	input:
		r1 = "{sample}_R1.fastq.gz",
		r2 = "{sample}_R2.fastq.gz"
	output:
		sam = "{sample}.sam",
		logInfo = "logs/{sample}_bowtie2.txt"

	threads: 8
	params:
		refgenome = REFGENOMEPATH,
		module = bowtie2Ver
	shell:
		"""
		module purge && module load {params.module}
		bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {params.refgenome} -1 {input.readOne} -2 {input.readTwo} -S {output.sam} 2> {output.logInfo}
		"""

rule bamsort:
	input:
		bam = "{sample}.bam"
	output:
		bam = "Bam/{sample}_{REFGENOME}_sorted.bam",
		idx = "Bam/{sample}_{REFGENOME}_sorted.bam.bai"
	threads: 4
	params: 
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools sort -@ {threads} -o {output.bam} {input} &&
		samtools index {output.bam}
		"""
# TODO:
# -- Continue here
#rule bamFilter:
#	input:
#		bam = "Bam/{sample}_{REFGENOME}_sorted.bam",
#		idx = "Bam/{sample}_{REFGENOME}_sorted.bam.bai"
#	output:
#		bam = "Bam/{sample}_{REFGENOME}_sorted_noYUHet.bam",
#		idx = "Bam/{sample}_{REFGENOME}_sorted_noYUHet.bam.bai"
#	params: module = samtoolsVer
#	threads: 4
#	shell: 
#		"""
#		module purge && module load {params.module}
#		"""
