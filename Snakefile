import re
import os

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

rule all:
	expand("Sam/{sample}.sam", sample = baseName),
	expand("Bam/{sample}.bam", sample = baseName),
	expand("Bed/{sample}_{size}.bed", sample = baseName, size = {"allFrags", "20to120", "150to170"}),
	expand("Bam/{sample}.bam", sample = baseName),
