from stat import *
import pandas as pd
import subprocess
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

bbmapVer = str('bbmap/38.22')

##############################

file_info_path = "ss.tsv"
basename_columns = ['sample','rep']
pool_basename_columns = ['sample']

# https://metagenomic-methods-for-microbial-ecologists.readthedocs.io/en/latest/day-1/
REFGENOMEPATH = '/proj/mckaylab/genomeFiles/dm3/RefGenome/dm3'
SPIKEGENOMEPATH = '/proj/seq/data/sacCer3_UCSC/Sequence/Bowtie2Index/genome'
controlDNAPath = '/proj/mckaylab/genomeFiles/dm3/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed'

REFGENOME = 'dm3'
SPIKEGENOME = 'sacCer3'

genomeSize = '121400000'
readLen = '75'

chromSize_Path = '/proj/mckaylab/genomeFiles/dm3/dm3.chrom.sizes'

speciesList  = [REFGENOME, SPIKEGENOME]
indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}
fragTypes    = ['allFrags', '20to120', '150to700']
normTypeList = ['', '_spikeNorm', '_rpgcNorm']

sampleSheet, pool_sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-")

sampleSheet['fastq_trim_r1'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1'])
sampleSheet['fastq_trim_r2'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['2'])
sampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam"})

for frag, norm in zip(fragTypes, normTypeList):
	# Add column per bedgraph
	bg_colName = 'bedgraph_{frag}{norm}'.format(frag = frag, norm = norm)
	sampleSheet[bg_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}.bg", sample = sampleSheet.baseName, species = REFGENOME, fragType = frag, normType = norm) 

	# Add column per bigwig
	bw_colName = 'bigwig_{frag}{norm}'.format(frag = frag, norm = norm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}.bw", sample = sampleSheet.baseName, species = REFGENOME, fragType = frag, normType = norm) 

	# Add column per peak call list
	peak_colName = 'peak_{frag}'.format(frag = frag)
	sampleSheet[peak_colName] = expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)
	
	bed_colName = 'bed_{frag}'.format(frag = frag)
	sampleSheet[bed_colName] = expand('Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bed', sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)

sampleSheet.to_csv('sampleSheet.tsv', sep = "\t", index = False)


####
# BEGIN PIPELINE:
####

rule all:
	input:
		expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1','2']),
		expand("Sam/{sample}_{species}_trim.sam", sample = sampleSheet.baseName, species = speciesList),
		expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = speciesList, ftype = {"bam", "bam.bai"}),
		expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes, normType = normTypeList, ftype = {"bw", "bg"}),
		expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes)

rule trim_adapter:
	input:
		r1 = "Fastq/{sample}_R1.fastq.gz",
		r2 = "Fastq/{sample}_R2.fastq.gz"
	output:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "Fastq/{sample}_R2_trim.fastq.gz"
	log:
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	params:
		module  = bbmapVer
	shell:
		"""
		module purge && module load {params.module}
		bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}
		"""
		#bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11

rule align:
	input:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "Fastq/{sample}_R2_trim.fastq.gz"
	output:
		sam = "Sam/{sample}_{species}_trim.sam",
		logInfo = "Logs/{sample}_{species}_bowtie2.txt"
	threads: 8
	params:
		refgenome = lambda wildcards: indexDict[wildcards.species],
		module = bowtie2Ver
	shell:
		"""
		module purge && module load {params.module}
		bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {params.refgenome} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
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
		module = str('picard/' + picardVer),
		picardPath = picardPath
	shell:
		"""
		module purge && module load {params.module}
		java -Xmx8g -jar {params.picardPath} SortSam INPUT= {input} OUTPUT= {output.sorted} SORT_ORDER=coordinate &&
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
		samtools view -@ 4 -bF 0x400 {input} > {output.bam} &&
		samtools index {output.bam} {output.index}
		"""

rule sortBam:
	input:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved.bam'
	output:
		bam = 'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam',
		idx = 'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam.bai'
	params:
		module = samtoolsVer
	threads: 4
	shell:
		"""
		module purge && module load {params.module}
		samtools sort -@ {threads} -o {output.bam} {input}
		samtools index {output.bam}	
		"""

#rule indexSortedBam:
#	input:
#		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam'
#	output:
#		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam.bai'
#	params:
#		module = samtoolsVer
#	shell:
#		"""
#		module purge && module load {params.module}
#		samtools index {input}	
#		"""

rule nameSortBam:
	input:
		'Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam'
	output:
		temp('Bam/{sample}_{species}_trim_q5_dupsRemoved_nameSorted.bam')
	params:
		module = samtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		samtools sort -n {input} -o {output}
		"""

rule convertBamToBed:
	input:
		bam = 'Bam/{sample}_{species}_trim_q5_dupsRemoved_nameSorted.bam',
	output:
		'Bed/{sample}_{species}_trim_q5_dupsRemoved.bed'
	params:
		module = bedtoolsVer
	shell:
		"""
		module purge && module load {params.module}
		bedtools bamtobed -bedpe -i {input.bam} | sort -k 1,1 -k 2,2n > {output}		
		"""

rule splitFragments:
	input:
		'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved.bed'
	output:
		allFrags = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed',
		smallFrags = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_20to120.bed',
		bigFrags = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_150to700.bed'
	shell:
		"""
		cut -f 1,2,6,7 {input} | awk -F '\t' '{{print $0, ($3-$2)}}' - > {output.allFrags}
		awk -v OFS='\t' '($5>20) && ($5<120) {{print $0}}' {output.allFrags} > {output.smallFrags}
		awk -v OFS='\t' '($5>150) && ($5<700) {{print $0}}' {output.allFrags} > {output.bigFrags}
		"""

rule makeFragmentBedGraphs:
	input:
		ref   = lambda wildcards : 'Bed/' + wildcards.sample + '_' + REFGENOME + '_trim_q5_dupsRemoved_' + wildcards.fragType + '.bed',
		spike = lambda wildcards : 'Bam/' + wildcards.sample + '_' + SPIKEGENOME + '_trim_q5_dupsRemoved.bam' 
	output:
		#unNorm    = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bg'),
		#spikeNorm = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_spikeNorm.bg'),
		#rpgcNorm  = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bg')
		unNorm    = 'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bg',
		spikeNorm = 'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_spikeNorm.bg',
		rpgcNorm  = 'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bg'
	params:
		genomeSize = genomeSize,
		chromSize_Path = chromSize_Path,
		sampleName = '{sample}',
		module = bedtoolsVer,
		readLen = readLen
	shell:
		"""
		module purge && module load {params.module}
		
		# Count reads in spike-in & inputs for normalization
		spikeCount=$(samtools view -c {input.spike})
		readCount=$(wc -l {input.ref} | sed -e 's/^  *//' -e 's/  */,/g' | cut -d , -f 1)
		spikeScale=$(echo "scale=5; 10000/${{spikeCount}}/" | bc)
		rpgcScale=$(echo "scale=5; {params.genomeSize}/(${{readCount}} * {params.readLen})" | bc)

		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} > {output.unNorm}
		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} -scale ${{spikeCount}} > {output.spikeNorm}
		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} -scale ${{rpgcScale}} > {output.rpgcNorm}
		"""

rule convertToBigWig:
	input:
		'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bg'
	output:
		'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bw'
	params:
		module = ucscVer,
		chromSize_Path = chromSize_Path
	shell:
		"""
		module purge && module load {params.module}
		wigToBigWig {input} {params.chromSize_Path} {output}
		"""

#rule zNormBigWig:
#	input:
#		'BigWig/{sample}_dm_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bw'
#	output:
#		zNorm = 'BigWig/{sample}_dm_trim_q5_dupsRemoved_{fragType}_rpgcNorm_zNorm.bw',
#		zStats = 'Logs/{sample}_dm_trim_q5_dupsRemoved_{fragType}.zNorm'
#	params:
#		module = rVer,
#		srcDirectory = srcDirectory
#	shell:
#		"""
#		module purge && module load {params.module}
#		Rscript --vanilla {params.srcDirectory}/zNorm.r {input} {output.zNorm} > {output.zStats}
#		"""

rule callPeaks:
	input:
		'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bed'
	output:
		'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak'
	params:
		module = macsVer,
		control = controlDNAPath,
		prefix = 'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}'
	shell:
		"""
		module purge && module load {params.module}
		macs2 callpeak -f BEDPE -c {params.control} -n {params.prefix} -g 121400000 -t {input}  --nomodel --seed 123 --extsize 125
		"""

#rule qcReport:
#	input:
#		expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = speciesList, ftype = ['bam', 'bai'])
#	output:
#		"multiqc_report.html"
#	params: moduleVer = python3Ver
#	shell:
#		"""
#		module purge && module load {params.moduleVer}
#		multiqc . -f -x *.out -x *.err 
#		"""


#rule makeFragmentSizePlots:
#	input:
#		bed = 'Bed/{sample}_dm_trim_q5_dupsRemoved_{fragType}.bed',
#		peaks = 'Peaks/{sample}_dm_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak'
#	output:
#		'Plots/{sample}_dm_trim_q5_{fragType}_cumulativeDistPlot.png'
#	params:
#		module = rVer,
#		srcDirectory = srcDirectory
#	shell:
#		"""
#		module purge && module load {params.module}
#		Rscript --vanilla {params.srcDirectory}/makeFragmentSizePlots.R {input.bed} {input.peaks} Plots/
#		"""
