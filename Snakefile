import pandas as pd
import preProcessSampleConfig as pre

configfile: 'config.json'

file_info_path = config['sampleInfo']
basename_columns = config['baseNameColumns']
pool_basename_columns = config['poolBaseNameColumns']

REFGENOME = config['refGenome']
SPIKEGENOME = config['spikeGenome']
REFGENOMEPATH = config['genome'][REFGENOME]['bowtie']
SPIKEGENOMEPATH = config['genome'][SPIKEGENOME]['bowtie']
controlDNAPath  = config['genome'][REFGENOME]['controlDNAPath']
chromSize_Path  = config['genome'][REFGENOME]['chrSize']

genomeSize = config['genome'][REFGENOME]['genomeSize']
readLen = config['readLen']

modules = config['module']
#########
# Validation 

if os.path.exists(file_info_path) == False:
	print('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name = file_info_path))

#########
# Generating sampleSheet outputs

speciesList  = [REFGENOME, SPIKEGENOME]
indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}

fragTypes    = ['allFrags', '20to120', '150to700']
normTypeList = ['', '_spikeNorm', '_rpgcNorm']

sampleInfo, sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-", fileDelimiter = config['sampleInfoDelimiter'])
poolSampleSheet = sampleSheet.copy()

sampleSheet['fastq_trim_r1'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1'])
sampleSheet['fastq_trim_r2'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['2'])
sampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam"})

for frag, norm in zip(fragTypes, normTypeList):
	# Add column per peak call list
	peak_colName = 'peak_{frag}'.format(frag = frag)
	sampleSheet[peak_colName] = expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)
	
	bed_colName = 'bed_{frag}'.format(frag = frag)
	sampleSheet[bed_colName] = expand('Bed/{sample}_{species}_trim_q5_dupsRemoved_{fragType}.bed', sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)
	

# Add columns for all combinations of fragType & normType
## Create column names:
fragTypeCombn = fragTypes * len(normTypeList)
fragTypeCombn.sort()

normTypeCombn = normTypeList * len(fragTypes)

fragNormCombn = [''.join(frag + norm) for frag, norm in zip(fragTypeCombn, normTypeCombn)]

## Add columns to sampleSheet
for fragNorm in fragNormCombn:

	# Add column per bigwig
	bw_colName = 'bigwig_{fragNorm}'.format(fragNorm = fragNorm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}.bw", sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm) 

	# Add column per zNorm bigwig
	bw_colName = 'zNorm_bigwig_{fragNorm}'.format(fragNorm = fragNorm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}_zNorm.bw", sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm) 

	# Threshold peakcalls:
	thresh_colName = 'threshold_peaks_{fragNorm}'.format(fragNorm = fragNorm)
	sampleSheet[thresh_colName] = expand('Threshold_PeakCalls/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}_thresholdPeaks.bed', sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm)

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
		expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes),
		expand('Threshold_PeakCalls/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}_thresholdPeaks.bed', sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes, normType = normTypeList),
		expand('FastQC/{sample}_R1_fastqc.html', sample = sampleSheet.baseName),
		expand('FastQC/{sample}_R1_trim_fastqc.html', sample = sampleSheet.baseName),
		expand('FQscreen/{sample}_R1_trim_screen.txt', sample = sampleSheet.baseName),
		expand('FQscreen/{sample}_R1_trim_screen.html', sample = sampleSheet.baseName),
		"multiqc_report.html",
		expand('Plots/FragDistInPeaks/{sample}_{REFGENOME}_trim_q5_allFrags_fragDistPlot.png', sample = sampleSheet.baseName, REFGENOME = REFGENOME),
		expand('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm_zNorm.bw', sample = sampleSheet.baseName, REFGENOME = REFGENOME, fragType = fragTypes)
		

rule combine_technical_reps:
	input:
		r1 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1,
		r2 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r2
	output:
		r1 = 'Fastq/{sample}_R1.fastq.gz',
		r2 = 'Fastq/{sample}_R2.fastq.gz'
	shell:
		"""
		cat {input.r1} > {output.r1} &&
		cat {input.r2} > {output.r2}
		"""

rule fastQC:
	input:
		'Fastq/{sample}_R1.fastq.gz',
	output:
		'FastQC/{sample}_R1_fastqc.html'
	params:
		module = config['module']['fastqcVer']
	shell:
		"""
		module purge && module load {params.module}

		fastqc -o ./FastQC/ -f fastq {input}
		"""

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
		module  = config['module']['bbmapVer']
	shell:
		"""
		module purge && module load {params.module}
		bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}
		"""
		#bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11

rule fastQC_trim:
	input:
		'Fastq/{sample}_R1_trim.fastq.gz',
	output:
		'FastQC/{sample}_R1_trim_fastqc.html'
	params:
		module = config['module']['fastqcVer']
	shell:
		"""
		module purge && module load {params.module}

		fastqc -o ./FastQC/ -f fastq {input}
		"""

rule fastqScreen:
	input:
		'Fastq/{sample}_R1_trim.fastq.gz'
	output:
		txt = 'FQscreen/{sample}_R1_trim_screen.txt',
		html = 'FQscreen/{sample}_R1_trim_screen.html'
	params:
		fqscreenPath = config['module']['fqscreenPath'],
		fqscreenConf = config['module']['fqscreenConf']
	threads: 4
	shell:
		"""
		{params.fqscreenPath} --threads {threads} --force --aligner bowtie2 -conf {params.fqscreenConf} {input} --outdir ./FQscreen/
		"""

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
		module = config['module']['bowtie2Ver']
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
		module = config['module']['samtoolsVer']
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
		module = modules['samtoolsVer']
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
		module = modules['picardVer'],
		picardPath = modules['picardPath']
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
		module = modules['samtoolsVer']
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
		module = modules['samtoolsVer']
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
		module = modules['samtoolsVer']
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
		module = modules['bedtoolsVer']
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
		unNorm    = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bg'),
		spikeNorm = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_spikeNorm.bg'),
		rpgcNorm  = temp('BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bg')
	params:
		genomeSize = genomeSize,
		chromSize_Path = chromSize_Path,
		sampleName = '{sample}',
		module = modules['bedtoolsVer'],
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
		module = modules['ucscVer'],
		chromSize_Path = chromSize_Path
	shell:
		"""
		module purge && module load {params.module}
		wigToBigWig {input} {params.chromSize_Path} {output}
		"""

rule zNormBigWig:
	input:
		'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bw'
	output:
		zNorm = 'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm_zNorm.bw',
		zStats = 'Logs/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.zNorm'
	params:
		module = modules['rVer']
	shell:
		"""
		module purge && module load {params.module}
		Rscript --vanilla scripts/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

rule callThresholdPeaks:
	input:
		'BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bw'
	output:
		'Threshold_PeakCalls/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}_thresholdPeaks.bed'
	params:
		module = modules['rVer']
	shell:
		"""
		module load {params.module}
		Rscript --vanilla scripts/callThresholdPeaks.R {input} {output}
		"""
	
rule callPeaks:
	input:
		'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bed'
	output:
		'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak'
	params:
		module = modules['macsVer'],
		control = controlDNAPath,
		prefix = 'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}'
	shell:
		"""
		module purge && module load {params.module}
		macs2 callpeak -f BEDPE -c {params.control} -n {params.prefix} -g 121400000 -t {input}  --nomodel --seed 123
		"""

rule qcReport:
	input:
		expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = speciesList, ftype = ['bam', 'bam.bai']),
		expand('FQscreen/{sample}_R1_trim_screen.txt', sample = sampleSheet.baseName),
		expand('FastQC/{sample}_R1_trim_fastqc.html', sample = sampleSheet.baseName)
	output:
		"multiqc_report.html"
	params: moduleVer = modules['multiqcVer']
	shell:
		"""
		module purge && module load {params.moduleVer}
		multiqc . -f -x *.out -x *.err
		"""

rule makeFragmentSizePlots_inPeaks:
	input:
		bed = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed',
		peaks = 'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags_peaks.narrowPeak'
	output:
		'Plots/FragDistInPeaks/{sample}_{REFGENOME}_trim_q5_allFrags_fragDistPlot.png'
	params:
		module = modules['rVer']
	shell:
		"""
		module purge && module load {params.module}
		Rscript --vanilla scripts/makeFragsizePlot.R {input.bed} {input.peaks} {output}
		"""

#rule makeFragmentSizePlots:
#	input:
#		bed = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed'
#	output:
#		'Plots/FragDist/{sample}_{REFGENOME}_trim_q5_allFrags_cumulativeDistPlot.png'
#	params:
#		module = modules['rVer'],
#		srcDirectory = srcDirectory
#	shell:
#		"""
#		module purge && module load {params.module}
#		Rscript --vanilla scripts/plotFragsizeDist.R {input.bed}
#		"""
