# CUT&RUN Pipeline
## Author: Spencer Nystrom, Chris Uyehara, Jeanne-Marie McPherson

## Quick Start:

Clone pipeline (Current stable version: v1.7.2)
```
git clone https://github.com/snystrom/cutNrun-pipeline.git --branch v1.7.2 --depth 1 && cd cutNrun-pipeline/ && rm -rf .git
```

Create `sampleInfo.tsv` ([see below](#sampleInfo)) with descriptive columns of data.
**Note:** `fastq_r1` and `fastq_r2` columns are required.
```
sample	rep	fastq_r1	fastq_r2
mySample	Rep1	path/to/mySample_R1.fastq.gz	path/to/mySample_R2.fastq.gz
```

edit `config.json` and set `baseNameColumns` to each column of `sampleInfo.tsv` which describes individual replicates (pipeline will automatically pool technical replicates).

Set desired reference and spike-in genome in `config.json` ([see below](#config)). Set multiple spike-in genomes by passing as an array (e.g. `["sacCer3", "droYak2"]`).

```
{
	"sampleInfo" : "sampleInfo.tsv",
	"sampleInfoDelimiter" : "\t",
	"baseNameColumns" : ["sample", "rep"],
	"refGenome" : "dm6",
	"spikeGenome" : "sacCer3",
	"readLen" : 75
}
```

Edit `slurmConfig.json` to configure default parameters if necessary.

Deploy submission with `sh slurmSubmission.sh`

## Description of Pipeline Steps

This pipeline in implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow manager for handling job dependencies and submissions to HPC clusters. Snakemake will automate the job submission process for each sample individually and will therefore often run some of these steps seemingly out of order. Below is described the general flow of information through the pipeline. In some instances many jobs are handled by snakemake to accomplish a single defined "step" below.

The rules found in the [Snakefile](Snakefile) for this pipeline are written generally in the order in which they happen from top to bottom. See there for more detail.

1. combine technical replicates
 - Autodetects technical replicates and pools reads for alignment
2. adapter trimming using [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
 - soft-clips illumina adapters from reads
3. Build combined reference + spike-in genome for alignment
 - chromosomes for each species are prefixed as: "{species}_{chromosome}"
3. Alignment with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
4. Conversion to Bam format & removing reads with QUAL score < 5 ([samtools](http://www.htslib.org/))
5. Mark and remove PCR Duplicates using [picard](https://broadinstitute.github.io/picard/)
6. Split out reads from each species into individual bam files `scripts/get_bam_reads_prefix.sh`
 - chromosomes for each species are renamed to their original forms
 - Bam headers for each split file will still contain references to other spike-in chromosomes in form: "{species}_{chromosome}"
6. Sort Bam file & convert to bed format with [bedtools][bedtools]
7. Split small (20-120bp) and large (150-700bp) fragments into 2 bed files with [awk](https://www.geeksforgeeks.org/awk-command-unixlinux-examples/)
8. Make coverage files (bigwig format) from allFragments, small Fragments, large Fragments with [bedtools][bedtools]. Generates unnormalized, RPGC-normalized, and spike-in normalized files.
 - Spike-in normalization is performed separately for each spike-in genome. Files are listed as "{species}-spikeNorm.bw"
9. Z-score normalize bigwig files ([zNorm.R](scripts/zNorm.r))
10. Call traditional peaks with [MACS2](https://github.com/taoliu/MACS) and call peaks using threshold peak-calling with [callThresholdPeaks.R](scripts/callThresholdPeaks.R)
11. Make fragment-size distribution plots for each sample ([makeFragsizePlot.R](scripts/makeFragsizePlot.R))
12. Compute QC metrics for all samples using [FastQC][fastqc], and [multiqc][multiqc]

## <a name="sampleInfo"></a> Sample Info Requirements

The `sampleInfo` file's purpose is to describe each experiment in as much
detail as desired. An arbitrary number of sample descriptive columns can be
created. Helpful columns can include: replicate number, genotype, antibody,
library prep batch, developmental stage, etc. Consider adding columns for any
variable (technical or biological) that may be useful in downstream analysis,
as these will be appended to the resulting `sampleSheet.tsv` output of the
pipeline which acts as a useful input in downstream analyses.

A minimal sampleInfo file looks as follows (here in tab-separated format).
```
sample	rep	fastq_r1	fastq_r2
mySample	Rep1	path/to/mySample_R1.fastq.gz	path/to/mySample_R2.fastq.gz
```

The `fastq_r1` and `fastq_r2` columns point to the location of one pair of
fastq files for a sample described on that row. Technical replicates (i.e.
deeper sequencing of the same library) should be added as their own row with the same column values except for the fastq file paths.

**Note:** The current version of this pipeline (`v1.7.2` as of this writing) only works with the `.fastq.gz` file format. 

The sampleInfo file can be delimited with any delimiter (default is
tab-delimited). To tell the pipeline which delimiter is being used, set the
`sampleInfoDelimiter` flag in `config.json`. For example, to use csv format set
the following: `"sampleInfoDelimiter" : ","",`.


## <a name="config"></a> config.json Requirements


[config.json](config.json) describes pipeline-specific variables in [json format](https://www.tutorialspoint.com/json/json_overview.htm). 

The config file is broken up into 3 sections:
### run config
The variables in the top few lines of the file describe run-specific information.

**sampleInfo** the path to the `sampleInfo.tsv` file if located somewhere other than the pipeline directory.

**sampleInfoDelimiter** the delimiter used in the `sampleInfo` file. Default: `"\t"`

**baseNameColumns** a list describing the columns from `sampleInfo` that distinguish replicates from eachother. The listed columns will be combined together and used to name all downstream files from each sample. The combination of these columns will also be used to detect and pool reads from technical replicates. **WARNING:** this means that any samples which evaluate to the same combination of columns will be combined into a single sample. Be sure to use enough descriptive columns to uniquely identify each separate sample.

**refGenome** the genome to align samples to as the reference. Must be a value
from the [genome](#configGenome) section. **NOTE:** there is currently no error
checking to ensure this value is found in "genome" (see [#28][i28]).

**spikeGenome** genome to be used as a spike-in normalization. Must be found in
[genome](#configGenome). To use multiple genomes, pass as an array (e.g. ["sacCer3", "droYak2"]) **NOTE:** no error checking currently (see
[#28][i28]).

**readLen** read length. If using variable read-lengths set to the smallest
value. Used in normalizing coverage files. **Will update in future to allow
sample-specific read-length adjustment.** (see [#25][i25]).

### <a name="configGenome"></a> genome config 
This section is a nested structure which points to the absolute path to the
bowtie index, whole genome fasta file, UCSC-format chrom.sizes file, control
file for peak calling, and a value to set the genome size.

For **Reference genomes** the `bowtie`, `chrSize`, `controlDNAPath`, and
`genomeSize` parameters must be set.

**Note:** calling peaks with MACS2 against a controlDNA file is probably not
the best control for CUT&RUN data. This is a holdover from other FAIRE and ATAC
pipelines we kept to account for copynumber variation, etc. The
`controlDNAPath` setting is therefore not required *per se*. If you wanted to
leave this out for now, simply set the value to `""` then set the `-k` flag in
`slurmSubmission.sh` and accept that the peak calling step will fail for all
samples. Alternatively, edit the peak calling step. (See [#29][i29])

For **Spike-in genomes** only the `bowtie` and `fasta` parameters are required.

**Adding New Genomes**
To quickly add a genome, copy an existing entry and paste it below another after the
`},` line. Edit the settings from there. To add a new spike-in genome, all that
is needed is a `fasta` entry pointing to a whole genome fasta file where each entry is 1 chromosome.

### module config
This setting describes the software versions for each software package used in the pipeline. UNC uses the Lmod package to manage packages, so these strings represent the `{package}` string when calling `module load {package}`. For adapting to a new computing environment, these values will need to be changed accordingly.

#### Local software dependencies
We run a pipeline step using [FastQ
Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) to
detect common contaminants. We use a locally installed version and thus point to
its path and the path of its config file in `fqscreenPath` and `fqscreenConf`,
respectively.

We use a local install of [mlr](https://github.com/johnkerl/miller) to parse
certain log files. Download the appropriate binary from the github release page
& set the location in `mlrPath` in `config.json`.

## Troubleshooting
Consider setting `--rerun-incomplete` in the Snakemake call in `slurmSubmission.sh` when initially implementing the pipeline on a new compute environment or testing new parameters.

Set `-k` in the Snakemake call in `slurmSubmission.sh` if certian steps fail only for 1 sample as this will cause the entire run to fail even if all other samples are good quality. Often this happens when a sample has few reads resulting in the peak calling steps failing. Because these samples have few reads, jobs complete quickly and cause the run to end early for other good samples. 

[i25]: https://github.com/snystrom/cutNrun-pipeline/issues/25
[i28]: https://github.com/snystrom/cutNrun-pipeline/issues/28
[i29]: https://github.com/snystrom/cutNrun-pipeline/issues/29
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[fastqc]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[multiqc]: https://multiqc.info/


## Evaluating CUT&RUN Experimental Success

The following is a non-exhaustive list of things to check after the pipeline run is finished. Suggestions involving exact numbers are more rules of thumb than precise cutoffs, so interpret your data as such.

1. Alignment Stats -- \>80% of reads align
3. FastqScreen does not identify any significant contaminating genomes
4. Most reads are kept post filtering
	- a significant fraction of reads filtered out can indicate poor library quality
	- high % of duplicates indicates low diversity libary
5. Determine fragment size distribution 
        - Check `Plots/FragDistInPeaks/`, compare to Bioanalyzer/Tapestation trace for each library
6. Determine fraction of reads in and out of peaks 
  	- \>15% FrIP seen in previous CUT&RUN experiments
	- ~3-5% typical for ChIP-seq (ENCODE guideline for ChIP is FrIP > 1%)
7. Signal correlation between replicates (inside peaks)
8. Profit	
