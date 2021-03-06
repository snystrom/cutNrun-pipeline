# Release Notes
* 1.7.0
  * Spike-in methodology has been overhauled. Now, reads will be aligned to a combined genome of the reference + all spike-in genomes. Filtering & Duplicate marking are done on the combined alignment file which is subsequently split into separate bam files for each genome after filtering. Based on work by Jeanne-Marie McPherson demonstrating this strategy provides a more accurate estimate of experimental:spike-in reads.
    * The combined genome alignment strategy has a small side effect of adding extra chromosomes to each bam header corresponding to the spike-in genome chromosomes. These shouldn't cause downstream issues, but please file an issue report if this causes unintended breakage.
    * Default align step runtime increased to 1:30:00 from 1:00:00 to account for larger genomes
  * Users can now set `spikeGenome` to an array of genome entries to perform spike-in normalization to multiple genomes
    * To support these changes, spike-in normalized bigwig files are now denoted: "{species}-spikeNorm.bw" instead of merely "spikeNorm.bw"
    * The output sampleSheet columns referencing spikeNorm files now also contain the spike-in genome name i.e. "bigwig_allFrags_{species}_spikeNorm" vs old "bigwig_allFrags_spikeNorm"
  * No longer uses temp() anywhere, intermediate files will be kept at all steps. These steps are candidates for merging in the future to eliminate the need to keep intermediates.
  * Computes genome alignment stats in: `Logs/{sample}*_genomeStats.tsv`
    * added local software dependency: [mlr](https://github.com/johnkerl/miller), set in `mlrPath` in config.json
  * Pipeline now uses envmodules (no more manual calling module purge & load), and `--use-envmodules` is now set in `slurmSubmission.sh`. This enables easier adaptation of other module management systems like conda, docker & singularity.
  
  
* 1.6.4
  * Added computation of alignment stats
  * No longer uses temp() for some intermediate files.
  
