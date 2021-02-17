#!/usr/bin/env bash
set -euo pipefail

# Usage:
# sh get_bam_reads_prefix.sh <bam> <prefix> <output_bam>
# For use when splitting reads out from aligning to a combined genome,
# where each species chromosomes are prefixed by {species}_
# e.g. dm6_chr2L, droYak2_chr2L, sacCer_chrV
# in the above example the prefix should be passed as dm6_ to select dm6 chromosomes.
#
# takes input bam file & selects all reads matching the <prefix>,
# then rewrites the header to remove the <prefix> from the entries.
# this results in a bam file which still has the other {species}_chr entries in the header,
# but the target prefixes are back to their native representation in the reference genome
# NOTE: input bam file must be indexed


bam=$1
prefix=$2
output=$3

filtered_bam=$bam".filter"

chrs=$(samtools idxstats $bam | cut -f1 | grep $prefix)

samtools view -b $bam $chrs > $filtered_bam

samtools view -H $filtered_bam | sed -e "s/$prefix//" | samtools reheader - $filtered_bam > $output

# Cleanup tmp file
rm $filtered_bam
