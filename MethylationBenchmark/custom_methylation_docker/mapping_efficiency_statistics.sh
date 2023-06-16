#!/bin/bash

# $1 = original SAM/BAM file
# $2 = filtered SAM/BAM file, usually sorted at this stage

raw_total_sequences_original=$(samtools stats $1 | grep ^SN | cut -f 2- | awk 'NR==1 {print $4}')
raw_total_sequences_filtered=$(samtools stats $2 | grep ^SN | cut -f 2- | awk 'NR==1 {print $4}')

mapping_efficiency=$(echo "scale=3; $raw_total_sequences_filtered / $raw_total_sequences_original" | bc)

echo $mapping_efficiency
