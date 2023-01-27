#!/bin/bash

bedA=$1
bedB=$2

nameA=$(basename -s ".bed" $(basename -s ".gz" $bedA))
nameB=$(basename -s ".bed" $(basename -s ".gz" $bedB))

echo "Computing intersection of inputs..."
bedtools intersect -a $bedA -b $bedB -sortout > intersect.bed

echo "Computing Jaccard on original inputs:"
bedtools jaccard -a $bedA -b $bedB | column -t

# Compute proportion of A arising from overlap in B
echo "Computing % of ${nameA} in ${nameB}:"
bedtools jaccard -a $bedA -b intersect.bed | column -t 

# Compute proportion of B arising from overlap in A
echo "Computing % of ${nameB} in ${nameA}:"
bedtools jaccard -a intersect.bed -b $bedB | column -t

rm intersect.bed