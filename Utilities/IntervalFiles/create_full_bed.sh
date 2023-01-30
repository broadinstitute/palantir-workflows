#!/bin/bash

ref_dict=$1
output=$2

grep "@SQ" $ref_dict | cut -f 2 | cut -d":" -f 2 > contigs_full.txt
grep "@SQ" $ref_dict | cut -f 3 | cut -d":" -f 2 > contigs_full_len.txt

LENGTH=$(wc -l < contigs_full.txt)

yes 0 | head -n $LENGTH > zeros_full.txt

paste contigs_full.txt zeros_full.txt contigs_full_len.txt > $output

rm contigs_full.txt
rm contigs_full_len.txt
rm zeros_full.txt
