#!/bin/bash

### Bash script for subsampling of bam files.
# First set seed. Second give list of percentages separated by spaces.

mkdir subsampling

file=$1
name="$(basename $file .bam)"

read -p "Set seed: " seed
read -p "Please give list of percentages for subsampling as integers: " list

for perc in $list
do
	echo "Subsample with percentage of ${perc}"
	samtools view -s ${seed}.${perc} -b $1 > subsampling/${name}_${perc}p.bam
done

echo "Subsampling finished. Files saved in subsampling/"