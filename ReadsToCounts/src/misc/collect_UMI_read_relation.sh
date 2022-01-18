#!/usr/bin/env bash

# This is to run the check_UMIs script over a list of CBCs
# The list has to be provided as a text file with one CBC per line

while read fraction;
do
	echo $fraction
	java -jar $PICARD DownsampleSam I=gene_function_tagged.bam P=$fraction M=downsample_metrics.txt O=/dev/stdout | python $HOME/nas_1/Data/Rebekka/Splitseq_mapping/check_UMIs_2.py "-"
done < fractions.txt
