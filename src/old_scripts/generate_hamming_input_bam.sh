#!/bin/bash

# This tool takes R1 and R2 fastq files, filters out reads under a length of 94, 
# and generates a bam file from the two files.

prefix=$1
read1=$2
read2=$3
out_dir=`pwd`

r1_filtered=${out_dir}/${prefix}_R1_prefiltered.fastq
r2_filtered=${out_dir}/${prefix}_R2_prefiltered.fastq
output_bam=${out_dir}/${prefix}_unmapped.bam

# Filtering of short reads
cutadapt --nextseq-trim=20 --minimum-length 35:94 -j 20 -o ${r1_filtered} -p ${r2_filtered} ${read1} ${read2}

java -jar /home/hpc/meier/software/splitseq_toolbox/external_tools/Drop-seq_tools-2.1.0/3rdParty/picard/picard.jar FastqToSam \
F1=${r1_filtered} F2=${r2_filtered} O=${output_bam} SM=${prefix}

echo "Read 1 ${read1} and Read 2 ${read2} filtered and integrated into BAM file ${output_bam}"
