#!/bin/bash

# This tool takes R1 and R2 fastq files, filters out reads under a specified length, 
# pads the R2 sequences to a specified length, and generates a bam file from the two files.

prefix=$1
read1=$2
read2=$3
out_dir=`pwd`

#read1_name=`basename ${read1} .fastq.gz`
#read2_name=`basename ${read2} .fastq.gz`

r1_filtered=${out_dir}/${prefix}_R1_prefiltered.fastq
r2_filtered=${out_dir}/${prefix}_R2_prefiltered.fastq
r2_pad=${out_dir}/${prefix}_R2_prefiltered_pad.fastq
output_bam=${out_dir}/${prefix}_pad_unmapped.bam

# Filtering of short reads
cutadapt --nextseq-trim=20 --minimum-length 35:90 -j 20 -o ${r1_filtered} -p ${r2_filtered} ${read1} ${read2}

python /home/hpc/meier/software/splitseq_toolbox/src/padfastq.py -l 94 -q "E" -d "N" -o ${r2_pad} -s 100000 < ${r2_filtered}

java -jar /home/hpc/meier/software/splitseq_toolbox/external_tools/Drop-seq_tools-2.1.0/3rdParty/picard/picard.jar FastqToSam \
F1=${r1_filtered} F2=${r2_pad} O=${output_bam} SM=${prefix}

# to save disk space the unpadded filtered version of read 2 is deleted
rm ${r2_filtered}

echo "Read 1 ${read1} and Read 2 ${read2} filtered, padded and integrated into BAM file ${output_bam}"
