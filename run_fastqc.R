#!/usr/bin/env Rscript

########################################
# Script to run FastQC
#
# Input: Fastq files of Illumina sequencing run.
########################################

# Library
library("fastqcr")

# Fetch input directory
fq_dir = commandArgs(trailingOnly=TRUE)
qc_dir <- paste(fq_dir, "qc", sep = "/") 

# run fastqc
print("Run FastQC...")
print("HTML results are saved in 'qc' subfolder of sequencing directory.")

fastqc(fq.dir = fq_dir,
       qc.dir = qc_dir,
       threads = 4)

#############
# HTML results are saved in 'qc' subfolder of sequencing directory.
#############