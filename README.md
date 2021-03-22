# DbitX_toolbox
 
## Generate fastq files from .bcl files

```
# run bcl2fastq without splitting the lanes
nohup /usr/local/bin/bcl2fastq --runfolder-dir 201130_NB552024_0025_AHG3GMAFX2/ --no-lane-splitting -r 20 -p 20 &
```

## Quality control using FastQC

FastQC is an R package which performs a basic quality control on a sequencing run. Following code was used to run FastQC:

```
########################################
# Script to run FastQC
########################################

# set working directory for this notebook
setwd(<working_directory>)

# library
# install.packages("fastqcr")
library("fastqcr")

# How to install and update FastQC in case it's not installed yet
# fastqc_install()

## Quality Control with FastQC

fq_dir <- "<Fastq_directory>"
qc_dir <- paste(fq_dir, "qc", sep = "/") 

# run fastqc
fastqc(fq.dir = fq_dir,
       qc.dir = qc_dir,
       threads = 4)

#############
# HTML results are saved in output directory and can be viewed via filezilla
#############
```

For installation see: https://cran.r-project.org/web/packages/fastqcr/readme/README.html

## About Adapter Trimming

I thought about it and did some research on adapter trimming again and came up with following conclusions:

1. From https://thesequencingcenter.com/wp-content/uploads/2019/04/illumina-adapter-sequences.pdf:
	- "When read length exceeds DNA insert size, a run can sequence beyond the DNA insert and read bases fromthe sequencing adapter. To prevent these bases from appearing in FASTQ files, the adapter sequence is trimmed from the 3' ends of reads. Trimming the adapter sequence improves alignment accuracy andperformance in Illumina FASTQ generation pipelines."

2. From https://emea.support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html:
	- "if the sequencing extends beyond the length of the DNA insert, and into the adapter on the opposite end of the library fragment, that adapter sequence will be found on the 3’ end of the read. Therefore, reads require adapter trimming only on their 3’ ends."
	- Only trimming on 3' ends necessary.

In case of Split-seq/Dbit-seq Read 1 starts with the cDNA sequence and ends with the polyA sequence and the following barcodes etc. Therefore, it is necessary to do polyA trimming (as done in the Split-seq pipeline) and adapter trimming would not add any additional value.
Read 2 starts with the UMI and continues with the barcodes etc. From this read we extract anyway only the barcodes and don't need it for any alignment. So even if there were any adapters at the end they will not disturb the analysis.

### Final conclusion about Adapter Trimming:
- No adapter trimming necessary.
- Poly(A) trimming is enough.

However, since the FastQC for 37_28 gave 5-10 % adaptor content for "Illumina Universal Adapter" and "Nextera Transposase Sequence" I added adapter trimming to the cutadapt filtering step:
```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length ${min_lengths} -j 4 \
-o ${outdir}/R1_filtered.fastq.gz -p ${outdir}/R2_filtered.fastq.gz ${r1} ${r2} | tee ${tmpdir}/cutadapt_filter.out
```



