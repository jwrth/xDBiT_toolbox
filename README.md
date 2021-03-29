# DbitX toolbox

This toolbox includes an analysis pipeline to compute the digital gene expression matrix from Dbit-seq or DbitX experiments.
It is based on the Split-seq toolbox: https://github.com/RebekkaWegmann/splitseq_toolbox.

The pipeline uses a bash script, custom Python scripts, and many tools from the Drop-seq toolbox (Mc Caroll lab, Harvard Medical school) as well as Picard (Broad institute), which are all included in this toolbox.

## Introduction

### Barcoding layout
![layout](https://user-images.githubusercontent.com/76480183/112812833-d859af00-907d-11eb-8c9c-613ac9166099.png)

### Read structure

![read](https://user-images.githubusercontent.com/76480183/112814912-0e982e00-9080-11eb-8ffe-d2c17c4660d8.png)



## Requirements

### 1. Drop-seq toolbox

This toolbox is based on the Drop-seq toolbox 2.1.0 which can be also downloaded here:
https://github.com/broadinstitute/Drop-seq/releases/tag/v2.1.0

### 2. Samtools
Download and instructions from: https://www.htslib.org/download/

```
# download samtools
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2

# unzip
tar -xf samtools-1.12.tar.bz2

# install
cd samtools-1.x    # and similarly for bcftools and htslib
./configure --without-curses --prefix=/where/to/install
make
make install

# add permanently to path variable
nano ~/.bashrc

# add following command at the end of .bashrc file
export PATH=/where/to/install/bin:$PATH

# save with CTRL+O and exit with CTRL+X
```

Use samtools for example with following command to show the head of a .bam file:
```
samtools view file.bam | head
```

### 3. Cutadapt
```
# create python3 environment (cutadapt needs python 3 to use multiple cores)
conda create -n dbitx_toolbox python=3
conda activate dbitx_toolbox

# install cutadapt
python3 -m pip install --user --upgrade cutadapt

# to access parts of the pipeline in a Jupyter notebook install a kernel for this environment
conda install -c anaconda ipykernel
python3 -m ipykernel install --user --name=dbitx_toolbox_kernel
```

## Usage

### Create barcode legend file

#### 1. Distribution of barcodes in wells

The `barcode_legend_empty.csv` file links barcode sequences to the well positions and three columns `X`, `Y`, and (optionally) `Z` which correspond to the dimensions that were barcoded in the experiment. Mandatory columns here are `WellPosition`, `Barcode`, `X`, `Y`. `Z` only if three barcoding rounds were applied.

| Row | Column | WellPosition | Name      | Barcode  | X | Y | Z |
|-----|--------|--------------|-----------|----------|---|---|---|
| A   | 1      | A1           | Round1_01 | AACGTGAT |   |   |   |
| B   | 1      | B1           | Round1_02 | AAACATCG |   |   |   |
| C   | 1      | C1           | Round1_03 | ATGCCTAA |   |   |   |
| D   | 1      | D1           | Round1_04 | AGTGGTCA |   |   |   |

#### 2. Assignment of well position to spatial coordinate

The `well_coord_assignment.csv` file contains three pairs of columns which assign well positions to spatial coordinates as shown in the following table.

| X_Coord | X_WellPosition | Y_Coord | Y_WellPosition | Z_Coord | Z_WellPosition |
|---------|----------------|---------|----------------|---------|----------------|
| 49      |                | 0       |                | 0       | A1             |
| 48      | B1             | 1       | B1             | 0       | C4             |
| 47      | C1             | 2       | C1             | 0       | D6             |
| 46      | D1             | 3       | D1             | 0       | F6             |

#### 3. Filling of the barcode legend file

To fill the columns `X`, `Y` (, `Z`) run following python script. The output is saved as `barcodes_legend.csv` to the input folder.

```
# fill barcode legend file
python /path/to/script/fill_barcode_legend.py well_coord_assignment.csv barcodes_legend_empty.csv
```

### Run pipeline

```
# activate environment
conda activate dbitx_toolbox

# run pipeline
nohup bash /path/to/script/DbitX_pipeline.sh -g <GenomeDir> -r <ReferenceFasta> \
-b <BarcodeDir> -n <ExpectedNumberOfCells> -m <RunMode> -j <jobs> r1.fastq r2.fastq &
```

## Supplementary notes

### Quality control using FastQC

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

### Create test files

```
mkdir test_files

# create test files for read 1 and read 2
gzip -d -c 210319_NB552024_0035_AHCJ33AFX2/Data/Intensities/BaseCalls/37_28_Dbitseq_S1_R1_001.fastq.gz | head -n 400000 > test_files/r1_100k.fastq

gzip -d -c 210319_NB552024_0035_AHCJ33AFX2/Data/Intensities/BaseCalls/37_28_Dbitseq_S1_R2_001.fastq.gz | head -n 400000 > test_files/r2_100k.fastq

```

### About Adapter Trimming

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
