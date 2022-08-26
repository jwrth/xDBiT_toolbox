# ReadsToCounts pipeline

This pipeline transforms raw NGS sequencing reads from xDbit experiments into count matrices.

Input formats: `.fastq` or `.fastq.gz`

Output format: `txt.gz`

## Preparation of raw sequencing reads

### Demultiplexing

```
nohup bcl2fastq --runfolder-dir 220111_A00623_0449_BHN3G2DRXY/ --sample-sheet 220111_A00623_0449_BHN3G2DRXY/SampleSheet_HN3G2DRXY.csv --no-lane-splitting --minimum-trimmed-read-length=8 --mask-short-adapter-reads=8 --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -r 20 -w 20 -p 40 &> nohup_demult.out &
```

### Quality control

#### FastQC analysis

Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

```
cd /path/to/fastqs

# run FastQC in background
nohup fastqc -t 40 <fastq_folders>/*fastq.gz &> nohup_fastqc.out &
```

## Prepare barcode to coordinate assignment

The xDbit sequencing reads have following structure:
![read-structure](../graphics/xdbit_read-structure.png)

Each barcode on the read lies at a predetermined position relative to the beginning of the read and corresponds to a specific coordinate. In order to assign each read a specific X-, Y-, and if used Z-coordinate, the pipeline requires information about the barcode-to-coordinate assignents.

Two types of files are used to provide this information:

1. Well-to-coordinate information for each of the nine xDbit samples in `.csv` format
2. Barcode-to-well information in `barcode_legend_empty.csv` file.

Those two files are used to fill the empty barcode legend file with all information required by the pipeline.

The following section gives a more detailed overview of the file structure discussed above.

### 1. Well-to-coordinate assignment file

Examples of the files can be found under `./barcodes/well_coord_assignments/`

Each of the nine files has following structure:
| X_Coord | X_WellPosition | Y_Coord | Y_WellPosition | Z_Coord | Z_WellPosition |
|---------|----------------|---------|----------------|---------|----------------|
| 37      |                | 0       |                | 0       | A1             |
| 36      | B3             | 1       | B1             | 0       | C4             |
| 35      | C3             | 2       | C1             | 0       | D6             |
| 34      | D3             | 3       | D1             | 0       | F6             |
| 33      | A4             | 4       | A2             |         |                |
| …       | …              | …       | …              |         |                |

- ``X_Coord``: X-coordinate/column in the xDbit grid
- ``X_WellPosition``: Well position in barcoding well plate of the respective ligation round. Needs to correspond to the labels in the barcode legend file.
- ``Y_Coord``: Y-coordinate/row in xDbit grid
- ``Y_WellPosition``: Well position in barcoding well plate of the respective ligation round. needs to correspond to the labels in the barcode legend file.
- ``Z_Coord``: Sample ID.
- ``Z_WellPosition``: Well position of barcodes used for the RevT-indexing of this sample.

### 2. Barcode-to-well assignment file

Example can be found in `./barcodes/barcodes_legend_empty.csv`

Structure of the file:
| Row | Column | WellPosition | Name      | Barcode  | X | Y | Z | string_matching_algorithm | X_maxdist | Y_maxdist | Z_maxdist |
|-----|--------|--------------|-----------|----------|---|---|---|---------------------------|-----------|-----------|-----------|
| A   | 1      | A1           | Round1_01 | AACGTGAT |   |   |   | levenshtein               | 1         | 1         | 3         |
| B   | 1      | B1           | Round1_02 | AAACATCG |   |   |   |                           |           |           |           |
| C   | 1      | C1           | Round1_03 | ATGCCTAA |   |   |   |                           |           |           |           |
| D   | 1      | D1           | Round1_04 | AGTGGTCA |   |   |   |                           |           |           |           |
| E   | 1      | E1           | Round1_05 | ACCACTGT |   |   |   |                           |           |           |           |
| ... | ...    | ...          | ...       | ...      |   |   |   |                           |           |           |           |

- ``WellPosition``: Well position in barcoding well plates. Assumes that both ligation rounds use the same set of barcodes with the same positions in the well plates.
- ``Barcode``: Barcode in this well position.
- ``X/Y/Z``: Filled by `./src/fill_barcode_legend.py` script.
- ``string_matching_algorithm``: Method used to calculate the distance between read sequence and barcodes. Allows `hamming` or `levenshtein`.
- ``X/Y/Z_maxdist``: Maximum distance between read and barcode being allowed to be considered a match.


#### Summarize results with MultiQC

Documentation: https://multiqc.info/

```
# run multiqc after FastQC analysis
multiqc <fastq_folders>/
```

# Installation of tools needed for processing of raw reads and ReadsToCounts

## Quality control
### FastQC

`conda install -c bioconda fastqc`

### MultiQC

`pip install multiqc`

### Samtools
Download and instructions from: https://www.htslib.org/download/

```
# download samtools
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2

# unzip
tar -xf samtools-1.12.tar.bz2

cd samtools-1.x    # and similarly for bcftools and htslib
./configure --without-curses --prefix=/where/to/install
make
make install

# open.bashrc file using a text editor (here nano)
nano ~/.bashrc
# add following command at the end of .bashrc file
export PATH=/where/to/install/bin:$PATH

# save with CTRL+O and exit with CTRL+X
```

Use samtools for example with following command to show the head of a .bam file:
```
samtools view file.bam | head
```

### Cutadapt (gets already installed with the environment)
```
# install or upgrade cutadapt
conda activate dbitx_toolbox
python3 -m pip install --user --upgrade cutadapt
```

### STAR aligner

Manual on: https://github.com/alexdobin/STAR

#### Installation
```
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
tar -xzf 2.7.8a.tar.gz
cd STAR-2.7.8a

# Compile
cd STAR/source
make STAR

# Recommended: Add STAR to path variables
nano ~/.bashrc
# add following command at the end of .bashrc file
export PATH=$PATH:/path/to/software/STAR/STAR-2.7.Xa/bin/Linux_x86_64/

# On the HPC server it is installed and accessible for the Meier lab using following path variable: 
# export PATH=$PATH:/home/hpc/meier/software/STAR/STAR-2.7.4a/bin/Linux_x86_64/
```

### Drop-seq toolbox

This toolbox is based on the Drop-seq toolbox 2.1.0. The toolbox is included in this repository so it will be downloaded when cloning it.

Just the same it could be downloaded here:
https://github.com/broadinstitute/Drop-seq/releases/tag/v2.1.0

## Usage


### Create STAR metadata

Example for human genome:

```
# download and unzip reference genome and annotation file from ensembl
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
gzip -d Homo_sapiens.GRCh38.100.gtf.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


# run metadata script
./path/to/DbitX_toolbox/external_tools/Drop-seq_tools-2.1.0/create_Drop-seq_reference_metadata.sh \
-s Homo_sapiens -n Hs_metadata -r Homo_sapiens.GRCh38.dna.primary_assembly.fa -g Homo_sapiens.GRCh38.100.gtf \
-d ./path/to/DbitX_toolbox/external_tools/Drop-seq_tools-2.1.0
```

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

Examples for the .csv files can be find in the folder `/barcodes/`

### Run pipeline

```
# activate environment
conda activate dbitx_toolbox

# run pipeline
nohup bash /path/to/script/DbitX_pipeline.sh -g <GenomeDir> -r <ReferenceFasta> \
-b <BarcodeFile> -n <ExpectedNumberOfCells> -m <RunMode> -j <jobs> r1.fastq r2.fastq &
```

### Test run for use on HPC server in Meier lab

```
cd /path/to/repo/DbitX_toolbox/
cd test_files/

# create barcode .csv file
python ../src/fill_barcode_legend.py well_coord_assignment.csv barcodes_legend_empty.csv

# activate environment
conda activate dbitx_toolbox

# run test run
nohup bash ../DbitX_pipeline.sh -g /home/hpc/meier/genomes_STAR/mm_STARmeta/STAR/ -r /home/hpc/meier/genomes_STAR/mm_STARmeta/Mm_metadata.fasta -b ./barcodes_legend.csv -n 2300 -m DbitX -j 1 ./r1_100k.fastq ./r2_100k.fastq &
```

## Supplementary notes

### Generate fastq files from .bcl files

If the sequencer did not generate fastq files one can use bcl2fastq to generate the fastq files.
```
# run bcl2fastq without splitting the lanes
nohup /usr/local/bin/bcl2fastq --runfolder-dir 201130_NB552024_0025_AHG3GMAFX2/ --no-lane-splitting -r 20 -p 20 &
```

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
