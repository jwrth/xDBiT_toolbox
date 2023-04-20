# ReadsToCounts pipeline

This pipeline transforms raw NGS sequencing reads from xDbit experiments into count matrices.

Input formats: `.fastq` or `.fastq.gz`

Output format: `txt.gz`

## Contents

1. [Installation of pipeline](#installation-of-pipeline)
2. [Test run](#test-run)
3. [Preparations](#preparations)
4. [Run ReadsToCounts pipeline](#run-readstocounts-pipeline)
5. [Installation of external tools](#installation-of-external-tools)

# Installation of pipeline
## Requirements

`ReadsToCounts` was tested on a Linux machine and the genome alignment steps using [STAR](https://github.com/alexdobin/STAR) require >30 GB RAM.

Installation of external tools is described [here](#installation-of-external-tools).

### Create conda environment

```
# create environment from file
conda env create -f ./ReadsToCounts.yml

# activate environment
conda activate ReadsToCounts
```

# Test run

To test whether the installation worked and the ReadsToCounts pipeline works on your device, it can be run without the alignment steps as follows. Settings for the test run are read from `./batch_parameters.csv`.

```
# test run
conda activate ReadsToCounts
nohup python xDbit_run_batch.py --skip_align --batch_file batch_parameters.csv &

# OR
nohup python xDbit_run_batch.py --skip_align --batch_file batch_parameters.xlsx &
```

The test run creates a `nohup.out` file in the current directory and for each batch specified in `./batch_parameters.csv` a log file in `./test_files/pipeline_batch{n}_{date}.out`. To check the test files one can use `cat` or `tail -f` to follow the addition of output in real time.

If the log file does not show any errors and the pipeline runs to the end all necessary packages beside the STAR alignment software is installed and all files, except for the STAR alignment files, are correctly formatted.

# Preparations
## Preparation of raw sequencing reads

### Demultiplexing

It is important to use the `--no-lane-splitting` option here.

```
nohup bcl2fastq --runfolder-dir 220111_A00623_0449_BHN3G2DRXY/ --sample-sheet 220111_A00623_0449_BHN3G2DRXY/SampleSheet_HN3G2DRXY.csv --no-lane-splitting --minimum-trimmed-read-length=8 --mask-short-adapter-reads=8 --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -r 20 -w 20 -p 40 &> nohup_demult.out &
```

## Quality control

### FastQC analysis

Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

```
cd /path/to/fastqs

# run FastQC in background
nohup fastqc -t 40 <fastq_folders>/*fastq.gz &> nohup_fastqc.out &
```

### Summarize results with MultiQC

Documentation: https://multiqc.info/

```
# run multiqc after FastQC analysis
multiqc <fastq_folders>/
```

## Prepare alignment

### Download genome

Example for human genome:

```
# download and unzip reference genome and annotation file from ensembl
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
gzip -d Homo_sapiens.GRCh38.100.gtf.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Create STAR metadata
```
# run metadata script
./path/to/xDbit_toolbox/external_tools/Drop-seq_tools-2.1.0/create_Drop-seq_reference_metadata.sh \
-s Homo_sapiens -n Hs_metadata -r Homo_sapiens.GRCh38.dna.primary_assembly.fa -g Homo_sapiens.GRCh38.100.gtf \
-d ./path/to/xDbit_toolbox/external_tools/Drop-seq_tools-2.1.0
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

### Fill empty barcode legend file

In this step a script is used to merge well-to-coordinate and well-to-barcode information to get barcode-to-coordinate information, which is required for the pipeline.

The script uses all `.csv` files it finds under `path/to/well_coord_assignments`

```
# fill the empty barcodes_legend'
python ./src/fill_barcode_legend.py path/to/barcodes_legend_empty.csv path/to/well_coord_assignments/
```
The output is saved into the input folder.

# Run ReadsToCounts pipeline

## Prepare batch parameter file

The pipeline can process multiple datasets in parallel. All required information is summarized in a batch parameter file.

An example for the batch parameter file is provided here: `./batch_parameters.csv`. Each row of the spreadsheet represent one sample.

Description of all mandatory parameters:
- ``batch``: Number of batch in which this sample should be processed. Allows integers starting from 1. How many samples can be run in one batch depends mainly on the RAM provided.
- ``fastq_R1/fastq_R2``: Path to fastq file of read 1 or read 2, respectively.
- ``legend``: Path to filled barcode legend file.
- ``STAR_genome``: Path to directory with STAR genome.
- ``STAR_fasta``: Path to STAR metadata `.fasta`
- ``expected_n_spots``: Expected number of spots. This is important to accelerate the calculation of certain QC outputs during analysis.
- ``mode``: Allows `xDbit` and `Dbit-seq`.
	- `xDbit`: Pipeline takes the Z barcode into account.
	- `Dbit-seq`: Pipeline ignores Z barcode.
- ``pipeline_dir``: Path to `xDbit_pipeline.sh` file.


## Run pipeline

```
# activate environment
conda activate ReadsToCounts

# go to working directory
cd path/to/working-directory

# run batch process in background
nohup python path/to/xDbit_run_batch.py --batch_file path/to/batch_parameter.csv & # output is printed to nohup.out
```

All output and temporary files are saved into folders `out`/`tmp` in the base directory of the `fastq_R1` path that was provided in the batch parameters file. Output generated during the pipeline is writtein into the log file `pipeline_{}.out` in the same folder.

# Installation of external tools

## FastQC

`conda install -c bioconda fastqc`

## MultiQC

`pip install multiqc`

## Samtools
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

## Cutadapt
```
# install or upgrade cutadapt
conda activate ReadsToCounts
python3 -m pip install --user --upgrade cutadapt
```

## STAR aligner

Manual on: https://github.com/alexdobin/STAR

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

## Drop-seq toolbox

This toolbox is based on the Drop-seq toolbox 2.1.0. The toolbox is included in this repository so it will be downloaded when cloning it.

Alternatively, it could be downloaded from here:
https://github.com/broadinstitute/Drop-seq/releases/tag/v2.1.0

# Supplementary notes

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
