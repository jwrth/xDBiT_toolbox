#!/usr/bin/env bash
#
# This is a script to process xDbit reads. It is a modified version of \
# Drop-seq_alignment.sh provided alongside of  Dropseqtools from Steve McCarroll's lab.
# Also the code is based on work with a Copyright by Rebekka Wegmann (Snijderlab, ETH Zurich, 2019) which has been published in following GitHub repository: 
# https://github.com/RebekkaWegmann/splitseq_toolbox
#
# Author: Johannes Wirth, Meier Lab, Helmholtz Pioneer Campus, Helmholtz Zentrum Munich, 2021
# Original version Copyrights: 
#   - Copyright (c) 2018 Broad Institute
#   - Copyright (c) Rebekka Wegmann, 2019

# MIT License
#
# Copyright (c) Johannes Wirth, 2022
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


outdir=`pwd`
genomedir=
reference=
pipeline=0
clear=0
echo_prefix=
readstocounts_root=$(dirname $0)
dropseq_root=${readstocounts_root}/external_tools/Drop-seq_tools-2.1.0
star_executable=STAR
cutadapt_executable=cutadapt
estimated_num_cells=500
progname=`basename $0`
barcode_file=
jobs=1
mode=xDbit
feature_pipeline=1
shorten_summary=
exclude_align=0
spatial_only=0


function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform xDbit tagging, barcode filtering, alignment and digital expression matrix calculation for RNA and feature reads
-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: Subdirectory of the xDbit toolbox.
-o <outputdir>      : Where to write output bam.  Default: out.
-t <tmpdir>         : Where to write temporary files.  Default: tmp.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-c <cutadapt_path>  : Full path of cutadapt. Default: cutadapt is found via PATH environment variable.
-b <barcode_file>   : Full path to barcode legend.
-f <feature_legend> : Full path to feature legend.
-n <num_cells>      : Estimated number of cells in the library. Only affects visualization of barcode filtering results. Default: 500.
-p                  : Reduce file I/O by pipeline commands together.  Requires more memory and processing power.
-e                  : Echo commands instead of executing them.  Cannot use with -p.
-a                  : String matching algorithm (hamming or levenshtein). Default: hamming.
-j                  : Number of threads. Default: 1.
-l                  : Delete unnecessary files.
-m                  : Mode. "xDbit" (Searches for three barcodes) or "Dbit-seq" (Searches for two barcodes).
-u                  : Type of features: "antibody" or "interact".
-x <exclude_align>  : Excludes the alignment steps for testing purposes.
-y <spatial_only> : Run the pipeline only on R2 to analyze the spatial barcodes. Can be used for spillover analysis.
EOF
}

#Error functions
function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e # exit immediately if any of the commands in the script fail
set -o pipefail # Fail if any of the commands in a pipeline fails

#functions
function show_time () {
    num=$1
    min=0
    hour=0
    day=0
    if((num>59));then
        ((sec=num%60))
        ((num=num/60))
        if((num>59));then
            ((min=num%60))
            ((hour=num/60))
        else
            ((min=num))
        fi
    else
        ((sec=num))
    fi
    echo "$hour"h "$min"m "$sec"s
}

#getopts parses input options. Options followed by a : expect an input argument. The : at the very beginning prevents standard error messages.
while getopts ":d:t:o:n:b:f:pg:r:s:c:ej:lm:u:hxy" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    n ) estimated_num_cells=$OPTARG;;
    b ) barcode_file=$OPTARG;;
    f ) feature_file=$OPTARG;;
    p ) pipeline=1;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    c ) cutadapt_executable=$OPTARG;;
    e ) echo_prefix="echo";;
    j ) jobs=$OPTARG;;
    l ) clear=1;;
    m ) mode=$OPTARG;;
    u ) feat_mode=$OPTARG;;
    x ) exclude_align=1;;
    y ) spatial_only=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

#Checking inputs
if [[ "$pipeline" == 1 && -n "$echo_prefix" ]]
then error_exit "-p and -e cannot be used together"
fi

if [[ $spatial_only == 0 ]]
then
    check_set "$genomedir" "Genome directory" "-g"
    check_set "$reference" "Reference fasta"  "-r"
fi

# Check input arguments
if (( $# == 0 ))
then 
    error_exit "No input file given."
elif (($# == 1 ))
then
    if (( $spatial_only == 0 ))
    then error_exit "Only one input file and no `spatial_only` flag."
    fi
elif (( $# == 2 ))
then
	echo "Two input fastq files found. Only RNA pipeline will be run"
	feature_pipeline=0
elif (( $# == 4 ))
then 
    echo "Four input fastq files found. RNA and Feature pipeline will be run."
else
    error_exit "Other number of input files than 1, 2 or 4 given."
fi

if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

if (( ${jobs} > 1 ))
then multithreading="-m"
else 
    multithreading=""
fi

if [[ ${mode} == "xDbit" ]]
then
    if [[ $spatial_only == 1 ]]
    then
        min_lengths="94"
    else
        min_lengths="35:94"
    fi
    multiwell=1
    echo "${mode} mode: Expects three spatial barcodes (X, Y, Z)"
elif [[ ${mode} == "Dbit-seq" ]]
then
    if [[ $spatial_only == 1 ]]
    then
        min_lengths="56"
    else
        min_lengths="35:56"
    fi
    multiwell=0
    echo "${mode} mode: Expects two spatial barcodes (X, Y)"
else
    echo "ERROR: ${mode} is an invalid variable for mode. (Valid: 'xDbit'/'Dbit-seq')"
    exit 1
fi




# Check if all necessary python packages are loaded
python -c "import sys, os, timeit, h5py, pysam, itertools, numpy, pandas, matplotlib, argparse, tqdm, datetime, subprocess, glob, collections, multiprocessing, Levenshtein, json"
result="$?"
if [ "$result" -ne 0 ]; then
    echo "Not all Python packages installed. Script interrupted."
    exit 1
else
    echo "All Python packages installed."
fi

# create directories for output and temporary files if they do not exist
rna_tmpdir=${outdir}/rna_tmp
rna_outdir=${outdir}/rna_out

if [[ ! -d $rna_tmpdir ]]
then mkdir -p ${rna_tmpdir}
fi

if [[ ! -d $rna_outdir ]]
then mkdir -p ${rna_outdir}
fi

if [[ $feature_pipeline == 1 ]]
then
    feat_tmpdir=${outdir}/feature_tmp
    feat_outdir=${outdir}/feature_out

    if [[ ! -d $feat_tmpdir ]]
    then mkdir -p ${feat_tmpdir}
    fi

    if [[ ! -d $feat_outdir ]]
    then mkdir -p ${feat_outdir}
    fi
fi

if (( $spatial_only == 0 ))
then
    reference_suffix=$(echo $reference | sed s/.*\\./\\./) #reference can be .fa or .fasta
    reference_basename=$(basename $reference $reference_suffix)
    refflat=$(dirname $reference)/$reference_basename.refFlat
    gene_intervals=$(dirname $reference)/$reference_basename.genes.intervals
    exon_intervals=$(dirname $reference)/$reference_basename.exon.intervals
    rRNA_intervals=$(dirname $reference)/$reference_basename.rRNA.intervals
fi

picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

tagged_unmapped_bam=${rna_tmpdir}/unaligned_tagged_BC_filtered.bam
aligned_sam=${rna_tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${rna_tmpdir}/aligned.sorted.bam
files_to_delete="${aligned_sorted_bam} ${aligned_sam} ${tagged_unmapped_bam}"

### PART 1: RNA
echo "Part 1: RNA matrix generation"
abs_start_time=`date +%s`

# Read fastq files
r1=$1
r2=$2

echo "RNA input fastq file read 1: ${r1}"
echo "RNA input fastq file read 2: ${r2}"

## Stage 0: Filter .fastq files and generate .bam file

# filter fastq files for minimum length and generate bam file
rna_unmapped_bam="${rna_outdir}/unmapped.bam"
if [[ $spatial_only == 1 ]]
then
    filter_fastq="${cutadapt_executable} -a CTGTCTCTTATACACATCTGACGCTGCCGACGA \
    --minimum-length ${min_lengths} -j 4 \
    -o ${rna_outdir}/R2_filtered.fastq.gz ${r1}"

    # generate .bam file
    generate_bam="java -jar ${picard_jar} FastqToSam F1=${rna_outdir}/R2_filtered.fastq.gz \
    O=${rna_unmapped_bam} SM=xDbitpipe TMP_DIR=${rna_tmpdir}"

    # specify on which read the spatial coordinates are
    spatial_coord_read=1
else
    filter_fastq="${cutadapt_executable} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
    --minimum-length ${min_lengths} -j 4 \
    -o ${rna_outdir}/R1_filtered.fastq.gz -p ${rna_outdir}/R2_filtered.fastq.gz ${r1} ${r2}"

    # generate .bam file
    generate_bam="java -jar ${picard_jar} FastqToSam F1=${rna_outdir}/R1_filtered.fastq.gz F2=${rna_outdir}/R2_filtered.fastq.gz \
    O=${rna_unmapped_bam} SM=xDbitpipe TMP_DIR=${rna_tmpdir}"

    # specify on which read the spatial coordinates are
    spatial_coord_read=2
fi

## Stage 1: pre-alignment tag
# Extract UMI (Bases 1-10 on Read 2)
echo "RNA input .bam file: ${rna_unmapped_bam}"
tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Molecular.bam_summary.txt \
    BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=${spatial_coord_read} DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${rna_unmapped_bam}"

# Extract the spatial barcodes
tag_cells_well="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_1.bam_summary.txt \
BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=${spatial_coord_read} DISCARD_READ=false TAG_NAME=XZ NUM_BASES_BELOW_QUALITY=1"

tag_cells_y="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_2.bam_summary.txt \
BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=${spatial_coord_read} DISCARD_READ=false TAG_NAME=XY NUM_BASES_BELOW_QUALITY=1"

tag_cells_x="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_3.bam_summary.txt \
BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=${spatial_coord_read} DISCARD_READ=true TAG_NAME=XX NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

#discard all reads where any one of the barcode regions has at least 1 base with quality < 10
filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

## Stage 2: Trim reads
#trim anything that follows >= 6 contiguous As (assuming this is the polyA tail)
trim_poly_a="${dropseq_root}/PolyATrimmer OUTPUT_SUMMARY=${rna_outdir}/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"

## Stage 3: Filter barcodes
# split bam files for multiplexing
split_bam="bash ${readstocounts_root}/src/splitbam.sh ${rna_tmpdir}/tmp_split ${jobs}"

# filter each split file
filter_barcodes="python ${readstocounts_root}/xDbit_filtering.py \
--mode ${mode} -t ${rna_tmpdir} -d ${rna_outdir} -n ${estimated_num_cells} \
--stride 500000 -b ${barcode_file} ${multithreading}"

# add header and merge all bam files for alignment
merge_filtered_bam="bash ${readstocounts_root}/src/mergebam.sh ${rna_tmpdir}/tmp_split ${tagged_unmapped_bam}"

# Stage 4: alignment
sam_to_fastq="java -Xmx500m -jar ${picard_jar} SamToFastq INPUT=${tagged_unmapped_bam} TMP_DIR=${rna_tmpdir}"
star_align="$star_executable --genomeDir ${genomedir} --runThreadN 20 --quantMode GeneCounts --outFileNamePrefix ${rna_tmpdir}/star."

# Stage 5: Merge and tag BAM
# sort aligned reads in queryname order (STAR does not necessarily emit reads in the same order as the input)
sort_aligned="java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} \
SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${rna_tmpdir}"

# merge and tag aligned reads
merge_bam="java -Xmx4000m -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false TMP_DIR=${rna_tmpdir}"

# This one is a more flexible version of ta with gene exon, introduced in version 2.0.0 of Drop-seq tools
tag_with_gene_interval="${dropseq_root}/TagReadWithInterval INTERVALS=${gene_intervals} TAG=XG"
tag_with_gene_function="${dropseq_root}/TagReadWithGeneFunction O=${rna_outdir}/gene_function_tagged.bam ANNOTATIONS_FILE=${refflat}"


#### Start pipeline
start_time=`date +%s`

# Stage 0
$echo_prefix $filter_fastq | tee ${rna_outdir}/cutadapt.out
$echo_prefix $generate_bam | tee ${rna_tmpdir}/FastqToSam.out

# Stage 1
$echo_prefix $tag_molecules OUTPUT=${rna_tmpdir}/unaligned_tagged_Molecular.bam

if (( $multiwell == 1))
then
    $echo_prefix $tag_cells_well INPUT=$rna_tmpdir/unaligned_tagged_Molecular.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MW.bam
    $echo_prefix $tag_cells_y INPUT=$rna_tmpdir/unaligned_tagged_MW.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MWY.bam
    $echo_prefix $tag_cells_x INPUT=$rna_tmpdir/unaligned_tagged_MWY.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MWYX.bam
    $echo_prefix $filter_bam INPUT=$rna_tmpdir/unaligned_tagged_MWYX.bam OUTPUT=$rna_tmpdir/unaligned_tagged_filtered.bam

    files_to_delete="$files_to_delete $rna_tmpdir/unaligned_tagged_Molecular.bam \
                $rna_tmpdir/unaligned_tagged_MW.bam $rna_tmpdir/unaligned_tagged_MWY.bam \
                $rna_tmpdir/unaligned_tagged_MWYX.bam $rna_tmpdir/unaligned_tagged_filtered.bam"
else
    $echo_prefix $tag_cells_y INPUT=$rna_tmpdir/unaligned_tagged_Molecular.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MY.bam
    $echo_prefix $tag_cells_x INPUT=$rna_tmpdir/unaligned_tagged_MY.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MYX.bam
    $echo_prefix $filter_bam INPUT=$rna_tmpdir/unaligned_tagged_MYX.bam OUTPUT=$rna_tmpdir/unaligned_tagged_filtered.bam

    files_to_delete="$files_to_delete $rna_tmpdir/unaligned_tagged_Molecular.bam \
                $rna_tmpdir/unaligned_tagged_MY.bam $rna_tmpdir/unaligned_tagged_MYX.bam \
                $rna_tmpdir/unaligned_tagged_filtered.bam"
fi

    # Stage 2
    $echo_prefix $trim_poly_a INPUT=$rna_tmpdir/unaligned_tagged_filtered.bam OUTPUT=$rna_tmpdir/unaligned_mc_tagged_polyA_filtered.bam

    # Stage 3
    if [[ "${multithreading}" == "-m" ]]; then 
        if [[ -d "${tmpdir}/tmp_split" ]] && [[ ! -z "$(ls -A ${rna_tmpdir}/tmp_split)" ]]; then
            echo "tmp_split exists but not empty. Remove all files now."
            rm ${tmpdir}/tmp_split/*
            echo "All files in tmp_split removed."
        fi
        $echo_prefix $split_bam ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam
    fi
    
    echo "${filter_barcodes} -i ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam"
    $echo_prefix $filter_barcodes -i ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam

    if [[ "${multithreading}" == "-m" ]]
    then $echo_prefix $merge_filtered_bam
    fi

if [[ $spatial_only == 1 ]]
then
    sleep 2 # necessary to make the run time calculation valid if the pipeline runs in < 1 second until here.
    end_time=`date +%s`
    run_time=`expr $end_time - $start_time`
    total_time=`show_time $run_time`

    echo "Spillover analysis finished in ${total_time}."
else
    # Continue with analysis pipeline
    # Stage 4
    $echo_prefix $sam_to_fastq FASTQ=$rna_tmpdir/unaligned_tagged_BC_filtered.fastq

    if (($exclude_align == 1))
    then
        echo
        echo ">>>TEST MODE<<<"
        echo "Info: exclude_align flag added. All steps starting from the alignment are only printed to the output."
        echo
        echo_prefix="echo"
    fi

    $echo_prefix $star_align --readFilesIn $rna_tmpdir/unaligned_tagged_BC_filtered.fastq
    files_to_delete="$files_to_delete $rna_tmpdir/unaligned_tagged_BC_filtered.fastq"
    # Stage 5
    $echo_prefix $sort_aligned
    $echo_prefix $merge_bam OUTPUT=$rna_tmpdir/merged.bam
    $echo_prefix $tag_with_gene_interval I=$rna_tmpdir/merged.bam O=$rna_tmpdir/gene_tagged.bam TMP_DIR=${rna_tmpdir}
    $echo_prefix $tag_with_gene_function INPUT=$rna_tmpdir/merged.bam
    files_to_delete="$files_to_delete $rna_tmpdir/merged.bam $rna_tmpdir/gene_tagged.bam"


    ## Stage 6: create DGE matrix
    # counting exonic reads only
    dge="${dropseq_root}/DigitalExpression I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/DGE_matrix_min100.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_dir=${rna_tmpdir}"
    $echo_prefix $dge

    # counting both intronic and exonic reads
    dge_with_introns="${dropseq_root}/DigitalExpression I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/DGE_matrix_with_introns_min100.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 LOCUS_FUNCTION_LIST=INTRONIC TMP_dir=${rna_tmpdir}"
    $echo_prefix $dge_with_introns

    # collect RNAseq metrics with PICARD
    rnaseq_metrics="java -jar ${picard_jar} CollectRnaSeqMetrics I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/rnaseq_metrics.RNA_Metrics REF_FLAT=${refflat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=${rRNA_intervals}"
    $echo_prefix $rnaseq_metrics


    sleep 2 # necessary to make the run time calculation valid if the pipeline runs in < 1 second until here.
    end_time=`date +%s`
    run_time=`expr $end_time - $start_time`
    total_time=`show_time $run_time`

    echo "xDbit RNA pipeline finished in ${total_time}"


    if (( $feature_pipeline == 1 ))
    then
        ### PART 2: FEATURE
        start_time=`date +%s`
        rna_dge=${rna_outdir}/DGE_matrix_with_introns_min100.txt.gz
        echo "Part 2: Feature matrix generation"

        if [[ ${feat_mode} == "interact" ]]
        then
            r1_min="42"
        elif [[ ${feat_mode} == "antibody" ]]
        then
            r1_min="6"
        else
            echo "ERROR: ${feat_mode} is an invalid variable for feature_mode. (Valid: 'antibody'/'interact'"
        fi

        if [[ ${mode} == "xDbit" ]]
        then
            min_lengths="${r1_min}:94"
            multiwell=1
            echo "${mode} mode: Expects three spatial barcodes (X, Y, Z)"
        elif [[ ${mode} == "Dbit-seq" ]]
        then
            min_lengths="${r1_min}:56"
            multiwell=0
            echo "${mode} mode: Expects two spatial barcodes (X, Y)"
        else
            echo "ERROR: ${mode} is an invalid variable for mode. (Valid: 'xDbit'/'Dbit-seq')"
            exit 1
        fi

        # Stage 0: Read fastq files and generate bam file
        r3=$3
        r4=$4

        echo "Feature input read 1 file: ${r3}"
        echo "Feature input read 2 file: ${r4}"

        # filter fastq files for minimum length
        filter_fastq="${cutadapt_executable} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length ${min_lengths} -j 4 \
        -o ${feat_outdir}/R1_filtered.fastq.gz -p ${feat_outdir}/R2_filtered.fastq.gz ${r3} ${r4}"

        # generate .bam file
        feat_input_bam="${feat_outdir}/unmapped.bam"
        generate_bam="java -jar ${picard_jar} FastqToSam F1=${feat_outdir}/R1_filtered.fastq.gz F2=${feat_outdir}/R2_filtered.fastq.gz O=${feat_input_bam} SM=xDbitpipe TMP_DIR=${feat_tmpdir}"

        ## Stage 1: Extraction of UMI, cellular barcode and feature barcodes
        # Extract UMI (bases 1-10 of read2)
        echo "Feature input file: ${feat_input_bam}"

        tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Molecular.bam_summary.txt \
            BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${feat_input_bam}"

        # Extract the spatial barcodes
        tag_cells_well="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_1.bam_summary.txt \
        BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XZ NUM_BASES_BELOW_QUALITY=1"

        tag_cells_y="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_2.bam_summary.txt \
        BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XY NUM_BASES_BELOW_QUALITY=1"

        tag_cells_x="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_3.bam_summary.txt \
        BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XX NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

        tag_feature="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Features.bam_summary.txt \
        BASE_RANGE=1-6 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XG NUM_BASES_BELOW_QUALITY=1"

        tag_interact="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Interactions.bam_summary.txt \
        BASE_RANGE=37-42 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XH NUM_BASES_BELOW_QUALITY=1"

        # discard all reads where any one of the barcode regions has at least 1 base with quality < 10
        filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

        ## Stage 3: Filter barcodes

        # filter each split file
        filter_barcodes="python ${readstocounts_root}/xDbit_filtering.py --mode ${mode} -t ${feat_tmpdir} -d ${feat_outdir} \
        -n ${estimated_num_cells} -b ${barcode_file} ${multithreading} --stride 500000 \
        -f ${feature_file} -r ${rna_dge}"

        if [[ ${feat_mode} == "interact" ]]
        then
            filter_barcodes="${filter_barcodes} --interact"
        fi

        # Stage 0
        $echo_prefix $filter_fastq | tee ${feat_outdir}/cutadapt.out
        $echo_prefix $generate_bam | tee ${feat_tmpdir}/FastqToSam.out

        # # Stage 1
        $echo_prefix $tag_molecules OUTPUT=$feat_tmpdir/feat_tagged_Molecular.bam

        if (( $multiwell == 1))
        then
            
            $echo_prefix $tag_cells_well INPUT=$feat_tmpdir/feat_tagged_Molecular.bam OUTPUT=$feat_tmpdir/feat_tagged_MW.bam
            $echo_prefix $tag_cells_y INPUT=$feat_tmpdir/feat_tagged_MW.bam OUTPUT=$feat_tmpdir/feat_tagged_MWY.bam
            $echo_prefix $tag_cells_x INPUT=$feat_tmpdir/feat_tagged_MWY.bam OUTPUT=$feat_tmpdir/feat_tagged_MWYX.bam
            $echo_prefix $tag_feature INPUT=$feat_tmpdir/feat_tagged_MWYX.bam OUTPUT=$feat_tmpdir/feat_tagged_MWYXF.bam

            echo "Feature mode: ${feat_mode}"
            if [[ ${feat_mode} == "interact" ]]
            then
                $echo_prefix $tag_interact INPUT=$feat_tmpdir/feat_tagged_MWYXF.bam OUTPUT=$feat_tmpdir/feat_tagged_MWYXFI.bam
                $echo_prefix $filter_bam INPUT=$feat_tmpdir/feat_tagged_MWYXFI.bam OUTPUT=$feat_tmpdir/feat_tagged_filtered.bam
            else
                $echo_prefix $filter_bam INPUT=$feat_tmpdir/feat_tagged_MWYXF.bam OUTPUT=$feat_tmpdir/feat_tagged_filtered.bam
            fi

            files_to_delete="$files_to_delete $feat_tmpdir/feat_tagged_Molecular.bam \
                    $feat_tmpdir/feat_tagged_MW.bam $feat_tmpdir/feat_tagged_MWY.bam \
                    $feat_tmpdir/feat_tagged_MWYX.bam $feat_tmpdir/feat_tagged_MWYXF.bam"
        else
            $echo_prefix $tag_cells_y INPUT=$feat_tmpdir/feat_tagged_Molecular.bam OUTPUT=$feat_tmpdir/feat_tagged_MY.bam
            $echo_prefix $tag_cells_x INPUT=$feat_tmpdir/feat_tagged_MY.bam OUTPUT=$feat_tmpdir/feat_tagged_MYX.bam
            $echo_prefix $tag_feature INPUT=$feat_tmpdir/feat_tagged_MYX.bam OUTPUT=$feat_tmpdir/feat_tagged_MYXF.bam
            $echo_prefix $filter_bam INPUT=$feat_tmpdir/feat_tagged_MYXF.bam OUTPUT=$feat_tmpdir/feat_tagged_filtered.bam

            files_to_delete="$files_to_delete $feat_tmpdir/feat_tagged_Molecular.bam \
                        $feat_tmpdir/feat_tagged_MY.bam $feat_tmpdir/feat_tagged_MYX.bam \
                        $feat_tmpdir/feat_tagged_MYXF.bam"
        fi

        # Stage 2
        echo "Filtering command: ${filter_barcodes} -i ${feat_tmpdir}/feat_tagged_filtered.bam"
        $echo_prefix $filter_barcodes -i ${feat_tmpdir}/feat_tagged_filtered.bam

        end_time=`date +%s`
        run_time=`expr $end_time + 1 - $start_time`
        total_time=`show_time $run_time`
        echo "Feature part finished in ${total_time}"

        abs_end_time=`date +%s`
        abs_run_time=`expr ${abs_end_time} + 1 - ${abs_start_time}`
        abs_total_time=`show_time $abs_run_time`
        echo "xDbit feature pipeline finished in ${abs_total_time}"
    fi

    if (($clear == 1 ))
    then
        echo "Delete temporary files." 
        $echo_prefix rm $files_to_delete
    fi
fi