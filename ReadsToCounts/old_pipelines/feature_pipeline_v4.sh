#!/usr/bin/env bash

#### Feature pipeline
# This script is a modified version of the Drop-seq_alignment.sh provided from Steve McCarroll's lab.
# Further we implemented parts of the Split-seq pipeline provided by Rebekka Wegmann.

# The algorithm does following steps:
# 1. Extraction of UMI, cellular barcodes and feature barcode.
# 2. Filtering of the tags based the known feature barcodes.
# 2. Generation of DGE matrix.

# Input file: .bam file

tmpdir=`pwd`
outdir=`pwd`
pipeline=0
clear=0
create_dge=
shorten_summary=
echo_prefix=
dropseq_root=$(dirname $0)/external_tools/Drop-seq_tools-2.1.0
estimated_num_cells=500
progname=`basename $0`
splitseq_root=$(dirname $0)
barcode_dir=$splitseq_root/data/barcode_lists

dist_alg=hamming
jobs=1
multithreading=

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Split-seq tagging, barcode filtering, alignment and digital expression matrix calculation

-r <rna_dge>        : RNA DGE matrix. Required if DGE matrix shall be generated.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: Subdirectory of the splitseq toolbox.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: current directory.
-b <barcode_dir>    : Full path to directory where the list of expected barcodes is stored. Default: subdirectory of the splitseq toolbox. 
-n <num_cells>      : Estimated number of cells in the library. Only affects visualization of barcode filtering results. Default: 500.
-e                  : Echo commands instead of executing them.  Cannot use with -p.
-a                  : String matching algorithm (hamming or levenshtein). Default: hamming.
-j                  : Number of threads. Default: 1.
-c                  : Delete unnecessary files?
-x                  : Create DGE matrix?
-s                  : Shorten summary?
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
while getopts ":d:t:o:r:esn:b:a:j:cxs" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    n ) estimated_num_cells=$OPTARG;;
    b ) barcode_dir=$OPTARG;;
    r ) rna_dge=$OPTARG;;
    e ) echo_prefix="echo";;
    a ) dist_alg=$OPTARG;;
    j ) jobs=$OPTARG;;
    c ) clear=1;;
    x ) create_dge="-x";;
    s ) shorten_summary="-s";;
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

if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi

if (( ${jobs} > 1 ))
then multithreading="-m"
else multithreading=""
fi

#Create output directories if they do not exist
if [[ ! -d $outdir ]]
then mkdir $outdir
fi

if [[ ! -d $tmpdir ]]
then mkdir $tmpdir
fi

input_bam=$1

## Stage 1: Extraction of UMI, cellular barcode and feature barcodes
start_time=`date +%s`
# Extract UMI (bases 1-10 of read2)
echo $input_bam

tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/feat_tagged_Molecular.bam_summary.txt \
    BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${input_bam}"

# Extract the 3 cellular barcodes
tag_cells_1="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/feat_tagged_Cellular_1.bam_summary.txt \
BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XD NUM_BASES_BELOW_QUALITY=1"

tag_cells_2="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/feat_tagged_Cellular_2.bam_summary.txt \
BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XE NUM_BASES_BELOW_QUALITY=1"

tag_cells_3="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/feat_tagged_Cellular_3.bam_summary.txt \
BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XF NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

tag_feature="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/feat_tagged_Cellular_3.bam_summary.txt \
BASE_RANGE=1-6 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XG NUM_BASES_BELOW_QUALITY=1"

#discard all reads where any one of the barcode regions has at least 1 base with quality < 10
filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

## Stage 3: Filter barcodes

# filter each split file
filter_barcodes="python ${splitseq_root}/src/feature_filtering_v4.py -t ${tmpdir} -d ${outdir} \
-n ${estimated_num_cells} -b ${barcode_dir} -a ${dist_alg} ${multithreading} \
-r ${rna_dge} ${shorten_summary} ${create_dge}"


# Stage 1
#$echo_prefix $tag_molecules OUTPUT=$tmpdir/feat_tagged_Molecular.bam
#$echo_prefix $tag_cells_1 INPUT=$tmpdir/feat_tagged_Molecular.bam OUTPUT=$tmpdir/feat_tagged_MC1.bam
#$echo_prefix $tag_cells_2 INPUT=$tmpdir/feat_tagged_MC1.bam OUTPUT=$tmpdir/feat_tagged_MC1C2.bam
#$echo_prefix $tag_cells_3 INPUT=$tmpdir/feat_tagged_MC1C2.bam OUTPUT=$tmpdir/feat_tagged_MC1C2C3.bam
#$echo_prefix $tag_feature INPUT=$tmpdir/feat_tagged_MC1C2C3.bam OUTPUT=$tmpdir/feat_tagged_MC1C2C3F.bam

#$echo_prefix $filter_bam INPUT=$tmpdir/feat_tagged_MC1C2C3F.bam OUTPUT=$tmpdir/feat_tagged_filtered.bam

# Stage 3

echo Command: $filter_barcodes -i ${tmpdir}/feat_tagged_filtered.bam
$echo_prefix $filter_barcodes -i ${tmpdir}/feat_tagged_filtered.bam

end_time=`date +%s`
run_time=`expr $end_time - $start_time`
total_time=`show_time $run_time`
echo "Split-seq pipeline finished using ${dist_alg} algorithm in ${total_time}"