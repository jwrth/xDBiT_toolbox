#!/usr/bin/env bash
#
# This is a script to process Split-seq reads for both RNA and antibody reads. It is a modified version of \
# Drop-seq_alignment.sh provided alongside of  Dropseqtools from Steve McCarroll's lab
# and the Split-seq toolbox by Rebekka Wegmann from Snijder Lab at ETH Zurich.

#### Multimodal Split-seq pipeline
# This script is a modified version of the Drop-seq_alignment.sh provided from Steve McCarroll's lab.
# Further we implemented parts of the Split-seq pipeline provided by Rebekka Wegmann.

# The algorithm does following steps:
# 1. Extraction of UMI, cellular barcodes and feature barcode.
# 2. Filtering of the tags based the known feature barcodes.
# 2. Generation of DGE matrix.

# Input file: .bam file

#
# Author: Johannes Wirth, Meier lab, Helmholtz Zentrum Munich
# Original version Copyright 2017 Broad Institute


# MIT License
#
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
splitseq_root=$(dirname $0)/../
dropseq_root=${splitseq_root}/external_tools/Drop-seq_tools-2.1.0
star_executable=STAR
cutadapt_executable=cutadapt
estimated_num_cells=500
progname=`basename $0`
barcode_dir=${splitseq_root}/data/barcode_lists
jobs=1
feature_pipeline=1
shorten_summary=
dist_alg=levenshtein
multithreading=


function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Split-seq tagging, barcode filtering, alignment and digital expression matrix calculation for both RNA and feature reads.

-g <genomedir>      : Directory of STAR genome metadata. Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle. Required.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: Subdirectory of the splitseq toolbox.
-o <outputdir>      : Where to write the output folders.  Default: Current directory.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-b <barcode_dir>    : Full path to directory where the list of expected barcodes is stored. Default: Subdirectory of the splitseq toolbox. 
-n <num_cells>      : Estimated number of cells in the library. Only affects visualization of barcode filtering results. Default: 500.
-p                  : Reduce file I/O by pipeline commands of RNA pipeline together. Requires more memory and processing power.
-e                  : Echo commands instead of executing them. Cannot use with -p.
-a                  : String matching algorithm (hamming or levenshtein). Levenshtein is default.
-j                  : Number of threads for python filtering script in RNA pipeline. Default: 1.
-c                  : Delete unnecessary files.
-k                  : Shorten summary?

First input file: Input .bam file from RNA sample with both read 1 and read 2.
Second input file: Input .bam file from Feature sample with both read 1 and read 2.

If only one input file is given only the RNA pipeline is executed.

EOF
}

#getopts parses input options. Options followed by a : expect an input argument. The : at the very beginning prevents standard error messages.
while getopts ":d:o:pg:r:es:n:b:a:j:cxk" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    o ) outdir=$OPTARG;;
    n ) estimated_num_cells=$OPTARG;;
    b ) barcode_dir=$OPTARG;;
    p ) pipeline=1;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    e ) echo_prefix="echo";;
    a ) dist_alg=$OPTARG;;
    j ) jobs=$OPTARG;;
    c ) clear=1;;
	k ) shorten_summary="-s";;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))


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


### Checking inputs
if [[ "$pipeline" == 1 && -n "$echo_prefix" ]]
then error_exit "-p and -e cannot be used together"
fi

check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"

# Check input arguments
if (( $# == 0 ))
then error_exit "No input file given."
fi

if (( $# == 2 ))
then
	echo "Two input fastq files found. Only RNA pipeline will be run"
	feature_pipeline=0
elif (( $# == 4 ))
then 
    echo "Four input fastq files found. RNA and Feature pipeline will be run."
else
    error_exit "Other number of input files than 2 or 4 given."
fi

# Check location STAR executable
if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

# Check for multithreading
if (( ${jobs} > 1 ))
then multithreading="-m"
else multithreading=""
fi

# create directories for output and temporary files if they do not exist
rna_tmpdir=${outdir}/rna_tmp
rna_outdir=${outdir}/rna_out
feat_tmpdir=${outdir}/feature_tmp
feat_outdir=${outdir}/feature_out

if [[ ! -d $rna_tmpdir ]]
then mkdir ${rna_tmpdir}
fi

if [[ ! -d $rna_outdir ]]
then mkdir ${rna_outdir}
fi

if [[ ! -d $feat_tmpdir ]]
then mkdir ${feat_tmpdir}
fi

if [[ ! -d $feat_outdir ]]
then mkdir ${feat_outdir}
fi

# set up variables for files created during the pipeline
reference_suffix=$(echo $reference | sed s/.*\\./\\./) #reference can be .fa or .fasta
reference_basename=$(basename $reference $reference_suffix)
refflat=$(dirname $reference)/$reference_basename.refFlat
gene_intervals=$(dirname $reference)/$reference_basename.genes.intervals
exon_intervals=$(dirname $reference)/$reference_basename.exon.intervals
rRNA_intervals=$(dirname $reference)/$reference_basename.rRNA.intervals
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

tagged_unmapped_bam=${rna_tmpdir}/unaligned_tagged_BC_filtered.bam
aligned_sam=${rna_tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${rna_tmpdir}/aligned.sorted.bam
files_to_delete="${aligned_sorted_bam} ${aligned_sam} ${tagged_unmapped_bam}"


### PART 1: RNA
echo "Part 1: RNA matrix generation"
abs_start_time=`date +%s`

# Stage 0: Read fastq files and generate bam file
r1=$1
r2=$2

echo "RNA input read 1 file: ${r1}"
echo "RNA input read 2 file: ${r2}"

# filter fastq files for minimum length
min_lengths="35:94"
filter_fastq="${cutadapt_executable} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length ${min_lengths} -j 4 \
-o ${rna_outdir}/R1_filtered.fastq.gz -p ${rna_outdir}/R2_filtered.fastq.gz ${r1} ${r2}"

# generate .bam file
rna_unmapped_bam="${rna_outdir}/unmapped.bam"
generate_bam="java -jar ${picard_jar} FastqToSam F1=${rna_outdir}/R1_filtered.fastq.gz F2=${rna_outdir}/R2_filtered.fastq.gz O=${rna_unmapped_bam} SM=IntActpip TMP_DIR=${rna_tmpdir}"

## Stage 1: pre-alignment tag
# Extract UMI (bases 1-10 of read 2)
echo "RNA input file: ${rna_unmapped_bam}"

tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Molecular.bam_summary.txt \
    BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${rna_unmapped_bam}"

# Extract the 3 cellular barcodes
tag_cells_1="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_1.bam_summary.txt \
BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XD NUM_BASES_BELOW_QUALITY=1"

tag_cells_2="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_2.bam_summary.txt \
BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XE NUM_BASES_BELOW_QUALITY=1"

tag_cells_3="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${rna_outdir}/unaligned_tagged_Cellular_3.bam_summary.txt \
BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XF NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

# correct for the Ns that were added during padding
#correct_xq="python ${splitseq_root}/src/correct_xq.py"

#discard all reads where any one of the barcode regions has at least 1 base with quality < 10
filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

## Stage 2: Trim reads
#trim away adapter (template switching oligo)
trim_starting_sequence="${dropseq_root}/TrimStartingSequence OUTPUT_SUMMARY=${rna_outdir}/adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"

#trim anything that follows >= 6 contiguous As (assuming this is the polyA tail)
trim_poly_a="${dropseq_root}/PolyATrimmer OUTPUT_SUMMARY=${rna_outdir}/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"

## Stage 3: Filter barcodes
# split bam files for multiplexing
split_bam="bash ${splitseq_root}/src/splitbam.sh ${rna_tmpdir}/tmp_split ${jobs}"

# filter each split file
filter_barcodes="python ${splitseq_root}/SplitSeq/Splitseq_barcode_filtering_flex.py -t ${rna_tmpdir} -d ${rna_outdir} -n ${estimated_num_cells} -b ${barcode_dir} -a ${dist_alg} ${multithreading}"

# add header and merge all bam files for alignment
merge_filtered_bam="bash ${splitseq_root}/src/mergebam.sh ${rna_tmpdir}/tmp_split ${tagged_unmapped_bam}"

# Stage 4: alignment
sam_to_fastq="java -Xmx500m -jar ${picard_jar} SamToFastq INPUT=${tagged_unmapped_bam}"
star_align="$star_executable --genomeDir ${genomedir} --runThreadN 20 --quantMode GeneCounts --outFileNamePrefix ${rna_tmpdir}/star."

# Stage 5: Merge and tag BAM
# sort aligned reads in queryname order (STAR does not necessarily emit reads in the same order as the input)
sort_aligned="java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} \
SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${rna_tmpdir}"

# merge and tag aligned reads
merge_bam="java -Xmx4000m -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false"

# This one is a more flexible version of ta with gene exon, introduced in version 2.0.0 of Drop-seq tools
tag_with_gene_interval="${dropseq_root}/TagReadWithInterval INTERVALS=${gene_intervals} TAG=XG"
tag_with_gene_function="${dropseq_root}/TagReadWithGeneFunction O=${rna_outdir}/gene_function_tagged.bam ANNOTATIONS_FILE=${refflat}"


#### Start pipeline
start_time=`date +%s`

# Stage 0
#$echo_prefix $filter_fastq | tee ${rna_outdir}/cutadapt.out
# $echo_prefix $generate_bam | tee ${rna_tmpdir}/FastqToSam.out

# # Stage 1
# $echo_prefix $tag_molecules OUTPUT=$rna_tmpdir/unaligned_tagged_Molecular.bam
# $echo_prefix $tag_cells_1 INPUT=$rna_tmpdir/unaligned_tagged_Molecular.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MC1.bam
# $echo_prefix $tag_cells_2 INPUT=$rna_tmpdir/unaligned_tagged_MC1.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MC1C2.bam
# $echo_prefix $tag_cells_3 INPUT=$rna_tmpdir/unaligned_tagged_MC1C2.bam OUTPUT=$rna_tmpdir/unaligned_tagged_MC1C2C3.bam
# ####$echo_prefix $correct_xq -i ${rna_tmpdir}/unaligned_tagged_MC1C2C3.bam -o ${rna_tmpdir}/unaligned_tagged_MC1C2C3_corrected.bam
# $echo_prefix $filter_bam INPUT=$rna_tmpdir/unaligned_tagged_MC1C2C3.bam OUTPUT=$rna_tmpdir/unaligned_tagged_filtered.bam

# # Stage 2
# $echo_prefix $trim_starting_sequence INPUT=$rna_tmpdir/unaligned_tagged_filtered.bam OUTPUT=$rna_tmpdir/unaligned_tagged_trimmed_smart.bam
# $echo_prefix $trim_poly_a INPUT=$rna_tmpdir/unaligned_tagged_trimmed_smart.bam OUTPUT=$rna_tmpdir/unaligned_mc_tagged_polyA_filtered.bam

# files_to_delete="$files_to_delete $rna_tmpdir/unaligned_tagged_Molecular.bam $rna_tmpdir/unaligned_tagged_MC1.bam $rna_tmpdir/unaligned_tagged_MC1C2.bam $rna_tmpdir/unaligned_tagged_MC1C2C3.bam \
#                 $rna_tmpdir/unaligned_tagged_filtered.bam $rna_tmpdir/unaligned_tagged_trimmed_smart.bam"

# # Stage 3
# if [[ "${multithreading}" == "-m" ]]; then 
#     if [[ -d "${rna_tmpdir}/tmp_split" ]] && [[ ! -z "$(ls -A ${rna_tmpdir}/tmp_split)" ]]; then
#         echo "tmp_split exists but not empty. Remove all files now."
#         rm ${rna_tmpdir}/tmp_split/*
#         echo "All files in tmp_split removed."
#     fi

#     $echo_prefix $split_bam ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam
# fi

# echo "${filter_barcodes} -i ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam"
# $echo_prefix $filter_barcodes -i ${rna_tmpdir}/unaligned_mc_tagged_polyA_filtered.bam

# if [[ "${multithreading}" == "-m" ]]
# then $echo_prefix $merge_filtered_bam
# fi

# # exit 0 #break here, for testing only

# # Stage 4
# $echo_prefix $sam_to_fastq FASTQ=$rna_tmpdir/unaligned_tagged_BC_filtered.fastq
# $echo_prefix $star_align --readFilesIn $rna_tmpdir/unaligned_tagged_BC_filtered.fastq
# files_to_delete="$files_to_delete $rna_tmpdir/unaligned_tagged_BC_filtered.fastq"

# # Stage 5
# $echo_prefix $sort_aligned
# $echo_prefix $merge_bam OUTPUT=$rna_tmpdir/merged.bam
# $echo_prefix $tag_with_gene_interval I=$rna_tmpdir/merged.bam O=$rna_tmpdir/gene_tagged.bam TMP_DIR=${rna_tmpdir}
# $echo_prefix $tag_with_gene_function INPUT=$rna_tmpdir/merged.bam

# files_to_delete="$files_to_delete $rna_tmpdir/merged.bam $rna_tmpdir/gene_tagged.bam"

# ## Stage 6: create DGE matrix

# # counting exonic reads only
# dge="${dropseq_root}/DigitalExpression I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/DGE_matrix_min100.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_dir=${rna_tmpdir}"
# $echo_prefix $dge

# # counting both intronic and exonic reads
# dge_with_introns="${dropseq_root}/DigitalExpression I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/DGE_matrix_with_introns_min100.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 LOCUS_FUNCTION_LIST=INTRONIC TMP_dir=${rna_tmpdir}"
# $echo_prefix $dge_with_introns

# # collect RNAseq metrics with PICARD
# rnaseq_metrics="java -jar ${picard_jar} CollectRnaSeqMetrics I=${rna_outdir}/gene_function_tagged.bam O=${rna_outdir}/rnaseq_metrics.RNA_Metrics REF_FLAT=${refflat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=${rRNA_intervals}"
# $echo_prefix $rnaseq_metrics

# end_time=`date +%s`
# run_time=`expr $end_time + 1 - $start_time`
# total_time=`show_time $run_time`
# echo "RNA part finished using ${dist_alg} algorithm in ${total_time}"

if (( $feature_pipeline == 1 ))
then
	### PART 2: FEATURE
    start_time=`date +%s`
	rna_dge=${rna_outdir}/DGE_matrix_with_introns_min100.txt.gz
	echo "Part 2: Feature matrix generation"

    # Stage 0: Read fastq files and generate bam file
    r3=$3
    r4=$4

    echo "Feature input read 1 file: ${r3}"
    echo "Feature input read 2 file: ${r4}"

    # filter fastq files for minimum length
    min_lengths="42:94"
    filter_fastq="${cutadapt_executable} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length ${min_lengths} -j 4 \
    -o ${feat_outdir}/R1_filtered.fastq.gz -p ${feat_outdir}/R2_filtered.fastq.gz ${r3} ${r4}"

    # generate .bam file
    feat_input_bam="${feat_outdir}/unmapped.bam"
    generate_bam="java -jar ${picard_jar} FastqToSam F1=${feat_outdir}/R1_filtered.fastq.gz F2=${feat_outdir}/R2_filtered.fastq.gz O=${feat_input_bam} SM=DbitXpipe TMP_DIR=${feat_tmpdir}"

	## Stage 1: Extraction of UMI, cellular barcode and feature barcodes
	# Extract UMI (bases 1-10 of read2)
	echo "Feature input file: ${feat_input_bam}"

	tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Molecular.bam_summary.txt \
	    BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${feat_input_bam}"

	# Extract the 3 cellular barcodes
	tag_cells_1="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_1.bam_summary.txt \
	BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XD NUM_BASES_BELOW_QUALITY=1"

	tag_cells_2="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_2.bam_summary.txt \
	BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XE NUM_BASES_BELOW_QUALITY=1"

	tag_cells_3="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_3.bam_summary.txt \
	BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XF NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

	tag_feature1="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_3.bam_summary.txt \
	BASE_RANGE=1-6 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XG NUM_BASES_BELOW_QUALITY=1"

    tag_feature2="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${feat_outdir}/feat_tagged_Cellular_3.bam_summary.txt \
	BASE_RANGE=37-42 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XH NUM_BASES_BELOW_QUALITY=1"

	#discard all reads where any one of the barcode regions has at least 1 base with quality < 10
	filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

	## Stage 3: Filter barcodes

	# filter each split file
	filter_barcodes="python ${splitseq_root}/IntAct/IntAct_feature_filtering.py -t ${feat_tmpdir} -d ${feat_outdir} \
	-n ${estimated_num_cells} -b ${barcode_dir} -a ${dist_alg} ${multithreading} \
	-r ${rna_dge} ${shorten_summary} -x"

    # Stage 0
    # $echo_prefix $filter_fastq | tee ${feat_outdir}/cutadapt.out
    # $echo_prefix $generate_bam | tee ${feat_tmpdir}/FastqToSam.out

	# # Stage 1
	# $echo_prefix $tag_molecules OUTPUT=$feat_tmpdir/feat_tagged_Molecular.bam
	# $echo_prefix $tag_cells_1 INPUT=$feat_tmpdir/feat_tagged_Molecular.bam OUTPUT=$feat_tmpdir/feat_tagged_MC1.bam
	# $echo_prefix $tag_cells_2 INPUT=$feat_tmpdir/feat_tagged_MC1.bam OUTPUT=$feat_tmpdir/feat_tagged_MC1C2.bam
	# $echo_prefix $tag_cells_3 INPUT=$feat_tmpdir/feat_tagged_MC1C2.bam OUTPUT=$feat_tmpdir/feat_tagged_MC1C2C3.bam
	# $echo_prefix $tag_feature1 INPUT=$feat_tmpdir/feat_tagged_MC1C2C3.bam OUTPUT=$feat_tmpdir/feat_tagged_MC1C2C3F1.bam
    # $echo_prefix $tag_feature2 INPUT=$feat_tmpdir/feat_tagged_MC1C2C3F1.bam OUTPUT=$feat_tmpdir/feat_tagged_MC1C2C3F2.bam

	# $echo_prefix $filter_bam INPUT=$feat_tmpdir/feat_tagged_MC1C2C3F2.bam OUTPUT=$feat_tmpdir/feat_tagged_filtered.bam

	files_to_delete="$files_to_delete $feat_tmpdir/feat_tagged_Molecular.bam $feat_tmpdir/feat_tagged_MC1.bam $feat_tmpdir/feat_tagged_MC1C2.bam $feat_tmpdir/feat_tagged_MC1C2C3.bam \
	$feat_tmpdir/feat_tagged_MC1C2C3F1.bam $feat_tmpdir/feat_tagged_MC1C2C3F2.bam"

	# Stage 2

	echo "Filtering command: ${filter_barcodes} -i ${feat_tmpdir}/feat_tagged_filtered.bam"
	$echo_prefix $filter_barcodes -i ${feat_tmpdir}/feat_tagged_filtered.bam

	end_time=`date +%s`
	run_time=`expr $end_time + 1 - $start_time`
	total_time=`show_time $run_time`
	echo "Feature part finished using ${dist_alg} algorithm in ${total_time}"

	if (($clear == 1 ))
	then $echo_prefix rm $files_to_delete
	fi

	abs_end_time=`date +%s`
	abs_run_time=`expr ${abs_end_time} + 1 - ${abs_start_time}`
	abs_total_time=`show_time $abs_run_time`
	echo "Multimodal Split-seq pipeline finished using ${dist_alg} algorithm in ${abs_total_time}"
fi