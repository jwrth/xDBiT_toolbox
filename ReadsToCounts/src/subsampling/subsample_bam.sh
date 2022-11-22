#!/bin/bash

### Bash script for subsampling of bam files.
# First set seed. Second give list of percentages separated by spaces.

manual=0
overwrite=0
dge=0
subsampling_root=$(dirname $0)
dropseq_root=${subsampling_root}/../../external_tools/Drop-seq_tools-2.1.0

function print_usage () {
    cat >&2 <<EOF
Usage: Bash script for subsampling of bam files.
-m <manual>      : Set seed and subsample list manually.
-o <overwrite>	 : Overwrite existing output folders and files.
-d <dge>		 : Create DGE matrix from subsampled bams.
EOF
}



while getopts "mod" flag; do
  case $flag in
    m ) manual=1;;
	o ) overwrite=1;;
	d ) dge=1;;
	h ) print_usage
		exit 1 ;;
    * ) print_usage
		exit 1 ;;
  esac
done
shift $(($OPTIND - 1))

file=$1
filedir=$(dirname $file)
outdir=${filedir}/subsampling
name="$(basename ${file} .bam)"

check if folder exists
if [[ -d $outdir ]]
then
	if [[ $overwrite == 1 ]]
	then
		rm -r $outdir
		mkdir $outdir
	else
		echo "Output folder subsampling exists. Exit script."
		exit 1
	fi
else
	mkdir $outdir
fi

if [[ $manual == 1 ]]
then
	read -p "Set seed: " seed
	read -p "Please give list of percentages for subsampling as integers: " list
else
	seed="0"
	list=(10 25 50 75)
fi

echo "Subsample for following values:" ${list[@]}
for perc in ${list[@]}
do
	echo "Subsample with percentage of ${perc}"
	samtools view -s ${seed}.${perc} -b $1 > $outdir/${name}_${perc}p.bam
done

echo "Subsampling finished. Files saved in subsampling/"

if [[ $dge == 1 ]]
then
	files=${outdir}/${name}_*.bam
	tmpdir=${outdir}/tmp

	# create tmp directory if it does not exist yet
	if [[ ! -d $tmpdir ]]
	then mkdir -p ${tmpdir}
	fi

	for file in $files
	do
		echo "Processing of ${file} started..."
		filename="$(basename ${file} .bam)"
		perc=`echo ${filename} | awk -F "_" '{print $NF}'`
		
		# counting exonic reads only
		dge="${dropseq_root}/DigitalExpression I=${file} O=${outdir}/DGE_matrix_${perc}.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_dir=${tmpdir}"
		$dge

		# counting both intronic and exonic reads
		dge_with_introns="${dropseq_root}/DigitalExpression I=${file} O=${outdir}/DGE_matrix_with_introns_${perc}.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 LOCUS_FUNCTION_LIST=INTRONIC TMP_dir=${tmpdir}"
		$dge_with_introns

		# collect RNAseq metrics with PICARD
		#rnaseq_metrics="java -jar ${picard_jar} CollectRnaSeqMetrics I=${outdir}/${file} O=${outdir}/rnaseq_metrics_${perc}.RNA_Metrics REF_FLAT=${refflat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=${rRNA_intervals}"
		#$echo_prefix $rnaseq_metrics

		echo "Processing of ${file} finished."
	done

	echo "Generation of DGE matrix of all subsampled files finished."
fi