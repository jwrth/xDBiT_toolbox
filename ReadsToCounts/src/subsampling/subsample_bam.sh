#!/bin/bash

### Bash script for subsampling of bam files.
# First set seed. Second give list of percentages separated by spaces.

manual=0
overwrite=0

function print_usage () {
    cat >&2 <<EOF
Usage: Bash script for subsampling of bam files.
-m <manual>      : Set seed and subsample list manually.
-o <overwrite>	 : Overwrite existing output folders and files.
EOF
}



while getopts "mo" flag; do
  case $flag in
    m ) manual=1;;
	o ) overwrite=1;;
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

# check if folder exists
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