#!/bin/bash

### Bash script for subsampling of bam files.
# First set seed. Second give list of percentages separated by spaces.

manual=false
overwrite=false

function print_usage () {
    cat >&2 <<EOF
Usage: Bash script for subsampling of bam files.
-m <manual>      : Set seed and subsample list manually.
EOF
}

while getopts 'mo' flag; do
  case "${flag}" in
    m) manual=true ;;
	o) overwrite=true;;
	h) print_usage
		exit 1 ;;
    *) print_usage
		exit 1 ;;
  esac
done
shift $(($OPTIND - 1))

# check if folder exists
if [ -d "subsampling" ]
then
	if [ "$overwrite" = true ]
	then
		rm -r subsampling
		mkdir subsampling
	else
		echo "Output folder subsampling exists."
		exit 1
	fi
else
	mkdir subsampling
fi

file=$1
name="$(basename $file .bam)"

if [ "$manual" = true ]
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
	samtools view -s ${seed}.${perc} -b $1 > subsampling/${name}_${perc}p.bam
done

echo "Subsampling finished. Files saved in subsampling/"