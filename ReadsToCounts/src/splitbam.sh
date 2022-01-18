#!/bin/bash

# This tool splits input bam files into a specified number of smaller bam files and extracts the header of the original bam file

# functions
timestamp() {
  date +"%Y-%m-%d %T"
}

show_time () {
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

# read input variables
split_dir=$1
nfiles=$2
file=$3

mkdir ${split_dir}

# count number of reads of bam file
size=$(samtools view ${file} | wc -l)

# calculate size of file chunks
chunksize=$(((${size}+(${nfiles}-1))/${nfiles}))

start_time=`date +%s`
echo "Splitting of bam files initiated..."
echo "`timestamp` - Creating header file"
# create header file
samtools view -H ${file} > ${split_dir}/header.sam
echo "`timestamp` - Header file created"

echo "`timestamp` - Splitting of files initiated..."
# split bam files into temp directory as sam files
samtools view ${file} | (cd ${split_dir}; split -l $chunksize)
echo "`timestamp` - Splitting finished."

echo "`timestamp` - Conversion to BAM files initiated..."
# make bam file out of sam files
for f in ${split_dir}/x*
do
 samtools view -b ${f} > ${f}.bam
done
echo "`timestamp` - Conversion to BAM files finished."

end_time=`date +%s`
run_time=`expr $end_time - $start_time`
total_time=`show_time $run_time`
echo "`timestamp` - ${file} split into ${nfiles} and saved to ${split_dir} in ${total_time}"