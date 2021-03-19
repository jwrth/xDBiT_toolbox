#!/bin/bash

## Reheading and merging
# This script adds the header back to the split files and then merges the files. This results 
# in a very long header in the split files because all the changes are documented in the header
# and summed up. However first merging and then reheading gives and error and the RG tag is lost.

# functions
timestamp() {
  date +"%Y-%m-%d %T"
}

split_dir=$1
out_file=$2

# add a header to each of the split files
echo "`timestamp` - Start adding headers to the split files..."
filtered_files=${split_dir}/out_x*
for file in $filtered_files
do
	filename="$(basename -- ${file})"
	samtools reheader ${split_dir}/header.sam ${file} > ${split_dir}/reheaded_${filename}
done
echo "`timestamp` - Headers added to files."

echo "`timestamp` - Merging of BAM files initiated..."
# merge the reheaded files
if samtools merge ${out_file} ${split_dir}/reheaded_* ; then
	echo "`timestamp` - Merging successful - ready to delete the unnecessary files"
	reheaded_files=${split_dir}/reheaded_*
	split_files=${split_dir}/x*
	header_file=${split_dir}/header.sam

	files_to_be_deleted="${filtered_files} ${reheaded_files} ${split_files} ${header_file}"

	if rm $files_to_be_deleted ; then
		echo "`timestamp` - Temporary files deleted."
	else
		echo "`timestamp` - Error during removal of temporary files."
	fi

else
	echo "`timestamp` - Merging failed"
fi