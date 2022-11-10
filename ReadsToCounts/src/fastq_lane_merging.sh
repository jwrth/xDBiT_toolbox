#!/bin/bash

'''
Script to merge .fastq.gz files across lanes in multiple directories.

Use case: If the --no-lane-splitting option was not used in the bcl2fastq tool this script can be used to create fastq.gz files that are not split by lane.

Usage:
    bash fastq_lane_merging.sh dir1 dir2 dir3 dir4 ...
    or
    bash fastq_lane_merging.sh pattern_in_dir*


Source:
https://www.biostars.org/p/317385/
'''

dirs=$@

for dir in $dirs
do
    echo "Process $dir"
    for i in $(find $dir -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
    do 
        echo "Merging R1 in $dir"
        cat $dir/"$i"_L00*_R1_001.fastq.gz > $dir/"$i"_merged_R1_001.fastq.gz

        echo "Merging R2 in $dir"
        cat $dir/"$i"_L00*_R2_001.fastq.gz > $dir/"$i"_merged_R2_001.fastq.gz

    done
done
