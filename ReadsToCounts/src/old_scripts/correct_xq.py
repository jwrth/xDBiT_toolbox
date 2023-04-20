#!/usr/bin/env python

"""
This tool is to modify the 'XQ' quality score in the tagged bam files.

It subtracts the number of padded Ns at the end of the 'XD' tag from the 'XQ' tag.
"""

# Library
import pysam
from argparse import ArgumentParser
import subprocess
from datetime import datetime, timedelta

# functions
def run_in_commandline(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout

# Setup input parser
parser = ArgumentParser()
parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Input bam file. Default to stdin.")
parser.add_argument("-o" "--output_bam", action="store", dest="output_bam", default="out.bam", help="Output bam file. Default out.bam.")

# Parse input
args = parser.parse_args()

input_bam = args.input_bam
output_bam = args.output_bam


infile =pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
outfile = pysam.AlignmentFile(output_bam, 'wb', template=infile)

bam_length = run_in_commandline('samtools view ' + input_bam + ' | wc -l')
total_reads = int(bam_length.decode())

# start timing
start_time_filtering = datetime.now()
start_time = datetime.now()
stride = 100000
filename = input_bam.split('/')[-1]

print("Correction of XQ tags started...")
for idx, entry in enumerate(infile.fetch(until_eof=True)):
    if entry.has_tag('XQ'):
        xq = entry.get_tag('XQ')
        xd = entry.get_tag('XD')
        
        # calculate new XQ
        xd_N = len(xd) - len(xd.rstrip('N'))
        xq -= xd_N

        if xq == 0:
            entry.set_tag('XQ', None)

        else:
            entry.set_tag('XQ', xq)

    if (idx+1) % stride == 0:
        totaltime = datetime.now() - start_time_filtering
        stride_steptime = datetime.now() - start_time
        time_to_go = (total_reads - (idx+1))/stride * stride_steptime
        print("File " + filename + " - Reads " + str(idx+1) + "/" + str(total_reads) + " processed. Time for last " + str(stride) + ": " + str(stride_steptime) + ", Total time: " + str(totaltime) + ". Time remaining: " + str(time_to_go))
        start_time = datetime.now()

    outfile.write(entry)
    
infile.close()
outfile.close()
print("Finished.")
