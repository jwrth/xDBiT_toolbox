#!/usr/bin/env python

"""
This is a tool to pad sequences from fasta files with a certain ladder to create reads of a homogeneous length.

It is used for the Split-seq pipeline to allow the analysis of reads that are shorter than the 
expected length of 94 nts.
"""
import itertools, string, sys
from argparse import ArgumentParser
from datetime import datetime, timedelta

# Setup input parser
parser = ArgumentParser()
parser.add_argument("-l" "--final_length", action="store", dest="final_length", default="-", help="Specify the final length for the padding.")
parser.add_argument("-q" "--quality_pad", action="store", dest="quality_pad", default="E", help="Specify the quality of the padded Ns.")
parser.add_argument("-d" "--letter_pad", action="store", dest="letter_pad", default="N", help="Letter to be used for the padding.")
parser.add_argument("-o" "--output_dir", action="store", dest="output_dir", default=".", help="Output directory. Default current directory.")
parser.add_argument("-s" "--stride", action="store", dest="stride", default=1000000, help="Output directory.")

# Parse input
args = parser.parse_args()

final_length = args.final_length
quality_pad = args.quality_pad
letter_pad = args.letter_pad
output_dir = args.output_dir
stride = args.stride

# strip newlines from the standard input
stream = map(lambda a: a.strip(), sys.stdin)

# create output file
output = open(output_dir, 'w')

print("Processing...")
# this works only if sequences span only one line
n = 0
for head in stream:
    
    # advance the stream
    seq  = next(stream)
    tmp  = next(stream)
    qual = next(stream)
    n+=1

    # this is how much to pad
    size = int(final_length) - len(seq)

    # pad each line
    seq  = seq  + str(letter_pad) * size
    qual = qual + str(quality_pad) * size

    # print the joined elements
    # print("\n".join( (head, seq, tmp, qual) ))

    # write to output file
    output.write("\n".join((head, seq, tmp, qual)) + "\n")

    if n % int(stride) == 0:
        print(str(n) + " reads process. Time: " + str(datetime.now()))

print("Finished")
