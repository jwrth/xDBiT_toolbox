#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is to check discarded reads after  barcode filtering

Copyright: Rebekka Wegmann, Snijderlab, ETH Zurich, 2019

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
"""

#%% libraries
import sys, os, timeit, h5py
import pysam
import itertools
import numpy as np
import pandas as pd
import matplotlib

#%% Functions
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

#%% Read in SAM file
#print(sys.argv[1])
infile = pysam.AlignmentFile(sys.argv[1], 'r', check_sq=False, check_header=False)
        
#%% Check read 2 for correct spacer sequences

spacer1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
spacer2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"



n=0
n1=0
n2=0
n3=0
n4=0
n5 = 0

for entry in infile.fetch(until_eof=True):
    if entry.flag == 77:
        continue
    n+=1
    if n <= 10:
        print(entry.query_name)
        print(entry.flag)
        print(entry.get_forward_sequence())

    spacer_seq_1 = entry.get_forward_sequence()[18:48]
    spacer_seq_2 = entry.get_forward_sequence()[56:86]

    if(hamming(spacer1, spacer_seq_1)<=3):
        n1+=1
    if(hamming(spacer2, spacer_seq_2)<=3):
        n2+=1

    # check if there are "off by one" errors, due to the read being dephased / shifted by one base
    if(hamming(spacer1, entry.get_forward_sequence()[17:47]) <= 3 or hamming(spacer1, entry.get_forward_sequence()[19:49]) <= 3):
        n3 +=1
    if(hamming(spacer2, entry.get_forward_sequence()[55:85]) <= 3 or hamming(spacer2, entry.get_forward_sequence()[57:87]) <= 3):
        n4 +=1  

    # check if there are cases where only the first part of spacer 2 is present (meaning that probably a BC3-BC2 dimer acted as a random primer?)
    if(hamming(spacer2[0:15], entry.get_forward_sequence()[56:71]) <= 3 and not hamming(spacer2[15:30], entry.get_forward_sequence()[71:86]) <= 3):
        n5 +=1  

print('Read %d read pairs' % n)
print('Found %d [%.2f%%] correct spacer 1 sequences' % (n1, n1/float(n)*100))
print('Found %d [%.2f%%] correct spacer 2 sequences' % (n2, n2/float(n)*100))
print('Found %d [%.2f%%] spacer 1 that are shifted +- 1 basepair' % (n3, n3/float(n)*100))
print('Found %d [%.2f%%] spacer 2 that are shifted +- 1 basepair' % (n4, n4/float(n)*100))
print('Found %d [%.2f%%] first half of spacer 2' % (n5, n5/float(n)*100))

infile.close()
