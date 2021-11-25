#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is to process Split-seq reads after tagging the bam file with cellular and molecular barcodes.

Input: BAM file with reads that have been filtered, polyA and SMART adapter trimmed, and tagged with the following:
    - XD: cell barcode 1
    - XE: cell barcode 2
    - XF: cell barcode 3
    - XM: molecular barcode (UMI)

Output:
    - Filtered bam file containing only reads with expected barcodes. Barcodes that are within 1 hamming distance of an expected
    barcode sequence are corrected. An additional tag XC is added, which corresponds to the concatenated barcode sequences 
    strating with the RT barcode, then round 1 and round 2 ligation barcodes.

Copyright: Rebekka Wegmann, Snijderlab, ETH Zurich, 2019

# MIT License
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
matplotlib.use('pdf') #prevents matplotlib from trying to open a figure window, which will fail on systems that do not have graphics
from matplotlib import pyplot as plt
from argparse import ArgumentParser
from datetime import datetime, timedelta
import subprocess

#%% Functions
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def well2ind(well):
    """Convert well positions to (zero-based) indices"""
    d = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7}
    row = d[well[0]]
    col = int(well[1:])-1
    return [row,col]

def make_plate_overview_mod(well_count):
    # bring the counts per well into a plate layout
    out = np.zeros([8, 12])
    for wellpos in well_count.items():
        xy = well2ind(wellpos[0])
        out[xy[0], xy[1]] = int(wellpos[1])
    return out

def get_bam_length(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout

def average_time(time_list):
    average_timedelta = sum(time_list, timedelta(0)) / len(time_list)
    return average_timedelta


#%%Setup input parser
parser = ArgumentParser()
parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Specify the input bam file. Defaults to stdin.")
parser.add_argument("-o" "--output_bam", action="store", dest="output_bam", default="-", help="Specify the output bam file. Defaults to stdout.")
parser.add_argument("-n" "--est_num_cells", action="store", dest="est_num_cells", default=2500, help="Estimated number of cells. Defaults to 2500.",type=int)
parser.add_argument("-d" "--out_dir", action="store", dest="out_dir", default=".", help="Directory to store logfiles and output plots. Defaults to the current directory.")
parser.add_argument("-b" "--bc_dir", action="store", dest="bc_dir", default=".", help="Directory where the expected barcode files are stored. Defaults to the directory this script is in.")
parser.add_argument("--debug_flag",action="store_true",help="Turn on debug flag. This will produce some additional output which might be helpful.")
parser.add_argument("--store_discarded",action="store_true",help="Store names of discarded reads?")

#%% Parse input
args = parser.parse_args()

debug_flag = args.debug_flag
store_discarded = args.store_discarded
input_bam = args.input_bam #use "-" for stdin, set flag to rb
output_bam = args.output_bam
est_num_cells = args.est_num_cells
out_dir = args.out_dir
bc_dir = args.bc_dir

#%% Determine number of reads in .bam file
bam_length = get_bam_length('samtools view ' + input_bam + ' | wc -l')
total_reads = int(bam_length.decode())
print("Number of reads to be filtered: " + str(total_reads))


if bc_dir==".":
    bc_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
    
#%% Write parameters to logfile
print('Splitseq barcode filtering log - based on hamming distance\n---------------------------------------\nParameters:', file = open(os.path.join(out_dir,'barcode_filtering_log.txt'), 'w'))
print('Input bam: %s' % input_bam, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Output bam: %s' % output_bam, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Estimated number of cells: %d' % est_num_cells, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Output directory: %s' % out_dir, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Barcode directory: %s' % bc_dir, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))

if store_discarded:
    if os.path.isfile(os.path.join(out_dir,'discarded_reads.txt')):
        os.remove(os.path.join(out_dir,'discarded_reads.txt'))
        print('Old version of discarded reads.txt deleted.')

#%% Start timing
#t_start = timeit.default_timer()
t_start = datetime.now()

#%% Get expected barcodes. 
bc1 = pd.read_csv(os.path.join(bc_dir,"32_19_barcodes_1.csv"))
bc2 = pd.read_csv(os.path.join(bc_dir,"32_19_barcodes_2.csv"))
bc3 = pd.read_csv(os.path.join(bc_dir,"32_19_barcodes_3.csv"))

#%% Generate dictionaries for barcode-well assignment
bc1_dict = {bc1.Barcode.values[x]:bc1.WellPosition.values[x][0]+bc1.WellPosition.values[x][1:].zfill(2) for x in range(len(bc1.Barcode.values))}
bc2_dict = {bc2.Barcode.values[x]:bc2.WellPosition.values[x][0]+bc2.WellPosition.values[x][1:].zfill(2) for x in range(len(bc2.Barcode.values))}
bc3_dict = {bc3.Barcode.values[x]:bc3.WellPosition.values[x][0]+bc3.WellPosition.values[x][1:].zfill(2) for x in range(len(bc3.Barcode.values))}

#%% Read in BAM file
infile =pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
outfile = pysam.AlignmentFile(output_bam, 'wb', template=infile)

#%% for debugging only
if debug_flag:
    for entry in itertools.islice(infile, 10 ):
        print(entry.query_name)
        print(entry.get_forward_sequence())
        print(entry.get_tag('XD'))
        print(entry.get_tag('XE'))
        print(entry.get_tag('XF'))
        print(entry.get_tag('XM'))
        
#%% Check and correct barcodes
n=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0

n_well1 = {x[0]+x[1:].zfill(2):0 for x in bc1.WellPosition.unique()}
n_well2 = {x[0]+x[1:].zfill(2):0 for x in bc2.WellPosition.unique()}
n_well3 = {x[0]+x[1:].zfill(2):0 for x in bc3.WellPosition.unique()}

all_bcs =  [None]*est_num_cells*100

# modified filtering algorithm
print('Filtering and correcting barcodes...')

start_time_filtering = datetime.now()
avg_steptime = []
for entry in infile.fetch(until_eof=True):
    start_time = datetime.now()
    n+=1
    keep=True
    xd = entry.get_tag('XD')
    xe = entry.get_tag('XE')
    xf = entry.get_tag('XF')
    
    if xd in bc1.Barcode.values and xe in bc2.Barcode.values and xf in bc3.Barcode.values:
        n1+=1
        n2+=1
        n3+=1
        n4+=1
        #n_bc1[xd]+=1
        #n_bc2[xe]+=1
        #n_bc3[xf]+=1
        
        # substitute barcode tag with well tag
        xd_well = bc1_dict[xd]
        xe_well = bc2_dict[xe]
        xf_well = bc3_dict[xf]
        
        entry.set_tag('XD', xd_well)
        entry.set_tag('XE', xe_well)
        entry.set_tag('XF', xf_well)
        
        n_well1[xd_well]+=1
        n_well2[xe_well]+=1
        n_well3[xf_well]+=1
        
    else:
        if xd in bc1.Barcode.values:
            n1+=1
            #n_bc1[xd]+=1
            xd_well = bc1_dict[xd]
            entry.set_tag('XD', xd_well)
            n_well1[xd_well]+=1
        else:
            d = [hamming(xd,bc) for bc in bc1.Barcode.values]
            idx = [i for i,e in enumerate(d) if e==1]
            if len(idx)==1:
                n5+=1
                xd = bc1.Barcode.values[idx[0]]
                xd_well = bc1_dict[xd]
                entry.set_tag('XD', xd_well)
                #entry.set_tag('XD',xd)
                #n_bc1[xd]+=1
                n_well1[xd_well]+=1
            else:
                keep=False
           
        if xe in bc2.Barcode.values:
            n2+=1
            #n_bc2[xe]+=1
            xe_well = bc2_dict[xe]
            entry.set_tag('XE', xe_well)
            n_well2[xe_well]+=1
        else:
            d = [hamming(xe,bc) for bc in bc2.Barcode.values]
            idx = [i for i,e in enumerate(d) if e==1]
            if len(idx)==1:
                n6+=1
                xe = bc2.Barcode.values[idx[0]]
                xe_well = bc2_dict[xe]
                entry.set_tag('XE', xe_well)
                #entry.set_tag('XE',xe)
                #n_bc2[xe]+=1
                n_well2[xe_well]+=1
            else:
                keep=False
           
        if xf in bc3.Barcode.values:
            n3+=1
            #n_bc3[xf]+=1
            xf_well = bc3_dict[xf]
            entry.set_tag('XF', xf_well)
            n_well3[xf_well]+=1
        else:
            d = [hamming(xf,bc) for bc in bc3.Barcode.values]
            idx = [i for i,e in enumerate(d) if e==1]
            if len(idx)==1:
                n7+=1
                xf = bc3.Barcode.values[idx[0]]
                xf_well = bc3_dict[xf]
                entry.set_tag('XF', xf_well)
                #entry.set_tag('XF', xf)
                #n_bc3[xf]+=1
                n_well3[xf_well]+=1
            else:
                keep=False
    if keep:
        xc_well = xd_well + xe_well + xf_well
        #xc = xd+xe+xf
        entry.set_tag('XC', xc_well)
        outfile.write(entry)
        if n8<est_num_cells*100:
            all_bcs[n8]=xc_well
        n8+=1
    else:
        if store_discarded:
            print(entry.query_name, file = open(os.path.join(out_dir,'discarded_reads.txt'),'a'))

    steptime = datetime.now() - start_time
    avg_steptime.append(steptime)
            
    if n % 500000 == 0:
        totaltime = datetime.now() - start_time_filtering
        avg_steptime = average_time(avg_steptime)
        time_to_go = (total_reads - n) * avg_steptime
        print("Reads " + str(n) + "/" + str(total_reads) + " processed. Step time: " + str(avg_steptime) + ", Total time: " + str(totaltime) + ". Time remaining: " + str(time_to_go))
        avg_steptime = []


all_bcs = all_bcs[:n8]
n_all_bcs = len(all_bcs)

filter_stop = datetime.now()
filter_elapsed = filter_stop - t_start

print('Read %d entries' % n, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] complete barcodes' % (n4, n4/float(n)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected BC1, corrected %d [%.2f%%] BC1' % (n1 , n1/float(n)*100, n5, n5/float(n)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected BC2, corrected %d [%.2f%%] BC2' % (n2 , n2/float(n)*100, n6, n6/float(n)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected BC3, corrected %d [%.2f%%] BC3' % (n3 , n3/float(n)*100, n7, n7/float(n)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Retained %d [%.2f%%] reads after correction and filtering' % (n8, n8/float(n)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Elapsed time for filtering: ' + str(filter_elapsed), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))


infile.close()
outfile.close()

#t_stop = timeit.default_timer()
t_stop = datetime.now()
t_elapsed = t_stop-t_start
print("Barcode filtering and correcting finished. Elapsed time: " + str(t_elapsed))
t_start=t_stop

#%% Create summary plots
print("Creating summary plots. This step might take a while...")
bc1_matrix = make_plate_overview_mod(n_well1)
bc2_matrix = make_plate_overview_mod(n_well2)
bc3_matrix = make_plate_overview_mod(n_well3)

# this part is super slow if you have a large number of cells, maybe omit
all_bc_counts = {i:all_bcs.count(i) for i in list(set(all_bcs))}
all_bc_cumsum = np.cumsum(sorted(list(all_bc_counts.values()), reverse=True))/n_all_bcs

# plotting
fig, axes = plt.subplots(2,2)
p1 = axes[0,0].imshow(np.log10(bc1_matrix+1));  axes[0,0].set_title('Number of reads per RT barcode')
clb = fig.colorbar(p1, ax=axes[0,0]); clb.set_label('No. BCs, log10')
p2 = axes[0,1].imshow(np.log10(bc2_matrix+1));  axes[0,1].set_title('Number of reads per round 2 barcode')
clb = fig.colorbar(p2, ax=axes[0,1]); clb.set_label('No. BCs, log10')
p3 = axes[1,0].imshow(np.log10(bc3_matrix+1));  axes[1,0].set_title('Number of reads per round 3 barcode')
clb = fig.colorbar(p3, ax=axes[1,0]); clb.set_label('No. BCs, log10')
p4 = axes[1,1].plot(all_bc_cumsum); axes[1,1].set_title('Cumulative fraction of reads per barcode'); axes[1,1].set_xlim(0, 5000)

fig.set_size_inches(12,7)
fig.savefig(os.path.join(out_dir,'barcode_filtering_summary.pdf'),bbox_inches='tight')

#%% Stop timing
t_stop = datetime.now()
t_elapsed = t_stop-t_start

#%% Save the QC output
f = h5py.File(os.path.join(out_dir,'splitseq_filtering_QC_data.hdf5'), 'w')

f.create_dataset('BC1_plate_overview', data = bc1_matrix)
f.create_dataset('BC2_plate_overview', data = bc2_matrix)
f.create_dataset('BC3_plate_overview', data = bc3_matrix)
f.create_dataset('reads_per_BC', data=list(all_bc_counts.values()))
f.create_dataset('labels_reads_per_BC', data=np.string_(list(all_bc_counts.keys())))

f.close()

#%% print info to stdout
print("Generating summary finished. Elapsed time: " + str(t_elapsed))
