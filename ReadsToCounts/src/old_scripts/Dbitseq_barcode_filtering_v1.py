#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is to process Split-seq reads after tagging the bam file with cellular and molecular barcodes.

Input: BAM file with reads that have been filtered, polyA and SMART adapter trimmed, and tagged with the following:
    - XD: cell barcode 1
    - XE: cell barcode 2
    - XF: cell barcode 3
    - XM: molecular barcode (UMI)

Algorithm:
    - Choice between levenshtein distance and hamming distance for correction

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
from matplotlib import pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm_notebook as tqdm
from datetime import datetime, timedelta
import subprocess
from fuzzywuzzy import fuzz, process
import glob
import collections
from multiprocessing import Pool
import Levenshtein as lv


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

def average_time(time_list):
    average_timedelta = sum(time_list, timedelta(0)) / len(time_list)
    return average_timedelta

def get_stdout(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout

def sum_dicts(dictlist):
    counter = collections.Counter() 
    for d in dictlist:  
        counter.update(d) 
      
    result = dict(counter)
    return result

def spatialfilter(in_bam):
    # count read length of input bam
    print("Count number of reads in input bam...", flush=True)
    bam_length = get_stdout('samtools view ' + in_bam + ' | wc -l')
    total_reads = int(bam_length.decode())

    # read input bam
    infile = pysam.AlignmentFile(in_bam, 'rb', check_sq=False)
    filename = in_bam.split('/')[-1]

    # generate output bam name and directory
    if multi:
        out_bam = '/'.join(in_bam.split('/')[:-1]) + '/out_' + in_bam.split('/')[-1]
    else:
        out_bam = os.path.join(tmp_dir, "unaligned_tagged_BC_filtered.bam")

    outfile = pysam.AlignmentFile(out_bam, 'wb', template=infile)

    # set up recording array and variables
    n = np.zeros(7, dtype=np.int32)

    n_wellx = {x[0]+x[1:].zfill(2):0 for x in bcx.WellPosition.unique()}
    n_welly = {x[0]+x[1:].zfill(2):0 for x in bcy.WellPosition.unique()}

    # set up list to collect the barcodes for later statistics
    all_bcs =  [None]*est_num_cells*100

    # start timing
    start_time_filtering = datetime.now()
    start_time = datetime.now()

    # start filtering
    print("Filtering started...", flush=True)

    for entry in infile.fetch(until_eof=True):
        n[0]+=1
        keep = True

        # extract barcodes from tags
        x = entry.get_tag('XX')
        y = entry.get_tag('XY')

        if x in bcx.Barcode.values and y in bcy.Barcode.values:
            # both barcodes directly found
            # record finding
            n[1]+=1
            n[2]+=1
            n[3]+=1

            # get well tag
            x_well = bcx_dict[x]
            y_well = bcy_dict[y]

            # get coordinate
            x_coord = x_assign_dict[x_well]
            y_coord = y_assign_dict[y_well]

            # save coordinate as tag
            entry.set_tag('XX', str(x_coord))
            entry.set_tag('XY', str(y_coord))

            n_wellx[x_well]+=1
            n_welly[y_well]+=1

        else:
            # checking for barcode 1
            if x in bcx.Barcode.values:
                # not all barcodes were found directly but bc1
                #print("Barcode 1 found directly")
                # determine well coordinate and set new tag with well coordinate
                x_well = bcx_dict[x]
                x_coord = x_assign_dict[x_well]
                entry.set_tag('XX', str(x_coord))

                # record
                n[1]+=1
                n_wellx[x_well]+=1

            else:
                # barcode 1 was not found directly. Fuzzy string matching is applied
                #print("Fuzzy string matching applied")
                d = [compute_dist(x,bc) for bc in bcx.Barcode.values]
                idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                if len(idx)==1:
                    x = bcx.Barcode.values[idx[0]]
                    x_well = bcx_dict[x]
                    x_coord = x_assign_dict[x_well]
                    entry.set_tag('XX', str(x_coord))

                    # record
                    n[4]+=1
                    n_wellx[x_well]+=1
                else:
                    keep=False

            # checking for barcode 2
            if y in bcy.Barcode.values:
                # not all barcodes were found directly but bc2
                #print("Barcode 2 found directly")
                # determine well coordinate and set new tag with well coordinate
                y_well = bcy_dict[y]
                y_coord = y_assign_dict[y_well]
                entry.set_tag('XY', str(y_coord))

                # record
                n[2]+=1
                n_welly[y_well]+=1

            else:
                # barcode 2 was not found directly. Fuzzy string matching is applied
                #print("Fuzzy string matching applied")

                d = [compute_dist(y,bc) for bc in bcy.Barcode.values]
                idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                if len(idx)==1:
                    y = bcy.Barcode.values[idx[0]]
                    y_well = bcy_dict[y]
                    y_coord = y_assign_dict[y_well]
                    entry.set_tag('XY', str(y_coord))

                    # record
                    n[5]+=1
                    n_welly[y_well]+=1
                else:
                    keep = False


        if keep:
            # concatenate coordinates
            xy_coord = str(x_coord) + "," + str(y_coord)

            # set coordinates as tag and write to output file
            entry.set_tag('XC', str(xy_coord))
            outfile.write(entry)
            if n[6]<est_num_cells*100/ncores:
                all_bcs[n[6]]=xy_coord
            n[6]+=1
        else:
            if store_discarded:
                print(entry.query_name, file = open(os.path.join(out_dir,'discarded_reads.txt'),'a'), flush=True)

        if n[0] % stride == 0:
            totaltime = datetime.now() - start_time_filtering
            stride_steptime = datetime.now() - start_time
            time_to_go = (total_reads - n[0])/stride * stride_steptime
            print("File " + filename + " - Reads " + str(n[0]) + "/" + str(total_reads) + " processed. Time for last " + str(stride) + ": " + str(stride_steptime) + ", Total time: " + str(totaltime) + ". Time remaining: " + str(time_to_go),
                flush=True)
            start_time = datetime.now()
            
    all_bcs = all_bcs[:n[6]]
    n_all_bcs = len(all_bcs)
    print("Filtering finished.", flush=True)
    
    infile.close()
    outfile.close()
    return [n, [n_wellx, n_welly], all_bcs, n_all_bcs]
    
#%%Setup input parser
parser = ArgumentParser()
parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Specify the input bam file. Defaults to stdin.")
parser.add_argument("-n" "--est_num_cells", action="store", dest="est_num_cells", default=2500, help="Estimated number of cells. Defaults to 2500.",type=int)
parser.add_argument("-d" "--out_dir", action="store", dest="out_dir", default=".", help="Directory to store logfiles and output plots. Defaults to the current directory.")
parser.add_argument("-t" "--tmp_dir", action="store", dest="tmp_dir", default=".", help="Temp directory")
parser.add_argument("-b" "--bc_dir", action="store", dest="bc_dir", default=".", help="Directory where the expected barcode files are stored. Defaults to the directory this script is in.")
parser.add_argument("--debug_flag",action="store_true",help="Turn on debug flag. This will produce some additional output which might be helpful.")
parser.add_argument("--store_discarded",action="store_true",help="Store names of discarded reads?")
parser.add_argument("-a" "--dist_alg", action="store", dest="dist_alg", default="hamming", help="Distance algorithm to be used: levenshtein or hamming")
parser.add_argument("-z" "--dist_threshold", action="store", dest="dist_threshold", default=1, help="Threshold to be used for levenshtein or hamming distance matching")
parser.add_argument('-m', action='store_true', help="Use multithreading?")


#%% Parse input
args = parser.parse_args()

debug_flag = args.debug_flag
store_discarded = args.store_discarded
input_bam = args.input_bam #use "-" for stdin, set flag to rb
est_num_cells = args.est_num_cells
out_dir = args.out_dir
tmp_dir = args.tmp_dir
bc_dir = args.bc_dir
dist_alg = args.dist_alg
dist_threshold = int(args.dist_threshold)
multi = args.m

split_dir = os.path.join(tmp_dir, 'tmp_split')

# define frequency of printed outputs during barcode filtering
stride = 500000

# retrieve information about algorithm
if dist_alg == "hamming":
    compute_dist = lv.hamming
    
if dist_alg == "levenshtein":
    compute_dist = lv.distance

if bc_dir==".":
    bc_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
    
#%% Write parameters to logfile
print('Splitseq barcode filtering log - based on %s algorithm\n---------------------------------------\nParameters:' % dist_alg, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'), 'w'))
print('Input bam: %s' % input_bam, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Path to output bam files: %s' % split_dir, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
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

#%% Import barcode and well assignment files
bc_dir = "../barcodes/"
bcx = pd.read_csv(glob.glob(os.path.join(bc_dir,"*barcodes_Lig*.csv"))[0])
bcy = pd.read_csv(glob.glob(os.path.join(bc_dir,"*barcodes_RT*.csv"))[0])
well_coord_ass = pd.read_csv(glob.glob(os.path.join(bc_dir,"*well_coord_assignment*.csv"))[0])

#%% Get expected barcodes and bridges
#bc_files = glob.glob(os.path.join(bc_dir,"*barcodes*.csv"))
#bc1 = pd.read_csv(bc_files[0])
#bc2 = pd.read_csv(bc_files[1])
#bc3 = pd.read_csv(bc_files[2])
#bridges = pd.read_csv(os.path.join(bc_dir, "32_19_bridges.csv"))

#%% Generate dictionaries for barcode-well-coordinate assignment
bcx_dict = {bcx.Barcode.values[x]:bcx.WellPosition.values[x][0]+bcx.WellPosition.values[x][1:].zfill(2) for x in range(len(bcx.Barcode.values))}
bcy_dict = {bcy.Barcode.values[x]:bcy.WellPosition.values[x][0]+bcy.WellPosition.values[x][1:].zfill(2) for x in range(len(bcy.Barcode.values))}

x_assign_dict = {well_coord_ass.X_Well[x][0]+well_coord_ass.X_Well[x][1:].zfill(2):well_coord_ass.X_Coord[x] for x in range(len(well_coord_ass)) if not pd.isnull(well_coord_ass.X_Well[x])}
y_assign_dict = {well_coord_ass.Y_Well[x][0]+well_coord_ass.Y_Well[x][1:].zfill(2):well_coord_ass.Y_Coord[x] for x in range(len(well_coord_ass)) if not pd.isnull(well_coord_ass.Y_Well[x])}

#%% Generate dictionaries for barcode-well assignment
#bc1_dict = {bc1.Barcode.values[x]:bc1.WellPosition.values[x][0]+bc1.WellPosition.values[x][1:].zfill(2) for x in range(len(bc1.Barcode.values))}
#bc2_dict = {bc2.Barcode.values[x]:bc2.WellPosition.values[x][0]+bc2.WellPosition.values[x][1:].zfill(2) for x in range(len(bc2.Barcode.values))}
#bc3_dict = {bc3.Barcode.values[x]:bc3.WellPosition.values[x][0]+bc3.WellPosition.values[x][1:].zfill(2) for x in range(len(bc3.Barcode.values))}

#%% for debugging only
if debug_flag:
    for entry in itertools.islice(infile, 10):
        print(entry.query_name)
        print(entry.get_forward_sequence())
        print(entry.get_tag('XD'))
        print(entry.get_tag('XE'))
        print(entry.get_tag('XF'))
        print(entry.get_tag('XM'))


if multi:
    #%% start multithreaded filtering   
    if __name__ == '__main__':
        files = glob.glob(os.path.join(split_dir, 'x*.bam'))
        ncores = len(files)
        multifilter = Pool(processes=ncores)
        results = multifilter.map(spatialfilter, files)

    # extract and sum up recording variables
    n = np.array([elem[0] for elem in results]).sum(axis=0)

    n_wellx = sum_dicts([elem[1][0] for elem in results])
    n_welly = sum_dicts([elem[1][1] for elem in results])

    all_bcs = [item for sublist in [elem[2] for elem in results] for item in sublist]
    n_all_bcs = np.array([elem[3] for elem in results]).sum()

else:
    # without multithreading
    ncores = 1
    results = np.array(spatialfilter(input_bam))

    n = results[0]

    n_wellx = results[1][0]
    n_welly = results[1][1]

    all_bcs = results[2]
    n_all_bcs = results[3]


filter_stop = datetime.now()
filter_elapsed = filter_stop - t_start

# print into log file
print('Read %d entries' % n[0], file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] complete barcodes' % (n[3], n[3]/float(n[0])*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected ligation barcode (x-coordinate), with %s matching %d [%.2f%%] (distance: %d)' % (n[1] , n[1]/float(n[0])*100, dist_alg, n[4], n[4]/float(n[0])*100, dist_threshold), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected RT barcode (y-coordinate), with %s matching %d [%.2f%%] (distance: %d)' % (n[2] , n[2]/float(n[0])*100, dist_alg, n[5], n[5]/float(n[0])*100, dist_threshold), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Retained %d [%.2f%%] reads after %s matching and filtering' % (n[6], n[6]/float(n[0])*100, dist_alg), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Elapsed time for filtering: %s' % str(filter_elapsed), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))


# create plate overview
print("Generate summary...", flush=True)
bcx_matrix = make_plate_overview_mod(n_wellx)
bcy_matrix = make_plate_overview_mod(n_welly)

# this part is super slow if you have a large number of cells, maybe omit
all_bc_counts = {i:all_bcs.count(i) for i in list(set(all_bcs))}

# calculate cumulative fraction of reads
all_bc_cumsum = np.cumsum(sorted(list(all_bc_counts.values()), reverse=True))/n_all_bcs

#%% plotting summary graphs
fig, axes = plt.subplots(2,2)
p1 = axes[0,0].imshow(np.log10(bcy_matrix+1));  axes[0,0].set_title('Number of reads per RT barcode (y-coordinate)')
clb = fig.colorbar(p1, ax=axes[0,0]); clb.set_label('No. BCs, log10')
p2 = axes[0,1].imshow(np.log10(bcx_matrix+1));  axes[0,1].set_title('Number of reads per ligation barcode (x-coordinate)')
clb = fig.colorbar(p2, ax=axes[0,1]); clb.set_label('No. BCs, log10')

# plotting the cumulative fraction of reads per barcode helps to determine the number of assayed cells (see Macosko, 2015)
p4 = axes[1,1].plot(all_bc_cumsum); axes[1,1].set_title('Cumulative fraction of reads per barcode'); axes[1,1].set_xlim(0, 5000)

fig.set_size_inches(12,7)
fig.savefig(os.path.join(out_dir,'barcode_filtering_summary.pdf'),bbox_inches='tight')

#%% Stop timing
t_stop = datetime.now()
t_elapsed = t_stop-t_start

#%% Save the QC output
f = h5py.File(os.path.join(out_dir,'splitseq_filtering_QC_data.hdf5'), 'w')

f.create_dataset('BCX_plate_overview', data = bcx_matrix)
f.create_dataset('BCY_plate_overview', data = bcy_matrix)
f.create_dataset('reads_per_BC', data=list(all_bc_counts.values()))
f.create_dataset('labels_reads_per_BC', data=np.string_(list(all_bc_counts.keys())))

f.close()

#%% print info to stdout
print("Summary generated. Elapsed time: " + str(t_elapsed), flush=True)