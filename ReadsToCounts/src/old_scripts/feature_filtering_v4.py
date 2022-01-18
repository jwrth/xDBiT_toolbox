#!/usr/bin/env python

"""
This is to process Split-seq reads after tagging the bam file with cellular and molecular barcodes.

Input: BAM file with reads that have been filtered, polyA and SMART adapter trimmed, and tagged with the following:
    - XD: cell barcode 1
    - XE: cell barcode 2
    - XF: cell barcode 3
    - XM: molecular barcode (UMI)
    - XG: Feature barcode

Algorithm:
    - Choice between levenshtein distance and hamming distance for correction

Output:
    - Filtered bam file containing only reads with expected barcodes. Barcodes that are within 1 hamming distance of an expected
    barcode sequence are corrected. An additional tag XC is added, which corresponds to the concatenated barcode sequences 
    strating with the RT barcode, then round 1 and round 2 ligation barcodes.
    - Further the feature barcodes are aligned to a reference feature barcode list and assigned to the correct antibody.

Copyright: Johannes Wirth, Meier Lab, Helmholtz Zentrum München, 2020

# The software is based on the "Split-seq toolbox" provided by Rebekka Wegmann,
# Snijderlab, ETH Zurich, 2019
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
import json
import gzip


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

def featurefilter(in_bam):
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
        out_bam = os.path.join(tmp_dir, "feat_tagged_BC_filtered.bam")

    outfile = pysam.AlignmentFile(out_bam, 'wb', template=infile)

    if create_dge:
        # read RNA cell list from RNA DGE matrix.
        if rna_dge_file.endswith('.txt.gz'):
            with gzip.open(rna_dge_file, 'rt') as f:
                first_line = f.readline()
                rna_cell_list = first_line.rstrip('\n').split('\t')[1:]
        else:
            if rna_dge_file.endswith('.txt'):
                with open(rna_dge_file, 'rt') as f:
                    first_line = f.readline()
                    rna_cell_list = first_line.rstrip('\n').split('\t')[1:]
            else:
                print('Wrong file format for RNA cell (Neither .txt not .txt.gz)', flush=True)


        # create dictionary to collect all UMIs per cell
        #umi_dict = {key:[] for key in rna_cell_list}
        umi_dict = {cell:{gene:[] for gene in fbc.Feature.values} for cell in rna_cell_list}
        #total_umi_dict = {key:[] for key in rna_cell_list}
        total_umi_dict = {cell:{gene:[] for gene in fbc.Feature.values} for cell in rna_cell_list}

        # create empty gene expression matrix
        dge = pd.DataFrame(0, index=fbc.Feature, columns=rna_cell_list,dtype=int)
    
    # set up recording array and variables
    n = np.zeros(15, dtype=np.int32)

    n_well1 = {x[0]+x[1:].zfill(2):0 for x in bc1.WellPosition.unique()}
    n_well2 = {x[0]+x[1:].zfill(2):0 for x in bc2.WellPosition.unique()}
    n_well3 = {x[0]+x[1:].zfill(2):0 for x in bc3.WellPosition.unique()}

    n_feature = {x:0 for x in fbc.Feature}

    # set up list to collect the barcodes for later statistics
    all_bcs = [None]*est_num_cells*100

    # record UMI collapsing over time
    record_collapse = []
    record_adds = []
    collapse_count = 0
    add_count = 0

    
    # start timing
    start_time_filtering = datetime.now()
    start_time = datetime.now()
    
    # start filtering
    print("Filtering started...", flush=True)
    for entry in infile.fetch(until_eof=True):
        n[0]+=1
        keep = True
        add = True

        # extract barcodes from tags
        xd = entry.get_tag('XD').rstrip('N')
        xe = entry.get_tag('XE')
        xf = entry.get_tag('XF')
        xg = entry.get_tag('XG') # feature barcode
        umi = entry.get_tag('XM')

        if len(xd) >= 6:    
            if xg in fbc.Barcode.values:

                # get feature name
                featurename = fbc_dict[xg]

                entry.set_tag('gn', featurename, value_type = 'Z')

                # record
                n[9]+=1
                n_feature[featurename]+=1

            else:
                # If barcode is not found in feature barcode list: Check for mismatches in feature barcode
                d = [compute_dist(xg,bc) for bc in fbc.Barcode.values]
                idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                if len(idx)==1:
                    xg = fbc.Barcode.values[idx[0]]
                    featurename = fbc_dict[xg]
                    entry.set_tag('gn', featurename, value_type = 'Z')

                    # record
                    n[10]+=1
                    n_feature[featurename]+=1
                else:
                    keep=False

            # Cellular barcode filtering
            if xd in bc1.Barcode.values and xe in bc2.Barcode.values and xf in bc3.Barcode.values:
                # all three barcodes directly found as in original split-seq protocol
                # record finding
                #print("All barcodes found directly")
                n[1]+=1
                n[2]+=1
                n[3]+=1
                n[4]+=1

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
                # checking for barcode 1
                if xd in bc1.Barcode.values:
                    # not all barcodes were found directly but bc1
                    #print("Barcode 1 found directly")
                    # determine well coordinate and set new tag with well coordinate
                    xd_well = bc1_dict[xd]
                    entry.set_tag('XD', xd_well)

                    # record
                    n[1]+=1
                    n_well1[xd_well]+=1

                else:
                    # barcode 1 was not found directly. Fuzzy string matching is applied
                    #print("Fuzzy string matching applied")
                    d = [compute_dist(xd,bc) for bc in bc1.Barcode.values]
                    idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                    if len(idx)==1:
                        xd = bc1.Barcode.values[idx[0]]
                        xd_well = bc1_dict[xd]
                        entry.set_tag('XD', xd_well)

                        # record
                        n[5]+=1
                        n_well1[xd_well]+=1
                    else:
                        keep=False

                # checking for barcode 2
                if xe in bc2.Barcode.values:
                    # not all barcodes were found directly but bc2
                    #print("Barcode 2 found directly")
                    # determine well coordinate and set new tag with well coordinate
                    xe_well = bc2_dict[xe]
                    entry.set_tag('XE', xe_well)

                    # record
                    n[2]+=1
                    n_well2[xe_well]+=1

                else:
                    # barcode 2 was not found directly. Fuzzy string matching is applied
                    #print("Fuzzy string matching applied")
                    
                    d = [compute_dist(xe,bc) for bc in bc2.Barcode.values]
                    idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                    if len(idx)==1:
                        xe = bc2.Barcode.values[idx[0]]
                        xe_well = bc2_dict[xe]
                        entry.set_tag('XE', xe_well)

                        # record
                        n[6]+=1
                        n_well2[xe_well]+=1
                    else:
                        keep = False

                # checking for barcode 3
                if xf in bc3.Barcode.values:
                    # not all barcodes were found directly but bc3
                    #print("Barcode 3 found directly")
                    # determine well coordinate and set new tag with well coordinate
                    xf_well = bc3_dict[xf]
                    entry.set_tag('XF', xf_well)

                    # record
                    n[3]+=1
                    n_well3[xf_well]+=1

                else:
                    # barcode 3 was not found directly. Fuzzy string matching is applied
                    #print("Fuzzy string matching applied")
                    
                    d = [compute_dist(xf,bc) for bc in bc3.Barcode.values]
                    idx = [i for i,e in enumerate(d) if e<=dist_threshold]

                    if len(idx)==1:
                        xf = bc3.Barcode.values[idx[0]]
                        xf_well = bc3_dict[xf]
                        entry.set_tag('XF', xf_well)

                        # record
                        n[7]+=1
                        n_well3[xf_well]+=1
                        #print("Keep read")
                    else:
                        #print("Not able to assign barcode uniquely. Do not keep read")
                        keep = False

        if keep:
            # concatenate well coordinates to cellular barcode
            xc_well = xd_well + xe_well + xf_well

            # set cellular barcode as tag and write to output file
            entry.set_tag('XC', xc_well)

            # add xc to cell list
            #if xc_well in cell_dict:
            #    cell_dict[xc_well]+=1
            #else:
            #    cell_dict[xc_well]=1

            # Write to output file
            outfile.write(entry)

            # record statistics
            if n[8]<est_num_cells*100/ncores:
                all_bcs[n[8]]=xc_well
            n[8]+=1

            if create_dge:
                # Generation of DGE matrix
                # If full cellular barcode was found: Check if cellular barcode of current read is in RNA cell list
                if xc_well in rna_cell_list:

                    current_umi_list = umi_dict[xc_well][featurename]

                    # record
                    n[11]+=1
                    total_umi_dict[xc_well][featurename].append(umi)


                    # Check if UMI was already used for this cell
                    if umi in current_umi_list:
                        # if UMI of current read was found in the UMI list of the respective cell: Discard read.
                        add = False

                        # record
                        n[12]+=1
                        collapse_count+=1
                        
                    else:
                        d = [lv.hamming(umi, u) for u in current_umi_list]
                        idx = [i for i,e in enumerate(d) if e<=1] # check if there is a UMI in the list with hamming distance smaller than 1

                        if len(idx) >= 1:
                            add = False

                            # record
                            n[13]+=1
                            collapse_count+=1

                else:
                    add = False

                if add:
                    # If cell barcode in RNA cell list and UMI was not used before: Count +1 in DGE.
                    dge.loc[featurename, xc_well]+=1
                    # And add UMI to UMI list
                    umi_dict[xc_well][featurename].append(umi)

                    # record
                    n[14]+=1
                    add_count+=1

            else:
                dge="No DGE generated."

        else:
            if store_discarded:
                print(entry.query_name, file = open(os.path.join(out_dir,'discarded_reads.txt'),'a'), flush=True)

        if n[0] % 100000 == 0:
            # record UMI collapsing over time
            record_collapse.append([int(n[0]), int(collapse_count)])
            record_adds.append([int(n[0]), int(add_count)])
            #collapse_count = 0
            #add_count = 0

        if n[0] % stride == 0:
            totaltime = datetime.now() - start_time_filtering
            stride_steptime = datetime.now() - start_time
            time_to_go = (total_reads - n[0])/stride * stride_steptime
            print("File " + filename + " - Reads " + str(n[0]) + "/" + str(total_reads) + " processed. Time for last " + str(stride) + ": " + str(stride_steptime) + ", Total time: " + str(totaltime) + ". Time remaining: " + str(time_to_go),
                flush=True)
            start_time = datetime.now()


    all_bcs = all_bcs[:n[8]]
    n_all_bcs = len(all_bcs)
    print("Filtering finished.", flush=True)
    
    infile.close()
    outfile.close()
    return [n, [n_well1, n_well2, n_well3], all_bcs, n_all_bcs, dge, umi_dict, record_collapse, record_adds, total_umi_dict]
    
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
parser.add_argument("-h" "--dist_threshold", action="store", dest="dist_threshold", default=1, help="Threshold to be used for levenshtein or hamming distance matching")
parser.add_argument('-m', action='store_true', help="Use multithreading?")
parser.add_argument('-s', action='store_true', help="Shorten summary?")
parser.add_argument('-r', "--rna_dge_file", action='store', dest="rna_dge_file", default="-", help="Specify RNA DGE matrix.")
parser.add_argument('-x', action='store_true', help="Create DGE matrix?")

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
dist_threshold = args.dist_threshold
multi = args.m
shorten_summary = args.s
rna_dge_file = args.rna_dge_file
create_dge = args.x

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

if multi and create_dge:
    sys.exit("Multithreading and DGE generation not compatible!")
    
#%% Write parameters to logfile
print('Splitseq barcode filtering log - based on %s algorithm\n---------------------------------------\nParameters:' % dist_alg, file = open(os.path.join(out_dir,'feature_filtering_log.txt'), 'w'))
print('Input bam: %s' % input_bam, file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Path to output bam files: %s' % split_dir, file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Estimated number of cells: %d' % est_num_cells, file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Output directory: %s' % out_dir, file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Barcode directory: %s' % bc_dir, file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))

if store_discarded:
    if os.path.isfile(os.path.join(out_dir,'discarded_reads.txt')):
        os.remove(os.path.join(out_dir,'discarded_reads.txt'))
        print('Old version of discarded reads.txt deleted.')

#%% Start timing
#t_start = timeit.default_timer()
t_start = datetime.now()

#%% Get expected barcodes and bridges
cellbc_files = glob.glob(os.path.join(bc_dir,"*cellbarcodes_*.csv"))
bc1 = pd.read_csv(cellbc_files[0])
bc2 = pd.read_csv(cellbc_files[1])
bc3 = pd.read_csv(cellbc_files[2])

# Get expected feature barcodes
fbc_file = glob.glob(os.path.join(bc_dir,"*featurebarcodes*.csv"))
fbc = pd.read_csv(fbc_file[0])


#%% Generate dictionaries for barcode-well and barcode-feature assignment
bc1_dict = {bc1.Barcode.values[x]:bc1.WellPosition.values[x][0]+bc1.WellPosition.values[x][1:].zfill(2) for x in range(len(bc1))}
bc2_dict = {bc2.Barcode.values[x]:bc2.WellPosition.values[x][0]+bc2.WellPosition.values[x][1:].zfill(2) for x in range(len(bc2))}
bc3_dict = {bc3.Barcode.values[x]:bc3.WellPosition.values[x][0]+bc3.WellPosition.values[x][1:].zfill(2) for x in range(len(bc3))}

fbc_dict = {fbc.Barcode.values[x]:fbc.Feature.values[x] for x in range(len(fbc))}

# Generate list with all cellular barcodes 
#cell_dict = dict()


#%% for debugging only
if debug_flag:
    for entry in itertools.islice(infile, 10 ):
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
        results = multifilter.map(featurefilter, files)

    # extract and sum up recording variables
    n = np.array([elem[0] for elem in results]).sum(axis=0)

    n_well1 = sum_dicts([elem[1][0] for elem in results])
    n_well2 = sum_dicts([elem[1][1] for elem in results])
    n_well3 = sum_dicts([elem[1][2] for elem in results])

    all_bcs = [item for sublist in [elem[2] for elem in results] for item in sublist]
    n_all_bcs = np.array([elem[3] for elem in results]).sum()

else:
    # without multithreading
    ncores = 1
    results = np.array(featurefilter(input_bam))

    n = results[0]

    n_well1 = results[1][0]
    n_well2 = results[1][1]
    n_well3 = results[1][2]

    all_bcs = results[2]
    n_all_bcs = results[3]
    dge = results[4]
    umi_dict = results[5]
    record_collapse = results[6]
    record_adds = results[7]
    total_umi_dict = results[8]

filter_stop = datetime.now()
filter_elapsed = filter_stop - t_start

if create_dge:
    # save DGE matrix file as .txt
    print('Write DGE matrix...', flush=True)
    dge.to_csv(os.path.join(out_dir, "DGE_matrix_features.txt.gz"), sep = '\t', compression= 'gzip', index=True, header=True)

# print into log file
print('Read %d entries' % n[0], file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] complete cellular barcodes' % (n[4], n[4]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected cellular BC1, with fuzzy matching %d [%.2f%%] BC1' % (n[1] , n[1]/float(n[0])*100, n[5], n[5]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected cellular BC2, with fuzzy matching %d [%.2f%%] BC2' % (n[2] , n[2]/float(n[0])*100, n[6], n[6]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected cellular BC3, with fuzzy matching %d [%.2f%%] BC3' % (n[3] , n[3]/float(n[0])*100, n[7], n[7]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] expected feature barcodes, with fuzzy matching %d [%.2f%%] feature barcodes' % (n[9], n[9]/float(n[0])*100, n[10], n[10]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Retained %d [%.2f%%] reads after fuzzy matching and filtering' % (n[8], n[8]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
print('Elapsed time for filtering: ' + str(filter_elapsed), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))

if create_dge:
    print('', file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('Statistics about DGE matrix generation:', file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('RNA cell list: ' + str(rna_dge_file), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('Found %d [%.2f%%] complete cellular barcodes in RNA cell list' % (n[11], n[11]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('%d [%.2f%%] cases where UMI was found directly in UMI list and was not added to DGE matrix.' % (n[12], n[12]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('%d [%.2f%%] cases where UMI was found with hamming distance of 1 in UMI list and was not added to DGE matrix.' % (n[13], n[13]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
    print('Found %d [%.2f%%] unique UMIs with correct cell and feature barcode that were added to the DGE matrix' % (n[14], n[14]/float(n[0])*100), file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))
else:
    print('No DGE matrix generated.', file = open(os.path.join(out_dir,'feature_filtering_log.txt'),'a'))




# create plate overview
print("Generate summary...", flush=True)
bc1_matrix = make_plate_overview_mod(n_well1)
bc2_matrix = make_plate_overview_mod(n_well2)
bc3_matrix = make_plate_overview_mod(n_well3)

# this part is super slow if you have a large number of cells, maybe omit
if not shorten_summary:
    all_bc_counts = {i:all_bcs.count(i) for i in list(set(all_bcs))}

    # calculate cumulative fraction of reads
    all_bc_cumsum = np.cumsum(sorted(list(all_bc_counts.values()), reverse=True))/n_all_bcs

#%% plotting summary graphs
fig, axes = plt.subplots(2,2)
p1 = axes[0,0].imshow(np.log10(bc1_matrix+1));  axes[0,0].set_title('Number of reads per RT barcode')
clb = fig.colorbar(p1, ax=axes[0,0]); clb.set_label('No. BCs, log10')
p2 = axes[0,1].imshow(np.log10(bc2_matrix+1));  axes[0,1].set_title('Number of reads per round 2 barcode')
clb = fig.colorbar(p2, ax=axes[0,1]); clb.set_label('No. BCs, log10')
p3 = axes[1,0].imshow(np.log10(bc3_matrix+1));  axes[1,0].set_title('Number of reads per round 3 barcode')
clb = fig.colorbar(p3, ax=axes[1,0]); clb.set_label('No. BCs, log10')

if not shorten_summary:
    # plotting the cumulative fraction of reads per barcode helps to determine the number of assayed cells (see Macosko, 2015)
    p4 = axes[1,1].plot(all_bc_cumsum); axes[1,1].set_title('Cumulative fraction of reads per barcode'); axes[1,1].set_xlim(0, 5000)

fig.set_size_inches(12,7)
fig.savefig(os.path.join(out_dir,'feature_filtering_summary.pdf'),bbox_inches='tight')

#%% Stop timing
t_stop = datetime.now()
t_elapsed = t_stop-t_start

#%% Save the QC output
f = h5py.File(os.path.join(out_dir,'splitseq_feature_filtering_QC_data.hdf5'), 'w')

f.create_dataset('BC1_plate_overview', data = bc1_matrix)
f.create_dataset('BC2_plate_overview', data = bc2_matrix)
f.create_dataset('BC3_plate_overview', data = bc3_matrix)
if not shorten_summary:
    f.create_dataset('reads_per_BC', data=list(all_bc_counts.values()))
    f.create_dataset('labels_reads_per_BC', data=np.string_(list(all_bc_counts.keys())))

f.close()

# Save the UMI dictionary as .json
print('Saving UMI dictionary...', flush=True)
umi_dict_file = open(os.path.join(tmp_dir, 'umi_dictionary.json'), 'w')
json.dump(umi_dict, umi_dict_file)
umi_dict_file.close()

# Save the recorded UMI collapsing events as .json
print('Saving UMI records...', flush=True)
record_umi_file = open(os.path.join(tmp_dir, 'umi_record.json'), 'w')
json.dump(record_collapse, record_umi_file)
record_umi_file.close()

# Save the recorded addition events as .json
print('Saving addition events...', flush=True)
record_adds_file = open(os.path.join(tmp_dir, 'addition_record.json'), 'w')
json.dump(record_adds, record_adds_file)
record_adds_file.close()

# Save the total umi list
print('Saving total UMI list...', flush=True)
total_umi_file = open(os.path.join(tmp_dir, 'total_umi_dict.json'), 'w')
json.dump(total_umi_dict, total_umi_file)
total_umi_file.close()


#%% print info to stdout
print("Summary generated. Elapsed time: " + str(t_elapsed), flush=True)