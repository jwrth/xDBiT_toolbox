#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is to process AbxDbit reads after tagging the bam file with cellular and molecular barcodes.

Input: BAM file with reads that have been filtered, polyA and SMART adapter trimmed, and tagged with the following:
    - XX: x-coordinate
    - XY: y-coordinate
    - (OPTIONAL) XZ: z-coordinate (well coordinate in xDbit)
    - XM: molecular barcode (UMI)

Algorithm:
    - Choice between levenshtein distance and hamming distance for correction

Output:
    - Filtered bam file containing only reads with expected barcodes. Barcodes that are within 1 hamming distance of an expected
    barcode sequence are corrected. An additional tag XC is added, which consists of the concatenated coordinates separated by an 'x' (XxYxZ).


Copyright: Johannes Wirth, Meier Lab, Helmholtz Pioneer Campus, Helmholtz Zentrum Munich, 2021


The code is based on work with a Copyright by Rebekka Wegmann (Snijderlab, ETH Zurich, 2019) which has been published in following GitHub repository: 
https://github.com/RebekkaWegmann/splitseq_toolbox

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

## libraries
import sys, os, h5py
import pysam
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from argparse import ArgumentParser
from datetime import datetime, timedelta
import subprocess
import glob
import collections
from multiprocessing import Pool
import Levenshtein as lv
import json
import gzip
import string
from itertools import combinations


## Functions

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

class xDbit_filter:
    def __init__(self, args):
        ## Start timing
        self.t_start = datetime.now()
        
        # retrieve arguments
        self.debug_flag = args.debug_flag
        self.store_discarded = args.store_discarded
        self.input_bam = args.input_bam #use "-" for stdin, set flag to rb
        self.est_num_cells = args.est_num_cells
        self.out_dir = args.out_dir
        self.tmp_dir = args.tmp_dir
        self.bc_file = args.bc_file
        self.multi = args.m
        self.mode = args.mode
        self.interact_mode = args.interact
        self.feat_file = args.feature_file
        self.rna_dge_file = args.rna_dge_file
        self.stride = int(args.stride)
        
        # for multithreading define temp directory for split files
        self.split_dir = os.path.join(self.tmp_dir, 'tmp_split')
        
        # generate different objects from input
        # check out mode
        if self.mode == 'Dbit-seq':
            self.coord_names = ['X', 'Y']
        elif self.mode == 'xDbit':
            self.coord_names = ['X', 'Y', 'Z']
        else:
            # exit script
            sys.exit('{} is no valid mode ["xDbit", "Dbit-seq"]'.format(self.mode))
            
        # Import barcode legend and extract information
        self.barcode_legend = pd.read_csv(self.bc_file)
        
        # generate barcode dictionary
        self.barcode_dicts = {name: self.create_barcode_dict(self.barcode_legend, name) for name in self.coord_names}
        
        # retrieve string matching settings
        self.dist_alg = self.barcode_legend['string_matching_algorithm'][0]
        self.dist_thresholds = {name: int(self.barcode_legend[name + "_maxdist"][0]) for name in self.coord_names}
        
        # check which distance algorithm to use
        if self.dist_alg == "hamming":
            self.compute_dist = lv.hamming
        elif self.dist_alg == "levenshtein":
            self.compute_dist = lv.distance
        else:
            print("Unknown string matching alogrithm. Used levenshtein as default.", flush=True)
            self.compute_dist = lv.distance
            
            # Check if feature legend file is given
        if self.feat_file is not None:
            # prepare everything for feature extraction if necessary
            # import feature legend
            self.feature_legend = pd.read_csv(self.feat_file)
            
            # read distance algorithm and distance thresholds for features
            self.feat_dist_alg = self.feature_legend['string_matching_algorithm'][0]
            self.feat_dist_threshold = self.feature_legend['maxdist'][0]

            # create feature barcode dict
            self.feature_dict = {self.feature_legend.Barcode.values[x]:self.feature_legend.Feature.values[x] for x in range(len(self.feature_legend))}
            self.feature_barcodes = self.feature_legend.Barcode.values

            # read RNA cell list from RNA DGE matrix.
            if self.rna_dge_file.endswith('.txt.gz'):
                with gzip.open(self.rna_dge_file, 'rt') as f:
                    first_line = f.readline()
                    self.rna_spot_list = first_line.rstrip('\n').split('\t')[1:]
            elif self.rna_dge_file.endswith('.txt'):
                with open(self.rna_dge_file, 'rt') as f:
                    first_line = f.readline()
                    self.rna_spot_list = first_line.rstrip('\n').split('\t')[1:]
            else:
                sys.exit('Wrong file format for RNA cell (Neither .txt not .txt.gz)')

            # create dictionary to collect all UMIs per cell
            self.umi_dict = {cell:{gene:[] for gene in self.feature_legend.Feature.values} for cell in self.rna_spot_list}
            self.total_umi_dict = {cell:{gene:[] for gene in self.feature_legend.Feature.values} for cell in self.rna_spot_list}
            
            if self.interact_mode:
                # add also all combinations of antibodies as entries in umi dictionaries
                for dic in [self.umi_dict, self.total_umi_dict]:
                    for cell in dic:
                        for v in self.feature_legend.Feature.values:
                            for w in self.feature_legend.Feature.values:
                                dic[cell][v + "+" + w] = []
                            
            # create empty gene expression matrix
            self.dge = pd.DataFrame(0, index=self.feature_legend.Feature, columns=self.rna_spot_list, dtype=int)
            
        ## Write parameters to logfile
        print('AbxDbit barcode filtering log - based on %s algorithm\n---------------------------------------\nParameters:' % self.dist_alg, file = open(os.path.join(self.out_dir,'filtering_log.txt'), 'w'))
        print('Input bam: %s' % self.input_bam, file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        print('Path to output bam files: %s' % self.split_dir, file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        print('Estimated number of cells: %d' % self.est_num_cells, file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        print('Output directory: %s' % self.out_dir, file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        #print('Barcode directory: %s' % bc_dir, file = open(os.path.join(out_dir,'filtering_log.txt'),'a'))
        print('Barcode legend file: %s' % self.bc_file, file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        
        if self.store_discarded:
            if os.path.isfile(os.path.join(self.out_dir,'discarded_reads.txt')):
                os.remove(os.path.join(self.out_dir,'discarded_reads.txt'))
                print('Old version of discarded reads.txt deleted.')
                
    def run_filtering(self):
        ### Run Filtering
        if self.multi:
            ## start multithreaded filtering   
            
            files = glob.glob(os.path.join(self.split_dir, 'x*.bam'))
            self.ncores = len(files)
            multifilter = Pool(processes=self.ncores)
            self.results = multifilter.map(self.readfilter, files)

            # extract recording dictionaries
            self.record_dicts = [elem[0] for elem in self.results]

            # Sum up the recording dictionaries
            first_level = ['total_count', 'total_count_kept', 'all_direct']
            self.record_dict = sum_dicts([{k: d[k] for k in first_level} for d in self.record_dicts])

            second_level = ['direct', 'corrected']

            for name in self.coord_names:
                dicts_per_name = [d[name] for d in self.record_dicts]
                self.record_dict[name] = sum_dicts([{k: d[k] for k in second_level} for d in dicts_per_name])
        
                self.record_dict[name]['well_counts'] = sum_dicts([d[name]['well_counts'] for d in self.record_dicts])

            # extract recording values outside of recording dictionary
            self.all_bcs = [item for sublist in [elem[1] for elem in self.results] for item in sublist]
            self.n_all_bcs = np.array([elem[2] for elem in self.results]).sum()

        else:
            # without multithreading
            self.ncores = 1
            self.record_dict, self.all_bcs, self.n_all_bcs = self.readfilter(self.input_bam)            
        # take time
        self.filter_stop = datetime.now()
        self.filter_elapsed = self.filter_stop - self.t_start
            
    def save_files(self):
        # Store recording dictionary as .json
        with open(os.path.join(self.out_dir, 'recording_dictionary.json'), 'w') as fp:
            json.dump(self.record_dict, fp)
            
        ## Saving feature results
        if self.feat_file is not None:
            if not np.any([elem.startswith("feat_") for elem in self.dge.index]):
                # add feat_ prefix to feature names
                self.dge.rename(index=lambda s: "feat_" + s, inplace=True)

            # save feature DGE matrix file as .txt
            self.dge_out_file = os.path.join(self.out_dir, "DGE_matrix_features.txt.gz")
            print('Write feature DGE matrix to {}...'.format(self.dge_out_file), flush=True)
            self.dge.to_csv(self.dge_out_file, sep = '\t', compression= 'gzip', index=True, header=True)

            # combining RNA and feature reads
            self.rna = pd.read_csv(self.rna_dge_file, compression='gzip', sep='\t', index_col=0)
            self.combined = pd.concat([self.rna, self.dge], axis=0)
            self.combined_out_file = os.path.join(self.out_dir, "DGE_matrix_rna_with_features.txt.gz")
            print('Write combined DGE matrix to {}...'.format(self.combined_out_file), flush=True)
            self.combined.to_csv(self.combined_out_file, sep = '\t', compression= 'gzip', index=True, header=True)


    def well2ind(self, well):
        """Convert well positions to (zero-based) indices"""
        d = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7}
        row = d[well[0]]
        col = int(well[1:])-1
        return [row,col]

    def make_plate_overview(self, well_count):
        # bring the counts per well into a plate layout
        out = np.zeros([8, 12])
        for wellpos in well_count.items():
            xy = self.well2ind(wellpos[0])
            out[xy[0], xy[1]] = int(wellpos[1])
        return out

    def plot_plate_overview(self, bc_matrices, all_bc_cumsum, coord_names, est_num_cells, out_dir):
        ## plotting summary graphs
        fig, axs = plt.subplots(2,2)
        axs = axs.ravel()

        for i, name in enumerate(coord_names):
            # plot heat map
            matrix = bc_matrices[name]
            p = axs[i].imshow(np.log10(matrix+1))
            
            # set title
            axs[i].set_title('Number of reads per {}-coordinate'.format(name))
            
            # set x and y ticks
            axs[i].set_xticks(list(range(0, 12)), list(range(1, 13)))
            axs[i].set_yticks(list(range(0, 8)), list(string.ascii_uppercase[:8]))
            axs[i].tick_params(axis='both', which='both',length=0)
            
            # set color bar
            clb = fig.colorbar(p, ax=axs[i])
            clb.set_label('No. BCs, log10')

        # plotting the cumulative fraction of reads per barcode helps to determine the number of assayed cells (see Macosko, 2015)
        axs[3].plot(all_bc_cumsum); axs[3].set_title('Cumulative fraction of reads per barcode')
        axs[3].set_xlim(0, est_num_cells * 1.2)
        axs[3].set_xlabel("Barcodes sorted by number of reads")
        axs[3].set_ylabel("Cumulative fraction of reads")

        fig.set_size_inches(12,7)
        fig.savefig(os.path.join(out_dir,'filtering_summary.pdf'),bbox_inches='tight')

    def create_barcode_dict(self, barcode_legend, coordinate_name):
        '''
        Function to generate barcode dictionary.
        Expects following inputs:
            1. 'barcode_legend': data from .csv file with columns: [WellPosition, Name, Barcode, X, Y, Z]
            2. 'coordinate_name': name of coordinate.
        '''
        d = {barcode: [well[0] + well[1:].zfill(2), int(coord)] for (barcode, well, coord) in zip(barcode_legend["Barcode"],
                                                            barcode_legend["WellPosition"],
                                                            barcode_legend[coordinate_name]) if not pd.isnull(coord)}
        return d

    def check_barcode_and_add_coord_to_tag(self, entry, barcodes, coord_name, 
                                           barcode_dictionaries, record_dictionary, 
                                           dist_threshold, compute_dist_function, 
                                           keep
                                           ):

        '''
        Function to check a barcode against a list and add a coordinate tag to a pysam entry.
        Procedure:
            1. Check barcode against list.
            2. If it does not match exactly a distance computation function is applied.
            3. If unique match is found the barcode is transformed to the spatial coordinate using the barcode dictionary.
            4. Pysam entry is tagged with coordinate.
            5. If no match is found the keep variable is set to False.
        '''
        
        barcode = barcodes[coord_name]
        barcode_dict = barcode_dictionaries[coord_name] # lookup in dictionaries is faster
        barcode_list = list(barcode_dict.keys()) # iteration through lists is slightly faster

        if barcode in barcode_dict:
            # Not all barcodes were found directly but bc1
            # Determine well coordinate and set new tag with well coordinate
            well = barcode_dict[barcode][0]
            coord = barcode_dict[barcode][1]

            entry.set_tag("X" + coord_name, str(coord))

            # record
            record_dictionary[coord_name]['direct']+=1
            record_dictionary[coord_name]['well_counts'][well]+=1

        else:
            # Barcode 1 was not found directly. Barcode with certain distance is searched.
            d = [compute_dist_function(barcode,bc) for bc in barcode_list]
            idx = [i for i,e in enumerate(d) if e <= dist_threshold]

            if len(idx)==1:
                barcode = barcode_list[idx[0]]

                well = barcode_dict[barcode][0]
                coord = barcode_dict[barcode][1]

                entry.set_tag("X" + coord_name, str(coord))

                # record
                record_dictionary[coord_name]['corrected']+=1
                record_dictionary[coord_name]['well_counts'][well]+=1

            else:
                keep = False
                coord = None
                well = None

        return entry, record_dictionary, coord, well, keep

    def readfilter(self, in_bam):
        # count read length of input bam
        print("Count number of reads in input bam...", flush=True)
        bam_length = get_stdout('samtools view ' + in_bam + ' | wc -l')
        total_reads = int(bam_length.decode())

        # read input bam
        infile = pysam.AlignmentFile(in_bam, 'rb', check_sq=False)
        filename = in_bam.split('/')[-1]

        # generate output bam name and directory
        if self.multi:
            out_bam = '/'.join(in_bam.split('/')[:-1]) + '/out_' + in_bam.split('/')[-1]
        else:
            out_bam = os.path.join(self.tmp_dir, "unaligned_tagged_BC_filtered.bam")

        outfile = pysam.AlignmentFile(out_bam, 'wb', template=infile)

        # setup recording dictionary
        self.record_dict = {
        'total_count': 0, 
        'total_count_kept': 0,
        'all_direct': 0
        }

        for name in self.coord_names:
            self.record_dict[name] = {}
            self.record_dict[name]['direct'] = 0
            self.record_dict[name]['corrected'] = 0
            self.record_dict[name]['well_counts'] = {e[0][0]+e[0][1:].zfill(2):0 for e in self.barcode_dicts[name].values()}

        # set up list to collect the barcodes for later statistics
        all_bcs =  [None]*self.est_num_cells*100

        if self.feat_file is not None:
            # add features category to record dictionary
            self.record_dict['features'] = {}
            self.record_dict['features']['direct'] = 0
            self.record_dict['features']['corrected'] = 0
            self.record_dict['features']['spot_found'] = 0
            self.record_dict['features']['umi_found'] = 0
            self.record_dict['features']['added'] = 0
            
            for v in self.feature_legend.Feature.values:
                self.record_dict['features'][v] = 0
            
            # for IntAct-seq
            self.record_dict['interact'] = {}
            self.record_dict['interact']['direct'] = 0
            self.record_dict['interact']['corrected'] = 0
            self.record_dict['interact']['spot_found'] = 0
            self.record_dict['interact']['umi_found'] = 0
            self.record_dict['interact']['added'] = 0
            
            interaction_combinations = ["+".join(map(str, comb)) for comb in combinations(self.feature_legend.Feature.values, 2)]
            for c in interaction_combinations:
                self.record_dict['interact'][c] = 0
            
        # start timing
        start_time_filtering = datetime.now()
        start_time = datetime.now()

        # start filtering
        print("Filtering started...", flush=True)

        for entry in infile.fetch(until_eof=True):
            self.record_dict['total_count']+=1
            keep = True
            add = True
            interact_signal = False

            # extract barcodes from entry
            barcodes = {name: entry.get_tag("X" + name) for name in self.coord_names}

            # check if all barcodes are in the respective barcode sets
            if np.all([barcodes[name] in self.barcode_dicts[name] for name in self.coord_names]):
                # if all barcodes directly found:
                # get well coordinate
                wells = {name: self.barcode_dicts[name][bc][0] for (name, bc) in barcodes.items()}

                # get coordinate
                coords = {name: self.barcode_dicts[name][bc][1] for (name, bc) in barcodes.items()}

                # save coordinate as tag
                for name in self.coord_names:
                    entry.set_tag('X' + name, str(coords[name]))


                # record findings
                for name in self.coord_names:
                    self.record_dict[name]['direct']+=1

                    well = wells[name]
                    self.record_dict[name]['well_counts'][well]+=1

                self.record_dict['all_direct']+=1            

            else:
                coords = {}
                wells = {}
                for name in self.coord_names:
                    
                    # check barcode and add coordinate to tag if matching barcode found.
                    entry, self.record_dict, coord, well, keep = self.check_barcode_and_add_coord_to_tag(
                        entry=entry, barcodes=barcodes, coord_name=name, 
                        barcode_dictionaries=self.barcode_dicts, record_dictionary=self.record_dict, 
                        dist_threshold=self.dist_thresholds[name], compute_dist_function=self.compute_dist,
                        keep=keep)

                    # collect the coordinates
                    coords[name] = coord
                    wells[name] = well

            if keep:
                # concatenate coordinates
                coords_list = [str(coords[name]) for name in self.coord_names]
                coord_concat = "x".join(coords_list)

                # check if features have to be extracted too
                if self.feat_file is not None:
                    xg = entry.get_tag('XG')
                    
                    # check if there is an IntAct read in the entry
                    if entry.has_tag('XH'):
                        xh = entry.get_tag('XH')
                        interact_signal = True
                    else:
                        xh = None
                    
                    umi = entry.get_tag('XM')

                    # check first possible feature tag
                    if xg in self.feature_dict:
                        # get feature name
                        featurename = self.feature_dict[xg]
                        entry.set_tag('gn', featurename, value_type='Z')

                        # record
                        self.record_dict['features']['direct']+=1
                        self.record_dict['features'][featurename]+=1

                    else:
                        # If barcode is not found in feature barcode list: Check for mismatches in feature barcode
                        d = [self.compute_dist(xg,bc) for bc in self.feature_barcodes]
                        idx = [i for i,e in enumerate(d) if e<=self.feat_dist_threshold]

                        if len(idx)==1:
                            xg = self.feature_barcodes[idx[0]]
                            featurename = self.feature_dict[xg]
                            entry.set_tag('gn', featurename, value_type = 'Z')

                            # record
                            self.record_dict['features']['corrected']+=1
                            self.record_dict['features'][featurename]+=1
                        else:
                            keep=False
                    
                    if keep and interact_signal:
                        # check first possible feature tag
                        if xh in self.feature_dict:
                            # add feature name to first featurename to mark an interact signal
                            featurename = featurename + "+" + self.feature_dict[xh]
                            entry.set_tag('gn', featurename, value_type='Z')

                            # record
                            self.record_dict['interact']['direct']+=1
                            
                            if featurename in self.record_dict['interact'].keys():
                                self.record_dict['interact'][featurename]+=1
                            else:
                                self.record_dict['interact'][featurename]=1

                        else:
                            # If barcode is not found in feature barcode list: Check for mismatches in feature barcode
                            d = [self.compute_dist(xh,bc) for bc in self.feature_barcodes]
                            idx = [i for i,e in enumerate(d) if e<=self.feat_dist_threshold]

                            if len(idx)==1:
                                xh = self.feature_barcodes[idx[0]]
                                featurename = featurename + "+" + self.feature_dict[xh]
                                entry.set_tag('gn', featurename, value_type = 'Z')

                                # record
                                self.record_dict['interact']['corrected']+=1
                                
                                if featurename in self.record_dict['interact'].keys():
                                    self.record_dict['interact'][featurename]+=1
                                else:
                                    self.record_dict['interact'][featurename]=1
                                
                            else:
                                keep=False
                        

            if keep:
                # set coordinates as tag and write to output file
                entry.set_tag('XC', str(coord_concat))
                outfile.write(entry)

                # record statistics
                if self.record_dict['total_count_kept'] < self.est_num_cells * 100 / self.ncores:
                    all_bcs[self.record_dict['total_count_kept']] = coord_concat
                self.record_dict['total_count_kept']+=1

                if self.feat_file is not None:
                    if coord_concat in self.rna_spot_list:
                        current_umi_list = self.umi_dict[coord_concat][featurename]

                        # record
                        if interact_signal:
                            self.record_dict['interact']['spot_found']+=1
                        else:
                            self.record_dict['features']['spot_found']+=1
                        
                        if featurename in self.total_umi_dict[coord_concat].keys():
                            self.total_umi_dict[coord_concat][featurename].append(umi)
                        else:
                            self.total_umi_dict[coord_concat][featurename] = [umi]

                        # Check if UMI was already used for this spot
                        if umi in current_umi_list:
                            # if UMI of current read was found in the UMI list of the respective spot: Discard read.
                            add = False

                            # record
                            if interact_signal:
                                self.record_dict['interact']['umi_found']+=1
                            else:
                                self.record_dict['features']['umi_found']+=1
                        
                    else:
                        add = False
                    
                    if add:
                        # If cell barcode in RNA cell list and UMI was not used before: Count +1 in DGE.
                        if featurename in self.dge.index:
                            self.dge.loc[featurename, coord_concat]+=1
                        else:
                            # if interaction featurename is not yet in dge dataframe add it here and add count
                            self.dge.loc[featurename] = [0] * self.dge.shape[1]
                            self.dge.loc[featurename, coord_concat]+=1
                        # And add UMI to UMI list
                        if featurename in self.umi_dict[coord_concat].keys():
                            self.umi_dict[coord_concat][featurename].append(umi)
                        else:
                            self.umi_dict[coord_concat][featurename] = [umi]

                        # record
                        if interact_signal:
                            self.record_dict['interact']['added']+=1
                        else:
                            self.record_dict['features']['added']+=1

            elif self.store_discarded:
                    print(entry.query_name, file = open(os.path.join(self.out_dir,'discarded_reads.txt'),'a'), flush=True)

            total_count = self.record_dict['total_count']
            if total_count % self.stride == 0:
                totaltime = datetime.now() - start_time_filtering
                stride_steptime = datetime.now() - start_time
                time_to_go = (total_reads - total_count) / self.stride * stride_steptime
                print("File " + filename + " - Reads " + str(total_count) + "/" + str(total_reads) + " processed. Time for last " + str(self.stride) + ": " + str(stride_steptime) + ", Total time: " + str(totaltime) + ". Time remaining: " + str(time_to_go),
                    flush=True)
                start_time = datetime.now()
                
        all_bcs = all_bcs[:self.record_dict['total_count_kept']]
        n_all_bcs = len(all_bcs)
        print("Filtering finished.", flush=True)
        
        infile.close()
        outfile.close()
        #return [n, [n_wellx, n_welly], all_bcs, n_all_bcs]
        return [self.record_dict, all_bcs, n_all_bcs]
    
    def print_to_log(self, ):
        # Print to log file
        print('Read %d entries' % self.record_dict['total_count'], file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        print('Found %d [%.2f%%] complete spot barcodes' % (self.record_dict['all_direct'], self.record_dict['all_direct']/float(self.record_dict['total_count'])*100), file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        for name in self.coord_names:
                print('Found %d [%.2f%%] expected %s-coordinates, with %s matching %d [%.2f%%] (distance: %d)' % (self.record_dict[name]['direct'] , 
                    self.record_dict[name]['direct'] / float(self.record_dict['total_count']) * 100, name, self.dist_alg, 
                    self.record_dict[name]['corrected'], self.record_dict[name]['corrected']/float(self.record_dict['total_count'])*100, self.dist_thresholds[name]),
                file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        
        print('Retained %d [%.2f%%] reads after %s matching and filtering' % (self.record_dict['total_count_kept'], self.record_dict['total_count_kept']/self.record_dict['total_count']*100, 
            self.dist_alg), file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))

        if self.feat_file is not None:
            print('', file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Statistics about DGE matrix generation for feature reads:', file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('RNA cell list retrieved from: ' + str(self.rna_dge_file), file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Found %d [%.2f%%] expected feature barcodes' % (self.record_dict['features']['direct'], self.record_dict['features']['direct']/self.record_dict['total_count']*100), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Found %d [%.2f%%] expected feature barcodes with distance %d' % (self.record_dict['features']['corrected'], self.record_dict['features']['corrected']/self.record_dict['total_count']*100, self.feat_dist_threshold), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Found %d [%.2f%%] complete cellular barcodes in RNA cell list' % (self.record_dict['features']['spot_found'], self.record_dict['features']['spot_found']/self.record_dict['total_count']*100), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Found %d [%.2f%%] cases where UMI was found directly in UMI list and was not added to DGE matrix.' % (self.record_dict['features']['umi_found'], self.record_dict['features']['umi_found']/self.record_dict['total_count']*100), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Found %d [%.2f%%] unique UMIs with correct cell and feature barcode that were added to the DGE matrix' % (self.record_dict['features']['added'], self.record_dict['features']['added']/self.record_dict['total_count']*100), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            
            # write log for antibody features
            for v, b in zip(self.feature_legend.Feature.values, self.feature_legend.Barcode.values):
                print('Antibody feature %s (%s) was found %d times' % (v, b, self.record_dict['features'][v]), 
                        file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            
            if self.interact_mode:
                print('-------------------------------', file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Interaction reads', file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Found %d [%.2f%%] expected feature barcodes' % (self.record_dict['interact']['direct'], self.record_dict['interact']['direct']/self.record_dict['total_count']*100), 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Found %d [%.2f%%] expected feature barcodes with distance %d' % (self.record_dict['interact']['corrected'], self.record_dict['interact']['corrected']/self.record_dict['total_count']*100, self.feat_dist_threshold), 
                        file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Found %d [%.2f%%] complete cellular barcodes in RNA cell list' % (self.record_dict['interact']['spot_found'], self.record_dict['interact']['spot_found']/self.record_dict['total_count']*100), 
                        file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Found %d [%.2f%%] cases where UMI was found directly in UMI list and was not added to DGE matrix.' % (self.record_dict['interact']['umi_found'], self.record_dict['interact']['umi_found']/self.record_dict['total_count']*100), 
                        file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                print('Found %d [%.2f%%] unique UMIs with correct cell and feature barcode that were added to the DGE matrix' % (self.record_dict['interact']['added'], self.record_dict['interact']['added']/self.record_dict['total_count']*100), 
                        file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                # write log for interaction features
                for c in self.record_dict['interact'].keys():
                    if c not in ["direct", "corrected", "spot_found", "umi_found", "added"]:
                        print('Interaction feature %s was found %d times' % (c, self.record_dict['interact'][c]), 
                                file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
                    
            
            print('Feature spot-count matrix was saved into %s.' % self.dge_out_file, 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
            print('Combined spot-count matrix was saved into %s.' % self.combined_out_file, 
                    file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))

        print('Elapsed time for filtering: %s' % str(self.filter_elapsed), file = open(os.path.join(self.out_dir,'filtering_log.txt'),'a'))
        
    def save_summary_and_qc(self):
        # Create plate overview
        print("Generate summary...", flush=True)
        bc_matrices = {name: self.make_plate_overview(self.record_dict[name]['well_counts']) for name in self.coord_names}

        # This part is super slow if you have a large number of cells, maybe omit
        all_bc_counts = {i:self.all_bcs.count(i) for i in list(set(self.all_bcs))}

        # Calculate cumulative fraction of reads
        all_bc_cumsum = np.cumsum(sorted(list(all_bc_counts.values()), reverse=True))/self.n_all_bcs
        
        # plot overview of plate
        self.plot_plate_overview(bc_matrices, all_bc_cumsum, self.coord_names, self.est_num_cells, self.out_dir)
        
        ## Save the QC output
        f = h5py.File(os.path.join(self.out_dir,'abxdbit_filtering_QC_data.hdf5'), 'w')

        for name in self.coord_names:
            f.create_dataset('BC{}_plate_overview'.format(name), data=bc_matrices[name])

        f.create_dataset('reads_per_BC', data=list(all_bc_counts.values()))
        f.create_dataset('labels_reads_per_BC', data=np.string_(list(all_bc_counts.keys())))

        f.close()

        # save dictionaries
        if self.feat_file is not None:
            # Save the UMI dictionary as .json
            print('Saving UMI dictionary...', flush=True)
            umi_dict_file = open(os.path.join(self.tmp_dir, 'umi_dictionary.json'), 'w')
            json.dump(self.umi_dict, umi_dict_file)
            umi_dict_file.close()
            
        # save record dict
        print('Saving record dictionary...', flush=True)
        record_dict_file = open(os.path.join(self.tmp_dir, 'record_dictionary.json'), 'w')
        json.dump(self.record_dict, record_dict_file)
        record_dict_file.close()
        
        print("Summary and QC files saved.", flush=True)
        
    def finish(self):
        ## Stop timing
        self.t_stop = datetime.now()
        self.t_elapsed = self.t_stop-self.t_start
        print("-------------------------------------------", flush=True)
        print("Elapsed time: " + str(self.t_elapsed), flush=True)
        print("-------------------------------------------", flush=True)
        
if __name__ == '__main__':
    ## Run Script
    print("--------------------------")
    print("--xDbit filtering script--")
    print("--------------------------")
    
    # Setup input parser
    parser = ArgumentParser()
    parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Specify the input bam file. Defaults to stdin.")
    parser.add_argument("-n" "--est_num_cells", action="store", dest="est_num_cells", default=2500, help="Estimated number of cells. Defaults to 2500.",type=int)
    parser.add_argument("-d" "--out_dir", action="store", dest="out_dir", default=".", help="Directory to store logfiles and output plots. Defaults to the current directory.")
    parser.add_argument("-t" "--tmp_dir", action="store", dest="tmp_dir", default=".", help="Temp directory")
    parser.add_argument("-b" "--bc_file", action="store", dest="bc_file", default=None, help="Path to the barcode legend. Defaults to a file barcodes_legend.csv in the current directory.")
    parser.add_argument("--debug_flag",action="store_true",help="Turn on debug flag. This will produce some additional output which might be helpful.")
    parser.add_argument("--store_discarded",action="store_true",help="Store names of discarded reads?")
    parser.add_argument("-m", action='store_true', help="Use multithreading?")
    parser.add_argument("--mode", action='store', default='xDbit', help="Which mode? xDbit or Dbit-seq?")
    parser.add_argument("-f" "--feature_file", action="store", dest="feature_file", default=None, help="Path to the feature legend. Defaults to a file feature_legend.csv in the current directory.")
    parser.add_argument("-r" "--rna_dge_file", action="store", dest="rna_dge_file", default='-', help="Specify RNA DGE matrix")
    parser.add_argument("--interact", action='store_true', help="Use flag if there are interaction reads in the data.")
    parser.add_argument("--stride", action="store", default=500000, help="Define after how many analysed reads the outputs should be printed during barcode filtering.")

    ## Parse input
    args = parser.parse_args()
    
    # initialize filter object
    filter = xDbit_filter(args)
    
    # run filtering
    filter.run_filtering()
    
    # Save files and print to log file
    filter.save_files()
    filter.print_to_log()
    filter.save_summary_and_qc()
    
    # closing statement
    filter.finish()
