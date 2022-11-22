#!/usr/bin/env python

'''
Script to plot the mean genes per cell for the subsampled matrices created with subsample_bam.sh and create_dge_from_subsampled_bams.sh

'''

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys
from pathlib import Path
import argparse

# # get inputs
# full_mtx_path = Path(sys.argv[1])
# subsampled_path = Path(sys.argv[2])

# Parse arguments
print("Parse arguments")
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--full_mtx_path", help="Path to full DGE matrix.")
parser.add_argument("-s", "--subsampled_path", help="Path to subsampled DGe matrices.")
parser.add_argument("-n", "--name", default="seq", help="Name to add as prefix to output file.")
args = parser.parse_args()

full_mtx_path = Path(args.full_mtx_path)
subsampled_path = Path(args.subsampled_path)
name = args.name

# create output file
output_file = subsampled_path / "{}_saturation.png".format(name)

# load full matrix and subsampled matrices
print("Find files...", flush=True)
adata_100 = sc.read_text(os.path.join(full_mtx_path, "DGE_matrix_with_introns_min100.txt.gz")).transpose()
#adata_100 = sc.read_text("out/DGE_matrix_with_introns_min100.txt.gz").transpose()

files = sorted(glob.glob(os.path.join(subsampled_path, 'DGE_matrix_with_introns_*')))
#files = glob.glob(os.path.join('out/subsampling/DGE_matrix_with_introns_*'))
print("Import files...", flush=True)
adatas = [sc.read_text(elem).transpose() for elem in files]

# append full adata set to subsampled sets
adatas.append(adata_100)

# give fractions in order of files
fractions = [int(elem.split("_")[-1].split("p")[0]) / 100 for elem in files] + [1]
#fractions = [0.1, 0.25, 0.5, 0.75, 1]

# sorting of adatas by fractions
adatas_sorted = [x for _,x in sorted(zip(fractions, adatas))]
fractions = sorted(fractions)

# filter all adata sets
print("Filter adata objects...", flush=True)
for adata in adatas:
    #sc.plotting.highest_expr_genes(adata, n_top = 20)
    sc.preprocessing.filter_cells(adata, min_genes=100)
    sc.preprocessing.filter_genes(adata, min_cells=3)

for adata in adatas:
    sc.preprocessing.filter_cells(adata, min_genes=100)
    sc.preprocessing.filter_genes(adata, min_cells=3)

# add parameters to adata sets
print("Add mito percentage and counts to adata object...", flush=True)
for adata in adatas:
    mito_genes = adata.var_names.str.contains('MT-', case = False)

    # sums up all mito genes per cell and normalizes to the sum of all genes per cell
    adata.obs['percent_mito'] = np.array(adata[:, mito_genes].X.sum(axis=1, dtype=int)) / np.array(adata.X.sum(axis=1, dtype=int)) * 100 # A1 to make dense array out of sparse matrix which was derived from X funtion

    # add the total counts per cell as observation
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1))

    # plot the total counts and mito counts in violin plot
    #sc.plotting.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
    #                  jitter=0.4, multi_panel=True,
    #                  save="_genes_counts_mito.png"
    #                  )

# Calculate mean of genes and UMIs per cell
print("Calculate mean of genes and UMIs...", flush=True)
genes_means = [np.mean(adata.obs['n_genes']) for adata in adatas]
umi_means = [np.mean(adata.obs['n_counts']) for adata in adatas]

# Plotting
print("Plotting...", flush=True)
fig, ax = plt.subplots(1, 2, figsize=(15, 5))

plt.suptitle('Analysis of sequencing saturation by subsampling of reads')

ax[0].plot(fractions, umi_means)
ax[0].set_xlabel('Fractions of reads')
ax[0].set_ylabel('Mean UMIs per cell')

ax[1].plot(fractions, genes_means)
ax[1].set_xlabel('Fractions of reads')
ax[1].set_ylabel('Mean genes per cell')

plt.savefig(output_file, bbox_inches='tight', dpi=100)

print("Plot saved in {}".format(output_file), flush=True)