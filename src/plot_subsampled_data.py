#!/usr/bin/env python

'''
Script to plot the mean genes per cell for the subsampled matrices created with subsample_bam.sh and create_dge_from_subsampled_bams.sh

'''

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# get inputs
full_mtx_path = sys.argv[1]
subsampled_path = sys.argv[2]

# load full matrix and subsampled matrices
adata_100 = sc.read_text(os.path.join(full_mtx_path, "DGE_matrix_with_introns_min100.txt.gz")).transpose()
#adata_100 = sc.read_text("out/DGE_matrix_with_introns_min100.txt.gz").transpose()

files = glob.glob(os.path.join(subsampled_path, 'DGE_matrix_with_introns_*'))
#files = glob.glob(os.path.join('out/subsampling/DGE_matrix_with_introns_*'))
adatas = [sc.read_text(elem).transpose() for elem in files]

# append full adata set to subsampled sets
adatas.append(adata_100)

# give fractions in order of files
fractions = [int(elem.split("_")[-1].split("p")[0]) / 100 for elem in files]
#fractions = [0.1, 0.25, 0.5, 0.75, 1]

# sorting of adatas by fractions
adatas_sorted = [x for _,x in sorted(zip(fractions, adatas))]
fractions = sorted(fractions)

# filter all adata sets
for adata in adatas:
    #sc.plotting.highest_expr_genes(adata, n_top = 20)
    sc.preprocessing.filter_cells(adata, min_genes=100)
    sc.preprocessing.filter_genes(adata, min_cells=3)

for adata in adatas:
    sc.preprocessing.filter_cells(adata, min_genes=100)
    sc.preprocessing.filter_genes(adata, min_cells=3)

# add parameters to adata sets
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
genes_means = [np.mean(adata.obs['n_genes']) for adata in adatas]
umi_means = [np.mean(adata.obs['n_counts']) for adata in adatas]

# Plotting
fig, ax = plt.subplots(1, 2, figsize=(15, 5))

plt.suptitle('Analysis of sequencing saturation by subsampling of reads')

ax[0].plot(fractions, umi_means)
ax[0].set_xlabel('Fractions of reads')
ax[0].set_ylabel('Mean UMIs per cell')

ax[1].plot(fractions, genes_means)
ax[1].set_xlabel('Fractions of reads')
ax[1].set_ylabel('Mean genes per cell')

plt.savefig('./seq_saturation.png', bbox_inches='tight', dpi=100)
