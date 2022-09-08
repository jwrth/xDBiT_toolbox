import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
from ..tools import get_nrows_maxcols, read_rna_metrics
from ..readwrite import save_and_show_figure



def rna_metrics_pie(metrics, max_cols=4,
    savepath=None, save_only=False, dpi_save=300,
    categories = ['PCT_RIBOSOMAL_BASES', 'PCT_CODING_BASES', 'PCT_UTR_BASES',
            'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES'],
    suptitle_fontsize=32, title_fontsize=28, label_fontsize=18
    ):

    '''
    Plot Alignment metrics from STAR alignment run. Input: Output dataframe from `db.tl.read_rna_metrics`
    '''

    # extract categories to be plotted
    metrics = metrics[metrics.metric.isin(categories)]

    # extract datanames
    data_names = list(metrics.index.unique(level=0))

    # Plotting
    n_plots, n_rows, n_cols = get_nrows_maxcols(data_names, max_cols=max_cols)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(12*n_cols, 8*n_rows))
    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, name in enumerate(metrics.index.unique(level=0)):
        pct_values = metrics.xs(name).value
        pct_labels = metrics.xs(name).metric
        
        axs[i].pie(pct_values, labels=pct_labels, autopct='%1.1f%%', 
                textprops={'fontsize': label_fontsize})
        axs[i].set_title(r"$\bf{" + name + "}$", fontsize=title_fontsize)

    plt.suptitle("RNA alignment metrics", fontsize=suptitle_fontsize, y=1.005)

    fig.tight_layout()
    save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)

