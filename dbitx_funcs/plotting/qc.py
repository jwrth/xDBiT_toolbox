import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
from ..tools import get_nrows_maxcols
from ..readwrite import save_and_show_figure

def qc_alignment(rnametrics_files, data_names=None, max_cols=4,
    savepath=None, save_only=False, show=True, dpi_save=300):

    '''
    Plot Alignment metrics from STAR alignment run.
    '''

    # make inputs to list
    rnametrics_files = [rnametrics_files] if isinstance(rnametrics_files, str) else list(rnametrics_files)

    if data_names is None:
        data_names = ["Sample\_{}".format(i) for i in range(len(rnametrics_files))]

    # read metrics file
    metrics_dict = {}

    for n, f in zip(data_names, rnametrics_files):
        metrics = []
        for line in open(f, "r"):
            if not line.startswith("#"):
                metrics.append(line.strip())
        metrics_dict[n] = metrics

    # extract metric names and values and collect in dataframe
    names_dict = {key: metrics_dict[key][1].split("\t")[:-3] for key in metrics_dict}
    values_dict = {key: metrics_dict[key][2].split("\t") for key in metrics_dict}

    df_dict = {key: pd.DataFrame({'metric': names_dict[key], 'value': values_dict[key]}) for key in names_dict}
    df = pd.concat(df_dict)

    # Plotting
    n_plots, n_rows, n_cols = get_nrows_maxcols(rnametrics_files, max_cols=max_cols)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_cols, 8*n_rows))
    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, name in enumerate(df.index.unique(level=0)):
        pct_values = df.xs(name).value[15:20]
        pct_labels = df.xs(name).metric[15:20]
        
        axs[i].pie(pct_values, labels=pct_labels, autopct='%1.1f%%')
        axs[i].set_title(r"$\bf{" + name + "}$")

    plt.suptitle("RNA alignment metrics", size=16)

    fig.tight_layout()
    save_and_show_figure(savepath=savepath, save_only=save_only, show=show, dpi_save=dpi_save)