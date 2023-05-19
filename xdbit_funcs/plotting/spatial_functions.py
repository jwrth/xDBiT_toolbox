from lib2to3.pytree import convert
from matplotlib import pyplot as plt
from matplotlib import patches, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
from ..calculations._calc import dist_points, minDistance
from sklearn.preprocessing import MinMaxScaler
from ..tools import extract_groups, check_raw, create_color_dict, get_nrows_maxcols, get_crange
from ..calculations import smooth_fit
from ..readwrite import save_and_show_figure
from ..images import set_histogram
from tqdm import tqdm
import warnings
from warnings import warn
from scipy.stats import zscore
from typing import Literal
from anndata import AnnData
import matplotlib.ticker as ticker

# ignore future warnings (suppresses annoying pandas warning)
warnings.simplefilter(action='ignore', category=FutureWarning)

def radial_expression(adata, coordinates, radius, min_dist=50, keys=None, mode='genes', raw=False, log_transform=False, show=True, axis=None, fig=None,
                      dpi_display=80, percent=False, legend_x=1.0, scaling=False):

    if raw:
        adata_var_names = adata.raw.var_names
        adata_X = adata.raw.X
    else:
        adata_var_names = adata.var_names
        adata_X = adata.X

    if axis is None:
        if isinstance(keys, str):
            n_plots = 1
            keys = [keys]
        elif isinstance(keys, list):
            n_plots = len(keys)
        else:
            print('Keys have unknown type')

    else:
        n_plots = 1
        if isinstance(keys, str):
            keys = [keys]

    coordinate_len = len(np.array(coordinates).shape)
    if coordinate_len == 1:
        dists = adata.obs.apply(lambda x: dist_points(
            [x.um_col, x.um_row], [coordinates[0], coordinates[1]]), axis=1)

    elif coordinate_len > 1:
        dists = adata.obs.apply(lambda x: np.min([minDistance(coordinates[i], coordinates[i + 1], 
            [x.um_col, x.um_row]) for i in range(len(coordinates) - 1)]), axis=1)

    else:
        print("Not able to determine shape of input coordinates.")
        return

    # filter for spots inside radius
    mask = list(dists < radius)
    dists_inradius = dists[mask]

    if mode == 'genes' and keys is None:
        print("Keys are missing.")
        return

    elif mode == 'genes':
        # extract the expression of key gene
        X_inradius = adata_X[mask, :]
        gene_expr = np.array(
            [X_inradius[:, adata_var_names.get_loc(key)] for key in keys]).T

        if log_transform:
            # transform log back to absolute values
            gene_expr = np.array([math.e**elem for elem in gene_expr])
            ylabel = "Expression"
        else:
            ylabel = "Log expression"

        # Make dataframe
        if len(keys) > 1 and scaling:
            scaler = MinMaxScaler()
            df = pd.DataFrame(scaler.fit_transform(gene_expr),
                              columns=keys, index=dists_inradius)
        else:
            df = pd.DataFrame(gene_expr, columns=keys, index=dists_inradius)

    elif mode == 'celltypes':
        if keys is None:
            df = adata.uns['stereoscope'].loc[mask, :].copy()
        else:
            df = adata.uns['stereoscope'].loc[mask, keys].copy()
        df.index = dists_inradius
        df = df.groupby(df.index).mean()

        if percent:
            df *= 100

    else:
        print("Mode not specified.")
        return

    if axis is None:
        fig, ax = plt.subplots(1, 1, dpi=dpi_display)
    else:
        ax = axis

    # Plotting
    if mode == 'genes':
        sns.lineplot(
            data=df,
            err_style="band", ci=95,
            ax=ax
        )
        ax.legend(bbox_to_anchor=(legend_x, 1), loc='upper left')
        ax.set_xlabel("Distance [µm]")
        ax.set_ylabel(ylabel)

    elif mode == 'celltypes':
        ax.stackplot(df.index, df.T, labels=df.columns)
        ax.legend(bbox_to_anchor=(legend_x, 1), loc='upper left')

        ax.set_xlim(left=min_dist)
        ax.set_xlabel("Distance [µm]")
        if percent:
            ax.set_ylabel("%")
        else:
            ax.set_ylabel("Fraction")

    if show:
        return plt.show()
    else:
        return fig, ax


def spatial_clusters(adata, save=False, savepath="figures/clusters.png", cluster_key='leiden', dpi=300, axis=None, show=True):

    if axis is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    else:
        ax = axis

    sns.scatterplot('um_col', 'um_row', data=adata.obs,
                    hue=cluster_key, marker='s', s=18, ax=ax, edgecolor="none")
    ax.set_xlabel('µm')
    ax.set_ylabel('µm')
    ax.set_title('Clusters')
    ax.set_facecolor('k')
    ax.invert_yaxis()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    ax.grid(False)

    if save:
        plt.savefig(savepath, dpi=dpi, bbox_inches='tight')

    if show:
        return plt.show()
    else:
        return ax
    
def expression_along_observation_value(*args, **kwargs):
    warn('This function is deprecated. Use instead `db.pl.expr_along_obs_val()`', DeprecationWarning, stacklevel=2)
    expr_along_obs_val(*args, **kwargs)

def expr_along_obs_val(adata: AnnData, 
                       keys: str, 
                       x_category: str, 
                       groupby: str, 
                       splitby: str = None, 
                       hue: str = None, 
                       method: Literal["lowess", "loess"] = 'loess',
                       stderr: bool = False, 
                       range_min=None, 
                       range_max=None, 
                       cmap="tab10",
                       linewidth=8,
                       extra_cats=None,
                       normalize=False,
                       nsteps=100,
                       show_progress=False,
                       use_raw=False,
                       max_cols=4,
                       xlabel=None,ylabel=None,
                       vline=None, hline=None, vlinewidth=4,
                       #values_into_title=None, title_suffix='', 
                       custom_titles=None,
                       legend_fontsize=24, 
                       plot_legend=True,
                       xlabel_fontsize=28, ylabel_fontsize=28, title_fontsize=20, tick_fontsize=24,
                       figsize=(8,6),
                       savepath=None, save_only=False, show=True, axis=None, return_data=False, fig=None,
                       dpi_save=300,
                       smooth=True, 
                       **kwargs
                       ):

    '''
    Plot the expression of a gene as a function of an observation value (e.g. the automatic expression histology value 
    given by SpatialDE).
    Grouping by one other observation is possible.

    Future ideas:
        -   Include the possibility of plotting categorical values like leiden plot (stacked line plot as done with 
            the radial expression and different cell types)

    '''

    # check type of input
    if isinstance(keys, dict):
        if custom_titles is not None:
            print("Attention: `custom_titles` was not None and `keys` was dictionary. Titles were retrieved from dictionary.")
        custom_titles = list(keys.keys())
        keys = list(keys.values())

    # make inputs to lists
    keys = [keys] if isinstance(keys, str) else list(keys)
    
    if hue is not None:
        hue_cats = list(adata.obs[hue].unique())
        cmap_colors = plt.get_cmap(cmap)
        color_dict = {a: cmap_colors(i) for i, a in enumerate(hue_cats)}

        if extra_cats is None:
            extra_cats = [hue]
        else:
            extra_cats.append(hue)

        

    #if show:
    if not return_data:
        # prepare plotting
        if axis is None:
            n_plots, n_rows, max_cols = get_nrows_maxcols(keys, max_cols)
            fig, axs = plt.subplots(n_rows,max_cols, figsize=(figsize[0]*max_cols, figsize[1]*n_rows))

        else:            
            axs = axis
            #fig = None
            n_plots = 1
            show = False # otherwise plotting into given axes wouldn't work
            
        if n_plots > 1:
            axs = axs.ravel()
        else:
            axs = [axs]

    data_collection = {}
    for i, key in (enumerate(tqdm(keys)) if show_progress else enumerate(keys)):
        # check if the keys are also grouped
        keys_grouped = isinstance(key, list)
        # select data per group
        groups = adata.obs[groupby].unique()
        added_to_legend = []

        # check if plotting raw data
        X, var, var_names = check_raw(adata, use_raw=use_raw)
        
        group_collection = {}
        for group in groups:
            #partial = extract_groups(adata, groupby=groupby, groups=group)
            group_mask = adata.obs[groupby] == group
            group_obs = adata.obs.loc[group_mask, :].copy()
            if hue is not None:
                _hue = adata.obs.loc[group_mask, hue][0]

            # hue_data = adata.obs.loc[group_mask, hue].copy()
            # print(hue_data)

            # select only group values from matrix
            group_X = X[group_mask, :]

            if splitby is None:
                # select x value
                x = group_obs.loc[:, x_category].values

                if keys_grouped:
                    # extract expression values of all keys in the group
                    idx = var.index.get_indexer(key)
                    dd = pd.DataFrame(group_X[:, idx], index=x)

                    if normalize:
                        #dd = dd.apply(minmax_scale, axis=0)
                        dd = dd.apply(zscore, axis=0)

                    dd = dd.reset_index().melt(id_vars="index") # reshape to get long list of x values
                    x = dd["index"].values
                    y = dd["value"].values
                
                elif key in var_names:
                    # extract expression values as y
                    idx = var.index.get_loc(key)
                    y = group_X[:, idx].copy()

                    if normalize:
                        #y = minmax_scale(y)
                        y = zscore(y)

                elif key in group_obs.columns:
                    y = group_obs.loc[:, key].values.copy()
                else:
                    print("Key '{}' not found.".format(key))                    

                if smooth:
                    # do smooth fitting
                    df = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps, method=method, stderr=stderr, **kwargs)
                else:
                    # set up dataframe without smooth fitting
                    df = pd.DataFrame({"x": x, "y_pred": y})

                if extra_cats is not None:
                    df = df.join(adata.obs.loc[group_mask, extra_cats].reset_index(drop=True))

            else:
                splits = group_obs[splitby].unique()
                df_collection = {}

                # get min and max values for x values
                x = group_obs[x_category].values
                range_min = x.min()
                range_max = x.max()
                for split in splits:
                    # extract x values                                            
                    split_mask = group_obs[splitby] == split
                    x = group_obs.loc[split_mask, x_category].values

                    # extract expression values as y
                    idx = var.index.get_loc(key)
                    y = group_X[split_mask, idx].copy()

                    # do smooth fitting
                    if smooth:
                        df_split = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps, method=method, stderr=stderr, **kwargs)
                    else:
                        # set up dataframe without smooth fitting
                        df = pd.DataFrame({"x": x, "y_pred": y})

                    # collect data
                    df_collection[split] = df_split
                
                df_collection = pd.concat(df_collection)
                
                # calculate mean and std
                df = df_collection[['x', 'y_pred']].groupby('x').mean()
                df['std'] = df_collection[['x', 'y_pred']].groupby('x').std()
                df['conf_lower'] = [a-b for a,b in zip(df['y_pred'], df['std'])]
                df['conf_upper'] = [a+b for a,b in zip(df['y_pred'], df['std'])]
                df.reset_index(inplace=True)
            if return_data:
                group_collection[group] = df
            else:
                # sort by x-value
                df.sort_values('x', inplace=True)

                # plotting
                cols = df.columns
                if 'conf_lower' in cols and 'conf_upper' in cols:
                    axs[i].fill_between(df['x'], 
                                    df['conf_lower'],
                                    df['conf_upper'],
                                    alpha = 0.2,
                                    color = 'grey')

                # determine label variable
                if hue is not None:
                    label = _hue if _hue not in added_to_legend else ""
                    color = color_dict[_hue]
                else:
                    label = group
                    color = None
                
                axs[i].plot(df['x'], 
                    df['y_pred'], 
                    label=label, 
                    color=color,
                    linewidth=linewidth)

                if hue is not None and _hue not in added_to_legend:
                    added_to_legend.append(_hue)

        # optionally add vertical or horizontal lines to plot
        if vline is not None:
            if isinstance(vline, dict):
                linecolors = list(vline.keys())
                vline = list(vline.values())
            else:
                vline = [vline] if isinstance(vline, int) or isinstance(vline, float) else list(vline)
                linecolors = ['k'] * len(vline)
            
            for c, v in zip(linecolors, vline):
                axs[i].axvline(x=v, ymin=0, ymax=1, c=c, linewidth=vlinewidth, linestyle='dashed')

        if hline is not None:
            if isinstance(hline, dict):
                linecolors = list(hline.keys())
                hline = list(hline.values())
            else:
                hline = [hline] if isinstance(hline, int) or isinstance(hline, float) else list(hline)
                linecolors = ['k'] * len(hline)
            
            for c, h in zip(linecolors, hline):
                axs[i].axhline(y=h, xmin=0, xmax=1, c=c, linewidth=4, linestyle='dashed')

        if not return_data:
            if xlabel is None:
                xlabel = x_category
            if ylabel is None:
                ylabel = "Gene expression"

            axs[i].set_xlabel(xlabel, fontsize=xlabel_fontsize)
            axs[i].set_ylabel(ylabel, fontsize=ylabel_fontsize)
            axs[i].tick_params(axis='both', which='major', labelsize=tick_fontsize)
            axs[i].set_xlim(0, 1)
            #axs[i].xaxis.set_major_locator(ticker.FixedLocator([0.1, 0.9]))

            if custom_titles is None:
                axs[i].set_title(key, fontsize=title_fontsize)
            else:
                assert len(custom_titles) == len(keys), "List of title values has not the same length as list of keys."
                axs[i].set_title(str(custom_titles[i]), fontsize=title_fontsize)

            if plot_legend:
                axs[i].legend(fontsize=legend_fontsize, 
                #bbox_to_anchor=(legend_x, 1), 
                loc='best'
                )
            else:
                axs[i].legend().remove()

            # if values_into_title is None:
            #     axs[i].set_title("{}{}".format(key, title_suffix), fontsize=title_fontsize)
            # else:
            #     assert len(values_into_title) == len(keys), "List of title values has not the same length as list of keys."
            #     axs[i].set_title("{}{}{}".format(key, title_suffix, round(values_into_title[i], 5)), fontsize=title_fontsize)

        if return_data:
            # collect data
            group_collection = pd.concat(group_collection)
            data_collection[key] = group_collection


    
    if return_data:
        # close plot
        plt.close()

        # return data
        data_collection = pd.concat(data_collection)
        data_collection.index.names = ['key', groupby, None]
        return data_collection

    else: 
        if n_plots > 1:

            # check if there are empty plots remaining
            while i < n_rows * max_cols - 1:
                i+=1
                # remove empty plots
                axs[i].set_axis_off()
        if show:
            #fig.tight_layout()
            save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save, tight=True)
        else:
            return fig, axs


def linear_expression_grouped(adata, splitby=None, keys=None, x_category=None, groupby=None, max_cols=4, n_bins=20, use_raw=False,
x_limit_labels=None, savepath=None, save_only=False, show=True, dpi_save=300, **kwargs):
    
    n_plots, n_rows, max_cols = get_nrows_maxcols(keys, max_cols)
    
    fig, axs = plt.subplots(n_rows,max_cols, figsize=(8*max_cols, 6*n_rows))

    if n_plots > 1:
        axs = axs.ravel()

    for i, key in enumerate(keys):
        expression_along_observation_value(adata, splitby=splitby, key=key, x_category=x_category, groupby=groupby, n_bins=n_bins, use_raw=use_raw,
        x_limit_labels=x_limit_labels, savepath=savepath, save_only=save_only, show=show, dpi_save=dpi_save, ax=axs[i])

    plt.tight_layout()
    plt.show()
