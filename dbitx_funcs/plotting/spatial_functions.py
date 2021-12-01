from matplotlib import pyplot as plt
from matplotlib import patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
import math
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
from ..calculations._calc import dist_points, minDistance
from sklearn.preprocessing import MinMaxScaler
from ..tools import extract_groups, check_raw, create_color_dict, get_nrows_maxcols, smooth_fit, get_crange
from ..readwrite import save_and_show_figure
from ..images import set_histogram
from tqdm import tqdm


def spatial(adata, keys, groupby=None, group=None, max_cols=4, pd_dataframe=None,
            header=None, headersize=18, header_names=None, raw=False, percent=False,
            dpi_save=300,
            obsm_key = 'spatial',
            spot_size=50, spot_type='s',
            axis=None, fig=None, show=True,
            patch_style=None, patch_xy=(1000, 1000), patch_radius=1000, patch_color='r',
            xlim=None, ylim=None, oversize=1, dpi_display=80, figsize=(8.2,6),
            #image=None, image_metadata=None, 
            image_key=None, lowres=True, histogram_setting=None,
            alpha=1, palette="tab10", color_dict=None,
            header_x=0.5, header_y=0.98, header_fontsize=12,
            save_only=False, savepath=None, save_background=None, crange=None,
            verbose=True):

    if isinstance(pd_dataframe, pd.DataFrame):
        assert len(pd_dataframe) == len(adata), "Given dataframe ({}) and anndata object ({}) do not have same length.".format(pd_dataframe.shape, adata.shape)
        data_in_dataframe = True
    else:
        data_in_dataframe = False

    # check if certain group is selected
    if groupby is not None:
        if group is not None:
            adata, mask = extract_groups(adata, groupby, group, extract_uns=False, return_mask=True)
            if data_in_dataframe:
                # slice pd_dataframe
                pd_dataframe = pd_dataframe[mask].copy()

        else:
            print("Subset variable `group` missing.")
            #return
    else:
        if isinstance(pd_dataframe, pd.DataFrame):
            assert len(pd_dataframe) == len(adata), "Given dataframe ({}) and anndata object ({}) do not have same length.".format(pd_dataframe.shape, adata.shape)

            data_in_dataframe=True
        data_in_dataframe=False

    # convert parameters to lists
    keys = [keys] if isinstance(keys, str) else list(keys)
    if xlim is not None:
        xlim = [xlim] if isinstance(xlim, str) else list(xlim)
    if xlim is not None:
        ylim = [ylim] if isinstance(ylim, str) else list(ylim)

    # check if plotting raw data
    adata_X, adata_var, adata_var_names = check_raw(adata, use_raw=raw)


    if header_names is not None:
        assert len(header_names) == len(keys)

    # get image data
    if image_key is not None:
        # check if image_key gives unique match
        image_key_list = adata.uns['spatial'].keys()
        image_key_matches = [k for k in image_key_list if image_key in k]
        if group is not None:
            # make sure that group name is also in image_key
            image_key_matches = [k for k in image_key_matches if group in k]

        if len(image_key_matches) == 0:
            print("No match for img_key [{}] found in list of image_keys: ".format(image_key, image_key_list))
            return
        elif len(image_key_matches) > 1:
            print("More than one possible match for img_key [{}] found: {}".format(image_key, image_key_matches))
            return
        else:
            image_key = image_key_matches[0]
        
        image_metadata = adata.uns['spatial'][image_key]['scalefactors']

        pixel_per_um = image_metadata["pixel_per_um"]
        if lowres:
            if 'lowres' in adata.uns['spatial'][image_key]['images'].keys():
                image = adata.uns['spatial'][image_key]['images']['lowres']
                scale_factor = image_metadata['tissue_lowres_scalef']
            else:
                print('`hires` image displayed because no `lowres` images was found.')
                image = adata.uns['spatial'][image_key]['images']['hires']
                scale_factor = image_metadata['tissue_hires_scalef']
        else:
            image = adata.uns['spatial'][image_key]['images']['hires']
            scale_factor = image_metadata['tissue_hires_scalef']

        if histogram_setting is not None:
            bit_type = np.uint8 if image.max() < 256 else np.uint16
            image = set_histogram(image, lower=histogram_setting[0], upper=histogram_setting[1], bit_type=bit_type)
    else:
        image = None
        scale_factor = 1
        pixel_per_um = 1
        image_metadata = None
        # search for pixel_per_um scalefactor
        if 'spatial' in adata.uns.keys():
            # get first entry of dictionary as image_metadata
            first_entry = next(iter(adata.uns['spatial'].values()))
            
            # extract image metadata if possible
            if 'scalefactors' in first_entry:
                image_metadata = first_entry['scalefactors']
                pixel_per_um = image_metadata["pixel_per_um"]
                scale_factor = image_metadata['tissue_hires_scalef']                  
            else:
                print("pixel_per_um scalefactor not found. Plotted pixel coordinates instead.")
        else:
            print("No key `spatial` in adata.uns. Therefore pixel_per_um scalefactor could not be found. Plotted pixel coordinates instead.") if verbose else None

    # extract x and y pixel coordinates and convert to micrometer
    x_coord = adata.obsm[obsm_key][:, 0] / pixel_per_um
    y_coord = adata.obsm[obsm_key][:, 1] / pixel_per_um

    # shift coordinates that they start at (0,0)
    x_offset = x_coord.min()
    y_offset = y_coord.min()
    x_coord -= x_offset
    y_coord -= y_offset

    if xlim is None:
        xmin = x_coord.min()
        xmax = x_coord.max()

        xlim = (xmin - spot_size, xmax + spot_size)
    else:
        xlim[0] -= spot_size
        xlim[1] += spot_size

    if ylim is None:
        ymin = y_coord.min()
        ymax = y_coord.max()       

        ylim = (ymin - spot_size, ymax + spot_size)
    else:
        ylim[0] -= spot_size
        ylim[1] += spot_size


    if axis is None:
        n_plots = len(keys)

    else:
        n_plots = 1

    if max_cols is None:
        max_cols = n_plots
        n_rows = 1
    else:
        if n_plots > max_cols:
            n_rows = math.ceil(n_plots / max_cols)
        else:
            n_rows = 1
            max_cols = n_plots

    # Plotting section
    # summarize data for seaborn plotting
    data = adata.obs.copy()
    data['x_coord'] = x_coord
    data['y_coord'] = y_coord

    # check plotting parameters
    if axis is None:
        fig, axs = plt.subplots(n_rows, max_cols, figsize=(
            figsize[0] * max_cols, figsize[1] * n_rows), dpi=dpi_display)
    else:
        axs = axis
    if n_plots > 1:
        axs = axs.ravel()
    for i, key in enumerate(keys):
        if n_plots > 1:
            ax = axs[i]
        else:
            ax = axs

        
        # locate gene in matrix and extract values
        if key is not None:
            if data_in_dataframe:
                color = pd_dataframe[key].copy()
                if color.dtype.name == 'category':
                    categorical = True
                else:
                    categorical = False
            else:
                if key in adata_var_names:
                    idx = adata_var.index.get_loc(key)
                    color = adata_X[:, idx].copy()
                    categorical = False
                    
                elif key in adata.obs.columns:
                    if adata.obs[key].dtype.name == 'category':
                        categorical = True
                    else:
                        color = adata.obs[key].values
                        categorical = False
                else:
                    print("Key '{}' not found.".format(key))
                    return

            if percent:
                color *= 100
        
        # set axis
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        if image_metadata is None:
            ax.set_xlabel('pixels')
            ax.set_ylabel('pixels')
        else:
            ax.set_xlabel('µm')
            ax.set_ylabel('µm')

        if header_names is None:
            ax.set_title(key, fontsize=header_fontsize, fontweight='bold')
        else:
            plot_name = header_names[i]
            ax.set_title(plot_name, fontsize=header_fontsize, fontweight='bold')
        
        ax.invert_yaxis()
        ax.grid(False)
        ax.set_aspect(1)
        ax.set_facecolor('k')

        # calculate marker size
        pixels_per_unit = ax.transData.transform(
            [(0, 1), (1, 0)]) - ax.transData.transform((0, 0))
        x_ppu = pixels_per_unit[1, 0]
        y_ppu = pixels_per_unit[0, 1]
        pxs = y_ppu * spot_size * oversize
        size = (72. / fig.dpi * pxs)**2

        if color_dict is None:
            color_dict = palette

        # plot
        if image is None and key is None:
            print("Nothing to plot.")
            return
        else:
            if key is not None:
                # plot transcriptomic data
                if not categorical:
                    s = ax.scatter(x_coord, y_coord, c=color, marker=spot_type,
                                   s=size, alpha=alpha, edgecolors=None)
                else:
                    sns.scatterplot(x='x_coord', y='y_coord', data=data,
                                    hue=key, marker=spot_type, s=size*1.6, ax=ax, edgecolor="none", palette=color_dict, alpha=alpha)
            if image is not None:
                ax.imshow(image, extent=(
                    -0.5 - x_offset, image.shape[1] / pixel_per_um / scale_factor - 0.5 - x_offset, 
                    image.shape[0] / pixel_per_um / scale_factor - 0.5 - y_offset, -0.5 - y_offset
                ), origin='upper', cmap='gray')

        if patch_style is not None:
            pats = []
            if patch_style == "spot":
                # Create a Rectangle patch
                pat = patches.Rectangle((patch_xy[0] - 25, patch_xy[1] - 25), 50, 50, linewidth=2, fill=True,
                                        edgecolor=patch_color, facecolor=patch_color)
                pats.append(pat)
                # Create a circular patch
                pat = patches.Circle(patch_xy, patch_radius, linewidth=2,
                                     edgecolor=patch_color, facecolor='none')
                pats.append(pat)
            if patch_style == "line":
                # Create a line patch
                pat = patches.Polygon(patch_xy, closed=False, linewidth=2,
                                      edgecolor=patch_color, facecolor='none')
                pats.append(pat)
            if patch_style == "polygon":
                # Create a polygon patch
                pat = patches.Polygon(patch_xy, closed=True, linewidth=2,
                                      edgecolor=patch_color, facecolor='none')
                pats.append(pat)

            # Add the patch to the Axes
            for pat in pats:
                ax.add_patch(pat)
        if key is not None:
            if categorical:
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
            else:
                divider = make_axes_locatable(ax)
                #cax = divider.append_axes("right", size="0%", pad=1)
                clb = fig.colorbar(s, ax=ax, shrink=0.8)
                #clb = plt.colorbar(s, ax=cax, shrink=1)
                #cax.axis('off')
                if crange is not None:
                    clb.mappable.set_clim(crange[0], crange[1])
                if percent:
                    clb.ax.set_title('%')

        
    if n_plots > 1:
        # check if there are empty plots remaining
        while i < n_rows * max_cols - 1:
            i+=1
            # remove empty plots
            axs[i].set_axis_off()

    if header is not None:
        plt.suptitle(header, fontsize=headersize, x=header_x, y=header_y)
    if savepath is not None:
        fig.tight_layout()
        print("Saving figure to file " + savepath)
        plt.savefig(savepath, dpi=dpi_save,
                    facecolor=save_background, bbox_inches='tight')
        print("Saved.")

    if save_only:
        plt.close(fig)
    if show:
        fig.tight_layout()
        return plt.show()
    else:
        return fig, ax


def spatial_grouped(adata, keys, groupby='well_name', groups=None, raw=False, max_cols=4, 
    spot_size=50, prefix_groups='', palette="tab10", groupheader_fontsize=20,
    savepath=None, dpi_save=300, show=True, save_only=False, pd_dataframe=None, normalize_crange_not_for=[], dpi_display=80, header_names=None,
    xlim=None, ylim=None,
    **kwargs):
    
    '''
    Creates spatial plot that is grouped by a certain parameter in `adata.obs`.
    The resulting plots has the groups as rows and columns as genes.
    '''
    
    # check keys and groups
    keys = [keys] if isinstance(keys, str) else list(keys)
    multikeys = False
    if len(keys) > 1:
        multikeys = True

    if groups is None:
        groups = list(adata.obs[groupby].unique())
    else:
        groups = [groups] if isinstance(groups, str) else list(groups)

    if header_names is not None:
        assert len(header_names) == len(keys)

    # check if dataframe is given
    if isinstance(pd_dataframe, pd.DataFrame):
        data_in_dataframe=True
    else:
        data_in_dataframe=False

    # determine x and y limits
    if xlim is None:
        xlim = [adata.obs["um_col"].min(), adata.obs["um_col"].max()]
    if ylim is None:
        ylim = [adata.obs["um_row"].min(), adata.obs["um_row"].max()]

    # determine the color range for each key
    crange_per_key_dict = {key: get_crange(adata, groupby, groups, key, 
                use_raw=raw, data_in_dataframe=data_in_dataframe, pd_dataframe=pd_dataframe) if key not in normalize_crange_not_for else None for key in keys}
    
    if multikeys:
        n_rows = len(groups)
        max_cols = len(keys)
        n_plots = n_rows * max_cols
        fig, axs = plt.subplots(n_rows, max_cols, figsize=(7.6 * max_cols, 6 * n_rows), dpi=dpi_display)

        i = 0
        for row, group in enumerate(groups):
            for col, key in enumerate(keys):
                # counter
                i+=1
                
                if header_names is not None:
                    header_name = [header_names[col]]
                else:
                    header_name = None

                # determine x and y limits
                if xlim is None:    
                    xlim = [adata.obs["um_col"].min(), adata.obs["um_col"].max()]
                if ylim is None:
                    ylim = [adata.obs["um_row"].min(), adata.obs["um_row"].max()]

                # create color dictionary if key is categorical
                color_dict = create_color_dict(adata, key, palette)

                spatial(adata, key, raw=raw, groupby=groupby,
                        group=group, fig=fig, axis=axs[row, col], show=False,
                        xlim=xlim, ylim=ylim, spot_size=spot_size, crange=crange_per_key_dict[key], 
                        palette=palette, color_dict=color_dict, pd_dataframe=pd_dataframe, header_names=header_name, **kwargs)

        for ax, row in zip(axs[:, 0], groups):
            ax.annotate(prefix_groups + row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size=groupheader_fontsize, ha='right', va='center', weight='bold')

    else:
        n_plots = len(groups)
        if n_plots > max_cols:
            n_rows = math.ceil(n_plots / max_cols)
        else:
            n_rows = 1
            max_cols = n_plots

        fig, axs = plt.subplots(n_rows, max_cols, figsize=(7.6 * max_cols, 6 * n_rows), dpi=dpi_display)

        if n_plots > 1:
            axs = axs.ravel()
        else:
            axs = [axs]

        # for k in range(len(axs)):
        #     print(axs[k].transData.transform([(0, 1), (1, 0)]) - axs[k].transData.transform((0, 0)))
        # print(axs[1].transData.transform([(0, 1), (1, 0)]) - axs[1].transData.transform((0, 0)))
        # print(axs[2].transData.transform([(0, 1), (1, 0)]) - axs[2].transData.transform((0, 0)))
        # print(axs[3].transData.transform([(0, 1), (1, 0)]) - axs[3].transData.transform((0, 0)))

        
        for i, group in enumerate(groups):
            key = keys[0]

            # determine x and y limits
            if xlim is None:    
                xlim = [adata.obs["um_col"].min(), adata.obs["um_col"].max()]
            if ylim is None:
                ylim = [adata.obs["um_row"].min(), adata.obs["um_row"].max()]

            # create color dictionary if key is categorical
            color_dict = create_color_dict(adata, key, palette)

            spatial(adata, key, raw=raw, groupby=groupby,
                    group=group, fig=fig, axis=axs[i], show=False,
                    xlim=xlim, ylim=ylim, spot_size=spot_size, crange=crange_per_key_dict[key],
                    palette=palette, color_dict=color_dict, pd_dataframe=pd_dataframe, **kwargs)

            axs[i].set_title("{} - {}{}".format(key, prefix_groups, group))
            

    if n_plots > 1:
        # check if there are empty plots remaining
        while i < n_rows * max_cols - 1:
            i+=1
            # remove empty plots
            axs[i].set_axis_off()
    
    fig.tight_layout()
    if savepath is not None:
        print("Saving figure to file " + savepath)
        plt.savefig(savepath, dpi=dpi_save, bbox_inches='tight')
        print("Saved.")
    if save_only:
        plt.close(fig)
    elif show:
        return plt.show()
    else:
        return fig, ax




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

def expression_along_observation_value(adata, keys, x_category, groupby, splitby=None,
    range_min=-10, range_max=10, stepsize=0.01, show_progress=False,
    n_bins=20, use_raw=False,
    max_cols=4,
    x_limit_labels=None, xlabel='PC1 of normalized AEH values', ylabel='Gene expression', values_into_title=None, title_suffix='',
    #ax=None, 
    savepath=None, save_only=False, show=True, 
    dpi_save=300, return_fig_axis=False, 
    loess=True, **kwargs):

    '''
    Plot the expression of a gene as a function of an observation value (e.g. the automatic expression histology value 
    given by SpatialDE).
    Grouping by one other observation is possible.

    Future ideas:
        -   Include the possibility of plotting categorical values like leiden plot (stacked line plot as done with 
            the radial expression and different cell types)

    '''

    # make inputs to lists
    keys = [keys] if isinstance(keys, str) else list(keys)

    if show:
        n_plots, n_rows, max_cols = get_nrows_maxcols(keys, max_cols)
        fig, axs = plt.subplots(n_rows,max_cols, figsize=(8*max_cols, 6*n_rows))

        if n_plots > 1:
            axs = axs.ravel()
        else:
            axs = [axs]
    else:
        n_plots = 1

    data_collection = {}
    for i, key in (enumerate(tqdm(keys)) if show_progress else enumerate(keys)):
        # data = {}
        # for part in sorted(adata.obs[splitby].unique()):
        #     partial = extract_groups(adata, groupby=splitby, groups=part)

        #     # check if plotting raw data
        #     X, var, var_names = check_raw(partial, use_raw=use_raw)

        #     # extract groups to groupby
        #     group = partial.obs[groupby]

        #     idx = var.index.get_loc(key)
        #     y = X[:, idx].copy()
        #     x = partial.obs[x_category]

        #     # bin x and calculate mean for bins
        #     lower = x.min()
        #     upper = x.max()

        #     bins = np.linspace(lower, upper, n_bins)
        #     digitized = np.digitize(x, bins, right=True)

        #     partial_data = pd.DataFrame({
        #         groupby: group,
        #         'x': x,
        #         'bin': digitized,
        #         'expression': y
        #     })
            
        #     data[part] = partial_data

        # data = pd.concat(data)
        # #return data
        # data.index.names = [splitby, 'id']

        ## Plotting
        if loess:
            groups = adata.obs[groupby].unique()

            # check if plotting raw data
            X, var, var_names = check_raw(adata, use_raw=use_raw)
            
            group_collection = {}
            for group in groups:
                #partial = extract_groups(adata, groupby=groupby, groups=group)
                group_mask = adata.obs[groupby] == group
                group_obs = adata.obs.loc[group_mask, :].copy()

                # select only group values from matrix
                group_X = X[group_mask, :]

                if splitby is None:
                    # select x value
                    #x = partial.obs.loc[:, x_category].values
                    x = adata.obs.loc[group_mask, x_category].values

                    # extract expression values as y
                    idx = var.index.get_loc(key)
                    #y = X[:, idx].copy()
                    y = group_X[:, idx].copy()

                    # x = data.query('{} == "{}"'.format(groupby, group)).x.values
                    # y = data.query('{} == "{}"'.format(groupby, group)).expression.values
                    
                    # determine the x-values to fit on
                    #x_to_fit_on = np.arange(x.min(), x.max() + stepsize/10, stepsize)
                    x_to_fit_on = np.arange(range_min, range_max, stepsize)

                    # fit curve
                    df = smooth_fit(x, y, x_to_fit_on=x_to_fit_on)
                else:
                    #splits = partial.obs[splitby].unique()
                    #splits = adata.obs.loc[group_mask, splitby].unique()
                    splits = group_obs[splitby].unique()
                    df_collection = {}
                    for split in splits:
                        # extract x values
                        #obs_mask = partial.obs[splitby] == split
                        #x = partial.obs.loc[obs_mask, x_category].values
                                            
                        #split_mask = group_mask & (adata.obs.loc[group_mask, splitby] == split)
                        split_mask = group_obs[splitby] == split
                        x = group_obs.loc[split_mask, x_category].values
                        #print(np.sum(split_mask))

                        # extract expression values as y
                        idx = var.index.get_loc(key)
                        #y = X[obs_mask, idx].copy()
                        #print(split_mask.shape)
                        #print(X.shape)
                        y = group_X[split_mask, idx].copy()

                        # determine the x-values to fit on
                        #x_to_fit_on = np.arange(x.min(), x.max() + stepsize/10, stepsize)
                        x_to_fit_on = np.arange(-10, 10, stepsize)

                        # fit curve
                        df_split = smooth_fit(x, y, x_to_fit_on=x_to_fit_on)

                        # collect data
                        df_collection[split] = df_split
                    
                    df_collection = pd.concat(df_collection)
                    
                    # calculate mean and std
                    df = df_collection[['x', 'y_pred']].groupby('x').mean()
                    df['std'] = df_collection[['x', 'y_pred']].groupby('x').std()
                    df.reset_index(inplace=True)
                    #return df_collection

                if show:
                    axs[i].fill_between(df['x'], 
                                    df['y_pred'] - df['std'],
                                    df['y_pred'] + df['std'],
                                    alpha = 0.2,
                                    color = 'grey')
                    axs[i].plot(df['x'], df['y_pred'], label=group, linewidth=8)

                else:
                    group_collection[group] = df

            if show:
                axs[i].set_title(key)
                axs[i].legend(fontsize=24)
                axs[i].set_xlabel(xlabel, fontsize=28)
                axs[i].set_ylabel(ylabel, fontsize=28)
                axs[i].tick_params(axis='both', which='major', labelsize=24)

                if values_into_title is None:
                    axs[i].set_title("{}{}".format(key, title_suffix), fontsize=32)
                else:
                    assert len(values_into_title) == len(keys), "List of title values has not the same length as list of keys."
                    axs[i].set_title("{}{}{}".format(key, title_suffix, round(values_into_title[i], 5)), fontsize=32)
        else:
            # plot binned values
            # sort out expression values of 0 since they disturbed the result in case of well C2 (bin 7 and 8 showed expression of 0)
            print("Plotting the binned version of the plot is deprecated and currently not activated.")
            # data = data.query('expression > 0')

            # data = data.groupby([groupby, 'bin', splitby]).mean()
            # data.reset_index(inplace=True)

            # # Plotting
            # # if ax is None:
            # #     fig, ax = plt.subplots(1,1)

            # sns.lineplot(data=data, x="bin", y="expression", hue=groupby, ax=axs[i])

            # xticks = [int(elem) for elem in sorted(data.bin.unique())]
            # xticks = xticks[::1] # in case one wants to skip ticks this can be done here

            # axs[i].xaxis.set_major_locator(mticker.FixedLocator(xticks))

            # if x_limit_labels is not None:
            #     xlabels = [str(elem) for elem in xticks]
            #     xlabels[0] = x_limit_labels[0]
            #     xlabels[-1] = x_limit_labels[1]
            #     axs[i].set_xticklabels(xlabels)

            # if xlabel is not None:
            #     axs[i].set_xlabel(xlabel)

            # axs[i].set_title(key)

        if not show:
            # collect data
            group_collection = pd.concat(group_collection)
            data_collection[key] = group_collection

    if n_plots > 1:

        # check if there are empty plots remaining
        while i < n_rows * max_cols - 1:
            i+=1
            # remove empty plots
            axs[i].set_axis_off()
    
    if show:
        plt.tight_layout()
        save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)
    else:
        data_collection = pd.concat(data_collection)
        data_collection.index.names = ['key', groupby, 'id']
        return data_collection


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