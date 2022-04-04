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
from ..tools import extract_groups, check_raw, create_color_dict, get_nrows_maxcols, get_crange
from ..calculations import smooth_fit
from ..readwrite import save_and_show_figure
from ..images import set_histogram
from tqdm import tqdm


def spatial_single(adata, keys, groupby=None, group=None, max_cols=4, pd_dataframe=None,
            header=None, headersize=18, header_names=None, raw=False, percent=False,
            dpi_save=300,
            obsm_key = 'spatial', plot_pixel=False,
            spot_size=50, spot_type='s', clb_pad=-0.05,
            axis=None, fig=None, show=True,
            patch_style=None, patch_xy=(1000, 1000), patch_radius=1000, patch_color='r',
            xlim=None, ylim=None, oversize=1, dpi_display=80, figsize=(8.2,6),
            #image=None, image_metadata=None, 
            image_key=None, lowres=True, histogram_setting=None,
            alpha=1, palette="tab10", color_dict=None,
            header_x=0.5, header_y=0.98, header_fontsize=20,
            save_only=False, savepath=None, save_background=None, crange=None, colorbar=True,
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
    if header_names is not None:
        header_names = [header_names] if isinstance(header_names, str) else list(header_names)

    # check if plotting raw data
    adata_X, adata_var, adata_var_names = check_raw(adata, use_raw=raw)


    if header_names is not None:
        assert len(header_names) == len(keys)

    # get image data and or metadata
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
            if histogram_setting == 'minmax':
                print('bin hier')
                image = set_histogram(image, lower=image.min(), upper=image.max(), bit_type=bit_type, )
            elif isinstance(histogram_setting, tuple):
                image = set_histogram(image, lower=histogram_setting[0], upper=histogram_setting[1], bit_type=bit_type)
            elif (histogram_setting > 0) & (histogram_setting < 1):
                image = set_histogram(image, lower=image.min(), upper=int(image.max() * histogram_setting), bit_type=bit_type)

            else:
                print('Unknown format of `histogram_setting`. Must be either "minmax" or (minval, maxval) or value between 0 and 1', flush=True)
    else:
        image = None
        scale_factor = 1
        pixel_per_um = 1
        image_metadata = None
        # search for pixel_per_um scalefactor
        if 'spatial' in adata.uns.keys():
            image_key_list = adata.uns['spatial'].keys()
            if group is not None:
                # make sure that group name is also in image_key and select first option
                image_key = [k for k in image_key_list if group in k][0]
                first_entry = adata.uns['spatial'][image_key]
            else:
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
    x_pixelcoord = adata.obsm[obsm_key][:, 0]
    y_pixelcoord = adata.obsm[obsm_key][:, 1]
    x_coord = x_pixelcoord / pixel_per_um
    y_coord = y_pixelcoord / pixel_per_um


    # shift coordinates that they start at (0,0)
    x_offset = x_coord.min()
    y_offset = y_coord.min()
    x_coord -= x_offset
    y_coord -= y_offset

    if xlim is None:
        if plot_pixel:
            xmin = x_pixelcoord.min()
            xmax = x_pixelcoord.max()
        else:
            xmin = x_coord.min()
            xmax = x_coord.max()

        xlim = (xmin - spot_size, xmax + spot_size)
    else:
        xlim[0] -= spot_size
        xlim[1] += spot_size

    if ylim is None:
        if plot_pixel:
            ymin = y_pixelcoord.min()
            ymax = y_pixelcoord.max()
        else:
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
        show = False # otherwise plotting into given axes wouldn't work

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
                    ax.set_axis_off()
                    return

            if percent:
                color *= 100
        
        # set axis
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        if image_metadata is None or plot_pixel:
            ax.set_xlabel('pixels', fontsize=14)
            ax.set_ylabel('pixels', fontsize=14)
        else:
            ax.set_xlabel('µm', fontsize=14)
            ax.set_ylabel('µm', fontsize=14)

        if header_names is None:
            ax.set_title(key, fontsize=header_fontsize, fontweight='bold')
        else:
            plot_name = header_names[i]
            ax.set_title(plot_name, fontsize=header_fontsize, fontweight='bold')
        
        ax.invert_yaxis()
        ax.grid(False)
        ax.set_aspect(1)
        ax.set_facecolor('k')
        ax.tick_params(labelsize=12)

        # calculate marker size
        pixels_per_unit = ax.transData.transform(
            [(0, 1), (1, 0)]) - ax.transData.transform((0, 0))
        #x_ppu = pixels_per_unit[1, 0]
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
            if plot_pixel:
                # plot transcriptomic data
                s = ax.scatter(x_pixelcoord, y_pixelcoord, c=color, marker=spot_type,
                                s=size, alpha=alpha, edgecolors=None)
                if image is not None:
                    ax.imshow(image, origin='upper', cmap='gray')
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
                if colorbar:
                    clb = fig.colorbar(s, ax=ax, shrink=1, 
                        pad=clb_pad, 
                        aspect=40)
                    clb.ax.tick_params(labelsize=14)
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


def spatial(adata, keys, groupby='well_name', groups=None, raw=False, max_cols=4, 
    spot_size=50, prefix_groups='', palette="tab10", groupheader_fontsize=20,
    savepath=None, dpi_save=300, show=True, save_only=False, pd_dataframe=None, normalize_crange_not_for=[], 
    dpi_display=80, header_names=None,
    xlim=None, ylim=None,
    **kwargs):
    
    '''
    Creates spatial plot that is grouped by a certain parameter in `adata.obs`.
    The resulting plots has the groups as rows and columns as genes.
    '''
    
    # check keys and groups
    keys = [keys] if isinstance(keys, str) else list(keys)
    multikeys = False
    multigroups = False
    if len(keys) > 1:
        multikeys = True

    if groups is None:
        groups = list(adata.obs[groupby].unique())
    else:
        groups = [groups] if isinstance(groups, str) else list(groups)
    if len(groups) > 1:
        multigroups = True

    if header_names is not None:
        assert len(header_names) == len(keys)

    # check if dataframe is given
    if isinstance(pd_dataframe, pd.DataFrame):
        data_in_dataframe=True
    else:
        data_in_dataframe=False

    # determine the color range for each key
    crange_per_key_dict = {key: get_crange(adata, groupby, groups, key, 
                use_raw=raw, data_in_dataframe=data_in_dataframe, pd_dataframe=pd_dataframe) if key not in normalize_crange_not_for else None for key in keys}
    if multigroups:
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

                    # create color dictionary if key is categorical
                    color_dict = create_color_dict(adata, key, palette)

                    spatial_single(adata, key, raw=raw, groupby=groupby,
                            group=group, fig=fig, axis=axs[row, col], show=False,
                            xlim=xlim, ylim=ylim, 
                            spot_size=spot_size, crange=crange_per_key_dict[key], 
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
            
            for i, group in enumerate(groups):
                key = keys[0]

                # create color dictionary if key is categorical
                color_dict = create_color_dict(adata, key, palette)

                spatial_single(adata, key, raw=raw, groupby=groupby,
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
    else:
        # if there is only one group use `spatial_single` function
        spatial_single(adata, keys, raw=raw, groupby=groupby, group=groups[0], show=True,
                        xlim=xlim, ylim=ylim, spot_size=spot_size, max_cols=max_cols,
                        palette=palette, pd_dataframe=pd_dataframe, 
                        savepath=savepath, save_only=save_only,
                        **kwargs
                        )


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
    range_min=None, range_max=None, 
    #stepsize=0.01, 
    nsteps=100,
    show_progress=False,
    n_bins=20, use_raw=False,
    max_cols=4,
    x_limit_labels=None, 
    xlabel=None,ylabel=None,
    values_into_title=None, title_suffix='',
    #ax=None, 
    legend_fontsize=24, xlabel_fontsize=28, ylabel_fontsize=28, title_fontsize=20, tick_fontsize=24,
    savepath=None, save_only=False, show=True, axis=None, return_data=False, fig=None,
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

    #if show:
    if not return_data:
        # prepare plotting
        if axis is None:
            n_plots, n_rows, max_cols = get_nrows_maxcols(keys, max_cols)
            fig, axs = plt.subplots(n_rows,max_cols, figsize=(8*max_cols, 6*n_rows))

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
                    x = adata.obs.loc[group_mask, x_category].values

                    # extract expression values as y
                    idx = var.index.get_loc(key)
                    #y = X[:, idx].copy()
                    y = group_X[:, idx].copy()
                    
                    # do smooth fitting
                    df = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps)

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
                        df_split = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps)

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
                    axs[i].fill_between(df['x'], 
                                    #df['y_pred'] - df['std'],
                                    #df['y_pred'] + df['std'],
                                    df['conf_lower'],
                                    df['conf_upper'],
                                    alpha = 0.2,
                                    color = 'grey')
                    axs[i].plot(df['x'], df['y_pred'], label=group, linewidth=8)



            if not return_data:
                if xlabel is None:
                    xlabel = x_category
                if ylabel is None:
                    ylabel = "Gene expression"


                axs[i].set_title(key, fontsize=title_fontsize)
                axs[i].legend(fontsize=legend_fontsize)
                axs[i].set_xlabel(xlabel, fontsize=xlabel_fontsize)
                axs[i].set_ylabel(ylabel, fontsize=ylabel_fontsize)
                axs[i].tick_params(axis='both', which='major', labelsize=tick_fontsize)

                if values_into_title is None:
                    axs[i].set_title("{}{}".format(key, title_suffix), fontsize=title_fontsize)
                else:
                    assert len(values_into_title) == len(keys), "List of title values has not the same length as list of keys."
                    axs[i].set_title("{}{}{}".format(key, title_suffix, round(values_into_title[i], 5)), fontsize=title_fontsize)
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

        if return_data:
            # collect data
            group_collection = pd.concat(group_collection)
            data_collection[key] = group_collection


    
    if return_data:
        # close plot
        plt.close()

        # return data
        data_collection = pd.concat(data_collection)
        data_collection.index.names = ['key', groupby, 'id']
        return data_collection

    else: 
        if n_plots > 1:

            # check if there are empty plots remaining
            while i < n_rows * max_cols - 1:
                i+=1
                # remove empty plots
                axs[i].set_axis_off()
        if show:
            fig.tight_layout()
            save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)
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
