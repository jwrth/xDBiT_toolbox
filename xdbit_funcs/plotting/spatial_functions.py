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
from scipy.stats import zscore

# ignore future warnings (suppresses annoying pandas warning)
warnings.simplefilter(action='ignore', category=FutureWarning)


def spatial_single(adata, keys, groupby=None, group=None, max_cols=4, pd_dataframe=None,
            header=None, headersize=18, header_names=None, 
            raw=False, layer=None,
            percent=False,
            dpi_save=300,
            obsm_key = 'spatial', plot_pixel=False,
            spot_size_unit=50, spot_type='s', clb_pad=0.05,
            axis=None, fig=None, show=True,
            patch_style=None, patch_xy=(1000, 1000), patch_radius=1000, patch_color='r',
            xlim=None, ylim=None, oversize=1, dpi_display=80, figsize=(8.2,6),
            #image=None, image_metadata=None, 
            image_key=None, lowres=True, histogram_setting=None,
            alpha=1, palette="tab10", color_dict=None, cmap='viridis', 
            header_x=0.5, header_y=0.98, header_fontsize=20,
            crange=None, crange_type='minmax', colorbar=True, clb_title=None,
            cmap_center=None,
            origin_zero=True, # whether to start axes ticks at 0
            margin=True, # whether to leave margin of one spot width around the plot
            save_only=False, savepath=None, save_background=None,
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
    adata_X, adata_var, adata_var_names = check_raw(adata, use_raw=raw, layer=layer)

    if header_names is not None:
        assert len(header_names) == len(keys)

    # check if cmap is supposed to be centered
    if cmap_center is None:
        normalize=None
    else:
        normalize = colors.CenteredNorm(vcenter=cmap_center)
        
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

        pixel_per_um = image_metadata["pixel_per_um_real"]
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

        bit_type = np.uint8 if image.max() < 256 else np.uint16
        if histogram_setting is None:
            # do min max scaling
            image = set_histogram(image, lower=image.min(), upper=np.percentile(image, 99), bit_type=bit_type)
        else:
            if histogram_setting == 'minmax':
                image = set_histogram(image, lower=image.min(), upper=image.max(), bit_type=bit_type)
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
                print(image_metadata)
                pixel_per_um = image_metadata["pixel_per_um_real"]
                scale_factor = image_metadata['tissue_hires_scalef']
            else:
                print("pixel_per_um scalefactor not found. Plotted pixel coordinates instead.")
        else:
            print("No key `spatial` in adata.uns. Therefore pixel_per_um scalefactor could not be found. Plotted pixel coordinates instead.") if verbose else None

    # extract x and y pixel coordinates and convert to micrometer
    x_pixelcoord = adata.obsm[obsm_key][:, 0].copy()
    y_pixelcoord = adata.obsm[obsm_key][:, 1].copy()
    x_coord = x_pixelcoord / pixel_per_um
    y_coord = y_pixelcoord / pixel_per_um


    # shift coordinates that they start at (0,0)
    if origin_zero:
        x_offset = x_coord.min()
        y_offset = y_coord.min()
        x_coord -= x_offset
        y_coord -= y_offset
    else:
        x_offset = y_offset = 0
    
    if xlim is None:
        if plot_pixel:
            xmin = x_pixelcoord.min() * scale_factor
            xmax = x_pixelcoord.max() * scale_factor
        else:
            # xmin = x_coord.min()
            # xmax = x_coord.max()
            xmin = np.min([x_coord.min(), y_coord.min()]) # make sure that result is always a square
            xmax = np.max([x_coord.max(), y_coord.max()])

        xlim = (xmin - spot_size_unit, xmax + spot_size_unit)
    elif margin:
        xlim[0] -= spot_size_unit
        xlim[1] += spot_size_unit

    if ylim is None:
        if plot_pixel:
            ymin = y_pixelcoord.min() * scale_factor
            ymax = y_pixelcoord.max() * scale_factor
        else:
            # ymin = y_coord.min()
            # ymax = y_coord.max() 
            ymin = np.min([x_coord.min(), y_coord.min()])
            ymax = np.max([x_coord.max(), y_coord.max()])

        ylim = (ymin - spot_size_unit, ymax + spot_size_unit)
    elif margin:
        ylim[0] -= spot_size_unit
        ylim[1] += spot_size_unit

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
        # x_ppu = pixels_per_unit[1, 0]
        y_ppu = pixels_per_unit[0, 1]
        pxs = y_ppu * spot_size_unit * oversize
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
                s = ax.scatter(x_pixelcoord * scale_factor, y_pixelcoord * scale_factor, c=color, marker=spot_type,
                                s=size, alpha=alpha, 
                                #edgecolors=None,
                                linewidths=0,
                                cmap=cmap
                                )
                if image is not None:
                    ax.imshow(image, origin='upper', cmap='gray')
            else:
                if key is not None:
                    # plot transcriptomic data
                    if not categorical:
                        s = ax.scatter(x_coord, y_coord, c=color, marker=spot_type,
                                    s=size, alpha=alpha, 
                                    #edgecolors=None, 
                                    linewidths=0,
                                    cmap=cmap, norm=normalize)
                    else:
                        sns.scatterplot(x='x_coord', y='y_coord', data=data,
                                        hue=key, marker=spot_type, s=size, 
                                        #edgecolor="none", 
                                        linewidth=0, 
                                        palette=color_dict, alpha=alpha,
                                        ax=ax,
                                        )
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
                if colorbar:
                    # divide axis to fit colorbar
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="4%", pad=0.1)
                    clb = fig.colorbar(s, cax=cax)
                    
                    # set colorbar
                    clb.ax.tick_params(labelsize=14)
                    
                    if clb_title is not None:
                        clb.ax.set_ylabel(clb_title, rotation=270, fontdict={"fontsize": 14}, labelpad=20)

                    if crange is not None:
                            clb.mappable.set_clim(crange[0], crange[1])
                    else:
                        if crange_type == 'percentile':
                            clb.mappable.set_clim(0, np.percentile(color, 99))
                        
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

    save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, show=show, dpi_save=dpi_save, save_background=save_background)


def spatial(adata, keys, groupby='id', groups=None, raw=False, layer=None, max_cols=4, 
    spot_size=50, prefix_groups='', palette="tab10", groupheader_fontsize=20,
    savepath=None, dpi_save=300, save_only=False, 
    pd_dataframe=None, normalize_crange_not_for=[], 
    dpi_display=80, header_names=None,
    crange=None, crange_type='minmax',
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

    if groupby is None:
        groups = [None]
    else:
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
                use_raw=raw, layer=layer, data_in_dataframe=data_in_dataframe, pd_dataframe=pd_dataframe, ctype=crange_type) if key not in normalize_crange_not_for else None for key in keys}

    if multigroups:
        if multikeys:
            n_rows = len(groups)
            max_cols = len(keys)
            n_plots = n_rows * max_cols
            fig, axs = plt.subplots(n_rows, max_cols, figsize=(7.6 * max_cols, 6 * n_rows), dpi=dpi_display)
            fig.tight_layout() # helps to equalize size of subplots. Without the subplots change parameters during plotting which results in differently sized spots.

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
                
                    if crange is None:
                        crange_ = crange_per_key_dict[key]
                    else:
                        crange_ = crange

                    spatial_single(adata, key, raw=raw, layer=layer, groupby=groupby,
                            group=group, fig=fig, axis=axs[row, col], show=False,
                            xlim=xlim, ylim=ylim, 
                            spot_size_unit=spot_size, crange=crange_, 
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
            fig.tight_layout() # helps to equalize size of subplots. Without the subplots change parameters during plotting which results in differently sized spots.

            if n_plots > 1:
                axs = axs.ravel()
            else:
                axs = [axs]
            
            for i, group in enumerate(groups):
                key = keys[0]

                # create color dictionary if key is categorical
                color_dict = create_color_dict(adata, key, palette)

                if crange is None:
                    crange_ = crange_per_key_dict[key]
                else:
                    crange_ = crange

                spatial_single(adata, key, raw=raw, layer=layer, groupby=groupby,
                        group=group, fig=fig, axis=axs[i], show=False, header_names=header_names,
                        xlim=xlim, ylim=ylim, spot_size_unit=spot_size, crange=crange_, crange_type=crange_type,
                        palette=palette, color_dict=color_dict, pd_dataframe=pd_dataframe, **kwargs)

                axs[i].set_title("{} - {}{}".format(key, prefix_groups, group))
                

        if n_plots > 1:
            # check if there are empty plots remaining
            while i < n_rows * max_cols - 1:
                i+=1
                # remove empty plots
                axs[i].set_axis_off()

        save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save)
    else:
        # if there is only one group use `spatial_single` function
        spatial_single(adata, keys, raw=raw, layer=layer, groupby=groupby, group=groups[0], show=True,
                        xlim=xlim, ylim=ylim, spot_size_unit=spot_size, max_cols=max_cols, header_names=header_names,
                        palette=palette, pd_dataframe=pd_dataframe, crange=crange,
                        savepath=savepath, save_only=save_only, crange_type=crange_type,
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

def expression_along_observation_value(adata, keys, x_category, groupby, splitby=None, hue=None, 
    range_min=None, range_max=None, cmap="tab10",
    extra_cats=None,
    normalize=False,
    nsteps=100,
    show_progress=False,
    use_raw=False,
    max_cols=4,
    xlabel=None,ylabel=None,
    vline=None, hline=None,
    #values_into_title=None, title_suffix='', 
    custom_titles=None,
    legend_fontsize=24, 
    plot_legend=True,
    xlabel_fontsize=28, ylabel_fontsize=28, title_fontsize=20, tick_fontsize=24,
    savepath=None, save_only=False, show=True, axis=None, return_data=False, fig=None,
    dpi_save=300,
    loess=True, 
    custom_lowess=False,
    stderr=True, 
    **kwargs):

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

                if loess:
                    # do smooth fitting
                    df = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps, custom_lowess=custom_lowess, stderr=stderr, **kwargs)
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
                    if loess:
                        df_split = smooth_fit(x, y, 
                                min=range_min, max=range_max,
                                nsteps=nsteps)
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
                    linewidth=8)

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
                axs[i].axvline(x=v, ymin=0, ymax=1, c=c, linewidth=4, linestyle='dashed')

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
            fig.tight_layout()
            save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save)
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
