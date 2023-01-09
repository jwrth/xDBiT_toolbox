import math
from typing import Any, Dict, List, Optional, Tuple, Union, Type, Literal
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from anndata import AnnData
from ..images import set_histogram, resize_image
from ..tools import check_raw, create_color_dict, extract_groups, get_crange, get_nrows_maxcols, deactivate_empty_plots
from ..readwrite import save_and_show_figure
import gc
from pathlib import Path
import dask.array as da
import pickle


class ImageData:
    '''
    Object extracting image data from anndata object.
    '''
    def __init__(
        self,
        adata: AnnData,
        image_key: Optional[str] = None,
        group: Optional[str] = None,
        lowres: bool = True,
        histogram_setting: Optional[Tuple[int, int]] = None,
        lowres_sf: float = 0.1,
        zarr_key: str = "zarr_directories"
        ):
        
        self.image_key = image_key
        self.group = group
        self.lowres = lowres
        self.histogram_setting = histogram_setting
        self.lowres_sf = lowres_sf
        self.zarr_key = zarr_key
        
        # get image data and or metadata
        if self.image_key is not None:
            self.check_source(adata)
            if self.load_from_zarr:
                self.load_image_from_zarr()
                self.process_image()
            else:
                self.load_image_from_uns(adata)
                self.process_image()
        else:
            self.no_image(adata)
    
    def check_source(self, adata):
        '''
        Check source from which the image should be loaded.
        '''
        self.load_from_zarr = False
        if self.zarr_key in adata.uns.keys():
            zarr_dirs = adata.uns[self.zarr_key][self.group]
            zarr_dirs = list(set(zarr_dirs)) # make list unique
            
            # select zarr_dirs that are a directory
            zarr_dirs = [Path(d) for d in zarr_dirs if Path(d).is_dir()]
            
            # unique_dir = False
            # dir_found = False
            # if len(zarr_dirs) == 1:
            #     # a unique file directory was found and will be used
            #     unique_dir = True
            #     dir_found = True
            # elif len(zarr_dirs) > 1:
            #     dir_found = True
            # else:
            #     dir_found = False
            
            image_found = False
            found_dirs = []
            for zd in zarr_dirs:
                channel_dirs = [elem for elem in list(zd.glob("*")) if elem.is_dir()]
                
                # check if image key is one of the channels
                chdir = [d for d in channel_dirs if self.image_key == d.stem]
                
                if len(chdir) == 1:
                    found_dirs.append(chdir[0])
                    image_found = True
            
            if image_found:
                if len(found_dirs) == 1:
                    self.zarr_path = found_dirs[0]
                    self.load_from_zarr = True
                elif len(found_dirs) > 1:
                    print("Warning: More than one possible zarr path found: {}. The first one was used for plotting.".format(found_dirs), flush=True)
                    self.zarr_path = found_dirs[0]
                    self.load_from_zarr = True
            else:
                self.load_from_zarr = False
                raise ValueError("No image path found for image key {} in zarr directories {}.".format(self.image_key, zarr_dirs))
                
        else:
            # check if image_key gives unique match
            image_key_list = adata.uns['spatial'].keys()
            image_key_matches = [k for k in image_key_list if self.image_key in k]
            if self.group is not None:
                # make sure that group name is also in image_key
                image_key_matches = [k for k in image_key_matches if self.group in k]

            if len(image_key_matches) == 0:
                raise ValueError("No match for img_key [{}] found in list of image_keys: {}".format(
                    self.image_key, image_key_list))

            elif len(image_key_matches) > 1:
                raise ValueError("More than one possible match for img_key [{}] found: {}".format(
                    self.image_key, image_key_matches
                    ))
                
            else:
                self.image_key = image_key_matches[0]
            
    def load_image_from_uns(self, adata):
        self.image_metadata = adata.uns['spatial'][self.image_key]['scalefactors']
        self.pixel_per_um = self.image_metadata["pixel_per_um_real"]
        
        if self.lowres:
            if 'lowres' in adata.uns['spatial'][self.image_key]['images'].keys():
                self.image = adata.uns['spatial'][self.image_key]['images']['lowres']
                self.scale_factor = self.image_metadata['tissue_lowres_scalef']
            else:
                print('`hires` image displayed because no `lowres` images was found.')
                self.image = adata.uns['spatial'][self.image_key]['images']['hires']
                self.scale_factor = self.image_metadata['tissue_hires_scalef']
        else:
            self.image = adata.uns['spatial'][self.image_key]['images']['hires']
            self.scale_factor = self.image_metadata['tissue_hires_scalef']
            
    def load_image_from_zarr(self):
        # load image
        self.zarr_path = Path(self.zarr_path)
        self.image = da.from_zarr(self.zarr_path)
        self.image = self.image.compute() # convert to numpy array
        
        if self.lowres:
            self.image = resize_image(self.image, scale_factor=self.lowres_sf)
            self.scale_factor = self.lowres_sf
        else:
            self.scale_factor = 1
            self.scale_factor = self.image_metadata['tissue_hires_scalef']

        # load and extract metadata
        meta_file = self.zarr_path.parent / "metadata.pkl"
        with open(meta_file, "rb") as f:
            self.image_metadata = pickle.load(f)
            
        self.pixel_per_um = self.image_metadata["pixel_per_um_real"]

    def process_image(self):
        self.bit_type = np.uint8 if self.image.max() < 256 else np.uint16
        if self.histogram_setting is None:
            # do min max scaling
            self.image = set_histogram(self.image, 
                                        lower=self.image.min(), upper=np.percentile(self.image, 99), 
                                        bit_type=self.bit_type)
        else:
            if self.histogram_setting == 'minmax':
                self.image = set_histogram(self.image, 
                                        lower=self.image.min(), upper=self.image.max(), 
                                        bit_type=self.bit_type)
            elif isinstance(self.histogram_setting, tuple):
                self.image = set_histogram(self.image, 
                                        lower=self.histogram_setting[0], upper=self.histogram_setting[1], 
                                        bit_type=self.bit_type)
            elif (self.histogram_setting > 0) & (self.histogram_setting < 1):
                self.image = set_histogram(self.image, lower=self.image.min(), upper=int(self.image.max() * self.histogram_setting), 
                                        bit_type=self.bit_type)

            else:
                raise ValueError('Unknown format of `histogram_setting`. Must be either "minmax" or (minval, maxval) or value between 0 and 1')
        
    def no_image(self, adata):
        self.image = None
        self.scale_factor = 1
        self.pixel_per_um = 1
        self.image_metadata = None
        # search for pixel_per_um scalefactor
        if 'spatial' in adata.uns.keys():
            image_key_list = adata.uns['spatial'].keys()
            if self.group is not None:
                # make sure that group name is also in image_key and select first option
                self.image_key = [k for k in image_key_list if self.group in k][0]
                first_entry = adata.uns['spatial'][self.image_key]
            else:
                # get first entry of dictionary as image_metadata
                first_entry = next(iter(adata.uns['spatial'].values()))
            
            # extract image metadata if possible
            if 'scalefactors' in first_entry:
                self.image_metadata = first_entry['scalefactors']
                self.pixel_per_um = self.image_metadata["pixel_per_um_real"]
                self.scale_factor = self.image_metadata['tissue_hires_scalef']
            else:
                print("pixel_per_um scalefactor not found. Plotted pixel coordinates instead.")
        else:
            print("No key `spatial` in adata.uns. Therefore pixel_per_um scalefactor could not be found. Plotted pixel coordinates instead.")

class ExpressionData:
    '''
    Object extracting spatial coordinates and expression data from anndata object.
    '''
    def __init__(
        self, 
        adata: AnnData,
        key: List[str],
        ImageDataObject: Type[ImageData],
        raw: bool = False,
        layer: Optional[str] = None,
        obsm_key: str = 'spatial',
        origin_zero: bool = True, # whether to start axes ticks at 0
        plot_pixel: bool = False,
        xlim: Optional[Tuple[int, int]] = None,
        ylim: Optional[Tuple[int, int]] = None,
        spot_size_unit: float = 50, # whether to leave margin of one spot width around the plot
        margin: bool = True
        ):
        
        # add arguments to object
        self.key = key
        self.raw = raw
        self.layer = layer
        self.xlim = xlim
        self.ylim = ylim
                
        ## Extract coordinates    
        # extract parameters from ImageDataObject
        pixel_per_um = ImageDataObject.pixel_per_um
        scale_factor = ImageDataObject.scale_factor
        
        # extract x and y pixel coordinates and convert to micrometer
        self.x_pixelcoord = adata.obsm[obsm_key][:, 0].copy()
        self.y_pixelcoord = adata.obsm[obsm_key][:, 1].copy()
        self.x_coord = self.x_pixelcoord / pixel_per_um
        self.y_coord = self.y_pixelcoord / pixel_per_um

        # shift coordinates that they start at (0,0)
        if origin_zero:
            self.x_offset = self.x_coord.min()
            self.y_offset = self.y_coord.min()
            self.x_coord -= self.x_offset
            self.y_coord -= self.y_offset
        else:
            self.x_offset = self.y_offset = 0
        
        if self.xlim is None:
            if plot_pixel:
                xmin = self.x_pixelcoord.min() * scale_factor
                xmax = self.x_pixelcoord.max() * scale_factor
            else:
                xmin = np.min([self.x_coord.min(), self.y_coord.min()]) # make sure that result is always a square
                xmax = np.max([self.x_coord.max(), self.y_coord.max()])

            self.xlim = (xmin - spot_size_unit, xmax + spot_size_unit)
        elif margin:
            self.xlim[0] -= spot_size_unit
            self.xlim[1] += spot_size_unit

        if self.ylim is None:
            if plot_pixel:
                ymin = self.y_pixelcoord.min() * scale_factor
                ymax = self.y_pixelcoord.max() * scale_factor
            else:
                # ymin = y_coord.min()
                # ymax = y_coord.max() 
                ymin = np.min([self.x_coord.min(), self.y_coord.min()])
                ymax = np.max([self.x_coord.max(), self.y_coord.max()])

            self.ylim = (ymin - spot_size_unit, ymax + spot_size_unit)
        elif margin:
            self.ylim[0] -= spot_size_unit
            self.ylim[1] += spot_size_unit
            
        ## Extract expression data
        # check if plotting raw data
        adata_X, adata_var, adata_var_names = check_raw(adata, 
                                                        use_raw=self.raw, 
                                                        layer=self.layer)
        
        self.data = adata.obs.copy()
        self.data['x_coord'] = self.x_coord
        self.data['y_coord'] = self.y_coord
        
        # locate gene in matrix and extract values
        if key in adata_var_names:
            idx = adata_var.index.get_loc(key)
            self.color = adata_X[:, idx].copy()
            self.categorical = False
            self.exists = True
            
        elif key in adata.obs.columns:
            self.exists = True
            if adata.obs[key].dtype.name == 'category':
                self.categorical = True
            else:
                self.color = adata.obs[key].values
                self.categorical = False
        else:
            self.exists = False
        
class MultiSpatialPlot:
    '''
    Class to plot xDbit spatial transcriptomics data.
    '''
    def __init__(self, 
                 adata: AnnData, 
                 keys: Union[str, List[str]], 
                 groupby: str = 'id',
                 groups: Optional[List[str]] = None,
                 raw: bool = False,
                 layer: Optional[str] = None,
                 max_cols: int = 4,
                 xlim: Optional[Tuple[float, float]] = None,
                 ylim: Optional[Tuple[float, float]] = None,
                 normalize_crange_not_for: List = [],
                 crange: Optional[List[int]] = None, 
                 crange_type: Literal['minmax', 'percentile'] = 'minmax',
                 palette: str = 'tab10',
                 cmap_center: Optional[float] = None,
                 dpi_display: int = 80,
                 obsm_key: str = 'spatial',
                 origin_zero: bool = True,
                 plot_pixel: bool = False,
                 spot_size_unit: float = 50,
                 margin: bool = True,
                 spot_type: str = 's',
                 cmap: str = 'viridis',
                 alpha: float = 1,
                 colorbar: bool = True,
                 clb_title: Optional[str] = None,
                 header: Optional[str] = None,
                 
                 # image stuff
                 image_key: Optional[str] = None,
                 histogram_setting: Optional[Tuple[int, int]] = None,
                 lowres: bool = True,
                 
                 # saving
                 savepath: Optional[str] = None,
                 save_only: bool = False,
                 dpi_save: int = 300,
                 
                 # less important
                 prefix_groups: str = '',
                 groupheader_fontsize: int = 20
                 ):
        
        self.adata = adata
        self.keys = keys
        self.groupby = groupby
        self.groups = groups
        self.raw = raw
        self.layer = layer
        self.max_cols = max_cols
        self.xlim = xlim
        self.ylim = ylim
        self.normalize_crange_not_for = normalize_crange_not_for
        self.crange = crange
        self.crange_type = crange_type
        self.palette = palette
        self.cmap_center = cmap_center
        self.dpi_display = dpi_display
        self.obsm_key = obsm_key
        self.origin_zero = origin_zero,
        self.plot_pixel = plot_pixel
        self.spot_size_unit = spot_size_unit
        self.margin = margin
        self.spot_type = spot_type
        self.cmap = cmap
        self.alpha = alpha
        self.colorbar = colorbar
        self.clb_title = clb_title
        self.header = header
        self.savepath = savepath
        self.save_only = save_only
        self.dpi_save = dpi_save
        
        # image stuff
        self.image_key = image_key
        self.histogram_setting = histogram_setting
        self.lowres = lowres
        
        # less important
        self.prefix_groups = prefix_groups
        self.groupheader_fontsize = groupheader_fontsize
        
        # check arguments
        self.check_arguments()
        
        # plotting
        self.setup_subplots()
        self.plot_to_subplots()
                
        save_and_show_figure(
            savepath=self.savepath,
            fig=self.fig,
            save_only=self.save_only,
            dpi_save=self.dpi_save
            )

        gc.collect()
        
    def check_arguments(self):
        # convert arguments to lists
        self.keys = [self.keys] if isinstance(self.keys, str) else list(self.keys)
        if self.xlim is not None:
            self.xlim = [self.xlim] if isinstance(self.xlim, str) else list(self.xlim)
        if self.xlim is not None:
            self.ylim = [self.ylim] if isinstance(self.ylim, str) else list(self.ylim)
            
        # check if cmap is supposed to be centered
        if self.cmap_center is None:
            self.normalize=None
        else:
            self.normalize = colors.CenteredNorm(vcenter=self.cmap_center)
            
        # set multiplot variables
        self.multikeys = False
        self.multigroups = False
        if len(self.keys) > 1:
            self.multikeys = True
            
        if self.groupby is None:
            self.groups = [None]
        else:
            if self.groups is None:
                self.groups = list(self.adata.obs[self.groupby].unique())
            else:
                self.groups = [self.groups] if isinstance(self.groups, str) else list(self.groups)
            if len(self.groups) > 1:
                self.multigroups = True

        # determine the color range for each key
        self.crange_per_key_dict = {
            key: get_crange(
                self.adata, self.groupby, self.groups, key, 
                use_raw=self.raw, layer=self.layer, 
                ctype=self.crange_type
                ) if key not in self.normalize_crange_not_for else None for key in self.keys
            }
        
    def setup_subplots(self):
        if self.multigroups:
            if self.multikeys:
                n_rows = len(self.groups)
                self.max_cols = len(self.keys)
                n_plots = n_rows * self.max_cols
                self.fig, self.axs = plt.subplots(n_rows, self.max_cols, 
                                                  figsize=(7.6 * self.max_cols, 6 * n_rows), 
                                                  dpi=self.dpi_display)
                self.fig.tight_layout() # helps to equalize size of subplots. Without the subplots change parameters during plotting which results in differently sized spots.
                
            else:
                n_plots = len(self.groups)
                if n_plots > self.max_cols:
                    n_rows = math.ceil(n_plots / self.max_cols)
                else:
                    n_rows = 1
                    self.max_cols = n_plots

                self.fig, self.axs = plt.subplots(n_rows, self.max_cols, 
                                        figsize=(7.6 * self.max_cols, 6 * n_rows), 
                                        dpi=self.dpi_display)
                self.fig.tight_layout() # helps to equalize size of subplots. Without the subplots change parameters during plotting which results in differently sized spots.

                if n_plots > 1:
                    self.axs = self.axs.ravel()
                else:
                    self.axs = [self.axs]
                    
                deactivate_empty_plots(
                    n_plots=n_plots, 
                    nrows=n_rows,
                    ncols=self.max_cols,
                    axis=self.axs
                    )
        
        else:
            n_plots = len(self.keys)
            if self.max_cols is None:
                self.max_cols = n_plots
                n_rows = 1
            else:
                if n_plots > self.max_cols:
                    n_rows = math.ceil(n_plots / self.max_cols)
                else:
                    n_rows = 1
                    self.max_cols = n_plots
                    
            self.fig, self.axs = plt.subplots(n_rows, self.max_cols, 
                                              figsize=(7.6 * self.max_cols, 6 * n_rows), 
                                              dpi=self.dpi_display)
            
            if n_plots > 1:
                self.axs = self.axs.ravel()
            else:
                self.axs = np.array([self.axs])
                
            # remove axes from empty plots
            deactivate_empty_plots(
                    n_plots=n_plots, 
                    nrows=n_rows,
                    ncols=self.max_cols,
                    axis=self.axs
                    )
            
        if self.header is not None:
            plt.suptitle(self.header, fontsize=18, x=0.5, y=0.98)
            
    def plot_to_subplots(self):
        i = 0
        for row, group in enumerate(self.groups):
            # extract this group
            ad, mask = extract_groups(self.adata, self.groupby, group, 
                                            extract_uns=False, return_mask=True)
            # get image data for this group
            ImgD = ImageData(
                    adata=ad,
                    image_key=self.image_key,
                    group=group,
                    lowres=self.lowres,
                    histogram_setting=self.histogram_setting,
                    )
            
            for col, key in enumerate(self.keys):
                # create color dictionary if key is categorical
                color_dict = create_color_dict(ad, key, self.palette)
                
                if self.crange is None:
                    crange_ = self.crange_per_key_dict[key]
                else:
                    crange_ = self.crange
                    
                # get axis to plot
                if len(self.axs.shape) == 2:
                    ax = self.axs[row, col]
                elif len(self.axs.shape) == 1:
                    if self.multikeys:
                        ax = self.axs[col]
                    else:
                        ax = self.axs[i]
                else:
                    raise ValueError("`len(self.axs.shape)` has wrong shape {}. Requires 1 or 2.".format(len(self.axs.shape)))
                
                # counter for axis
                i+=1
                
                # get data
                ExpD = ExpressionData(
                    adata=ad,
                    key=key,
                    ImageDataObject=ImgD,
                    raw=self.raw,
                    layer=self.layer,
                    obsm_key=self.obsm_key,
                    origin_zero=self.origin_zero,
                    plot_pixel=self.plot_pixel,
                    xlim=self.xlim,
                    ylim=self.ylim,
                    spot_size_unit=self.spot_size_unit,
                    margin=self.margin
                )
                
                if ExpD.exists:
                    # set axis
                    ax.set_xlim(ExpD.xlim[0], ExpD.xlim[1])
                    ax.set_ylim(ExpD.ylim[0], ExpD.ylim[1])

                    if ImgD.image_metadata is None or self.plot_pixel:
                        ax.set_xlabel('pixels', fontsize=14)
                        ax.set_ylabel('pixels', fontsize=14)
                    else:
                        ax.set_xlabel('µm', fontsize=14)
                        ax.set_ylabel('µm', fontsize=14)
                    
                    ax.invert_yaxis()
                    ax.grid(False)
                    ax.set_aspect(1)
                    ax.set_facecolor('k')
                    ax.tick_params(labelsize=12)
                    
                    if self.multigroups and not self.multikeys:
                        ax.set_title(group + "\n" + ExpD.key, fontsize=20, fontweight='bold')
                    else:
                        # set titles
                        ax.set_title(ExpD.key, fontsize=20, fontweight='bold')
                        
                        if col == 0:
                            ax.annotate(group, 
                                        xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                                        xycoords=ax.yaxis.label, textcoords='offset points',
                                        size=20, 
                                        ha='right', va='center', weight='bold')

                    # calculate marker size
                    pixels_per_unit = ax.transData.transform(
                        [(0, 1), (1, 0)]) - ax.transData.transform((0, 0))
                    # x_ppu = pixels_per_unit[1, 0]
                    y_ppu = pixels_per_unit[0, 1]
                    pxs = y_ppu * self.spot_size_unit
                    size = (72. / self.fig.dpi * pxs)**2

                    if color_dict is None:
                        color_dict = self.palette
                    
                    # plot single spatial plot in given axis
                    self.single_spatial(
                        ExpData=ExpD,
                        ImgData=ImgD,
                        size=size,
                        axis=ax,
                        color_dict=color_dict,
                        crange=crange_
                                    )
                else:
                    print("Key '{}' not found.".format(key), flush=True)
                    ax.set_axis_off()
                    
                # free RAM
                del ExpD
                gc.collect()
                
            # free RAM
            del ImgD
            gc.collect()
                    
    def single_spatial(
        self, 
        ExpData: Type[ExpressionData],
        ImgData: Type[ImageData],
        size: float, 
        axis: plt.Axes,
        color_dict: Dict,
        crange: Optional[Tuple[float, float]]
        ):
        
        if self.plot_pixel:
            # plot image data
            if ImgData.image is not None:
                axis.imshow(ImgData.image, origin='upper', cmap='gray')
                
            # plot transcriptomic data
            s = axis.scatter(
                ExpData.x_pixelcoord * ImgData.scale_factor, 
                ExpData.y_pixelcoord * ImgData.scale_factor, 
                c=ExpData.color, marker=self.spot_type,
                s=size, alpha=self.alpha,
                linewidths=0,
                cmap=self.cmap
                )
        else:
            # plot image data
            if ImgData.image is not None:
                axis.imshow(ImgData.image, 
                            extent=(
                                -0.5 - ExpData.x_offset, 
                                ImgData.image.shape[1] / ImgData.pixel_per_um / ImgData.scale_factor - 0.5 - ExpData.x_offset, 
                                ImgData.image.shape[0] / ImgData.pixel_per_um / ImgData.scale_factor - 0.5 - ExpData.y_offset, 
                                -0.5 - ExpData.y_offset
                                ), 
                            origin='upper', cmap='gray')
            # plot transcriptomic data
            if not ExpData.categorical:
                s = axis.scatter(ExpData.x_coord, ExpData.y_coord, c=ExpData.color, marker=self.spot_type,
                            s=size, alpha=self.alpha, 
                            linewidths=0,
                            cmap=self.cmap, norm=self.normalize)
            else:
                sns.scatterplot(x='x_coord', y='y_coord', data=ExpData.data,
                                hue=ExpData.key, marker=self.spot_type, s=size, 
                                #edgecolor="none", 
                                linewidth=0, 
                                palette=color_dict, alpha=self.alpha,
                                ax=axis,
                                )
                
        # plot legend
        if ExpData.categorical:
            axis.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
        else:
            if self.colorbar:
                # divide axis to fit colorbar
                divider = make_axes_locatable(axis)
                cax = divider.append_axes("right", size="4%", pad=0.1)
                clb = self.fig.colorbar(s, cax=cax)
                
                # set colorbar
                clb.ax.tick_params(labelsize=14)
                
                if self.clb_title is not None:
                    clb.ax.set_ylabel(self.clb_title, 
                                      rotation=270, 
                                      fontdict={"fontsize": 14}, 
                                      labelpad=20)

                if crange is not None:
                        clb.mappable.set_clim(crange[0], crange[1])
                else:
                    if self.crange_type == 'percentile':
                        clb.mappable.set_clim(0, np.percentile(ExpData.color, 99))
