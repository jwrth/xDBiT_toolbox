import math
from typing import Any, Dict, List, Optional, Tuple, Union, Type, Literal
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from anndata import AnnData
from ..tools import check_raw, create_color_dict, extract_groups, get_crange, get_nrows_maxcols, deactivate_empty_plots
from ..readwrite import save_and_show_figure
import gc
from ..datasets import ImageData, ExpressionData

        
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
                 fig: plt.figure = None,
                 ax: plt.Axes = None,
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
                 show: bool = True,
                 
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
        self.fig = fig
        self.ax = ax
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
        self.show = show
        
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
        if self.ax is None:
            self.setup_subplots()
        else:
            assert self.fig is not None, "If axis for plotting is given, also a figure object needs to be provided via `fig`"
            assert len(self.keys) == 1, "If single axis is given not more than one key is allowed."
            
        self.plot_to_subplots()
                
        save_and_show_figure(
            savepath=self.savepath,
            fig=self.fig,
            save_only=self.save_only,
            show=self.show,
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
                if self.ax is None:
                    if len(self.axs.shape) == 2:
                        ax = self.axs[row, col]
                    elif len(self.axs.shape) == 1:
                        if self.multikeys:
                            ax = self.axs[col]
                        else:
                            ax = self.axs[i]
                    else:
                        raise ValueError("`len(self.axs.shape)` has wrong shape {}. Requires 1 or 2.".format(len(self.axs.shape)))
                else:
                    ax = self.ax
                
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
