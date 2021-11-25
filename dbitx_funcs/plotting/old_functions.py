from matplotlib import pyplot as plt
from matplotlib import patches
import math
import pandas as pd

def genes(adata, keys, max_cols=None, pd_expr_data=None, save=False, header=None, headersize=18, raw=False,
                       plot_names = None, savepath="figures/spatial_genes.png", save_background=None, 
                       dpi_save=300, resolution=50, axis=None, fig=None, inline=True, 
                       patch_style=None, patch_xy=(1000,1000), patch_width=1000, patch_color='r',
                       xlim = None, ylim = None, oversize=1, dpi_display=80, image=None, image_metadata=None, alpha=1, 
                       header_x=0.5, header_y=0.98):

    if raw:
        adata_var = adata.raw.var
        adata_X = adata.raw.X
    else:
        adata_var = adata.var
        adata_X = adata.X
        
    if plot_names is not None:
        assert len(plot_names) == len(keys)
        
    if xlim is None:
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            x_origin = image_metadata["upper_left_point"][0]/pixel_per_um + resolution/2
        else:
            x_origin = 0
        xmin = adata.obs["um_col"].min()
        xmax = adata.obs["um_col"].max()
        xlim = (xmin - resolution, xmax + resolution)
    else:
        xlim[0] -= resolution
        xlim[1] += resolution
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            x_origin = image_metadata["upper_left_point"][0]/pixel_per_um + resolution/2
        else:
            x_origin = 0
        
    if ylim is None:
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            y_origin = image_metadata["upper_left_point"][1]/pixel_per_um + resolution/2
        else:
            y_origin = 0
        ymin = adata.obs["um_row"].min()
        ymax = adata.obs["um_row"].max()
        ylim = (ymin - resolution, ymax + resolution)
    else:
        ylim[0] -= resolution
        ylim[1] += resolution
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            y_origin = image_metadata["upper_left_point"][1]/pixel_per_um + resolution/2
        else:
            y_origin = 0
    
    if axis is not None:
        n_plots = 1
        if isinstance(keys, str):
            key = keys
        elif isinstance(keys, list):
            key = keys[0]
            
    else:
        if isinstance(keys, str):
            n_plots = 1
            key = keys
        elif isinstance(keys, list) and len(keys) > 1:
            n_plots = len(keys)
        elif isinstance(keys, list) and len(keys) == 1:
            n_plots = len(keys)
            key = keys[0]
        else:
            print('Keys have unknown type')
    
    if max_cols is not None:
        if n_plots > max_cols:
            n_rows = math.ceil(n_plots / max_cols)
        else:
            n_rows = 1
    else:
        max_cols = n_plots
        n_rows = 1
        
    # Plotting
    if axis is None:
        fig, ax = plt.subplots(n_rows, max_cols, figsize=(7.6*max_cols, 6*n_rows), dpi=dpi_display)
    else:
        ax = axis
    
    if n_plots > 1:
        ax = ax.ravel()
        for i, key in enumerate(keys):
            # locate gene in matrix and extract values
            if isinstance(pd_expr_data, pd.DataFrame):
                color = pd_expr_data[key]
                
            else:
                idx = adata_var.index.get_loc(key)
                color = adata_X[:, idx]
            
            
            # calculate marker size
            ax[i].set_xlim(xlim[0], xlim[1])
            ax[i].set_ylim(ylim[0], ylim[1])
            pixels_per_unit = ax[i].transData.transform([(0,1),(1,0)])-ax[i].transData.transform((0,0))
            x_ppu = pixels_per_unit[1, 0]
            y_ppu = pixels_per_unit[0, 1]
            pxs = y_ppu * resolution * oversize
            size = (72. / fig.dpi * pxs)**2

            # plot
            if image is None:
                s = ax[i].scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', s=size)
            else:
                s = ax[i].scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', 
                               s=size, alpha=alpha, edgecolors=None)
                ax[i].imshow(image, 
                          extent=[0 - x_origin, image.shape[1]/pixel_per_um - x_origin, 
                                  0 - y_origin, image.shape[0]/pixel_per_um - y_origin], 
                          origin='lower'
                         )
            fig.colorbar(s, ax=ax[i])
            ax[i].set_xlabel('µm')
            ax[i].set_ylabel('µm')
            if plot_names == None:
                ax[i].set_title(key)
            else:
                ax[i].set_title(plot_names[i])
            ax[i].set_facecolor('k')
            ax[i].invert_yaxis()
            ax[i].grid(False)
    else:
        # locate gene in matrix and extract values
        if isinstance(pd_expr_data, pd.DataFrame):
            color = pd_expr_data[key]
        else:
            idx = adata_var.index.get_loc(key)
            color = adata_X[:, idx]
        
        # calculate marker size
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        pixels_per_unit = ax.transData.transform([(0,1),(1,0)])-ax.transData.transform((0,0))
        x_ppu = pixels_per_unit[1, 0]
        y_ppu = pixels_per_unit[0, 1]
        pxs = y_ppu * resolution * oversize
        size = (72. / fig.dpi * pxs)**2

        # plot
        if image is None:
            s = ax.scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', s=size)
        else:            
            # plot data
            s = ax.scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', 
                           s=size, alpha=alpha, edgecolors=None)
            ax.imshow(image, 
                      extent=[0 - x_origin, image.shape[1]/pixel_per_um - x_origin, 
                              0 - y_origin, image.shape[0]/pixel_per_um - y_origin], 
                      origin='lower'
                     )
        if patch_style is not None:
            if patch_style is "rectangle":
                # Create a Rectangle patch
                rect = patches.Rectangle(patch_xy, patch_width, patch_width, linewidth=2, 
                                         edgecolor=patch_color,facecolor='none')
            if patch_style is "circle":
                # Create a Rectangle patch
                rect = patches.Circle(patch_xy, patch_width, linewidth=2, 
                                         edgecolor=patch_color,facecolor='none')
            if patch_style is "line":
                # Create a Rectangle patch
                rect = patches.Polygon(patch_xy, closed=False, linewidth=2, 
                                         edgecolor=patch_color,facecolor='none')
            if patch_style is "polygon":
                # Create a Rectangle patch
                rect = patches.Polygon(patch_xy, closed=True, linewidth=2, 
                                         edgecolor=patch_color,facecolor='none')

            # Add the patch to the Axes
            ax.add_patch(rect)
        
        fig.colorbar(s, ax=ax)
        ax.set_xlabel('µm')
        ax.set_ylabel('µm')
        ax.set_title(key)
        ax.set_facecolor('k')
        ax.invert_yaxis()
        ax.grid(False)
       
    if header is not None:
        plt.suptitle(header, fontsize=headersize, x=header_x, y=header_y)
    if save:
        plt.savefig(savepath, dpi=dpi_save, facecolor=save_background, bbox_inches='tight')
        
    if inline:
        return plt.show()
    else:
        return fig, ax

def celltypes(adata, keys, max_cols=None, df_celltypes=None, save=False, header=None, headersize=18, raw=False,
                       plot_names = None, savepath="figures/spatial_genes.png", save_background=None, 
                       dpi_save=300, resolution=50, 
                       xlim = None, ylim = None, oversize=1, dpi_display=80, image=None, image_metadata=None, alpha=1,
                       header_x=0.5, header_y=0.98):
        
    # Get input and parameters
        
    if df_celltypes is None:
        df_celltypes = adata.uns["stereoscope"]
        
    assert len(df_celltypes) == len(adata), "Celltype list and adata object have not the same length."
    
    if raw:
        adata_var = adata.raw.var
        adata_X = adata.raw.X
    else:
        adata_var = adata.var
        adata_X = adata.X
        
    if plot_names is not None:
        assert len(plot_names) == len(keys)
        
    if xlim is None:
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            x_origin = image_metadata["upper_left_point"][0]/pixel_per_um + resolution/2
        else:
            x_origin = 0
        xmin = adata.obs["um_col"].min()
        xmax = adata.obs["um_col"].max()
        xlim = (xmin - resolution, xmax + resolution)
    else:
        xlim[0] -= resolution
        xlim[1] += resolution
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            x_origin = image_metadata["upper_left_point"][0]/pixel_per_um + resolution/2
        else:
            x_origin = 0
        
    if ylim is None:
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            y_origin = image_metadata["upper_left_point"][1]/pixel_per_um + resolution/2
        else:
            y_origin = 0
        ymin = adata.obs["um_row"].min()
        ymax = adata.obs["um_row"].max()
        ylim = (ymin - resolution, ymax + resolution)
    else:
        ylim[0] -= resolution
        ylim[1] += resolution
        if image is not None:
            pixel_per_um = image_metadata["pixel_per_um"]
            y_origin = image_metadata["upper_left_point"][1]/pixel_per_um + resolution/2
        else:
            y_origin = 0
    
    if isinstance(keys, str):
        n_plots = 1
        key = keys
    elif isinstance(keys, list) and len(keys) > 1:
        n_plots = len(keys)
    elif isinstance(keys, list) and len(keys) == 1:
        n_plots = len(keys)
        key = keys[0]
    else:
        print('Keys have unknown type')
        
    if max_cols != None:
        if n_plots > max_cols:
            n_rows = math.ceil(n_plots / max_cols)
        else:
            n_rows = 1
    else:
        max_cols = n_plots
        n_rows = 1
        
    # Plotting
    fig, ax = plt.subplots(n_rows, max_cols, figsize=(7.6*max_cols, 6*n_rows), dpi=dpi_display)
    
    if n_plots > 1:
        ax = ax.ravel()
        for i, key in enumerate(keys):
            # locate cell type in dataframe and extract values
            color = df_celltypes.loc[:, key] * 100   
            
            # calculate marker size
            ax[i].set_xlim(xlim[0], xlim[1])
            ax[i].set_ylim(ylim[0], ylim[1])
            pixels_per_unit = ax[i].transData.transform([(0,1),(1,0)])-ax[i].transData.transform((0,0))
            x_ppu = pixels_per_unit[1, 0]
            y_ppu = pixels_per_unit[0, 1]
            pxs = y_ppu * resolution * oversize
            size = (72. / fig.dpi * pxs)**2

            # plot
            if image is None:
                s = ax[i].scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', s=size)
            else:
                s = ax[i].scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', 
                               s=size, alpha=alpha, edgecolors=None)
                ax[i].imshow(image, 
                          extent=[0 - x_origin, image.shape[1]/pixel_per_um - x_origin, 
                                  0 - y_origin, image.shape[0]/pixel_per_um - y_origin], 
                          origin='lower'
                         )
            clb = fig.colorbar(s, ax=ax[i])
            clb.ax.set_title('%')
            ax[i].set_xlabel('µm')
            ax[i].set_ylabel('µm')
            if plot_names == None:
                ax[i].set_title(key)
            else:
                ax[i].set_title(plot_names[i])
            ax[i].set_facecolor('k')
            ax[i].invert_yaxis()
            ax[i].grid(False)
    else:
        # locate cell type in dataframe and extract values
        color = df_celltypes.loc[:, key] * 100
        
        # calculate marker size
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        pixels_per_unit = ax.transData.transform([(0,1),(1,0)])-ax.transData.transform((0,0))
        x_ppu = pixels_per_unit[1, 0]
        y_ppu = pixels_per_unit[0, 1]
        pxs = y_ppu * resolution * oversize
        size = (72. / fig.dpi * pxs)**2

        # plot
        if image is None:
            s = ax.scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', s=size)
        else:            
            # plot data
            s = ax.scatter(adata.obs.um_col, adata.obs.um_row, c=color, marker='s', 
                           s=size, alpha=alpha, edgecolors=None)
            ax.imshow(image, 
                      extent=[0 - x_origin, image.shape[1]/pixel_per_um - x_origin, 
                              0 - y_origin, image.shape[0]/pixel_per_um - y_origin], 
                      origin='lower'
                     )
        clb = fig.colorbar(s, ax=ax)
        clb.ax.set_title('%')
        ax.set_xlabel('µm')
        ax.set_ylabel('µm')
        ax.set_title(keys)
        ax.set_facecolor('k')
        ax.invert_yaxis()
        ax.grid(False)
       
    if header is not None:
        plt.suptitle(header, fontsize=headersize, x=header_x, y=header_y)
    if save:
        plt.savefig(savepath, dpi=dpi_save, facecolor=save_background, bbox_inches='tight')
    
    return plt.show()