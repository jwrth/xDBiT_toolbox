import numpy as np
import matplotlib.pyplot as plt
from ..tools import get_nrows_maxcols, extract_groups
from .spatial_class import MultiSpatialPlot as spatial

def view_rgbs(results, max_cols=4, min=None, max=None):
    '''
    Function to plot results from `napari_to_rgb()`.
    '''
    # plotting
    nplots, nrows, ncols = get_nrows_maxcols(results, max_cols=max_cols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6*nrows))

    axs = [axs] if nplots == 1 else axs

    for i, im in enumerate(results):
        if min is not None or max is not None:
            im = np.clip(im, a_min=min, a_max=max)
        axs[i].imshow(im)
        axs[i].set_title(i)

    fig.tight_layout()
    plt.show()
    
def view_spots_rgb(adata, genes, regions, region_id, metadata, alpha=0.7):
    '''
    Function to view image with overlayed spots in region selected by napari.
    '''

    xlim = regions['xlims'][region_id]
    ylim = regions['ylims'][region_id]
    coords = adata.obsm['spatial']
    mask = ((coords[:, 0] > xlim[0]) & (coords[:, 0] < xlim[1])) & ((coords[:, 1] > ylim[0]) & (coords[:, 1] < ylim[1]))

    subset_crop = extract_groups(adata, groupby=None, strip=True)[mask, :].copy()

    # offset coordinates
    coords = subset_crop.obsm['spatial'] # extract coordinates from cropped adata
    coords[:, 0] -= xlim[0]
    coords[:, 1] -= ylim[0]

    # add image data to cropped adata
    im = regions['images'][region_id]
    imkey = "rgb"

    # add categories
    subset_crop.uns['spatial'] = {}
    subset_crop.uns['spatial'][imkey] = {}
    subset_crop.uns['spatial'][imkey]['images'] = {}
    subset_crop.uns['spatial'][imkey]['scalefactors'] = {}

    # add image
    subset_crop.uns['spatial'][imkey]['images']['hires'] = im

    # add scale and positioning information
    subset_crop.uns['spatial'][imkey]['scalefactors'] = metadata
    subset_crop.uns['spatial'][imkey]['scalefactors']['xlim'] = regions['xlims'][region_id]
    subset_crop.uns['spatial'][imkey]['scalefactors']['ylim'] = regions['ylims'][region_id]

    scalefactors = subset_crop.uns['spatial']['rgb']['scalefactors']
    ppm = scalefactors['pixel_per_um_real']

    spatial(subset_crop, 
                keys=genes, 
                groupby=None, 
                image_key="rgb", 
                lowres=False,
                plot_pixel=False,
                #xlim=(0, xlim[1]-xlim[0]), ylim=(0, ylim[1]-ylim[0])
                xlim=(0, (xlim[1]-xlim[0])/ppm), ylim=(0, (ylim[1]-ylim[0])/ppm),
                origin_zero=False, margin=False,
                alpha=alpha
                )