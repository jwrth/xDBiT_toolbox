import napari
import numpy as np
import matplotlib.pyplot as plt
from ..tools import get_nrows_maxcols

def interactive(adata, images, genes, channel_axis=2, channel_names=None, scalebar=True, spatial_key='spatial'):
    '''
    Interactive viewer for xDbit data using napari.
    '''
    # get information on the image scale
    keys = list(adata.uns['spatial'].keys())
    scalefactors = adata.uns['spatial'][keys[0]]['scalefactors']
    ppm = scalefactors['pixel_per_um_real']
    res = scalefactors['resolution']

    if channel_names is None:
        channel_names = ["ch" + str(i) for i in range(images.shape[channel_axis])] # set default channel names

    img_scale = tuple([1/ppm] * 2) if scalebar else (1,1)
    
    # create viewer
    # load multichannel image in one line
    viewer = napari.view_image(
        images, 
        channel_axis=channel_axis,
        name=channel_names,
        colormap=["gray", "blue", "green"], 
        rgb=False,
        scale=img_scale
    )

    # get position of points
    points = np.flip(adata.obsm[spatial_key].copy(), axis=1) # switch x and y

    if scalebar:
        points /= ppm
        unit = "um"
    else:
        res *= ppm
        unit = "pixel"    

    if genes is not None:
        # add points layer
        points_layer = {}
        for i, gene in enumerate(genes):
            # get expression values
            geneid = adata.var_names.get_loc(gene)
            expr = adata.X[:, geneid]

            point_properties = {
                #'good_point': np.array([True]*len(expr)),
                'confidence': expr,
            }

            visible = False
            if i == 0:
                visible = True

            points_layer[gene] = viewer.add_points(points, name=gene,
                                                   properties=point_properties,
                                                   symbol='s',
                                                   size=res,
                                                   face_color='confidence',
                                                   face_colormap='viridis',
                                                   opacity=0.7,
                                                   visible=visible, 
                                                   edge_width=0
                                            )

    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = unit

    napari.run()
    return viewer

def napari_to_rgb(viewer, shape_layer_name="Shapes", alpha=0.7):
    '''
    Convert shape selection in napari into RGB by additive blending.
    '''
    # get shape and image layers
    shapelays = {elem.name: elem for elem in viewer.layers if isinstance(elem, napari.layers.shapes.Shapes)}
    imlays = {elem.name: elem for elem in viewer.layers if isinstance(elem, napari.layers.image.Image)}

    blended_results = []
    for shape in shapelays[shape_layer_name].data:
        # get x and y limits
        ylim = (int(np.min(shape[:, 0])), int(np.max(shape[:, 0])))
        xlim = (int(np.min(shape[:, 1])), int(np.max(shape[:, 1])))

        # create empty blended image
        blended = np.zeros((ylim[1] - ylim[0], xlim[1] - xlim[0]) + (4,))

        # select image
        ch = 'dapi'
        imlay = imlays[ch]

        for ch, imlay in imlays.items():
            # see if this channel is visible
            if imlay.visible != 0:
                # selecct image
                im = imlay.data

                # do cropping
                cropped = im[ylim[0]:ylim[1], xlim[0]:xlim[1]]

                # normalize data by clims and map colors
                clims = imlay.contrast_limits
                normalized = (cropped - clims[0]) / (clims[1] - clims[0])
                colormapped = imlay.colormap.map(normalized.flatten())
                colormapped = colormapped.reshape(normalized.shape + (4,))

                # do additive blending
                blended = blended + colormapped

        # # clip to range [0...255]
        # blended = blended / blended.max() * 255
        # blended = blended.astype(np.uint8)

        blended[..., 3] = 1 # set alpha channel to 1
        #blended[:, 3] = 1 # set alpha channel to 1

        # collect results in list
        blended_results.append(blended)
    
    return blended_results

def view_rgbs(results, max_cols=4, **kwargs):
    '''
    Function to plot results from `napari_to_rgb()`.
    '''

    # plotting
    nplots, nrows, ncols = get_nrows_maxcols(results, max_cols=max_cols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6*nrows))

    axs = [axs] if nplots == 1 else axs

    for i, im in enumerate(results):
        axs[i].imshow(im, **kwargs)

    fig.tight_layout()
    plt.show()