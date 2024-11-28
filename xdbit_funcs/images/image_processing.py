import numpy as np
import cv2
from ..calculations._calc import dist_points
from ..tools import extract_groups
from ..calculations import dist_points, rotation_angle
from sklearn.preprocessing import minmax_scale
import numpy as np
from scipy import ndimage
from typing import Optional, Tuple, Union, List, Dict, Any, Literal, Callable
from PIL import Image
from skimage.color import rgb2hed, hed2rgb
from ..datasets import ImageData
from anndata import AnnData
from pathlib import Path
import pandas as pd
import warnings

def rotateImage(img, angle, pivot, imagetype="grayscale", PIL=True, order=3):
    """
    Rotate counterclockwise by given angle in degrees around pivot point (format [x,y]).
    """

    if PIL:
        # use Python Image Processing Library
        img = Image.fromarray(img)
        imgR = np.array(img.rotate(angle, center=tuple(pivot)))

    else:
        # use scipy function
        padX = [int(img.shape[0] - pivot[0]), int(pivot[0])]
        padY = [int(img.shape[1] - pivot[1]), int(pivot[1])]

        if imagetype == "grayscale":
            imgP = np.pad(img, [padY, padX], 'constant')
        elif imagetype == "rgb":
            imgP = np.pad(img, [padY, padX, [0, 0]], 'constant')
        else:
            print("Unknown image type.")
            return
        
        imgR = ndimage.rotate(imgP, angle, reshape=False, order=order)
        imgR = imgR[padY[0]: -padY[1], padX[0]: -padX[1]]

    return imgR


def set_histogram(data, lower=0, upper=None, bit_type=np.uint8, clip=True):
    '''
    Set histogram of image.
        data: image array
        lower: lower threshold for histogram
        upper: upper threshold for histogram
        bit_type: np.uint8 or np.uint16. Default is np.uint8.
        clip: if True histogram is clipped to upper and lower threshold. If not image is only transformed to new bit_type.
    '''

    if bit_type is np.uint8:
        max_int = 255
    elif bit_type is np.uint16:
        max_int = 65535
    else:
        print("Unknown bit type.")
        return

    if lower is None:
        lower = np.min(data)

    if upper is None:
        upper = np.max(data)

    norm = ((data - lower) / (upper - lower))

    if clip:
        norm = np.clip(norm, a_min=0,  a_max=1)

    norm *= max_int

    return norm.astype(bit_type)


def single_grayscale_to_rgb(image, bit_type="8bit", lower=None, upper=None):
    '''
    Function to transform single grayscale image into a rgb image.
    '''
    if bit_type == "8bit":
        bit_type = np.uint8
    elif bit_type == "16bit":
        bit_type = np.uint16
    else:
        print("Unknown bit type.")
        return

    shape = image.shape

    image = set_histogram(image, bit_type=bit_type,
                          lower=lower, upper=upper)

    rgb = cv2.merge((image, image, image))

    return rgb


def multi_grayscale_to_rgb(r=None, g=None, b=None, bit_type="8bit", lowers=[None] * 3, uppers=[None] * 3):
    '''
    Function to transform multiple grayscale images into a rgb image.
    '''
    if bit_type == "8bit":
        bit_type = np.uint8
    elif bit_type == "16bit":
        bit_type = np.uint16
    else:
        print("Unknown bit type.")
        return

    if r is not None:
        shape = r.shape
    elif g is not None:
        shape = g.shape
    elif b is not None:
        shape = b.shape
    else:
        print("All channels empty.")
        return

    if r is None:
        r = np.zeros(shape, dtype=bit_type)
    else:
        r = set_histogram(r, bit_type=bit_type,
                          lower=lowers[0], upper=uppers[0])

    if g is None:
        g = np.zeros(shape, dtype=bit_type)
    else:
        g = set_histogram(g, bit_type=bit_type,
                          lower=lowers[1], upper=uppers[1])

    if b is None:
        b = np.zeros(shape, dtype=bit_type)
    else:
        b = set_histogram(b, bit_type=bit_type,
                          lower=lowers[2], upper=uppers[2])

    rgb = cv2.merge((r, g, b))

    return rgb

def convert_to_8bit(img, save_mem=True):
    '''
    Convert numpy array image to 8bit.
    '''
    if save_mem:
        # for a 16-bit image at least int32 is necessary for signed integers because the value range is [-65535,...,0,...,65535]
        # or uint16 can be used as unsigned integer with only positive values
        img = np.uint16(img)
    img = (img / img.max()) * 255
    img = np.uint8(img)
    return img


def resize_image(img, dim=None, scale_factor=None):
    '''
    Resize image by scale_factor
    '''
    # make sure the image is np.uint8
    img = img.astype(np.uint8)
    
    if scale_factor is not None:
        width = int(img.shape[1] * scale_factor)
        height = int(img.shape[0] * scale_factor)
        dim = (width, height)

    return cv2.resize(img, dim)
        


def resize_images_in_adata(adata, scale_factor):
    '''
    Resizes images of all img_keys in an adata by a certain `scale_factor` and saves them in `adata.uns['spatial'][img_key]['images']['lowres']`.
    Key `tissue_lowres_scalef` is added to image metadata.
    '''

    img_keys = adata.uns['spatial'].keys()

    for key in img_keys:
        img = adata.uns['spatial'][key]['images']['hires']

        adata.uns['spatial'][key]['images']['lowres'] = resize_image(img, scale_factor=scale_factor)
        adata.uns['spatial'][key]['scalefactors']['tissue_lowres_scalef'] = scale_factor


def recalculate_scale(adata, groupby, group, ppm_given=None, spatial_key='spatial', 
    save_scale=True, return_angle_and_pivot=False):
    '''
    Recalculates the scale "pixel_per_um" of all images of a given `group` in an anndata object.
    Expects the images in `adata.uns[spatial_key]`.
    '''

    a = extract_groups(adata, groupby=groupby, groups=group, extract_uns=True)

    # find first row
    min_row = a.obs['array_row'].min()
    # find last row
    max_row = a.obs['array_row'].max()

    min_col = 0

    # find two orthogonal spots with maximal distance
    while True:
        # extract spot in first row
        spot1 = a.obs[(a.obs['array_row'] == min_row) &
                      (a.obs['array_col'] == min_col)]

        # extract spot in last row
        spot2 = a.obs[(a.obs['array_row'] == max_row) &
                      (a.obs['array_col'] == min_col)]

        # check if both spots exist
        if (len(spot1) == 1) & (len(spot2) == 1):
            # extract pixel and array coordinates in shape (x,y)
            upper_spot_px = spot1[['pixel_col', 'pixel_row']].values[0]
            lower_spot_px = spot2[['pixel_col', 'pixel_row']].values[0]

            upper_spot_ar = spot1[['array_col', 'array_row']].values[0]
            lower_spot_ar = spot2[['array_col', 'array_row']].values[0]

            if save_scale:
                # fetch images keys
                keys = a.uns[spatial_key].keys()

                for key in keys:
                    # get resolution of current image
                    res = a.uns[spatial_key][key]['scalefactors']['resolution']

                    # calculate distances in pixel and um
                    d_px = dist_points(upper_spot_px, lower_spot_px)
                    d_um = dist_points(upper_spot_ar, lower_spot_ar) * res * 2

                    pixel_per_um = d_px / d_um
                    
                    adata.uns[spatial_key][key]['scalefactors']['pixel_per_um'] = pixel_per_um
                    adata.uns[spatial_key][key]['scalefactors']['spot_diameter_real'] = res * pixel_per_um
                    adata.uns[spatial_key][key]['scalefactors']['pixel_per_um_real'] = ppm_given
            break
        else:
            # if not both spots exists switch into next column
            min_col += 1

    if return_angle_and_pivot:
        # calculate and return the rotation angle of the spots
        rot_angle = rotation_angle(upper_spot_px, lower_spot_px)
        return rot_angle, upper_spot_px

def calc_image_param_per_spot(adata, groupby='id', channel_pattern='dapi',
                              fun=np.mean, fun_descriptor='mean', lowres=False,
                              normalize=True, 
                              ):
    '''
    Function to apply a function spotwise to an image in a ST dataset. 
    The dataset is expected to be in anndata format.

    The image is expected in `adata.uns['spatial'][img_key]['images'][resolution_key]` where `img_key` consists
    of '{group}_{channel_pattern}'.
    In most cases the group is the well_name which is extracted from `adata.obs[groupby]` with groupby='well_name'.

    Spot coordinates are expected to be stored in adata.obsm['spatial'] and should match the order of `adata.obs`.

    The function to be applied can be specified with `fun`. Default is `np.mean`.

    Changes to anndata object are made in-place.

    '''

    if lowres:
        resolution_key = 'lowres'
    else:
        resolution_key = 'hires'

    results_list = []
    old_group = None
    for i, (index, row) in enumerate(adata.obs.iterrows()):

        # extract group
        group = row[groupby]

        if group != old_group:
            # search for correct img_key
            img_key = [elem for elem in adata.uns['spatial'].keys() if (
                channel_pattern in elem) & (group in elem)]

            if len(img_key) == 1:
                if lowres:
                    scalef = adata.uns['spatial'][img_key[0]
                                                ]['scalefactors']['tissue_lowres_scalef']
                else:
                    scalef = adata.uns['spatial'][img_key[0]
                                                ]['scalefactors']['tissue_hires_scalef']

                img = adata.uns['spatial'][img_key[0]
                                        ]['images'][resolution_key]
                px_dia = adata.uns['spatial'][img_key[0]
                                            ]['scalefactors']['spot_diameter_real'] * scalef

            else:
                print("No or no unique image key found for spot {}: {}".format(
                    index, img_key))
                img = None

        if img is not None:
            # extract spot
            spot = adata.obsm['spatial'][i] * scalef
            # determine region of spot
            region = img[int((spot[1] - px_dia / 2)): int((spot[1] + px_dia / 2)),
                         int((spot[0] - px_dia / 2)): int((spot[0] + px_dia / 2))]

            # apply function to region and record result
            if fun is None:
                results_list.append(region)
            else:
                results_list.append(fun(region))

        else:
            # if no unique image was found record NaN
            results_list.append(np.nan)

        # save group for next loop
        old_group = group    

    obshead = channel_pattern + '_' + fun_descriptor
    adata.obs[obshead] = results_list

    if normalize:
        adata.obs[obshead + '_norm'] = adata.obs[[groupby, obshead]].groupby(groupby).transform(minmax_scale)
        
def crop_center(image, center, width):
    '''
    Crop image from given center and width.
    '''
    return image[int(center[0] - width/2): int(center[0] + width/2), int(center[1] - width/2): int(center[1] + width/2)]

def deconvolve_he(
    img: np.ndarray,
    return_type: Literal["grayscale", "greyscale", "rgb"] = "grayscale"
    ) -> np.ndarray:
    '''
    Deconvolves H&E image to separately extract hematoxylin and eosin stainings.
    
    from: https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_ihc_color_separation.html
    
    --------------
    Returns:
    For return_type "grayscale": Numpy array with shape (h, w, 3) where the 3 channels correspond to 
        hematoxylin, eosin and a third channel
        
    For return_type "rgb": Three separate RGB images as numpy array. Order: Hematoxylin, eosin, and third green channel.
    '''
    # perform deconvolution
    ihc_hed = rgb2hed(img)
    
    if return_type in ["grayscale", "greyscale"]:
        return ihc_hed
    
    else:
        # Create an RGB image for each of the stains
        null = np.zeros_like(ihc_hed[:, :, 0])
        ihc_h = hed2rgb(np.stack((ihc_hed[:, :, 0], null, null), axis=-1))
        ihc_e = hed2rgb(np.stack((null, ihc_hed[:, :, 1], null), axis=-1))
        ihc_d = hed2rgb(np.stack((null, null, ihc_hed[:, :, 2]), axis=-1))

        return ihc_h, ihc_e, ihc_d
    
def calc_image_param_per_spot_zarr(
    adata: AnnData,
    groupby: str = "id",
    channel_pattern: str = "nephrin",
    lowres: bool = True,
    lowres_sf: float = 0.1,
    zarr_key: str = "zarr_directories",
    fun: Callable = np.mean, # function. If None the region image will be extracted.
    fun_descriptor: str = "mean",
    normalize: bool = True,
    overwrite: bool = False,
):
    '''
    Function to apply a function spotwise to an image in a ST dataset. 
    The dataset is expected to be in anndata format.

    The image is expected to be stored as zarr in adata.uns['zarr_directories'].

    Spot coordinates are expected to be stored in adata.obsm['spatial'] and should match the order of `adata.obs`.

    The function to be applied can be specified with `fun`. Default is `np.mean`.

    Changes to anndata object are made in-place.
    '''

    groups = adata.obs[groupby].unique()

    results = []
    indices = []
    for group in groups:
        # get mask of group
        mask = adata.obs[groupby] == group
        
        # select obs and coordinates
        obs = adata.obs[mask]
        coords = adata.obsm['spatial'][mask]
        
        # check if channel pattern is in current data of current group
        zarr_dir = Path(adata.uns[zarr_key][group][0])
        
        # get directory for given channel pattern
        img_dir = [elem for elem in zarr_dir.glob("*") if elem.is_dir() and channel_pattern in elem.stem] 
        
        # check if there was a unique channel pattern found
        if len(img_dir) == 1:
            img_dir = img_dir[0]
        elif len(img_dir) == 0:
            img_dir = None
        else:
            raise ValueError("More than one possible channel found for given channel pattern {}: {}".format(channel_pattern, img_dir))
        
        if img_dir is not None:
            # load image data
            ImgD = ImageData(adata=adata, image_key=channel_pattern, group=group, 
                                lowres=lowres, lowres_sf=lowres_sf, zarr_key=zarr_key
                                )
            img = ImgD.image
            scalef = ImgD.scale_factor
            px_dia = ImgD.image_metadata['spot_diameter_real'] * scalef
            
            # iterate through observations and do calculations per spot
            results_list = []
            indices_list = []
            for i, (index, row) in enumerate(obs.iterrows()):
                # extract spot
                spot = coords[i] * scalef
                # determine region of spot
                region = img[int((spot[1] - px_dia / 2)): int((spot[1] + px_dia / 2)),
                            int((spot[0] - px_dia / 2)): int((spot[0] + px_dia / 2))]
                
                # apply function to region and record result
                if fun is None:
                    results_list.append(region) # collect region image
                else:
                    results_list.append(fun(region)) # collect result of calculation
                    
                # collect index
                indices_list.append(index)
            
            # collect results and indices
            results = results + results_list
            indices = indices + indices_list

        else:
            # add NANs as results
            ImgD = None
            results = results + [np.nan] * len(obs)
            indices = indices + list(obs.index)

    # create name to be added to column
    obshead = "{}_{}".format(channel_pattern, fun_descriptor)

    # create series of results
    res = pd.Series(results, index=indices,
                    name=obshead
                )

    if obshead in adata.obs.columns:
        if overwrite:
            # drop column if it exists already
            adata.obs = adata.obs.drop(obshead, axis=1)
        else:
            raise AssertionError("Column {} exists already in `adata.obs`. If you still want to do the calculations use `overwrite=True`".format(obshead))

    # add results to adata.obs considering the indices
    adata.obs = pd.merge(left = adata.obs, right=res, left_index=True, right_index=True)

    if normalize:
        norm_obshead = obshead + "_norm"
        if norm_obshead in adata.obs.columns:
            if overwrite:
                # drop column if it exists already
                adata.obs = adata.obs.drop(norm_obshead, axis=1)
            else:
                raise AssertionError("Column {} exists already in `adata.obs`. If you still want to do the calculations use `overwrite=True`".format(norm_obshead))

        # due to the insertion of NaNs per group one group have only NaNs. Warnings are ignored.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            adata.obs[norm_obshead] = adata.obs[[groupby, obshead]].groupby(groupby).transform(minmax_scale)
        