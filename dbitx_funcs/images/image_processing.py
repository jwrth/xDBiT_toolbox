import numpy as np
import cv2
from ..calculations._calc import dist_points
from ..tools import extract_groups
from ..calculations import dist_points, rotation_angle
from sklearn.preprocessing import minmax_scale
import numpy as np
from scipy import ndimage
from typing import Optional, Tuple, Union, List, Dict, Any
from PIL import Image

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

def convert_to_8bit(img):
    '''
    Convert numpy array image to 8bit.
    '''
    img = (img / img.max()) * 255
    img = np.uint8(img)
    return img


def resize_image(img, dim=None, scale_factor=None):
    '''
    Resize image by scale_factor
    '''
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
                              normalize=True
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