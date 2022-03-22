import numpy as np
import cv2
from ..calculations._calc import order_points_clockwise, dist_points
#import imutils
import matplotlib.pyplot as plt
from datetime import datetime
from ..tools import extract_groups, rotatePoint
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
    return (img)


def align_image(image, vertices, frame: int = 100, return_grayscale=True):
    '''
    Function to align image based on four vertices using perspective transformation.
    '''

    # order vertices clockwise
    vertices = order_points_clockwise(vertices)

    # calculate the width of the alignment marker square
    pxl_width = int(dist_points(vertices[0], vertices[1]).astype(int))

    # x/y coordinates need to be flipped for opencv
    pts1 = np.flip(vertices, axis=1).astype(np.float32)
    pts2 = np.float32([[frame, frame], [pxl_width + frame, frame],
                       [pxl_width + frame, pxl_width + frame], [frame, pxl_width + frame]])

    # calculate transformation matrix
    M = cv2.getPerspectiveTransform(pts1, pts2)

    # align image
    aligned = cv2.warpPerspective(
        image, M, (int(pxl_width + 2 * frame), int(pxl_width + 2 * frame)))

    if return_grayscale:
        if len(aligned.shape) == 3:
            print("Convert image to grayscale...")
            aligned = cv2.cvtColor(aligned, cv2.COLOR_BGR2GRAY)

    return aligned


def align_to_dict(images, labels, vertices, resolution: int, n_channels: int,
                  frame: int = 100, bit_type='8bit'):

    # Part 1: Align image
    print("     Align images...")
    aligned = [align_image(img, vertices, frame) for img in images]

    # Transform to RGB
    aligned = [single_grayscale_to_rgb(
        img, bit_type=bit_type) for img in aligned]

    # Part 2: Create metadata
    print("     Create metadata...")
    # calculate the width of the alignment marker square
    pxl_width = dist_points(vertices[0], vertices[1]).astype(int)
    pts2 = np.float32([[frame, frame], [pxl_width + frame, frame],
                       [pxl_width + frame, pxl_width + frame], [frame, pxl_width + frame]])

    # calculate metadata
    um_width = ((n_channels - 1) * resolution * 2)
    pixel_per_um = pxl_width / um_width
    spot_diameter = pixel_per_um * resolution

    # summarize metadata in dictionary
    image_metadata = {
        # coordinates need to be converted to int
        "upper_left_spot_coord": list(map(int, pts2[0])),
        # this parameter is set to 60 since this is used by squidpy for plotting and gives best results.
        "spot_diameter_fullres": 60,
        "spot_diameter_real": spot_diameter,
        "pixel_per_um": pixel_per_um,
        "resolution": resolution,
        "tissue_hires_scalef": 1.0
    }

    # Part 3: Summarize aligned image and metadata
    print("     Summarized aligned images and metadata.")
    image_sum = {lab: {'images': {'hires': img}, 'scalefactors': image_metadata}
                 for (lab, img) in zip(labels, aligned)}

    return image_sum

def resize_image(img, scale_factor):
    '''
    Resize image by scale_factor
    '''
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

        adata.uns['spatial'][key]['images']['lowres'] = resize_image(img, scale_factor)
        adata.uns['spatial'][key]['scalefactors']['tissue_lowres_scalef'] = scale_factor


def register_image(image, template, maxFeatures=500, keepFraction=0.2, scale_factor=1,
                   debug=False, method="sift", ratio_test=True, flann=True, do_registration=True,
                   return_grayscale=True):

    if len(image.shape) == 3:
        print("Convert image to grayscale...")
        image_scaled = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    if len(template.shape) == 3:
        print("Convert template to grayscale...")
        template_scaled = cv2.cvtColor(template, cv2.COLOR_BGR2GRAY)

    if scale_factor < 1:
        print("Scale images before registration")
        image_scaled = resize_image(img=image, scale_factor=scale_factor)
        template_scaled = resize_image(img=template, scale_factor=scale_factor)
    else:
        image_scaled = image
        template_scaled = template

    print("{}: Get features...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    # Get features
    if method == "sift":
        print("     Method: SIFT...")
        # sift
        sift = cv2.SIFT_create()

        (kpsA, descsA) = sift.detectAndCompute(image_scaled, None)
        (kpsB, descsB) = sift.detectAndCompute(template_scaled, None)

    elif method == "surf":
        print("     Method: SURF...")
        surf = cv2.xfeatures2d.SURF_create(400)

        (kpsA, descsA) = surf.detectAndCompute(image_scaled, None)
        (kpsB, descsB) = surf.detectAndCompute(template_scaled, None)

    else:
        print("Unknown method. Aborted.")
        return

    if flann:
        print("{}: Compute matches...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        # FLANN parameters
        FLANN_INDEX_KDTREE = 1
        index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
        search_params = dict(checks=50)   # or pass empty dictionary

        # runn Flann matcher
        flann = cv2.FlannBasedMatcher(index_params, search_params)
        matches = flann.knnMatch(descsA, descsB, k=2)

    else:
        print("{}: Compute matches...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        # feature matching
        #bf = cv2.BFMatcher(cv2.NORM_L1, crossCheck=True)
        bf = cv2.BFMatcher()
        matches = bf.knnMatch(descsA, descsB, k=2)

    if ratio_test:
        print("{}: Filter matches...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        # store all the good matches as per Lowe's ratio test.
        good_matches = []
        for m, n in matches:
            if m.distance < 0.7*n.distance:
                good_matches.append(m)
    else:
        print("{}: Filter matches...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        # sort the matches by their distance (the smaller the distance, the "more similar" the features are)
        matches = sorted(matches, key=lambda x: x.distance)
        # keep only the top matches
        keep = int(len(matches) * keepFraction)
        good_matches = matches[:keep]

        print("Number of matches used: {}".format(len(good_matches)))

    # elif method=="orb":
    #     # use ORB to detect keypoints and extract (binary) local
    #     # invariant features
    #     orb = cv2.ORB_create(maxFeatures)
    #     (kpsA, descsA) = orb.detectAndCompute(image, None)
    #     (kpsB, descsB) = orb.detectAndCompute(template, None)

    #     # match the features
    #     method = cv2.DESCRIPTOR_MATCHER_BRUTEFORCE_HAMMING
    #     matcher = cv2.DescriptorMatcher_create(method)
    #     matches = matcher.match(descsA, descsB, None)

    #     # sort the matches by their distance (the smaller the distance,
    #     # the "more similar" the features are)
    #     matches = sorted(matches, key=lambda x:x.distance)
    #     # keep only the top matches
    #     keep = int(len(matches) * keepFraction)
    #     good_matches = matches[:keep]

    # check to see if we should visualize the matched keypoints
    if debug:
        print("{}: Debugging mode - Display matches...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        matchedVis = cv2.drawMatches(image_scaled, kpsA, template_scaled, kpsB,
                                     good_matches, None)
        #matchedVis = imutils.resize(matchedVis, width=1000)
        matchedVis = resize_image(matchedVis, scale_factor=0.1)
        plt.imshow(matchedVis)
        plt.show()
    
    # Compute homography matrix

    print("{}: Fetch keypoints...".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    # allocate memory for the keypoints (x, y)-coordinates from the
    # top matches -- we'll use these coordinates to compute our
    # homography matrix
    ptsA = np.zeros((len(good_matches), 2), dtype="float")
    ptsB = np.zeros((len(good_matches), 2), dtype="float")
    # loop over the top matches
    for (i, m) in enumerate(good_matches):
        # indicate that the two keypoints in the respective images
        # map to each other
        ptsA[i] = kpsA[m.queryIdx].pt
        ptsB[i] = kpsB[m.trainIdx].pt

    if debug:
        print("{}: Debugging mode - Display image and template with keypoints...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        # plot keypoints for image
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        axs[0].imshow(image_scaled)
        axs[0].scatter(x=ptsA[:, 0], y=ptsA[:, 1])
        axs[0].set_title('Image')

        axs[1].imshow(template_scaled)
        axs[1].scatter(x=ptsB[:, 0], y=ptsB[:, 1])
        axs[1].set_title('Template')

        plt.show()

    # compute the homography matrix between the two sets of matched
    # points
    print("{}: Compute homography matrix...".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}"))

    # apply scale_factor to points
    ptsA /= scale_factor
    ptsB /= scale_factor

    # determine homography matrix
    (H, mask) = cv2.findHomography(ptsA, ptsB, method=cv2.RANSAC)

    # use the homography matrix to register the images
    (h, w) = template.shape[:2]
    if do_registration:
        print("{}: Register image...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))

        # convert to 8-bit and scale
        image = convert_to_8bit(image)

        # warping
        registered = cv2.warpPerspective(image, H, (w, h))

        if return_grayscale:
            if len(registered.shape) == 3:
                print("Convert registered image to grayscale...")
                registered = cv2.cvtColor(registered, cv2.COLOR_BGR2GRAY)
    else:
        registered = None

    # return the registered image
    return registered, H


def recalculate_scale(adata, groupby, group, spatial_key='spatial', 
    save_scale=True, return_angle_and_pivot=False):
    '''
    Recalculates the scale `pixel_per_um` of all images of a given `group` in an anndata object.
    Expects the images in `adata.uns[spatial_key]`
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
            break
        else:
            # if not both spots exists switch into next column
            min_col += 1

    if return_angle_and_pivot:
        # calculate and return the rotation angle of the spots
        rot_angle = rotation_angle(upper_spot_px, lower_spot_px)
        return rot_angle, upper_spot_px


def register_adata_coords_to_new_images(adata_in, groupby, image_dir_dict, groups=None, reg_channel='dapi',
                                        spatial_key='spatial', 
                                        hires_key='hires', lowres_key='lowres', scale_factor_before_reg=None,
                                        rot_threshold=0, do_registration=False, scale_factor=0.1,
                                        keepFraction=0.2, method='sift', debug=False, in_place=False):
    '''
    Function to add new images to an adata object and register the coordinates accordingly.

    The images are expected to be in adata.uns[spatial_key][image_id].
    The image_id is expected to be in the format "{well}_{channel}".
    '''
    if not in_place:
        adata = adata_in.copy()

    if groups is None:
        groups = list(adata.obs[groupby].unique())

    groups = [groups] if isinstance(groups, str) else list(groups)

    for group in groups:
        print("{}: Process group {}...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", group))

        # extract group from adata
        adata_subset, obs_mask = extract_groups(
            adata, groupby=groupby, groups=group, extract_uns=True, return_mask=True)

        # fetch name of registration channel
        image_keys = list(adata_subset.uns[spatial_key].keys())
        reg_key = [k for k in image_keys if reg_channel in k]
        #other_keys = [k for k in image_keys if reg_channel not in k]
        if len(reg_key) == 1:
            reg_key = reg_key[0]
        else:
            raise AssertionError(
                "No unique registration channel found: {}".format(reg_key))

        # extract image from subset
        image_adata = adata_subset.uns[spatial_key][reg_key]['images'][hires_key]
        
        # extract image metadata
        lowres_metadata = adata_subset.uns[spatial_key][reg_key]['scalefactors']
        pixel_per_um = lowres_metadata['pixel_per_um']

        # load images and convert to grayscale if necessary
        hq_image_dict = {}
        #for key in image_keys:
        for key in image_dir_dict:
            print("{}: Load image for key {}...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}", key))
            hq_image_dict[key] = cv2.imread(image_dir_dict[key], 0)

            # convert to grayscale
            if len(hq_image_dict[key].shape) == 3:
                hq_image_dict[key] = cv2.cvtColor(hq_image_dict[key], cv2.COLOR_BGR2GRAY)

        # extract the registration image from the dict
        image_to_register = hq_image_dict[reg_key]

        if debug:
            # scale down the images
            image_adata = resize_image(image_adata, scale_factor=0.1)
            image_to_register = resize_image(image_to_register, scale_factor=0.1)

            # image_adata = imutils.resize(
            #     image_adata, width=int(image_adata.shape[1]*0.1))
            # image_to_register = imutils.resize(
            #     image_to_register, width=int(image_to_register.shape[1]*0.1))


        # register images and extract homography matrix
        print("{}: Register image {}...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", reg_key))
        registered_img, H = register_image(image_adata, image_to_register, do_registration=do_registration, scale_factor=scale_factor_before_reg,
                              keepFraction=keepFraction, method=method, debug=debug)

        if do_registration:
            # save registered image in adata
            print("Save registered image in adata...")

            # convert registered image to grayscale if necessary
            if len(registered_img.shape) == 3:
                    registered_img = cv2.cvtColor(registered_img, cv2.COLOR_BGR2GRAY)

            if 'registered' in adata.uns.keys():
                adata.uns['registered'][group] = registered_img
            else:
                adata.uns['registered'] = {}
                adata.uns['registered'][group] = registered_img

        # transform transcriptome coordinates using homography matrix
        # extract coordinates from subset
        coords = adata_subset.obsm['spatial']

        # reshape and transform
        print("{}: Perspective transformation of coordinates...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
            
        coords_reshaped = coords.reshape(-1, 1, 2)
        coords_trans = cv2.perspectiveTransform(coords_reshaped, H)
        coords_trans = coords_trans.reshape(-1, 2)

        # recalculate scale `pixel_per_um` in each image of this group and get rotation angle and pivot point
        print("{}: Calculate rotation angle...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        rot_angle, pivot_point = recalculate_scale(
            adata, groupby=groupby, group=group, save_scale=False, return_angle_and_pivot=True)

        # correct for rotational shifting of coordinates and image
        # image rotation
        if np.absolute(rot_angle) > rot_threshold:
            for key in hq_image_dict:
                print("{}: Rotate image {} by {} degrees...".format(
                    f"{datetime.now():%Y-%m-%d %H:%M:%S}", key, rot_angle))
                image_rotated = rotateImage(
                    hq_image_dict[key], angle=-rot_angle, pivot=pivot_point, PIL=True)
                hq_image_dict[key] = image_rotated

            # coordinate rotation
            print("{}: Rotate coordinates of group {} by {} degrees...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}", group, rot_angle))
            coords_transrot = np.array([rotatePoint(
                origin=pivot_point, point=elem, angle=-rot_angle) for elem in coords_trans])
        else:
            print("Rotation angle {} below threshold {}. Images and coordinates not rotated.".format(
                rot_angle, rot_threshold))
            coords_transrot = coords_trans

        # delete old images and store new images in adata
        print("{}: Store images and coordinates of group {} in anndata object...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", group))
        
        # get metadata from first key
        metadata = adata.uns['spatial'][image_keys[0]]['scalefactors'].copy()

        # delete old images
        for key in image_keys:
            adata.uns[spatial_key].pop(key, None) # remove key regardless of whether it is in the dict or not
            #del adata.uns[spatial_key][key]
            # res_keys = list(adata.uns[spatial_key][key]['images'].keys()).copy()
            # for res_key in res_keys:
            #     del adata.uns[spatial_key][key]['images'][res_key]
        
        # store new images in adata and add metadata
        for key in hq_image_dict:
            if key not in adata.uns[spatial_key]:
                adata.uns[spatial_key][key] = {}
                adata.uns[spatial_key][key]['images'] = {}
            
            # add image
            adata.uns[spatial_key][key]['images'][hires_key] = hq_image_dict[key]
            # add metadata
            adata.uns[spatial_key][key]['scalefactors'] = metadata

        # store transformed and rotated coordinates in adata
        adata.obsm['spatial'][obs_mask] = coords_transrot

        # substitute the obs coordinates with the .obsm coordinates
        adata.obs['pixel_row'] = adata.obsm['spatial'][:, 1]
        adata.obs['pixel_col'] = adata.obsm['spatial'][:, 0]
    
        # recalculate scale and save in `scalefactors`
        print("{}: Recalculate scale...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        recalculate_scale(adata, groupby=groupby, group=group, 
            save_scale=True, return_angle_and_pivot=False)

    # resize images and store as `lowres` for plotting
    print("{}: Resize all images and store as `lowres`...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    resize_images_in_adata(adata, scale_factor=scale_factor)

    if not in_place:
        return adata


def calc_image_param_per_spot(adata, groupby='well_name', channel_pattern='dapi',
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