from cv2 import resize
import numpy as np
import cv2
from ..calculations._calc import order_points_clockwise, dist_points
from datetime import datetime
from ..tools import extract_groups, rotatePoint
from ..calculations import dist_points
import numpy as np
from .image_processing import resize_image, recalculate_scale, rotateImage, resize_images_in_adata, convert_to_8bit

def register_image(image, template, maxFeatures=500, keepFraction=0.2, maxpx=None,
                   method="sift", ratio_test=True, flann=True, 
                   perspective_transform=False, 
                   do_registration=True,
                   return_grayscale=True):

    if len(image.shape) == 3:
        print("Convert image to grayscale...")
        image_scaled = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    if len(template.shape) == 3:
        print("Convert template to grayscale...")
        template_scaled = cv2.cvtColor(template, cv2.COLOR_BGR2GRAY)

    # scale_factor = 0.2
    # if scale_factor < 1:
    #     print("Scale images before registration by factor {}".format(scale_factor))
    #     image_scaled = resize_image(img=image, scale_factor=scale_factor)
    #     template_scaled = resize_image(img=template, scale_factor=scale_factor)
    # else:
    #     image_scaled = image
    #     template_scaled = template

    # dim = (4000,4000)
    if maxpx is not None:
        if np.max(image.shape) > maxpx:
            dim_image = tuple([int(elem / np.max(image.shape) * maxpx) for elem in image.shape])
        else:
            dim_image = image.shape
        
        if np.max(template.shape) > maxpx:        
            dim_template = tuple([int(elem / np.max(template.shape) * maxpx) for elem in template.shape])
        else:
            dim_template = template.shape
            
        print("Rescale image to following dimensions: {}".format(dim_image))
        print("Rescale template to following dimensions: {}".format(dim_template))
        image_scaled = resize_image(img=image, dim=dim_image)
        template_scaled = resize_image(img=template, dim=dim_template)
        print("Dim of image: {}".format(image_scaled.shape))
        print("Dim of template: {}".format(template_scaled.shape))
    else:
        image_scaled = image
        template_scaled = template

    # convert and normalize images to 8bit for registration
    print("Convert scaled images to 8 bit")
    image_scaled = convert_to_8bit(image_scaled)
    template_scaled = convert_to_8bit(template_scaled)

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

    # check to see if we should visualize the matched keypoints
    print("{}: Display matches...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    matchedVis = cv2.drawMatches(image_scaled, kpsA, template_scaled, kpsB,
                                    good_matches, None)
    
    # Get keypoints
    print("{}: Fetch keypoints...".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    # allocate memory for the keypoints (x, y)-coordinates of the top matches
    ptsA = np.zeros((len(good_matches), 2), dtype="float")
    ptsB = np.zeros((len(good_matches), 2), dtype="float")
    # loop over the top matches
    for (i, m) in enumerate(good_matches):
        # indicate that the two keypoints in the respective images map to each other
        ptsA[i] = kpsA[m.queryIdx].pt
        ptsB[i] = kpsB[m.trainIdx].pt

    # calculate scale factors for x and y dimension for image and template
    x_sf_image = dim_image[0] / image.shape[0]
    y_sf_image = dim_image[1] / image.shape[1]
    x_sf_template = dim_template[0] / template.shape[0]
    y_sf_template = dim_template[1] / template.shape[1]

    # apply scale factors to points - separately for each dimension
    ptsA[:, 0] = ptsA[:, 0] / x_sf_image
    ptsA[:, 1] = ptsA[:, 1] / y_sf_image
    ptsB[:, 0] = ptsB[:, 0] / x_sf_template
    ptsB[:, 1] = ptsB[:, 1] / y_sf_template

    # # apply scale_factor to points
    # ptsA /= scale_factor
    # ptsB /= scale_factor

    if perspective_transform:
        # compute the homography matrix between the two sets of matched
        # points
        print("{}: Compute homography matrix...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        (H, mask) = cv2.findHomography(ptsA, ptsB, method=cv2.RANSAC)
    else:
        print("{}: Estimate 2D affine transformation matrix...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
        (H, mask) = cv2.estimateAffine2D(ptsA, ptsB)

    # use the homography matrix to register the images
    (h, w) = template.shape[:2]
    if do_registration:
        

        if perspective_transform:
            # warping
            print("{}: Register image by perspective transformation...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
            
            image = convert_to_8bit(image)
            registered = cv2.warpPerspective(image, H, (w, h))
        else:
            print("{}: Register image by affine transformation...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
            
            image = convert_to_8bit(image)
            registered = cv2.warpAffine(image, H, (w, h))

        if return_grayscale:
            if len(registered.shape) == 3:
                print("Convert registered image to grayscale...")
                registered = cv2.cvtColor(registered, cv2.COLOR_BGR2GRAY)
    else:
        registered = None

    # return the registered image
    return registered, H, matchedVis

def register_adata_coords_to_new_images(adata_in, groupby, image_dir_dict, groups=None, reg_channel='dapi',
                                        spatial_key='spatial', to_8bit=False, maxpx_before_reg=None, ppmhq=None, 
                                        hires_key='hires', lowres_key='lowres', perspective_transform=False,
                                        rot_threshold=0, do_registration=False, lowres_factor=0.2,
                                        keepFraction=0.2, method='sift', in_place=False):
    '''
    Function to add new images to an adata object and register the coordinates accordingly.

    The images are expected to be in adata.uns[spatial_key][image_id].
    The image_id is expected to be in the format "{well}_{channel}".
    ------------
    Returns:
    registered_img, H, matchedVis
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
            #hq_image_dict[key] = cv2.imread(image_dir_dict[key], 0)
            hq_image_dict[key] = cv2.imread(image_dir_dict[key], -1) # load unchanged

            if to_8bit:
                hq_image_dict[key] = convert_to_8bit(hq_image_dict[key])

            # convert to grayscale
            if len(hq_image_dict[key].shape) == 3:
                hq_image_dict[key] = cv2.cvtColor(hq_image_dict[key], cv2.COLOR_BGR2GRAY)

        # extract the registration image from the dict
        image_to_register = hq_image_dict[reg_key]

        # register images and extract homography matrix
        print("{}: Register image {}...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", reg_key))
        registered_img, H, matchedVis = register_image(image_adata, image_to_register, do_registration=do_registration, 
                                maxpx=maxpx_before_reg, perspective_transform=perspective_transform,
                                keepFraction=keepFraction, method=method)

        if matchedVis is not None:
            if 'matchedVis' in adata.uns.keys():
                adata.uns['matchedVis'][group] = matchedVis
            else:
                adata.uns['matchedVis'] = {}
                adata.uns['matchedVis'][group] = matchedVis

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
        if perspective_transform:
            print("{}: Perspective transformation of coordinates...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
            
            coords_reshaped = coords.reshape(-1, 1, 2)
            coords_trans = cv2.perspectiveTransform(coords_reshaped, H)
            coords_trans = coords_trans.reshape(-1, 2)
        else:
            print("{}: Affine transformation of coordinates...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
            coords_reshaped = np.array([[p] for p in coords])
            coords_trans = cv2.transform(coords_reshaped, H)
            coords_trans = coords_trans.reshape(-1, 2)

        # get rotation angle and pivot point
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
        recalculate_scale(adata, groupby=groupby, group=group, ppm_given=ppmhq,
            save_scale=True, return_angle_and_pivot=False)

    # resize images and store as `lowres` for plotting
    print("{}: Resize all images and store as `lowres`...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"))
    resize_images_in_adata(adata, scale_factor=lowres_factor)

    if not in_place:
        return adata

def align_image(image, vertices, frame: int = 100, return_grayscale=True, perspective_transform=False):
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

    # calculate transformation matrix and align
    if perspective_transform:
        print("\t\tPerspective transformation...")
        # get transformation matrix
        M = cv2.getPerspectiveTransform(pts1, pts2)
        
        # align image
        aligned = cv2.warpPerspective(image, M, (int(pxl_width + 2 * frame), int(pxl_width + 2 * frame)))
    else:
        print("\t\tAffine transformation...")
        # get transformation matrix
        (M, mask) = cv2.estimateAffine2D(pts1, pts2)

        # align image
        aligned = cv2.warpAffine(image, M, (int(pxl_width + 2 * frame), int(pxl_width + 2 * frame)))
    
    if return_grayscale:
        if len(aligned.shape) == 3:
            print("\t\tConvert image to grayscale...")
            aligned = cv2.cvtColor(aligned, cv2.COLOR_BGR2GRAY)

    return aligned


def align_to_dict(images, labels, vertices, resolution: int, n_channels: int,
                    frame: int = 100, ppm_given=None, **kwargs):

    # Part 1: Align image
    print("\tAlign images...")
    aligned = []
    for i, img in enumerate(images):
        print("\t\tProcessing image {} of {}:".format(i+1, len(images)))
        aligned.append(align_image(img, vertices, frame, **kwargs))

    # Part 2: Create metadata
    print("\tCreate metadata...")
    # calculate the width of the alignment marker square
    pxl_width = dist_points(vertices[0], vertices[1]).astype(int)
    pts2 = np.float32([[frame, frame], [pxl_width + frame, frame],
                       [pxl_width + frame, pxl_width + frame], [frame, pxl_width + frame]])

    # calculate resolution from distance of vertices and resolution
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
        "tissue_hires_scalef": 1.0,
        "pixel_per_um_real": ppm_given
    }

    # Part 3: Summarize aligned image and metadata
    print("\tSummarized aligned images and metadata.")
    image_sum = {lab: {'images': {'hires': img}, 'scalefactors': image_metadata}
                 for (lab, img) in zip(labels, aligned)}

    return image_sum