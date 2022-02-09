import cv2
import numpy as np
from ..calculations.contours import circularity, centroid
import pandas as pd
import math
from tqdm import tqdm

def segment_vessels(img, threshold_ratio=0.01, exclude_edge=True, min_points=20, min_diameter=0,
    min_circularity=0.2, scale=None, return_raw=False):
    '''
    Segment vessels or other structures that are defined by no or background signal.
    '''

    # check if data should be filtered
    do_filtering = False
    if (min_diameter > 0) or (min_circularity > 0):
        do_filtering=True

    # thresholding of background
    t = round(img.max() * threshold_ratio)

    # thresh = img <= t
    # thresh = thresh.astype(np.uint8) * 255
    # med = cv2.medianBlur(thresh, 9, 0)
    
    # contours, hierarchy = cv2.findContours(med, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    # if exclude_edge:
    #     contours = filter_edge_contours(contours, img)

    if len(img.shape) == 3:
        print("Convert image to grayscale...")
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # segment section
    #signal = img >= t
    #signal = signal.astype(np.uint8) * 255
    signal = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY_INV,11,2)
    signal = cv2.medianBlur(signal, 19, 0)
    signal = cv2.morphologyEx(signal, cv2.MORPH_CLOSE, np.ones((9,9),np.uint8))
    sections, _ = cv2.findContours(signal, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    if len(sections) > 1:
        # remove sections that are too small
        section_areas = [cv2.contourArea(cnt) for cnt in sections]
        area_threshold = np.max(section_areas) * 0.05 # 5% of the biggest section
        sections = [s for s, a in zip(sections, section_areas) if a > area_threshold]

    # segment background structures
    #structures = img <= t
    #structures = structures.astype(np.uint8) * 255
    structures = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,11,2)
    structures = cv2.medianBlur(structures, 19, 0)
    structures, _ = cv2.findContours(structures, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

    # check which of the segmented structures are actually inside the sections and therefore potential vessels
    in_section = []
    for section in sections:
        in_section.append([np.all(np.array([cv2.pointPolygonTest(section, tuple(map(int, p[0])), False) for p in s]) > 0) for s in structures])

    # check which structures are inside any of the sections and filter
    in_section = np.array(in_section).any(axis=0)
    vessels = [v for v, i in zip(structures, in_section) if i]

    if min_points is not None:
        # filter out contours with less than three points
        vessels = [cnt for cnt in vessels if len(cnt) >= min_points]

    # calculate area, perimeter and circularity of contours
    areas = np.array([cv2.contourArea(cnt) for cnt in vessels])
    peris = np.array([cv2.arcLength(cnt, True) for cnt in vessels])
    circs = np.array([circularity(cnt, area=areas[i], peri=peris[i]) for i, cnt in enumerate(vessels)])
    avg_dia = np.array([np.sqrt(area / math.pi) * 2 for area in areas])
    centroids = np.array([centroid(cnt) for cnt in vessels])

    # save everything in dataframe
    df = pd.DataFrame({
        "contour": vessels,
        "area": areas,
        "perimeter": peris,
        "circularity": circs,
        "avg_dia": avg_dia,
        "centroid_x": centroids[:, 0],
        "centroid_y": centroids[:, 1]
    })

    if do_filtering:
        # save unfiltered version
        df_raw = df.copy()

        if scale is not None:
            # convert pixel into µm
            min_diameter *= scale
        
        # calculate area from diameter
        #min_area = math.pi * (min_diameter / 2)**2

        # filter dataframe by area
        #df.query('area >= {}'.format(min_area), inplace=True)
        df.query('avg_dia >= {}'.format(min_diameter), inplace=True)
        
        df.query('circularity >= {}'.format(min_circularity), inplace=True)
    else:
        df_raw = None
    
    img_with_contours = img.copy()
    for cnt in df['contour'].values:
        cv2.drawContours(img_with_contours,[cnt],0,255,4)
    
    if return_raw:
        return df, img_with_contours, df_raw
    else:
        return df, img_with_contours



def segment_vessels_from_adata(adata, key_added='vessels', img_key_pattern='dapi', uns_key='spatial', 
    lowres=False, threshold_ratio=0.01, exclude_edge=True, min_diameter=0, use_scale=True,
    min_circularity=0.2):
    '''
    Segment vessels from adata object.
    '''
    # find image keys that have image key pattern
    img_keys = [k for k in adata.uns[uns_key].keys() if img_key_pattern in k]

    if lowres:
        res_key = 'lowres'
    else:
        res_key = 'hires'

    # iterate through images and segment vessels
    dfs = {}
    df_raws = {}
    imgs = {}
    metadata = {}
    for img_key in tqdm(img_keys):
        img = adata.uns[uns_key][img_key]['images'][res_key].copy()

        if use_scale == True:
            if lowres:
                scale_factor = adata.uns[uns_key][img_key]['scalefactors']['tissue_lowres_scalef']
            else:
                scale_factor = 1

            ppm = adata.uns[uns_key][img_key]['scalefactors']['pixel_per_um'] * scale_factor
        else:
            ppm = None

        df, img_with_contours, df_raw = segment_vessels(img, threshold_ratio=threshold_ratio, exclude_edge=exclude_edge,
                                                    scale=ppm, min_diameter=min_diameter, min_circularity=min_circularity,
                                                    return_raw=True)

        # add to dict of dataframes
        dfs[img_key] = df
        df_raws[img_key] = df_raw
        imgs[img_key] = img_with_contours
        metadata[img_key] = {
            "resolution_key": res_key,
            "scalefactor": scale_factor,
            "pixel_per_um": adata.uns[uns_key][img_key]['scalefactors']['pixel_per_um'],
            "channel": img_key_pattern
            }

    # concatenate dataframes of all img_keys
    df = pd.concat(dfs)

    # add category to .uns to save the data
    adata.uns[key_added] = {}    
    adata.uns[key_added]['data'] = df
    adata.uns[key_added]['images'] = imgs
    adata.uns[key_added]['metadata'] = metadata

    if not np.all([elem is None for elem in list(df_raws.values())]):
        # add unfiltered data to anndata
        df_raw = pd.concat(df_raws)
        adata.uns[key_added]['raw'] = df_raw


def filter_edge_contours(contours, img):
        ## filter out contours that touch the edge of the image
        # determine dimensions of image
        w = img.shape[1]-1
        h = img.shape[0]-1

        # get a unique list of x and y values of each contour
        xs = [np.unique(cnt[:, :, 0].flatten()) for cnt in contours]
        ys = [np.unique(cnt[:, :, 1].flatten()) for cnt in contours]

        # test if any point of the contours touches the edges
        contours = [cnt for cnt, x, y in zip(contours, xs, ys) if not (0 in x or w in x or 0 in y or h in y)]

        return contours