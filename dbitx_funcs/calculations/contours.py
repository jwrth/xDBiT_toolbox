import cv2
import math
import numpy as np
from ..tools import extract_groups
from scipy import stats
import pandas as pd

def centroid(contour):
    '''
    Calculate the centroid of a contour.
    '''
    # calculate moments
    M = cv2.moments(contour)
    
    # calculate coordinates of centroid
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    return cx, cy

def circularity(contour, area=None, peri=None):
    '''
    Calculate circularity of contour using following equation:
    circularity = 4pi(area/perimeter^2)
    '''
    if area is None:
        area = cv2.contourArea(contour)
    if peri is None:
        peri = cv2.arcLength(contour, True)
    
    circularity = 4 * math.pi * (area / peri**2)
    return circularity

def distance_weight(distances, sigma=20):

    wts = [np.exp(-elem / sigma) for elem in distances]
    sum_wts = np.sum(wts)

    return [w / sum_wts for w in wts]

def calculate_neighborhood_score(adata, gene_set, distance_threshold, key_added, groupby='id',
                                 sigma=20, uns_key='vessels', contour_key='contour', return_results=False):
    '''
    Calculate a distance-weighted neighborhood score for structures in spatial transcriptomics dataset.
    Structures are expected as OpenCV contours saved in a pandas DataFrame in `adata.uns[uns_key]['data']`.
    Analysis can be grouped by a category in `adata.obs`. If no grouping is needed choose `groupby=None`.
    '''
    results = {}
    for idx in adata.obs[groupby].unique():
        # get subset of adata
        if groupby is not None:
            subset = extract_groups(adata, groupby=groupby, groups=idx, strip=True)
        else:
            subset = adata.copy()

        # extract data from subset
        gene_mask = subset.var_names.get_indexer(gene_set)
        subset = subset[:, gene_mask].copy()
        subset_obs = subset.obs
        subset_X = subset.X

        # standardize expression values of subset
        subset_X = stats.zscore(subset_X, axis=0)

        # find correct image keys in the structure data and extract structure data
        img_key = [elem for elem in adata.uns[uns_key]['data'].index.unique(level=0) if idx in elem][0]
        structures = adata.uns[uns_key]['data'].xs(img_key).copy()

        # retrieve metadata
        sf = adata.uns[uns_key]['metadata'][img_key]['scalefactor']
        ppm = adata.uns[uns_key]['metadata'][img_key]['pixel_per_um']
        img = adata.uns[uns_key]['images'][img_key]

        # iterate through vessels
        scores = []
        for _, v in structures.iterrows():
            # select vessel
            vc = v[contour_key]

            # iterate through spots
            dists = []
            for _, spot in subset_obs.iterrows():
                x = int(spot['pixel_col'] * sf)
                y = int(spot['pixel_row'] * sf)

                d = abs(cv2.pointPolygonTest(vc, (x, y), True)) / (ppm * sf) # calculate distance in µm
                dists.append(d)
            dists = np.array(dists)

            # calculate neighborhood
            inside_nh = dists <= distance_threshold # determine spots that are inside neighborhood
            N = subset_obs[inside_nh][["pixel_row", "pixel_col"]] # extract coordinates of neighborhood spots
            dists = dists[inside_nh]
            Nx = subset_X[inside_nh, :]

            # calculate distance weight
            wts = distance_weight(distances=dists, sigma=sigma)

            # calculate weighted expression for each spot and each gene
            weighted_expr = np.array([Nx[i, :] * wts[i] for i in range(len(wts))])

            # sum up weighted expression for the whole neighborhood
            weighted_expr = np.sum(weighted_expr, axis=0)

            # sum up across marker genes to get score and save
            score = np.sum(weighted_expr)
            scores.append(score)
        structures[key_added] = scores
        results[img_key] = structures
    
    results = pd.concat(results)
    
    adata.uns[uns_key]['data'] = results.copy()

    if return_results:
        return results

