import pandas as pd
import numpy as np
import math
from typing import Optional, Tuple, Union, List, Dict, Any


def decode_list(l):
    new_list = []
    for elem in l:
        try:
            new_elem = elem.decode()
            new_list.append(new_elem)
        except (UnicodeDecodeError, AttributeError):
            new_list.append(elem)
    
    return new_list


def rotatePoint(origin, point, angle, radians=False):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in degrees or radians (with `radians=True`).
    """
    if not radians:
        angle = -math.radians(angle)

    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return [qx, qy]


def get_nrows_maxcols(keys, max_cols):

    n_plots = len(keys)
    if n_plots > max_cols:
        n_rows = math.ceil(n_plots / max_cols)
    else:
        n_rows = 1
        max_cols = n_plots

    return n_plots, n_rows, max_cols

def check_list(List, list_to_compare):
    '''
    Compare two lists and return the elements that are in both lists.
    If not all elements are in both lists give message telling which are not.
    '''

    not_in = []
    List = [l if l in list_to_compare else not_in.append(l) for l in List]

    # remove None values
    List = [elem for elem in List if elem is not None]

    if len(not_in) > 0:
        print("Following elements not found: {}".format(", ".join(not_in)), flush=True)

    return List

def set_dict_difference(set_dict, set_key):
    '''
    Computes difference of one set with multiple other sets.
    Sets need to be part of a dictionary.
    `set_dict`: Dictionary containing the sets.
    `set_key`: Key for the set of which to compare the difference to all other sets in the dictionary.
    
    Returns set of unique elements.
    '''
    s = set_dict[set_key]
    keys = [key for key in set_dict.keys() if key != set_key]
    for key in keys:
        s = s.difference(set_dict[key])
    return s
    



def nth_repl_all(s, sub, repl, nth):

    '''
    Replaces a pattern `sub` with `repl` in a string `s` after every `nth` occurrence.
    '''

    find = s.find(sub)
    # loop util we find no match
    i = 1
    while find != -1:
        # if i  is equal to nth we found nth matches so replace
        if i == nth:
            s = s[:find]+repl+s[find + len(sub):]
            i = 0
        # find + len(sub) + 1 means we start after the last match
        find = s.find(sub, find + len(sub) + 1)
        i += 1
    return s

def save_multisheet_excel(dataframe, groupby, savefile):
    '''
    Saves pandas dataframe to xlsx file grouping the file into multiple sheets based on values in one column.
    '''

    writer = pd.ExcelWriter(savefile)

    for group, data in dataframe.groupby(groupby):
        data.to_excel(writer,group)
    writer.save()

    return print('Saved dataframe into {}'.format(savefile))


def select_keys(dictionary, keys):
    return {key: dictionary[key] for key in dictionary if key in keys}

def otsu_threshold(counts, bins_num = 256, is_normalized=False):
    '''
    Implementation of Otsu thresholding algorithm from:
    https://learnopencv.com/otsu-thresholding-with-opencv/

    Re-implementation was necessary since we couldn't find out how to use the thresholding algorithms
    from skimage and cv2 on counts instead of images.
    '''
    # Get the histogram
    hist, bin_edges = np.histogram(counts, bins=bins_num)

    # Get normalized histogram if it is required
    if is_normalized:
        hist = np.divide(hist.ravel(), hist.max())

    # Calculate centers of bins
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2.

    # Iterate over all thresholds (indices) and get the probabilities w1(t), w2(t)
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]

    # Get the class means mu0(t)
    mean1 = np.cumsum(hist * bin_mids) / weight1
    # Get the class means mu1(t)
    mean2 = (np.cumsum((hist * bin_mids)[::-1]) / weight2[::-1])[::-1]

    inter_class_variance = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    # Maximize the inter_class_variance function val
    index_of_max_val = np.argmax(inter_class_variance)

    threshold = bin_mids[:-1][index_of_max_val]
    return int(round(threshold))

def rle(inarray):
        """ run length encoding. Partial credit to R rle function. 
            Multi datatype arrays catered for including non Numpy
            returns: tuple (runlengths, startpositions, values) 
            from: https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi"""
        ia = np.asarray(inarray)                # force numpy
        n = len(ia)
        if n == 0: 
            return (None, None, None)
        else:
            y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
            i = np.append(np.where(y), n - 1)   # must include last element posi
            z = np.diff(np.append(-1, i))       # run lengths
            p = np.cumsum(np.append(0, z))[:-1] # positions

            d = {
                'run_lengths': z,
                'positions': p,
                'labels': ia[i]
            }
            return d

def read_rna_metrics(rnametrics_files, data_names=None):
    '''
    Read RNA Metrics from files
    '''
    
    # make inputs to list
    rnametrics_files = [rnametrics_files] if isinstance(rnametrics_files, str) else list(rnametrics_files)

    if data_names is None:
        data_names = ["Sample\_{}".format(i) for i in range(len(rnametrics_files))]

    # read metrics file
    metrics_dict = {}

    for n, f in zip(data_names, rnametrics_files):
        metrics = []
        for line in open(f, "r"):
            if not line.startswith("#"):
                metrics.append(line.strip())
        metrics_dict[n] = metrics

    # extract metric names and values and collect in dataframe
    names_dict = {key: metrics_dict[key][1].split("\t")[:-3] for key in metrics_dict}
    values_dict = {key: metrics_dict[key][2].split("\t") for key in metrics_dict}

    df_dict = {key: pd.DataFrame({'metric': names_dict[key], 'value': values_dict[key]}) for key in names_dict}
    df = pd.concat(df_dict)

    return df

    