from matplotlib.pyplot import install_repl_displayhook
import NaiveDE
import SpatialDE
import pandas as pd
from collections import OrderedDict
import numpy as np
from scipy import ndimage
import seaborn as sns
import math
from typing import Optional, Tuple
import anndata
from PIL import Image
import scanpy as sc
from .. import utils
from ..exceptions import ModuleNotFoundOnWindows

#from ..plotting import volcano_plot

# Some modules make problems on Windows. Therefore I imported them conditionally
# try:
#     from skmisc.loess import loess
# except ImportError:
#     print("Import issue with loess function detected. Probably due to Windows. Error ignored.")
#     pass

# try:
#     import bbknn
# except ModuleNotFoundError:
#     print("Package 'bbknn' not installed here. Error ignored.")
#     pass



def spatialde_run(adata, layer=None, run=True, normalize=True, output_name='spatialde', use_raw=False):

    print("Prepare data for SpatialDE analysis")
    if layer is not None:
        print("Use data from layer {}.".format(layer))
        counts = pd.DataFrame(adata.layers[layer], columns=adata.var_names)
    elif use_raw:
        print("Use raw genes.")
        counts = pd.DataFrame(adata.raw.X, columns=adata.raw.var_names)
    else:
        print("Use filtered genes.")
        counts = pd.DataFrame(adata.X, columns=adata.var_names)

    index_names = [str(col) + "x" + str(row)
                   for col, row in zip(adata.obs.um_col, adata.obs.um_row)]
    counts.index = index_names

    sample_info = adata.obs[['um_col', 'um_row', 'total_counts']]
    sample_info.columns = ['x', 'y', 'total_counts']
    sample_info.index = index_names

    spatialde_sum = OrderedDict()

    if normalize:
        print("Normalize counts...")
        norm_expr = NaiveDE.stabilize(counts.T).T
        resid_expr = NaiveDE.regress_out(
            sample_info, norm_expr.T, 'np.log(total_counts)').T

        spatialde_sum['norm_expr'] = norm_expr
        spatialde_sum['resid_expr'] = resid_expr
        spatialde_sum['counts'] = counts

    else:
        resid_expr = counts
        spatialde_sum['counts'] = counts

    # run all genes

    if run:
        print("Run SpatialDE...")
        X = sample_info[['x', 'y']]
        results = SpatialDE.run(X, resid_expr)

        spatialde_sum['results'] = results
        spatialde_sum['X'] = X

    adata.uns[output_name] = spatialde_sum

def decode_list(l):
    new_list = []
    for elem in l:
        try:
            new_elem = elem.decode()
            new_list.append(new_elem)
        except (UnicodeDecodeError, AttributeError):
            new_list.append(elem)
    
    return new_list

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

def extract_groups(adata, groupby, groups, extract_uns=False, uns_key='spatial', uns_exclusion_pattern=None, 
    return_mask=False, strip=False):

    '''
    Function to extract a group from a dataframe or anndata object. 
    If `anndata_object` is True it expects the dataframe in `adata.obs`    
    '''

    # convert groups into list
    groups = [groups] if isinstance(groups, str) else list(groups)

    if type(adata) == anndata._core.anndata.AnnData:
        anndata_object = True
    elif type(adata) == pd.core.frame.DataFrame:
        anndata_object = False
        extract_uns = False
    else:
        raise ValueError("Unknown type of input object.")
    
    # select dataframe on which to filter on
    if anndata_object:
        obs = adata.obs
    else:
        obs = adata

    # # check if we want to filter `.uns`
    # if uns_exclusion_pattern is not None:
    #     extract_uns = True

    if groupby in obs.columns:

        # create filtering mask for groups
        mask = obs[groupby].isin(groups).values

        # filter dataframe or anndata object
        if anndata_object:
            adata = adata[mask, :].copy()
        else:
            adata = adata.loc[mask, :].copy()

        if len(adata) == 0:
            print("Subset variables '{}' not in groupby '{}'. Object not returned.".format(groups, groupby))
            return
        else:
            # check if all groups are in groupby category
            groups_found = [group for group in groups if group in obs[groupby].unique()]
            groups_notfound = [group for group in groups if group not in groups_found]
            
            if len(groups_found) != len(groups):
                print("Following groups were not found in column {}: {}".format(groupby, groups_notfound))

            if extract_uns or uns_exclusion_pattern is not None:
                new_uns = {key:value for (key,value) in adata.uns[uns_key].items() if np.any([group in key for group in groups])}
                
                if uns_exclusion_pattern is not None:
                    new_uns = {key:value for (key,value) in new_uns.items() if uns_exclusion_pattern not in key}
                adata.uns[uns_key] = new_uns

            if strip:
                # remove all annotations but .var, .obs and .obsm
                del adata.uns
                del adata.varm
                del adata.layers
                del adata.obsp
            
            if return_mask:
                return adata, mask
            else:
                return adata
    
    else:
        print("Subset category '{}' not found".format(groupby))
        return
    
    # if not data_in_dataframe:
    #     # search for groupby category in adata.obs
    #     if groupby in adata.obs.columns:
    #         if group in adata.obs[groupby].values:
    #             mask = adata.obs[groupby] == group
    #             adata = adata[mask, :].copy()
    #         else:
    #             print("Subset variable '{}' not in groupby '{}'".format(
    #                 group, groupby))
    #             return

    #         if extract_uns:
    #             new_uns = {key:value for (key,value) in adata.uns[uns_key].items() if group in key}
    #             adata.uns[uns_key] = new_uns

    #     else:
    #         print("Subset category '{}' not found".format(groupby))
    #         return

    # else:
    #     # search for groupby category in columns of dataframe
    #     if groupby in adata.columns:
    #         if group in adata[groupby].values:
    #             mask = adata[groupby] == group
    #             adata = adata[mask].copy()
    #         else:
    #             print("Subset variable '{} not in groupby '{}'".format(group, groupby))
    #             return

    #     else:
    #         print("Subset cateogry '{}' not found".format(groupby))
    #         return

    # if return_mask:
    #     return adata, mask
    # else:
    #     return adata


def check_raw(adata, use_raw):
    # check if plotting raw data
    if use_raw:
        adata_X = adata.raw.X
        adata_var = adata.raw.var
        adata_var_names = adata.raw.var_names
    else:
        adata_X = adata.X
        adata_var = adata.var
        adata_var_names = adata.var_names

    return adata_X, adata_var, adata_var_names
    
def create_color_dict(adata, key, palette):
    if key in adata.obs.columns and adata.obs[key].dtype.name == 'category':
        # create color dictionary from color palette
        categories = adata.obs[key].values.categories
        color_dict = dict(zip(categories, sns.color_palette(palette=palette, n_colors=len(categories))))
    else:
        color_dict = None

    return color_dict

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
    
def smooth_fit(xs : np.ndarray, ys : np.ndarray, x_to_fit_on=None,
dist_thrs : Optional[float] = None,):

    """Smooth curve using loess

    will perform curve fitting using skmisc.loess,
    points above 'dist_thrs' will be excluded.
    Parameters:
    ----------
    xs : np.ndarray
        x values
    ys : np.ndarray
        y values
    dist_thrs : float
        exclude (x,y) tuples where x > dist_thrs
    Returns:
    -------
    A tuple with included x and y-values (xs',ys'), as well
    as fitted y-values (ys_hat) together with associated
    standard errors. The tuple is on the form
    (xs',ys',y_hat,std_err)

    From: https://github.com/almaan/ST-mLiver
    """

    try:
        from skmisc.loess import loess
    except ModuleNotFoundError as e:
        raise ModuleNotFoundOnWindows(e)

    srt = np.argsort(xs)
    xs = xs[srt]
    ys = ys[srt]

    if dist_thrs is None:
        dist_thrs = np.inf

    keep = np.abs(xs) < dist_thrs
    xs = xs[keep]
    ys = ys[keep]

    # generate loess class object
    ls = loess(xs,
               ys,
              )
    # fit loess class to data
    ls.fit()

    # check if there are specific x values to fit the data on
    if x_to_fit_on is not None:
        # make sure that all values to fit on lie inside the give x values
        x_to_fit_on = x_to_fit_on[(x_to_fit_on > xs.min()) & (x_to_fit_on < xs.max())]
        xs = np.asarray(x_to_fit_on)

    # predict on data
    pred =  ls.predict(xs,
                       stderror=True)
    # get predicted values
    ys_hat = pred.values
    # get standard error
    stderr = pred.stderr

    df = pd.DataFrame({
        'x': xs,
        #'y': ys,
        'y_pred': ys_hat,
        'std': stderr
    })

    #return (xs,ys,ys_hat,stderr)
    return df


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


def create_deg_df(adata, keys, groups=None, include_means=True, added_cats=None):
    '''
    Function to generate pandas dataframe from DEG analyze which used the `sc.tl.rank_genes_groups()` function.
    '''

    # transform keys to list
    keys = [keys] if isinstance(keys, str) else list(keys)

    # check if there were any additional categories added
    cats = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    if added_cats is not None:
        added_cats = [added_cats] if isinstance(added_cats, str) else list(added_cats)
        cats += added_cats

    final_dict = {}
    for key in keys:
        if groups is None:
            groups_ = list(adata.uns[key]['names'].dtype.names)
        else:
            groups_ = groups

        # transform groups to list
        groups_ = [groups_] if isinstance(groups_, str) else list(groups_)

        group_list = []
        for group in groups_:
            length = len(adata.uns[key][cats[0]][group])

            # create dataframe
            df = pd.DataFrame()
            df['group'] = [group] * length
            df['reference'] = [adata.uns[key]['params']['reference']] * length
            df['method'] = [adata.uns[key]['params']['method']] * length

            for cat in cats:
                df[cat] = adata.uns[key][cat][group]

            if include_means:
                df['names'] = decode_list(df['names']) # decode byte-coded strings if there are any
                df['means'] = [adata.var.loc[name]['means'] for name in df['names']]
            
            group_list.append(df)
        
        group_df = pd.concat(group_list)
        final_dict[key] = group_df
    
    final_df = pd.concat(final_dict)
    final_df.index.set_names(['key', 'idx'], inplace=True)

    return final_df


def collect_deg_data(adata, keys, groups=None, foldchanges_label='logfoldchanges', pvals_label='pvals_adj', gene_names='names', 
                fc_threshold=1, pval_threshold=0.05, **kwargs):

    '''
    Collect data of a DEG analysis.
    Data is expected in `adata.uns[key]` as output of `sc.tl.rank_genes_groups`.

    Returns:
    `groups, keys, key_data_dict, key_updown_dict, params_dict`
    Return dataframes for each key/group instead of a plot.
    '''

    # transform keys and groups to list
    keys = [keys] if isinstance(keys, str) else list(keys)

    groups_given = False if groups is None else True
    
    key_data_dict = {}
    key_updown_dict = {}
    params_dict = {}
    for key in keys:
        # extract information about DEG analysis
        extracted_groups = list(adata.uns[key]['names'].dtype.names)

        if groups_given:
            groups = check_list(groups, extracted_groups)

            if len(groups) == 0:
                return
        else:
            groups = extracted_groups

        # extract data and collect in dictionary
        datas = {g: create_deg_df(adata, keys=key, groups=g, **kwargs) for g in groups}

        data_dict = {}
        updown_dict = {}
        for group in groups:
            # select one group
            data = datas[group]

            # remove zeros from pvals to allow log transformation
            data[pvals_label] += data[pvals_label][data[pvals_label] > 0].min()

            # convert to -log10(foldchanges)
            data['-logpvals'] = -np.log10(data[pvals_label].values)
            logpvals = data['-logpvals'].copy()
            # extract parameters
            logfold = data[foldchanges_label].copy()
            
            # extract significantly up/down regulated genes
            up = data[(logfold > fc_threshold) & (logpvals > -np.log10(pval_threshold))].copy()
            up['combined'] = [a * b for a,b in zip(up[foldchanges_label], up['-logpvals'])]
            up = up.sort_values('combined', ascending=False)
            up['updown'] = 'up'

            down = data[(logfold < -fc_threshold) & (logpvals > -np.log10(pval_threshold))].copy()
            down['combined'] = [a * b for a,b in zip(down[foldchanges_label], down['-logpvals'])]
            down = down.sort_values('combined', ascending=True)
            down['updown'] = 'down'

            # summarize up down 
            updown = pd.concat([up, down], keys=['up', 'down'])

            # save data to dictionary for later plotting or return
            data_dict[group] = data
            updown_dict[group] = updown
        
        # collect data per key
        key_data_dict[key] = data_dict
        key_updown_dict[key] = updown_dict
        params_dict[key] = adata.uns[key]['params']
        
    return groups, keys, key_data_dict, key_updown_dict, params_dict

def select_keys(dictionary, keys):
    return {key: dictionary[key] for key in dictionary if key in keys}

def collect_spatialde_results(adata, uns_key):
    '''
    Collects the results of a grouped SpatialDE run. Returns dataframe.
    '''
    adata_ = adata.copy()

    results_df = {}
    for well in adata_.uns[uns_key].keys():
        results_df[well] = adata_.uns[uns_key][well]["results"]
        results_df[well]["well_name"] = well
        #results_df[well]["age"] = exp_conditions[well]

    results_df = pd.concat(results_df, ignore_index=True)
    #results_df = results_df[['age', 'well_name', 'g', 'l', 'qval']]
    return results_df

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

def scanorama(adata, batch, hvg=False, hvg_key='highly_variable', **kwargs):

    '''
    Function to perform Scanorama batch correction (https://github.com/brianhie/scanorama/).
    Code partially from: https://github.com/theislab/scib.
    '''

    import scanorama        

    utils.check_sanity(adata, batch, hvg, hvg_key)

    hvg_genes = list(adata.var.index[adata.var[hvg_key]])

    split, categories = utils.split_batches(adata.copy(), batch, hvg=hvg_genes, return_categories=True)
    corrected = scanorama.correct_scanpy(split, return_dimred=True, **kwargs)
    corrected = anndata.AnnData.concatenate(
        *corrected, batch_key=batch, batch_categories=categories, index_unique=None
    )
    corrected.obsm['X_emb'] = corrected.obsm['X_scanorama']
    # corrected.uns['emb']=True

    # add scanorama results to original adata - make sure to have correct order of obs
    X_scan = corrected.obsm['X_scanorama']
    orig_obs_names = list(adata.obs_names)
    cor_obs_names = list(corrected.obs_names)
    adata.obsm['X_scanorama'] = np.array([X_scan[orig_obs_names.index(o)] for o in cor_obs_names])
    adata.obsm['X_emb'] = adata.obsm['X_scanorama']

    #return corrected
    return adata

def standard_preprocessing(adata_in, batch_key_hvg=None, do_lognorm=True, regress_out=None, 
    filter_hvg=False, batch_correction_key=None, dim_reduction=True):

    '''
    Function to perform standard preprocessing on ST data. Adapted from Squidpy Napari Tutorial.
    '''

    try:
        import bbknn
    except ModuleNotFoundError as e:
        raise ModuleNotFoundOnWindows(e)

    adata = adata_in.copy()

    # if count matrix does not consist of raw counts abort preprocessing
    if not np.all(np.modf(adata.X)[0] == 0):
        print("`adata.X` does not contain raw counts. Preprocessing aborted.")
        return 

    if do_lognorm:
        # store raw counts in layer
        print("Store raw counts in adata.layers['counts']...")
        adata.layers['counts'] = adata.X.copy()

        # preprocessing according to napari tutorial in squidpy
        print("Normalization, log-transformation...")
        sc.pp.normalize_total(adata)
        adata.layers['norm_counts'] = adata.X.copy()
        sc.pp.log1p(adata)
        print("Calculate highly-variable genes...")
        
    sc.pp.highly_variable_genes(adata, batch_key=batch_key_hvg)

    if filter_hvg:
        print("Filter for highly-variable genes...")
        # add the normalized and log data to raw
        adata.raw = adata
        # Filter for highly variable genes
        adata = adata[:, adata.var.highly_variable].copy()

    if regress_out is not None:
        print("Regress out {}...".format(regress_out))
        sc.pp.regress_out(adata, regress_out)

    if dim_reduction:
        if batch_correction_key is None:
            # dimensionality reduction
            print("Dimensionality reduction...")
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.tsne(adata)

        else:
            # PCA
            sc.pp.pca(adata)

            neigh_uncorr_key = 'neighbors_uncorrected'
            sc.pp.neighbors(adata, key_added=neigh_uncorr_key)

            # dim reduction with uncorrected data
            #sc.tl.umap(adata, neighbors_key=neigh_uncorr_key)
            #sc.tl.tsne(adata)
            # clustering
            sc.tl.leiden(adata, neighbors_key=neigh_uncorr_key, key_added='leiden_uncorrected')  

            # batch correction
            print("Batch correction for {}...".format(batch_correction_key))
            bbknn.bbknn(adata, batch_key=batch_correction_key, metric='euclidean') # is used as alternative to sc.pp.neighbors

            # dim reduction with corrected data
            print("Dimensionality reduction with batch corrected data...")
            sc.tl.umap(adata)
            sc.tl.tsne(adata)

        # clustering
        print("Leiden clustering...")
        sc.tl.leiden(adata)

    # if not in_place:
    #     return adata
    return adata

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

def get_crange(adata, groupby, groups, key, use_raw, data_in_dataframe=False, pd_dataframe=None):

    if data_in_dataframe and pd_dataframe is not None:
        obs = pd_dataframe
    else:
        obs = adata.obs
    
    adata_X, adata_var, adata_var_names = check_raw(adata, use_raw=use_raw)

    if key in adata_var_names:
        #c = adata_X[[elem in groups for elem in adata.obs[groupby]], adata_var_names == key]
        c = adata_X[:, adata_var_names == key]
        cmin = c.min()
        cmax = c.max()
        crange = [cmin, cmax]
    elif key in obs.columns:
        if obs[key].dtype.name.startswith('float') or obs[key].dtype.name.startswith('int'):
            #c = obs[key][[elem in groups for elem in obs[groupby]]]
            c = obs[key]
            cmin = c.min()
            cmax = c.max()
            crange = [cmin, cmax]
        else:
            return
    else:
        print("Key not in var_names of adata object. Use raw?")
        return

    return crange
