from xml.etree.ElementInclude import include
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, Tuple, Union, List, Dict, Any
import anndata
from .tools import decode_list, check_list


def extract_groups(adata, groupby, groups=None, extract_uns=False, uns_key='spatial', uns_exclusion_pattern=None, 
    return_mask=False, strip=False):

    '''
    Function to extract a group from a dataframe or anndata object. 
    If `anndata_object` is True it expects the dataframe in `adata.obs`    
    '''

    # convert groups into list
    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)

    if type(adata) == anndata.AnnData:
        anndata_object = True
    elif type(adata) == pd.DataFrame:
        anndata_object = False
        extract_uns = False
        strip = False
    else:
        raise ValueError("Unknown type of input object.")
    
    # select dataframe on which to filter on
    if anndata_object:
        obs = adata.obs
    else:
        obs = adata

    # check if filtering is wanted
    filtering = True
    if groupby is None:
        filtering = False

    # # check if we want to filter `.uns`
    # if uns_exclusion_pattern is not None:
    #     extract_uns = True

    if groupby in obs.columns or not filtering:

        # create filtering mask for groups
        if groupby is None:
            mask = np.full(len(adata), True) # no filtering
        else:
            mask = obs[groupby].isin(groups).values

        # filter dataframe or anndata object
        if anndata_object:
            adata = adata[mask, :].copy()
        else:
            adata = adata.loc[mask, :].copy()

        if len(adata) == 0:
            print("Subset variables '{}' not in groupby '{}'. Object not returned.".format(groups, groupby))
            return
        elif filtering:
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
            # remove annotations in .uns and obsp
            stores = [adata.uns, adata.obsp]
            for s in stores:
                keys = list(s.keys())
                for k in keys:
                    del s[k]
        
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


def check_raw(adata, use_raw, layer=None):
    # check if plotting raw data
    if use_raw:
        adata_X = adata.raw.X
        adata_var = adata.raw.var
        adata_var_names = adata.raw.var_names
    else:
        if layer is None:
            adata_X = adata.X
        else:
            #adata_X = adata.layers[layer].toarray()
            adata_X = adata.layers[layer]
        adata_var = adata.var
        adata_var_names = adata.var_names
    
    return adata_X, adata_var, adata_var_names
    
def create_color_dict(adata, key, palette):
    if key in adata.obs.columns and adata.obs[key].dtype.name == 'category':
        # create color dictionary from color palette
        categories = adata.obs[key].values.categories
        color_dict = dict(zip(categories, sns.color_palette(palette=palette, n_colors=len(categories))))
    else:
        color_dict = palette

    return color_dict

def create_deg_df(adata, keys, groups=None, include_means=False, added_cats=None, use_raw=False):
    '''
    Function to generate pandas dataframe from DEG analyze which used the `sc.tl.rank_genes_groups()` function.
    '''

    # transform keys to list
    keys = [keys] if isinstance(keys, str) else list(keys)

    # check if .raw should be used
    _, adata_var, _ = check_raw(adata, use_raw=use_raw)

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
                df['means'] = [adata_var.loc[name]['means'] for name in df['names']]
            
            group_list.append(df)
        
        group_df = pd.concat(group_list)
        final_dict[key] = group_df
    
    final_df = pd.concat(final_dict)
    final_df.index.set_names(['key', 'idx'], inplace=True)

    return final_df

def collect_deg_data(adata, keys, groups=None, foldchanges_label='logfoldchanges', pvals_label='pvals_adj', gene_names='names', 
                fc_threshold=1, pval_threshold=0.05, include_means=False, **kwargs):

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
        datas = {g: create_deg_df(adata, keys=key, groups=g, include_means=include_means, **kwargs) for g in groups}

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

def get_crange(adata, groupby, groups, key, use_raw, 
    layer=None, data_in_dataframe=False, pd_dataframe=None, 
    ctype='minmax', cmin_at_zero=True
    ):

    if data_in_dataframe and pd_dataframe is not None:
        obs = pd_dataframe
    else:
        obs = adata.obs
    
    adata_X, adata_var, adata_var_names = check_raw(adata, use_raw=use_raw, layer=layer)

    if key in adata_var_names:
        #c = adata_X[[elem in groups for elem in adata.obs[groupby]], adata_var_names == key]
        c = adata_X[:, adata_var_names == key]
        if cmin_at_zero:
            cmin = 0
        else:
            cmin = c.min()

        if ctype == 'percentile':
            cmax = np.percentile(c, 95)
        else:
            cmax = c.max()
        crange = [cmin, cmax]
    elif key in obs.columns:
        if obs[key].dtype.name.startswith('float') or obs[key].dtype.name.startswith('int'):
            #c = obs[key][[elem in groups for elem in obs[groupby]]]
            c = obs[key]
            if cmin_at_zero:
                cmin = 0
            else:
                cmin = c.min()
            cmax = np.percentile(c, 95)
            crange = [cmin, cmax]
        else:
            return
    else:
        print("Key not in var_names of adata object. Use raw?")
        return

    return crange

def remove_images(adata, uns_key='spatial', reg_key='registered', other_keys=None,
    hires_only=False, hires_key='hires'):
    '''
    Remove all images from adata.
    Images are expected to be in `adata.uns[uns_key]`
    '''
    image_keys = list(adata.uns[uns_key].keys())
    #adata = adata.copy()
    for key in image_keys:
        res_keys = list(adata.uns[uns_key][key]['images'].keys())

        if hires_only:
            del adata.uns[uns_key][key]['images'][hires_key]
        else:
            for res_key in res_keys:
                del adata.uns[uns_key][key]['images'][res_key]

    if not hires_only:
        # remove registered images if available
        if reg_key in adata.uns:
            del adata.uns[reg_key]

    # remove further keys if given
    if other_keys is not None:
        other_keys = [other_keys] if isinstance(other_keys, str) else list(other_keys)
        for k in other_keys:
            if k in adata.uns:
                del adata.uns[k]

    #return adata

def replace_images(adata, image_adata, spatial_key="spatial", inplace=True, verbose=True):
    '''
    Replace images in `adata` with images from `image_adata`.
    Images in both adatas are expected to be stored in `.uns[spatial_key]` under a unique id that must
    match between the two adatas.
    '''
    
    if inplace:
        adata_ = adata
    else:
        adata_ = adata.copy()
    
    # get image names from adata
    img_keys = adata_.uns[spatial_key].keys()
    
    for img_key in img_keys:
        if img_key in image_adata.uns[spatial_key].keys():
            adata_.uns[spatial_key][img_key] = image_adata.uns[spatial_key][img_key].copy()
            
        elif verbose:
            print("Image key {} not found in `image_adata`. Was skipped.".format(img_key))
            
    if not inplace:
        return adata_

def get_images_from_adata(adata, groupby, group, hires=False):
    '''
    Extract images from adata.
    '''
    
    if hires:
        res_key = 'hires'
    else:
        res_key = 'lowres'
    
    if groupby is not None:
        # get subset
        adata = extract_groups(adata, groupby=groupby, groups=group, extract_uns=True)
    
    # get keys
    keys = list(adata.uns['spatial'].keys())
    
    imgs = {}
    for key in keys:
        # get image and channel name
        channel = key.split("-")[1]
        img = adata.uns['spatial'][key]['images'][res_key]
        
        # collect images in dictionary
        imgs[channel] = img
    return imgs

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