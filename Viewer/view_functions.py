import cv2
import numpy as np
import anndata
import pandas as pd

# own functions
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

def extract_groups(adata, groupby, groups, extract_uns=True, uns_key='spatial', uns_exclusion_pattern=None, return_mask=False):

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

    if groupby in obs.columns:

        # create filtering mask for groups
        mask = obs[groupby].isin(groups)

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

            if extract_uns:
                new_uns = {key:value for (key,value) in adata.uns[uns_key].items() if np.any([group in key for group in groups])}
                
                if uns_exclusion_pattern is not None:
                    new_uns = {key:value for (key,value) in new_uns.items() if uns_exclusion_pattern not in key}
                adata.uns[uns_key] = new_uns
            
            if return_mask:
                return adata, mask
            else:
                return adata
    
    else:
        print("Subset category '{}' not found".format(groupby))
        return