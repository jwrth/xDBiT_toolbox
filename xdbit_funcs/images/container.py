import squidpy as sq
from ..tools import extract_groups

def createImageContainer(adata, groupby, group, hires=False, return_adata=False):
    '''
    Create ImageContainer from adata.
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
    
    # create image container
    imgc = sq.im.ImageContainer()    
    
    for key in keys:
        # get image and channel name
        channel = key.split("-")[1]
        img = adata.uns['spatial'][key]['images'][res_key]
        
        # add to container
        imgc.add_img(img, layer=channel)
        
    # add library id
    lib_id_to_add = keys[0]
    imgc.library_ids = [lib_id_to_add]
    
    # get scaling settings
    ppm = adata.uns['spatial'][lib_id_to_add]['scalefactors']['pixel_per_um_real']
    res = adata.uns['spatial'][lib_id_to_add]['scalefactors']['resolution']
    
    if not hires:
        sf = adata.uns['spatial'][lib_id_to_add]['scalefactors']['tissue_lowres_scalef']
    else:
        sf = 1

    # change settings
    adata.uns['spatial'][lib_id_to_add]['scalefactors']['spot_diameter_fullres'] = res * ppm * sf # square width in pixel
        
    # scale coordinates
    adata.obsm['spatial'] *= sf
            
    if return_adata:
        return imgc, adata
    else:
        return imgc