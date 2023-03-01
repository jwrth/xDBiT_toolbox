from scipy.interpolate import UnivariateSpline
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale
from anndata import AnnData
import matplotlib.pyplot as plt
from ..preprocessing import standard_preprocessing
from ..readwrite import save_and_show_figure
import gc
from tqdm import tqdm
import numpy as np
import pandas as pd

def generate_axis(
    adata: AnnData,
    gene_to_sort: str,
    groupby: str,
    umap_key: str = 'X_umap',
    obs_key: str = 'umap_spline',
    max_cols: int = 6,
    k_spline: int = 3,
    s: int = 5000,
    plot: bool = True,
    savepath: bool = None,
    save_only: bool = False,
    **kwargs
):
    from ..tools import get_nrows_maxcols, extract_groups
    
    groups = adata.obs[groupby].unique()
    # check whether to calculate the UMAP per group or use existing UMAP calculation saved in adata.obsm[umap_key]
    if umap_key is None:
        # calculate UMAP separately per group
        umap_per_group = True
        
        obsm_umap = {}
        for group in tqdm(groups):
            #print("Processing group {}...".format(group))
            # do preprocessing for one group
            ad = extract_groups(adata, groupby=groupby, groups=group)
            ad.uns['log1p']['base'] = None
            adpp = standard_preprocessing(ad,
                                          hvg_n_top_genes=2000,
                                          do_lognorm=False,
                                          dim_reduction=True,
                                          umap=True, tsne=False,
                                          verbose=False)
            # collect results
            obsm_umap[group] = adpp.obsm['X_umap']
            
            # free RAM
            del adpp
            gc.collect()
    else:
        assert umap_key in adata.obsm, "`umap_key` not in `adata.obsm`."
        umap_per_group = False
        obsm_umap = None

    # decide whether to generate a plot or not
    if plot:
        n_plots, nrows, ncols = get_nrows_maxcols(groups, max_cols=max_cols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6*nrows))

        if len(axs.shape) > 1:
            axs = axs.ravel()

    # start processing
    intpol = {}
    for i, group in enumerate(groups):
        # generate selection mask
        mask = adata.obs[groupby] == group

        # extract points, indices and expression of gene
        if umap_per_group:
            pts = obsm_umap[group]
        else:
            pts = adata.obsm[umap_key][mask]
        
        idxs = adata.obs_names[mask]
        expr = adata.X[:, adata.var_names.get_loc(gene_to_sort)][mask]

        # perform PCA to rotate datapoints and always have the longest data axis as x
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(pts)

        # sort points and indices by x values
        sortmask = pcs[:, 0].argsort()
        pcs = pcs[sortmask]
        idxs = idxs[sortmask]
        expr = expr[sortmask]

        # extact coordinates
        x = pcs[:, 0]
        y = pcs[:, 1]

        # adjust x values according to expression of selected gene
        slope = np.polyfit(x, expr, 1)[0]
        slope = -1 if slope < 0 else 1
        x *= slope

        # sort again by x
        sortmask = np.argsort(x)
        y = y[sortmask]
        expr = expr[sortmask]
        idxs = idxs[sortmask]
        x.sort()

        # calculate spline
        spl = UnivariateSpline(x, y, k=k_spline, s=s)
        ys = spl(x)

        # create array of spline points
        Xs = np.array([x, ys]).T

        # collect results
        intpol[group] = pd.DataFrame(Xs, index=idxs, columns=['x', 'y'])

        if plot:
            axs[i].scatter(x, y, c=expr)
            axs[i].plot(x, ys, 'r', lw=3)
            axs[i].set_title(group)

    if plot:
        save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, **kwargs)
        
    x_new = {}
    for mid, Xs in intpol.items():
        # calculate distance of consecutive points
        d = np.diff(Xs.values, axis=0)
        segdists = np.hypot(d[:,0], d[:,1])

        # calculate cumulative sum of distances and add 0 at beginning
        cumsum = np.insert(np.cumsum(segdists), 0, 0)

        # min max scaling
        cumsum = minmax_scale(cumsum)

        x_new[mid] = pd.Series(cumsum, index=Xs.index)

    # concatenate results and reshape index
    x_new = pd.concat(x_new)

    x_new.index = x_new.index.droplevel(0)

    # add data to obs
    adata.obs[obs_key] = x_new
    
    