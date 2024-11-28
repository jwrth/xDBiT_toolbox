import pandas as pd
from collections import OrderedDict
import numpy as np
from typing import Optional, Tuple, Union, List, Dict, Any
import anndata
import scanpy as sc
from .. import utils
from ..exceptions import ModuleNotFoundOnWindows

def spatialde_run(adata, layer=None, run=True, normalize=True, output_name='spatialde', use_raw=False):

    import NaiveDE
    import SpatialDE

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

    spatialde_sum = dict()

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


def standard_preprocessing(adata_in, 
    hvg_batch_key=None, hvg_flavor='seurat', hvg_n_top_genes=None,
    do_lognorm=True, regress_out=None, filter_hvg=False, 
    dim_reduction=True, umap=True, tsne=True,
    batch_correction_key=None,  batch_correction_method="scanorama", verbose=True,
    tsne_lr=1000, tsne_jobs=8,
    random_state=0,
    **kwargs):
    '''
    Function to perform standard preprocessing on ST data. Adapted from Squidpy Napari Tutorial.
    '''

    available_batch_methods = ["bbknn", "scanorama"]
    assert batch_correction_method in available_batch_methods, \
        "`batch_correction_method` {} not available. Available methods are: {}".format(batch_correction_method, available_batch_methods)

    if hvg_flavor in ["seurat", "cell_ranger"]:
        hvg_layer = None
    elif hvg_flavor == "seurat_v3":
        hvg_layer = "counts" # seurat v3 method expects counts data

        # n top genes must be specified for this method
        if hvg_n_top_genes is None:
            print("HVG computation: For flavor {} `hvg_n_top_genes` is mandatory".format(hvg_flavor)) if verbose else None
            return
    else:
        print("Unknown value for `hvg_flavor`: {}. Possible values: {}".format(hvg_flavor, ["seurat", "cell_ranger", "seurat_v3"])) if verbose else None


    if batch_correction_method == "bbknn":
        try:
            import bbknn
        except ModuleNotFoundError as e:
            raise ModuleNotFoundOnWindows(e)

    adata = adata_in.copy()

    if do_lognorm:

        # if count matrix does not consist of raw counts abort preprocessing
        if not np.all(np.modf(adata.X)[0] == 0):
            print("`adata.X` does not contain raw counts. Preprocessing aborted.")
            return

        # store raw counts in layer
        print("Store raw counts in adata.layers['counts']...") if verbose else None
        adata.layers['counts'] = adata.X.copy()

        # preprocessing according to napari tutorial in squidpy
        print("Normalization, log-transformation...") if verbose else None
        sc.pp.normalize_total(adata)
        adata.layers['norm_counts'] = adata.X.copy()
        sc.pp.log1p(adata)

    if hvg_batch_key is None:
        print("Calculate highly-variable genes across all samples using {} flavor...".format(hvg_flavor)) if verbose else None
    else:
        print("Calculate highly-variable genes per batch key {} using {} flavor...".format(hvg_batch_key, hvg_flavor)) if verbose else None

    sc.pp.highly_variable_genes(adata, batch_key=hvg_batch_key, flavor=hvg_flavor, layer=hvg_layer, n_top_genes=hvg_n_top_genes)

    if filter_hvg:
        print("Filter for highly-variable genes...") if verbose else None
        # add the normalized and log data to raw
        adata.raw = adata
        # Filter for highly variable genes
        adata = adata[:, adata.var.highly_variable].copy()

    if regress_out is not None:
        print("Regress out {}...".format(regress_out)) if verbose else None
        sc.pp.regress_out(adata, regress_out)

    if dim_reduction:
        if batch_correction_key is None:
            # dimensionality reduction
            print("Dimensionality reduction...") if verbose else None
            sc.pp.pca(adata)
            if umap:
                sc.pp.neighbors(adata)
                sc.tl.umap(adata, random_state=random_state)
            if tsne:
                sc.tl.tsne(adata, n_jobs=tsne_jobs, learning_rate=tsne_lr, random_state=random_state)

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
            if batch_correction_method == "bbknn":
                print("Batch correction using {} for {}...".format(batch_correction_method, batch_correction_key)) if verbose else None
                bbknn.bbknn(adata, batch_key=batch_correction_key, metric='euclidean') # is used as alternative to sc.pp.neighbors

                # dim reduction with corrected data
                print("Dimensionality reduction with batch corrected data...") if verbose else None
                sc.tl.umap(adata, random_state=random_state)
                sc.tl.tsne(adata, random_state=random_state)
            elif batch_correction_method == "scanorama":
                print("Batch correction using {} for {}...".format(batch_correction_method, batch_correction_key)) if verbose else None
                hvgs = list(adata.var_names[adata.var['highly_variable']])
                adata = scanorama(adata, batch=batch_correction_key, hvg=hvgs, verbose=False, **kwargs)

                # find neighbors
                sc.pp.neighbors(adata, use_rep="X_scanorama")
                sc.tl.umap(adata, random_state=random_state)
                sc.tl.tsne(adata, use_rep="X_scanorama", random_state=random_state)

        # clustering
        print("Leiden clustering...") if verbose else None
        sc.tl.leiden(adata)
        
    # if not in_place:
    #     return adata
    return adata


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