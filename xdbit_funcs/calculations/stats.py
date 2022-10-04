import numpy as np
from statsmodels.stats.multitest import fdrcorrection

def benjamini_hochberg_nan(pvals, **kwargs):
    '''
    Version of Benjamini-Hochberg correction for multiple testing that ignores NaNs.
    '''

    # create empty array to exclude nans from analysis
    mask = np.isfinite(pvals)
    _pvals = np.empty(pvals.shape)
    _pvals.fill(np.nan)
    
    # do correction
    _pvals[mask] = fdrcorrection(pvals[mask])[1]
    
    return _pvals