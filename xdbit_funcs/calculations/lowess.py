import numpy as np
from scipy import stats
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess as  sm_lowess

class lowess_confidence_intervals:
    '''
    Class to store confidence intervals of lowess prediction.
    '''
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

class lowess_prediction:
    '''
    Class to store results from lowess prediction.
    '''
    def __init__(self, values, stderr, smooths):
        self.values = values
        self.stderr = stderr
        self.smooths = smooths
        
    def confidence(self, alpha=0.05, percentile_method=False):
        if percentile_method:
            # use 2.5 and 97.5% percentiles to calculate the 95% confidence interval
            # This approach is mentioned here: https://acclab.github.io/bootstrap-confidence-intervals.html
            # However I am not sure if it is also valid for low numbers of bootstrapping cycles.
            lower = np.nanpercentile(self.smooths, 2.5, axis=1) #2.5 percent
            upper = np.nanpercentile(self.smooths, 97.5, axis=1) # 97.5 percent
        else:            
            # calculate 95% CI use formula for confidence interval
            self.smooths_mean = np.nanmean(self.smooths, axis=1)
            lower, upper = stats.norm.interval(1-alpha, loc=self.smooths_mean, scale=self.stderr)
            
        return lowess_confidence_intervals(lower, upper)
        

class lowess:
    '''
    Function to perform LOWESS regression and optionally calculate the confidence intervals using bootstrapping.
    
    Adapted from: https://james-brennan.github.io/posts/lowess_conf/
    '''
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.fitted = False
        
    def predict(self, newdata, stderror=False, verbose=False, K=100, **kwargs):
        # make sure the fit() function was run before
        assert self.fitted, "Values have not been fitted yet. Run .fit() first."
        
        # regularly sample it onto the grid (statsmodel does not provide the solution on interpolated values)
        # for the statistics later it is important that the same grid of x values is used
        if not verbose:
            # suppress runtime error for division by nan
            np.seterr(invalid='ignore')
            
        values = scipy.interpolate.interp1d(self.pred_x, self.pred_y, 
                                            fill_value='extrapolate')(newdata)
        
        if stderror:
            self.calc_stderror(newdata, K=K, **kwargs)
        else:
            self.stderr = None
            self.bootstrap_result = None
            
        return lowess_prediction(values, self.stderr, self.bootstrap_result)
    
    def bootstrap(self, x, y, newdata, sample_frac=0.5):
        
        samples = np.random.choice(len(x), int(len(x)*sample_frac), replace=True)
        
        y_s = y[samples]
        x_s = x[samples]
        y_sm = sm_lowess(y_s, x_s, frac=self.frac, it=5,
                         return_sorted = False)
        # regularly sample it onto the grid (statsmodel does not provide the solution on interpolated values)
        # for the statistics later it is important that the same grid of x values is used
        return scipy.interpolate.interp1d(x_s, y_sm, fill_value='extrapolate')(newdata)
    
    def calc_stderror(self, newdata, sample_frac=0.5, K=100, **kwargs):
        '''
        Interesting input on topic from:
            - https://acclab.github.io/bootstrap-confidence-intervals.html
            - https://stackoverflow.com/questions/28242593/correct-way-to-obtain-confidence-interval-with-scipy
        '''        
        # calculate confidence interval using bootstrapping approach
        self.bootstrap_result = np.stack([self.bootstrap(self.x, self.y, newdata, sample_frac=sample_frac, **kwargs) for k in range(K)]).T
        
        # calc mean and stderr of smooths
        self.stderr = np.nanstd(self.bootstrap_result, axis=1, ddof=0) # OR std
        
    def fit(self, frac=0.3, **kwargs):
        self.pred_x, self.pred_y = sm_lowess(endog=self.y, exog=self.x, frac=frac, 
                               it=3, return_sorted=True, **kwargs).T
        
        # save that object was fitted
        self.fitted = True
        self.frac = frac # save frac setting for functions run later
