import numpy as np
import pandas as pd
from math import sqrt
from scipy.spatial import distance as dist
#from sympy import lowergamma
from ..exceptions import ModuleNotFoundOnWindows, ImportErrorLoess
from typing import Optional, Tuple, Union, List, Dict, Any, Literal
from .lowess import lowess
from warnings import warn


def minDistance(A, B, E) :
	# Python3 implementation of the approach from https://www.geeksforgeeks.org/minimum-distance-from-a-point-to-the-line-segment-using-vectors/   
	# Function to return the minimum distance  
	# between a line segment AB and a point E 

    # vector AB  
    AB = [None, None];  
    AB[0] = B[0] - A[0];  
    AB[1] = B[1] - A[1];  
  
    # vector BP  
    BE = [None, None]; 
    BE[0] = E[0] - B[0];  
    BE[1] = E[1] - B[1];  
  
    # vector AP  
    AE = [None, None]; 
    AE[0] = E[0] - A[0]; 
    AE[1] = E[1] - A[1];  
  
    # Variables to store dot product  
  
    # Calculating the dot product  
    AB_BE = AB[0] * BE[0] + AB[1] * BE[1];  
    AB_AE = AB[0] * AE[0] + AB[1] * AE[1];  
  
    # Minimum distance from  
    # point E to the line segment  
    reqAns = 0;  
  
    # Case 1  
    if (AB_BE > 0) : 
  
        # Finding the magnitude  
        y = E[1] - B[1];  
        x = E[0] - B[0];  
        reqAns = sqrt(x * x + y * y);  
  
    # Case 2  
    elif (AB_AE < 0) : 
        y = E[1] - A[1];  
        x = E[0] - A[0];  
        reqAns = sqrt(x * x + y * y);  
  
    # Case 3  
    else: 
  
        # Finding the perpendicular distance  
        x1 = AB[0];  
        y1 = AB[1];  
        x2 = AE[0];  
        y2 = AE[1];  
        mod = sqrt(x1 * x1 + y1 * y1);  
        reqAns = abs(x1 * y2 - y1 * x2) / mod;  
      
    return reqAns;

def rotation_angle(upper, lower):

    '''
    Returns counterclockwise rotation angle µ of two points (format [x,y]) with one point above the other:

        x upper
        |\
        | \
        |µ \
        |__/\
        |    x lower
    
    '''

    dx = lower[0] - upper[0]
    dy = lower[1] - upper[1]

    alpha = np.rad2deg(np.arctan(dx/dy))
    return alpha

#def Distance(p1, p2):
#    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def dist_points(p1, p2):
    dist = np.sqrt( (p2[1] - p1[1])**2 + (p2[0] - p1[0])**2 )
    return dist

# def rotateImage(img, angle, pivot, imagetype="grayscale"):
#     padX = [img.shape[1] - pivot[1], pivot[1]]
#     padY = [img.shape[0] - pivot[0], pivot[0]]
#     if imagetype == "grayscale":
#         imgP = np.pad(img, [padY, padX], 'constant')
#     elif imagetype == "rgb":
#         imgP = np.pad(img, [padY, padX, [0, 0]], 'constant')
#     else:
#         print("Unknown image type.")
#         return
#     imgR = ndimage.rotate(imgP, angle, reshape=False)
#     return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]

# function to compute the mid point of a square
def centroid_mean(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    centroid = [int(sum(x) / len(points)), int(sum(y) / len(points))]
    return centroid

# convert coordinates into µm
def coord_to_um(coord, resolution):
    x = coord * resolution * 2
    return x

# convert grid coordinates into pixel coordinates
def coord_to_pixel(coord, resolution, scale, offset):
    x = coord * resolution * 2 * scale + offset
    return x

# import the necessary packages
def order_points_clockwise(pts, mode='yx'):

    # check mode which specifies if the coordinates are in order xy (in opencv) or yx (python standard)
    if mode == 'yx':
        x = 1
        y = 0
    elif mode == 'xy':
        x = 0
        y = 1
    else:
        print(str(mode) + 'is unknown mode. Select either "xy" or "yx"')
        return

    # sort the points based on their x-coordinates
    xSorted = pts[np.argsort(pts[:, x]), :]
    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2, :]
    rightMost = xSorted[2:, :]
    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    leftMost = leftMost[np.argsort(leftMost[:, y]), :]
    (tl, bl) = leftMost
    # now that we have the top-left coordinate, use it as an
    # anchor to calculate the Euclidean distance between the
    # top-left and right-most points; by the Pythagorean
    # theorem, the point with the largest distance will be
    # our bottom-right point
    D = dist.cdist(tl[np.newaxis], rightMost, "euclidean")[0]
    (br, tr) = rightMost[np.argsort(D)[::-1], :]
    # return the coordinates in top-left, top-right,
    # bottom-right, and bottom-left order
    return np.array([tl, tr, br, bl], dtype="float32")


def smooth_fit(xs: np.ndarray, ys: np.ndarray, 
    #x_to_fit_on: Optional[List] = None,
    #dist_thrs: Optional[float] = None, 
    min: Optional[float] = None, 
    max: Optional[float] = None,
    #stepsize: Optional[float] = None,
    nsteps: Optional[float] = None,
    method: Literal["lowess", "loess"] = "lowess",
    stderr: bool = True,
    K: int = 100
    ):

    """Smooth curve using loess

    will perform curve fitting using skmisc.loess,
    points above 'dist_thrs' will be excluded.
    Parameters:
    ----------
    xs : np.ndarray
        x values
    ys : np.ndarray
        y values
    min: Optional[float]
        exclude (x,y) tuples were x < min
    max: Optional[float]
        exclude (x,y) tuples were x > max
    nsteps = Optional[float]
        Number of steps x is divided into for the prediction.
    method: str
        Which method to use for the smoothing. Options: "loess" or "lowess".
        "loess": Uses `skmisc.loess`. This is not implemented for Windows.
        "lowess": Uses `statsmodels.nonparametric.smoothers_lowess.sm_lowess` which is available on all tested platforms.
    stderr: bool
        Whether to calculate standard deviation of the prediction. Depends on the method used:
        "loess": Standard deviation returned by package is used.
        "lowess" Bootstrapping is used to estimate a confidence interval.
    K: int
        Only needed for `method="lowess"`. Determines number of bootstrapping cycles.
    
    Returns:
    -------
    A tuple with included x and y-values (xs',ys'), as well
    as fitted y-values (ys_hat) together with associated
    standard errors. The tuple is on the form
    (xs',ys',y_hat,std_err)

    From: https://github.com/almaan/ST-mLiver
    """
    
    # check method
    if method == "loess":
        loess = True
    elif method == "lowess":
        loess = False
    else:
        raise ValueError('Invalid `method`. Expected is one of: ["loess", "lowess"')

    if loess:
        try:
            from skmisc.loess import loess
            
        except ModuleNotFoundError as e:
            raise ModuleNotFoundOnWindows(e)

        except ImportError as e:
            #warn('{}`skmisc.loess` package not found. Used custom lowess impolementation instead.'.format(e))
            raise ImportErrorLoess(e)

    # sort x values
    srt = np.argsort(xs)
    xs = xs[srt]
    ys = ys[srt]

    # determine min and max x values and select x inside this range    
    if min is None:
        min = xs.min()
    if max is None:
        max = xs.max()    

    keep = (xs >= min) & (xs <= max)
    xs = xs[keep]
    ys = ys[keep]

    # generate loess class object
    if loess:
        ls = loess(xs, ys)
    else:
        ls = lowess(xs, ys)      

    # fit loess class to data
    ls.fit()

    # if stepsize is given determine xs to fit the data on
    if nsteps is not None:
        stepsize = (max - min) / nsteps
        #if stepsize is not None:
        xs_pred = np.asarray(np.arange(min, max+stepsize, stepsize))
        xs_pred = np.linspace(min, max, nsteps)
        xs_pred = xs_pred[(xs_pred < xs.max()) & (xs_pred > xs.min())]

    # predict on data
    if loess:
        pred =  ls.predict(xs_pred, stderror=stderr)
    else:
        pred =  ls.predict(xs_pred, stderror=stderr, K=K)
        
    # get predicted values
    ys_hat = pred.values
    
    if stderr:
        # get standard error
        stderr = pred.stderr
        conf = pred.confidence()
        lower = conf.lower
        upper = conf.upper
    else:
        lower = np.nan
        upper = np.nan
        
    df = pd.DataFrame({
        'x': xs_pred,
        #'y': ys,
        'y_pred': ys_hat,
        'std': stderr,
        'conf_lower': lower,
        'conf_upper': upper
    })

    #return (xs,ys,ys_hat,stderr)
    return df

def center_of_mass(values, weights):
    '''
    Calculate center of mass of an array of values with a corresponding array of weights.
    Can be also used to calculate the center of expression by using:
        - `values`: position along an axis
        - `weights`: expression value
    '''

    if np.sum(weights) > 0:
        com = np.average(values, weights=weights)
    else:
        com = np.nan
    return com

def cohens_d(a, b, paired=False, correct_small_sample_size=True):
    '''
    Function to calculate the Cohen's D measure of effect sizes.

    Function with correction was adapted from:
    https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/cohens-d/

    To allow measurement for different sample sizes following sources were used:
    https://stackoverflow.com/questions/21532471/how-to-calculate-cohens-d-in-python
    https://en.wikipedia.org/wiki/Effect_size#Cohen's_d

    For paired samples following websites were used:
    https://www.datanovia.com/en/lessons/t-test-effect-size-using-cohens-d-measure/
    The results were tested here: https://statistikguru.de/rechner/cohens-d-gepaarter-t-test.html
    '''
    if not paired:
        # calculate parameters
        mean1 = np.mean(a)
        mean2 = np.mean(b)
        std1 = np.std(a, ddof=1)
        std2 = np.std(b, ddof=1)
        n1 = len(a)
        n2 = len(b)
        dof = n1 + n2 - 2 # degrees of freedom
        SDpool = np.sqrt(((n1-1) * std1**2 + (n2 - 1) * std2**2) / dof) # pooled standard deviations

        if SDpool == 0:
            d = np.nan
        else:
            d = (mean1 - mean2) / SDpool

        n = np.min([n1, n2])
        if correct_small_sample_size and n < 50:
            # correct for small sample size
            corr_factor = ((n - 3) / (n-2.25)) * np.sqrt((n - 2) / n)
            d *= corr_factor

    else:
        assert len(a) == len(b), "For paired testing the size of both samples needs to be equal."
        diff = np.array(a) - np.array(b)
        d = np.mean(diff) / np.std(diff)
        
    return d

def jaccard_dist(a, b):
    intersect = set(a) & set(b)
    union = set(a) | set(b)
    j = 1 - (len(intersect) / len(union))
    return j