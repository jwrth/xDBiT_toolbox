import cv2
import math
import numpy as np
from scipy import stats
import pandas as pd

def centroid(contour):
    '''
    Calculate the centroid of a contour.
    '''
    # calculate moments
    M = cv2.moments(contour)
    
    # calculate coordinates of centroid
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    return cx, cy

def circularity(contour, area=None, peri=None):
    '''
    Calculate circularity of contour using following equation:
    circularity = 4pi(area/perimeter^2)
    '''
    if area is None:
        area = cv2.contourArea(contour)
    if peri is None:
        peri = cv2.arcLength(contour, True)
    
    circularity = 4 * math.pi * (area / peri**2)
    return circularity

def distance_weight(distances, sigma=20):

    wts = [np.exp(-elem / sigma) for elem in distances]
    sum_wts = np.sum(wts)

    return [w / sum_wts for w in wts]



