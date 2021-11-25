import numpy as np
from math import sqrt
from scipy import ndimage
from scipy.spatial import distance as dist
import cv2


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