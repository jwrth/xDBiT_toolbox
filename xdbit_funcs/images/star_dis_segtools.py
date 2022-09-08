from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.transform import rescale
from skimage import color
import numpy as np
import stardist
from stardist.models import StarDist2D 
from stardist import random_label_cmap, _draw_polygons, export_imagej_rois
from csbdeep.utils import Path, normalize
from skimage.measure import regionprops
import pandas as pd
from dbitx.calculations._calc import coord_to_pixel 
import matplotlib.pyplot as plt

np.random.seed(6)
lbl_cmap = random_label_cmap()


model_versatile = StarDist2D.from_pretrained('2D_versatile_fluo')
model = model_versatile


#Compute all pixel spot's coordinate in adata and return a panda dataframe with index and pixel coordinate and a list of array with each y,x pixel coordinate
def spot_pixelcoord(image_metadata, adata, resolution):
    scale = image_metadata["pixel_per_um"]
    offset_row =image_metadata[ 'pivot_spot'][0]
    offset_col = image_metadata[ 'pivot_spot'][1]
    adata.obs['pixel_spot_row']= np.array([coord_to_pixel(c, resolution, scale, offset_row) for c in adata.obs['array_row']])
    adata.obs['pixel_spot_col']= np.array([coord_to_pixel(c, resolution, scale, offset_col) for c in adata.obs['array_col']])

    
    pixel_spot_col = pd.DataFrame(adata.obs['pixel_spot_col'])
    pixel_spot_row = pd.DataFrame(adata.obs['pixel_spot_row'])
    try_ = pd.concat([pixel_spot_col,pixel_spot_row],axis=1 )
    try_= try_.drop_duplicates()
    pixel_spot_coord = try_
    pixel_spot_coord_list= []

    for i in range(len(pixel_spot_coord.index)):
        pixel_spot_coord_i = np.array(pixel_spot_coord.loc[pixel_spot_coord.index[i],['pixel_spot_col','pixel_spot_row']])
        pixel_spot_coord_list.append(pixel_spot_coord_i)
    return pixel_spot_coord, pixel_spot_coord_list


def image_preprocessing(img,layer):
    image_gray=color.rgb2gray(img[layer])
    return image_gray

def plot_segmentation(labels, img, show_dist,details):
    #Plotting the stardist segmentation result with this function
    plt.figure(figsize=(13,10))
    img_show = img if img.ndim==2 else img[...,0]
    coord, points, prob = details['coord'], details['points'], details['prob']
    plt.subplot(121); plt.imshow(img_show, cmap='gray'); plt.axis('off')
    a = plt.axis()
    _draw_polygons(coord, points, prob, show_dist=show_dist)
    plt.axis(a)
    plt.subplot(122); plt.imshow(img_show, cmap='gray'); plt.axis('off')
    plt.imshow(labels, cmap=lbl_cmap, alpha=0.5)
    plt.tight_layout()
    plt.show()

def crop_image_segmentation_nuclei(index,image_gray, pixel_spot_coord_list,image_metadata,pixel_spot_coord ,
                                   layer= 'dapi', model = model,
                                   returnlabel=False, figure = True, show_dist = True, Count = False):
# cropping the image and convert to grayscale according to a panda dataframe stored the pixel coordinates.
    x = int(pixel_spot_coord_list[index][0])
    y= int( pixel_spot_coord_list[index][1])
    dis= int(image_metadata['spot_diameter_fullres'])
    image_i= image_gray[y:y+dis,x:x+dis]
#then denoise images    
    patch_kw = dict(patch_size=10,      # 5x5 patches
                patch_distance=15,  # 13x13 search area
                multichannel=False)
    sigma_est = np.mean(estimate_sigma(image_i, multichannel=False))
    image_i = denoise_nl_means(image_i, h=0.6 * sigma_est, sigma=sigma_est,
                                 fast_mode=False, **patch_kw)
#Normalise value for stardist model    
    image_i_norm = normalize(image_i, 3,99)
#Predict value with pretrained model, return mask in labels and other infromatino in details    
    labels, details = model.predict_instances(image_i_norm)
#Count the number of objects in label to get the number of nuclei per iamge.    
    objects = regionprops(label_image = labels, intensity_image=image_i)
    print("The index {} has {} of labelled nuclei ".format(pixel_spot_coord.index[index],len(objects)))
#Specify return segmentation figure or the image mask
    if figure == True:
        plot_segmentation(labels = labels, img=image_i, show_dist = show_dist, details = details)
    if returnlabel == True:
        return labels
    if Count == True:
        return len(objects)

    
def validation(index, img, pixel_spot_coord,image_metadata, 
               pixel_spot_coord_list,coordinates=None, layer = 'dapi'):
    #This function is to segment individual spot with the index in obs or coordinate.
    if coordinates is not None:
        index = pixel_spot_coord.index.get_loc(coordinates)
    print(pixel_spot_coord.index[index])
    print(pixel_spot_coord_list[index])
    print(pixel_spot_coord.loc[pixel_spot_coord.index[index],['pixel_spot_col','pixel_spot_row']])
    image_gray = color.rgb2gray(img['dapi'])

    x = int(pixel_spot_coord_list[index][0])
    y= int( pixel_spot_coord_list[index][1])
    dis= int(image_metadata['spot_diameter_fullres'])
    image_i= image_gray[y:y+dis,x:x+dis]
    print("index location = {}".format(index))
    f, axarr = plt.subplots(1,2)
    axarr[0].imshow(image_i)
    axarr[1].imshow(img['bf'][y:y+dis,x:x+dis])
    crop_image_segmentation_nuclei(index = index, image_gray = image_gray, layer = 'dapi', pixel_spot_coord_list = pixel_spot_coord_list,pixel_spot_coord =pixel_spot_coord,image_metadata=image_metadata )
    

def bulk_get_nuclei_count_perspots(img, pixel_spot_coord_list,pixel_spot_coord, adata=None, image_metadata=None ):
    #Batch processing multiple spots to return the number of segmented nuclei from the pixel_spot_coord_list. 
    
    if pixel_spot_coord_list is None or pixel_spot_coord is None:
        pixel_spot_coord,pixel_spot_coord_list=spot_pixelcoord(
            image_metadata=image_metadata,adata = adata)
    
    image_gray = color.rgb2gray(img['dapi'])
    print('Total spots = {}'.format(len(pixel_spot_coord_list)))    
    count=[]
    i = 0
    total_task = len(pixel_spot_coord_list)
    for index in range(len(pixel_spot_coord_list)):
        print(index)
        count.append(
            crop_image_segmentation_nuclei(
                index = index,
                image_gray=image_gray,returnlabel=False,figure=False,
                Count=True,
                pixel_spot_coord_list= pixel_spot_coord_list,
                pixel_spot_coord =pixel_spot_coord,
                image_metadata = image_metadata))
        i=i+1
        if i == int(total_task*0.3):
            print("30% done")
        if i == int(total_task*0.5):
            print("50% done")
        if i ==int(total_task*0.7):
            print("70% done")
        if i == int(total_task*0.9):
            print("90% done")
    
    return count