
import sys
import os
module_path = os.path.abspath(os.path.join("D:/Users/johannes.wirth/Documents/GitHub/DbitX/"))
if module_path not in sys.path:
    sys.path.append(module_path)

import dbitx as db
import scanpy as sc
import napari
import cv2
import json

def main():

    # input parameters
    alignment_id = 1
    image_labels = ['bf', 'align', 'phalloidin', 'dapi']
    n_channels = 38
    spot_width = 50 # µm
    frame = 100 # pixel
    matrix_name = "DGE_matrix_with_introns_min100.txt.gz"

    # output name
    output_prefix = "37_30_adata_raw_with_images"

    # image path
    
    image_dir = "N:/01 HPC/03 Team Meier/10_Resources/08_Johannes Wirth/Nextcloud/DbitX/data/37_30/images/"
    #image_dir = "C:/Users/Johannes/Nextcloud/DbitX/data/37_30/images/"

    # transcriptome path
    #matrix_dir = "C:/Users/Johannes/Nextcloud/DbitX/data/37_30/matrices/"
    matrix_dir = "N:/01 HPC/03 Team Meier/10_Resources/08_Johannes Wirth/Nextcloud/DbitX/data/37_30/matrices/wells/"

    #output_dir = "C:/Users/Johannes/Nextcloud/DbitX/data/37_30/h5ad"
    output_dir = matrix_dir

    # get names of well folders
    image_wells = os.listdir(image_dir)
    matrix_wells = os.listdir(matrix_dir)

    # take only folders that are in both the matrix and the image path
    well_names = list(set(image_wells) & set(matrix_wells))
    well_names.sort()

    # iterate through wells and collect corner spots
    vertices_dict = {}
    for well_name in well_names:
        image_names = os.listdir(os.path.join(image_dir,well_name))

        # load images
        #images = [cv2.imread(os.path.join(image_dir, well_name, name), -1) for name in image_names]
        #alignment_image = images[alignment_id]

        alignment_image = cv2.imread(os.path.join(image_dir, well_name, image_names[alignment_id]), -1)

        ### Select corner spots in alignment image using napari viewer

        with napari.gui_qt():
            viewer = napari.view_image(alignment_image, title="Select corner spots in alignment image of well " + well_name)
            viewer.window._qt_window.raise_()

        # fetch vertices (center points of )
        corner_spots_center = viewer.layers["Points"].data.astype(int)

        # collect corner spot
        vertices_dict[well_name] = corner_spots_center

        # create and save adata objects for each well
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    #for well_name in well_names:
        images = [cv2.imread(os.path.join(image_dir, well_name, name), -1) for name in image_names]

        save_dir = os.path.join(output_dir, well_name)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        output_name = output_prefix + ".h5ad"
        savepath = os.path.join(save_dir, output_name)


        db.dbitseq_to_squidpy(matrix_path=os.path.join(matrix_dir, well_name, matrix_name), 
        images=images, labels=image_labels, vertices=vertices_dict[well_name], 
        resolution=spot_width, n_channels=n_channels, frame=frame, mode="Dbit-seq", savepath=savepath)

    # convert numpy array in list to make it json serializable
    for key in vertices_dict.keys():
        vertices_dict[key] = vertices_dict[key].tolist()
    with open(os.path.join(image_dir, 'vertices.json'), 'w') as fp:
        json.dump(vertices_dict, fp)

if __name__ == "__main__":
    main()