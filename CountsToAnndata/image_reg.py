#!/usr/bin/env python

'''
Script to add images from high-resolution imaging round to dataset.

Usage:

    Input:
        - .csv file giving following parameters:
            - 

'''
## library
print("Import packages...", flush=True)
import sys
import os
import pandas as pd
import numpy as np
from glob import glob
import gc
from datetime import datetime
from pathlib import Path
import cv2
import napari
import matplotlib.pyplot as plt

# functions
def registration_qc_plot(adata, adata_trans, output_dir, unique_id, reg_channel_label):
    '''
    Plot to check how well the registration worked.
    '''
    # create plot to check success of registration
    plotpath = os.path.join(output_dir, unique_id + "_check_registration.png")

    # parameters
    idxs = adata.obs['id'].unique()
    c_names = ["Low resolution image", "High-resolution image"]
    gene = "Actb"

    fig, axs = plt.subplots(len(idxs), 2, figsize=(8*2, 6*len(idxs)))

    if len(idxs) > 1:
        axs.ravel()

    for r, idx in enumerate(idxs):
        for c, ad in enumerate([adata, adata_trans]):
            db.pl.spatial(ad, keys=gene, groupby='id', group=idx, image_key=reg_channel_label, 
                        lowres=False,
                        xlim=(1800,2000), ylim=(1800,2000), alpha=0.5, 
                        axis=axs[r+c], fig=fig)
            if r == 0:
                axs[r+c].set_title(c_names[c] + "\n" + gene, fontsize=12, fontweight='bold')
                
    fig.tight_layout()
    plt.savefig(plotpath, dpi=200)
            
    plt.show()

#######
## Protocol start
#######

if __name__ == "__main__":
    # read file location
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # get path of dbitx module and import functions
    module_path = os.path.abspath(os.path.join(script_dir, "../.."))
    if module_path not in sys.path:
        sys.path.append(module_path)

    import dbitx_funcs as db

    ## Set basic parameters
    register_hq = False

    ## Read parameters
    print("Starting Alignment script...", flush=True)

    # read parameters file
    #settings_file = sys.argv[1]
    #settings_file = "/home/jwirth/projects/experiments/37_43/CountsToAnndata/37_43_CtoA_params+hq_withvertices.csv"
    settings_file = "C:\Users\Johannes\Documents\homeoffice\37_38\CountsToAnndata\37_38_CtoA_params+hq_wovertices_ho.csv"

    print("Reading batch parameters from {}".format(settings_file))
    lines = open(settings_file).readlines()
    param_start = [i for i, line in enumerate(lines) if line.startswith(">parameters")][0]
    dir_start = [i for i, line in enumerate(lines) if line.startswith(">directories")][0]

    # read settings file
    settings = pd.read_csv(settings_file, header=None)

    # read parameters
    parameters = pd.read_csv(settings_file,
                        nrows=dir_start-1).dropna(how='all', axis=0).dropna(how='all', axis=1)

    # check headers of parameters file and set category as index
    param_headers = ["category", "value"]
    assert np.all([elem in parameters.columns for elem in param_headers]), \
        "Parameters section does not have correct headers ({})".format(param_headers)                       

    parameters = parameters.set_index('category')

    # read directories
    n_headers_dir = len(lines[dir_start].split(",")) # get number of headers in directory line
    directories = pd.read_csv(settings_file,
                            skiprows=dir_start, usecols=range(1,n_headers_dir)).dropna(how='all', axis=0)

    ## Check if all necessary parameters are in the file
    param_cats = ["n_channels", "spot_width", "frame", 
        "align_images:align_channel", "align_images:dapi_channel", "hq_images:channel_names", "hq_images:channel_labels"]
    dir_cats = ["experiment_id", "unique_id", "main_dir", 
        "input_transcriptome", "align_images", "hq_images", "output", 
        "vertices_x", "vertices_y"]

    assert np.all([elem in parameters.index for elem in param_cats]), \
        "Not all required categories found in parameter section {}".format(param_cats)
    assert np.all([elem in directories.columns for elem in dir_cats]), \
        "Not all required column headers found in directory section {}".format(dir_cats)
    assert ~np.any([pd.isnull(parameters.loc[cat, "value"]) for cat in param_cats]), \
        "Not all required categories in parameter section have a value.\nRequired categories are: ({})".format(param_cats)

    # determine extra categories which are added later to adata.obs
    extra_cats_headers = [elem for elem in directories.columns if elem not in dir_cats]
    extra_cats_headers = ["experiment_id"] + extra_cats_headers

    ## extract image parameters
    # determine alignment channel
    image_cats = ["align_images:align_channel", "align_images:dapi_channel", 
        "hq_images:channel_names", "hq_images:channel_labels"]

    if not np.any([pd.isnull(parameters.loc[elem, "value"]) for elem in image_cats]):
        alignment_channel = parameters.loc["align_images:align_channel", "value"]
        dapi_channel = parameters.loc["align_images:dapi_channel", "value"]

        # get channel names of hq images
        hq_ch_names = str(parameters.loc["hq_images:channel_names", "value"]).split(" ")
        hq_ch_labels = str(parameters.loc["hq_images:channel_labels", "value"]).split(" ")

        # determine channel on which registration is done most probably the dapi
        reg_id = [i for i, elem in enumerate(hq_ch_names) if "*" in elem][0]
        #reg_channel_name = hq_ch_names[reg_id]
        reg_channel_label = hq_ch_labels[reg_id]
        #reg_channel_name = [elem for elem in hq_ch_names if "*" in elem][0]


        # check if channel names and labels have same length
        assert len(hq_ch_names) == len(hq_ch_labels), \
            "Number of channel_names ({}) and channel_labels ({}) differ.".format(len(hq_ch_names), len(hq_ch_labels))

    # get ids of vertices_x and vertices_y
    vertx_id = directories.columns.tolist().index('vertices_x')
    verty_id = directories.columns.tolist().index('vertices_y')

    # create full directories from main_dir and the relative paths
    for cat in ["input_transcriptome", "align_images", "hq_images"]:
        directories[cat] = [os.path.join(m, p) for m, p in zip(directories["main_dir", cat])]
    
    # create output file names from output_dir
    directories["output"] = [os.path.join(m, p, 
        "{}_{}_adata_raw_with_images.h5ad".format(e, u)) for e, u, m, p in zip(
            directories["experiment_id"],
            directories["unique_id"],
            directories["main_dir"],
            directories["output_dir"]
            )]

    # check if all input matrix files exist
    try:
        assert np.all([os.path.isfile(f) for f in directories["input_transcriptome"]]), \
            "Not all input transcriptome files exist."
    except AssertionError as e:
        # check which files are missing
        missing_files = [f for f in directories["input_transcriptome"] if not os.path.isfile(f)]
        print("{} Following transcriptome files are missing: {}".format(e, missing_files))
        #sys.exit()
        exit()

    # check if all input images exist
    assert np.all([os.path.isfile(img) for img_d in directories['align_images'] for img in glob(img_d)]), \
        "Not all alignment image input files exist."

    assert np.all([os.path.isfile(img) for img_d in directories['hq_images'] for img in glob(img_d)]), \
        "Not all hq image input files exist."

    # check if all output names are unique
    assert len(np.unique(directories["output"])) == len(directories["output"]), \
        "Output files are not unique. This would cause that one file is overwritten by another: {}".format(directories["output"])

    # assert that for each dataset there are either both hq_image and align_image or None of both
    assert np.all([pd.notnull(a) == pd.notnull(b) for a, b in zip(directories["align_images"], directories["hq_images"])]), \
        "For some datasets only alignment images OR hq images were given.Either both are given or None.\n" \
            "If you wish to use the alignment images as hq images put the alignment directories also into\n" \
                "the hq image section. If you wish to create an anndata without image data leave both sections empty."

    ## prepare image alignment
    n_datasets = len(directories)
    vertices_list = [None]*n_datasets

    # create path for new settings_file to include vertices
    settings_new = settings_file.rstrip(".csv") + "_" + f"{datetime.now():%Y%m%d}" + "_withvertices.csv"

    print("Processing {} datasets".format(n_datasets))
    #for i in range(0, n_datasets):
    for i, dirs in directories.iterrows():
        print("{} : Start processing dataset {} of {}".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, n_datasets), 
            flush=True)

        # extract directories for this dataset
        #dirs = directories.loc[i, :]

        if len(extra_cats_headers) > 0:
            extra_cats = dirs[extra_cats_headers]
        else:
            extra_cats = None

        # get parameters for this dataset
        matrix_file = dirs["input_transcriptome"]
        output_dir = dirs["output_dir"]
        output_file = dirs["output"]
        #output_dir = os.path.dirname(output_file)
        output_filename = os.path.basename(output_file)
        tmp_dir = os.path.join(output_dir, "tmp")
        unique_id = dirs["experiment_id"] + "_" + dirs["unique_id"]

        # check what images are given for this dataset
        images_given = pd.notnull(dirs["align_images"]) and pd.notnull(dirs["hq_images"])

        # create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if images_given:
            # check if the vertices are given in the settings file
            vertices_not_given = pd.isnull(dirs["vertices_x"]) or pd.isnull(dirs["vertices_y"])

            # get image directories
            align_images = glob(dirs["align_images"])
            hq_images = glob(dirs["hq_images"])

            # check if number of images matches number of channel names
            assert len(hq_images) == len(hq_ch_names), \
                "Number of detected hq images [{}] does not match number of channel names [{}] in parameters file.".format(
                    len(hq_images), len(hq_ch_names))

            # check whether paths to alignment images and hq images are identical
            if not len(set(align_images) & set(hq_images)) == len(align_images):
                register_hq = True
            
            # detect alignment marker image
            align_img = [i for i in align_images if alignment_channel in i][0]
            dapi_img = [i for i in align_images if dapi_channel in i][0]

            if vertices_not_given:
                # read alignment image
                alignment_image = cv2.imread(align_img)

                ### Select corner spots in alignment image using napari viewer
                # with napari.gui_qt():
                # https://napari.org/guides/stable/event_loop.html#intro-to-event-loop
                print("No vertices given. Select them from napari viewer.", flush=True)

                points_fetched = False
                while points_fetched is not True:
                    try:
                        viewer = napari.view_image(alignment_image, 
                            title="Select corner spots in alignment image {} of {} ".format(i+1, n_datasets))
                        napari.run()

                        # fetch vertices (center points at cross points of alignment channels)
                        assert "Points" in viewer.layers, "No Points selected. Select exactly 4 points as vertices."
                        corner_spots_center = viewer.layers["Points"].data.astype(int)

                        assert len(corner_spots_center) == 4, "Number of selected points is not correct. Select exactly 4 points as vertices."

                        points_fetched = True
                    except AssertionError as k:
                        print(k, flush=True)
                        #print("No or not the right number of points (4) was selected. Try again.")
                        pass

                # collect corner spot
                vertices_list[i] = corner_spots_center

                # save information about vertices in settings file
                # save y coordinates (row coordinates)
                settings.loc[dir_start+i+1, verty_id+1] = " ".join([str(elem[0]) for elem in vertices_list[i]])
                # save x coordinates (column coordinates)
                settings.loc[dir_start+i+1, vertx_id+1] = " ".join([str(elem[1]) for elem in vertices_list[i]])
                # save settings file with coordinates of vertices
                settings.to_csv(settings_new, index=None, header=None)
                print("Vertices selected and saved.", flush=True)

            else:
                # extract coordinates from directory input
                xs = [int(elem) for elem in directories.loc[i, "vertices_x"].split(" ")]
                ys = [int(elem) for elem in directories.loc[i, "vertices_y"].split(" ")]

                # make sure four vertices are given
                assert len(xs) == 4 and len(ys) == 4, \
                    "Not enough x- [{}] and/or y-values [{}] given.".format(len(xs), len(ys)) 

                # add extracted coordinates to list of vertices
                vertices_list[i] = np.array([[a, b] for a, b in zip(ys, xs)])

            # read images
            print("{} : Read images...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
                flush=True)
            
            if register_hq:
                image_paths = [dapi_img]
                channel_labels = [reg_channel_label]
            else:
                image_paths = hq_images # is equal to align_images in this case
                channel_labels = hq_ch_labels
            
            images = [cv2.imread(d, -1) for d in image_paths]
        else:
            images = None
            channel_labels = None

        ### Generation of squidpy formatted anndata object
        print("{} : Generate anndata object from matrix file and images...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
            flush=True)

        if register_hq:
            # create tmp directory
            Path(tmp_dir).mkdir(parents=True, exist_ok=True)

            # determine output file for lowres adata object
            lowres_outfile = os.path.join(tmp_dir, output_filename.replace(".h5ad", "_lowres.h5ad"))
            return_adata = True
        else:
            lowres_outfile = output_file
            return_adata = False


        adata = db.dbitseq_to_squidpy(matrix_path=matrix_file, images=images,
            resolution=int(parameters.loc["spot_width"]), 
            unique_id=unique_id, extra_categories=extra_cats,
            n_channels=int(parameters.loc["n_channels"]), 
            frame=int(parameters.loc["frame"]),
            #dbitx=False, 
            labels=channel_labels, vertices=vertices_list[i], 
            savepath=lowres_outfile, return_adata=return_adata)

        ### Registration of hiqh-quality imaging data
        if register_hq:
            # create image_dir_dict
            image_dir_dict = {}
            for i, n in enumerate(hq_ch_names):
                selected_image_path = [elem for elem in hq_images if n.strip("*") in os.path.basename(elem)][0]
                image_dir_dict[unique_id + "_" + str(hq_ch_labels[i])] = selected_image_path

            # start transformation
            adata_trans = db.im.register_adata_coords_to_new_images(adata, groupby='id', image_dir_dict=image_dir_dict, 
                                            reg_channel=reg_channel_label, in_place=False, debug=False, do_registration=True)
            
            adata_trans.write(output_file)

            # plot registration QC plot
            registration_qc_plot(adata, adata_trans, output_dir, unique_id, reg_channel_label)

    print("{} : Finished all datasets. Vertices saved into {}".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}", settings_new), 
        flush=True)
