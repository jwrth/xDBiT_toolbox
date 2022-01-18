#!/usr/bin/env python

'''
Script to align matrices and images from DbitX pipeline.

Usage:

    Input:
        - .csv file giving following parameters:
            - 

'''

print("Import packages...", flush=True)
import sys
import os
import gc

# read file location
script_dir = os.path.dirname(os.path.realpath(__file__))

# get path of dbitx module and import functions
module_path = os.path.abspath(os.path.join(script_dir, ".."))
if module_path not in sys.path:
    sys.path.append(module_path)

import dbitx_funcs as db
#import scanpy as sc
import napari
import cv2
#import json
import pandas as pd
import numpy as np
from glob import glob
from pathlib import Path
from datetime import datetime

## Read parameters
print("Starting Alignment script...", flush=True)

# read parameters file
print("Reading batch parameters from {}".format(sys.argv[1]))
settings_file = sys.argv[1]
lines = open(settings_file).readlines()
param_start = [i for i, line in enumerate(lines) if line.startswith(">parameters")][0]
dir_start = [i for i, line in enumerate(lines) if line.startswith(">directories")][0]

# read settings file
settings = pd.read_csv(settings_file, header=None)

# read parameters
parameters = pd.read_csv(settings_file, 
                       nrows=dir_start-1).dropna(how='all', axis=0).dropna(how='all', axis=1).set_index('category')

# read directories
n_headers_dir = len(lines[dir_start].split(",")) # get number of headers in directory line
directories = pd.read_csv(settings_file, 
                          skiprows=dir_start, usecols=range(1,n_headers_dir)).dropna(how='all', axis=0)

## Check if all necessary parameters are in the file
param_cats = ["channel_names", "channel_labels", "n_channels"]
dir_cats = ["experiment_id", "unique_id", "input_transcriptome", 
    "input_images", "output", "vertices_x", "vertices_y"]


assert np.all([elem in parameters.index for elem in param_cats]), \
    "Not all required categories found in parameter section {}".format(param_cats)
assert np.all([elem in directories.columns for elem in dir_cats]), \
    "Not all required column headers found in directory section {}".format(dir_cats)

# determine extra categories which are added later to adata.obs
extra_cats_headers = [elem for elem in directories.columns if elem not in dir_cats]
extra_cats_headers = ["experiment_id"] + extra_cats_headers

# extract parameters
if not (pd.isnull(parameters.loc["channel_names", "value"]) or pd.isnull(parameters.loc["channel_labels", "value"])):
    channel_names = str(parameters.loc["channel_names", "value"]).split(" ")
    channel_labels = str(parameters.loc["channel_labels", "value"]).split(" ")

    # determine alignment channel
    alignment_channel = [elem for elem in channel_names if "*" in elem][0]

    # check if channel names and labels have same length
    assert len(channel_names) == len(channel_labels), \
        "Numbers of channel_names and channel_labels differ."

# get ids of vertices_x and vertices_y
vertx_id = directories.columns.tolist().index('vertices_x')
verty_id = directories.columns.tolist().index('vertices_y')

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

# check if all input images are exist
assert np.all([os.path.isfile(img) for img_d in directories['input_images'] for img in glob(img_d)]), \
    "Not all input image files exist."

# check if all output names are unique
assert len(np.unique(directories["output"])) == len(directories["output"]), \
    "Output files are not unique. This would cause that one file is overwritten by another."

n_datasets = len(directories)
vertices_list = [None]*n_datasets

# create path for new settings_file to include vertices
settings_new = settings_file.rstrip(".csv") + "_" + f"{datetime.now():%Y%m%d}" + "withvertices.csv"

print("Processing {} datasets".format(n_datasets))
for i in range(0, n_datasets):
    print("{} : Start processing dataset {} of {}".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, n_datasets), 
        flush=True)

    # extract directories for this dataset
    dirs = directories.loc[i, :]

    if len(extra_cats_headers) > 0:
        extra_cats = dirs[extra_cats_headers]
    else:
        extra_cats = None

    # get parameters for this dataset
    matrix_file = dirs["input_transcriptome"]
    output_file = dirs["output"]
    output_dir = os.path.dirname(output_file)
    unique_id = dirs["experiment_id"] + "_" + dirs["unique_id"]

    # check if there are images given for this dataset
    images_given = pd.notnull(dirs["input_images"])

    # create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if images_given:
        # get image directories
        image_dirs = glob(dirs["input_images"])

        # check if the vertices are given in the settings file
        vertices_not_given = pd.isnull(directories.loc[i, "vertices_x"]) or pd.isnull(directories.loc[i, "vertices_y"])

        # check if number of images matches number of channel names
        assert len(image_dirs) == len(channel_names), \
            "Number of detected images [{}] does not match number of channel names [{}] in parameters file.".format(
                len(image_dirs), len(channel_names))

        # detect get alignment marker image
        alignimg_dir = [d for d in image_dirs if alignment_channel.strip("*") in d][0]

        # read alignment image
        alignment_image = cv2.imread(alignimg_dir)

        if vertices_not_given:
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

            # add extracted coordinates to list of vertices
            vertices_list[i] = np.array([[a, b] for a, b in zip(ys, xs)])

        # read all images
        print("{} : Read all images...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
            flush=True)
        
        images = [cv2.imread(d, -1) for d in image_dirs]
    else:
        images = None
        channel_labels = None

    # generate squidpy formatted anndata object
    print("{} : Generate anndata object from matrix file and images...".format(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
        flush=True)

    db.dbitseq_to_squidpy(matrix_path=matrix_file, images=images,
        resolution=int(parameters.loc["spot_width"]), 
        unique_id=unique_id, extra_categories=extra_cats,
        n_channels=int(parameters.loc["n_channels"]), 
        frame=int(parameters.loc["frame"]),
        #dbitx=False, 
        labels=channel_labels, vertices=vertices_list[i], 
        savepath=output_file)

# free memory
del images
del alignment_image
gc.collect()
    
print("{} : Finished all datasets. Vertices saved into {}".format(
    f"{datetime.now():%Y-%m-%d %H:%M:%S}", settings_new), 
    flush=True)