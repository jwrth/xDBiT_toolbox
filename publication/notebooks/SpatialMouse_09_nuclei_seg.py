#!/usr/bin/env python

#/lustre/groups/hpc/meier_lab/miniconda3/envs/spatial_gpu/bin/python

print("Start script...", flush=True)

import os
from cellpose import models
import torch
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
from glob import glob
import squidpy as sq
import gc
from datetime import datetime


def cellpose(img, min_size=15):
    model = models.Cellpose(model_type='nuclei', gpu=gpu)
    res, _, _, _ = model.eval(
        img,
        channels=[0, 0],
        diameter=None,
        min_size=min_size,
    )
    return res

# Check for GPU 
gpu = False
if torch.cuda.is_available():
    gpu = True

print(f'GPU available: {gpu}', flush=True)

## Load hq files
input_file = "/lustre/groups/hpc/meier_lab/projects/SpatialMouse/analysis/adata_filtered_unnormalized_hqimages.h5ad"

## Set save dir
savedir = "/lustre/groups/hpc/meier_lab/projects/SpatialMouse/analysis/segmentation/"

print("Input summary:\nInput file: {}\nSave directory: {}".format(input_file, savedir), flush=True)

print("{}: Load input file {}".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", input_file), flush=True)
adata = sc.read(input_file)

# get indices and set channel to process
indices = list(adata.obs['id'].unique())
channel = "dapi"

for i, idx in enumerate(indices):
    print("{}: Processing {}... ({}/{})".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", idx, i+1, len(indices)), flush=True)
    # generate img_key
    img_key = idx + "_" + channel

    # extract image and create container
    hq_img = adata.uns["spatial"][img_key]["images"]["hires"].copy()
    image = sq.im.ImageContainer(hq_img, lazy=True, chunks=200)

    ## Apply cellpose
    print("\t{}: Start segmentation".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    sq.im.segment(img=image, layer="image", channel=0, method=cellpose)
    print("\t{}: Segmentation finished.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    
    # Extract masks
    masks = np.where(image['segmented_custom']==0, 0, image['segmented_custom']).squeeze()

    ## Save masks
    savefile = os.path.join(savedir, "segmasks_{}.npy".format(img_key))
    np.save(savefile, masks)
    print("\t{}: Masks saved into {}".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", savefile), flush=True)

print("{}: All images processed.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)