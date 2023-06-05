#!/usr/bin/env python

#/lustre/groups/hpc/meier_lab/miniconda3/envs/spatial_gpu/bin/python

print("Starting script...", flush=True)

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
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from pathlib import Path


def cellpose(img, min_size=15):
    model = models.Cellpose(model_type='nuclei', gpu=gpu)
    res, _, _, _ = model.eval(
        img,
        channels=[0, 0],
        diameter=None,
        min_size=min_size,
    )
    return res

def stardist_2D_versatile_fluo(img, nms_thresh=None, prob_thresh=None):
    # Make sure to normalize the input image beforehand or supply a normalizer to the prediction function.
    # this is the default normalizer noted in StarDist examples.
    img = normalize(img, 1, 99.8, axis=(0,1))
    model = StarDist2D.from_pretrained('2D_versatile_fluo')
    model.config.use_gpu = gpu
    labels, _ = model.predict_instances(img, nms_thresh=nms_thresh, prob_thresh=prob_thresh)
    return labels

# Check for GPU 
gpu = False
if torch.cuda.is_available():
    gpu = True

print(f'GPU available: {gpu}', flush=True)

## Load hq files
input_file = "/lustre/groups/hpc/meier_lab/projects/SpatialMouse/analysis/adata_filtered_unnormalized_hqimages.h5ad"

## Set save dir
savedir = "/lustre/groups/hpc/meier_lab/projects/SpatialMouse/analysis/segmentation/"

print("Input summary:\n\tInput file: {}\n\tSave directory: {}".format(input_file, savedir), flush=True)

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
    print("\t{}: Start segmentation with cellpose...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    sq.im.segment(
        img=image, 
        layer="image", 
        channel=0,
        layer_added='segmented_cellpose',
        method=cellpose
        )
    print("\t{}: Segmentation finished.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    
    # Extract masks
    masks = np.where(image['segmented_cellpose']==0, 0, image['segmented_cellpose']).squeeze()
    
    # save masks
    cp_savedir = os.path.join(savedir, "cellpose")
    Path(cp_savedir).mkdir(parents=True, exist_ok=True)
    savefile = os.path.join(cp_savedir, "segmasks_cellpose_{}.npy".format(img_key))
    
    # save
    np.save(savefile, masks)
    print("\t{}: Cellpose masks saved into {}".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", savefile), flush=True)
    
    # free memory
    print("\tFree memory...", flush=True)
    del masks
    gc.collect()
    

    print("\t{}: Start segmentation with StarDist...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    sq.im.segment(
        img=image,
        layer="image",
        channel=0,
        method=stardist_2D_versatile_fluo,
        layer_added='segmented_stardist',
        nms_thresh=None,
        prob_thresh=None
    )
    print("\t{}: Segmentation finished.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    
    # extract masks
    masks = np.where(image['segmented_stardist']==0, 0, image['segmented_stardist']).squeeze()

    # Save masks
    sd_savedir = os.path.join(savedir, "stardist")
    Path(sd_savedir).mkdir(parents=True, exist_ok=True)
    savefile = os.path.join(sd_savedir, "segmasks_stardist_{}.npy".format(img_key))

    # save
    np.save(savefile, masks)
    print("\t{}: StarDist masks saved into {}".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", savefile), flush=True)
    
    # free memory
    print("\tFree memory...", flush=True)
    del masks
    gc.collect()

print("{}: All images processed.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)