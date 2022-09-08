#!/usr/bin/env python3

## Import the custom library
print("Import library...", flush=True)
import os
import sys

# add xDbit toolbox path to path variable
module_path = os.path.abspath("../../")
if module_path not in sys.path:
    sys.path.append(module_path)

import xdbit_funcs as db

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path
import seaborn as sns
from mofapy2.run.entry_point import entry_point
from datetime import datetime

## Load data per organ
### Set directories
data_dir = "out/"
model_dir = os.path.join(data_dir, "models")

# create model path if necessary
Path(model_dir).mkdir(parents=True, exist_ok=True)

### Read spatial data
input_file_pattern = "pp_wohires_deg"
input_files = [elem for elem in glob(os.path.join(data_dir, "*.h5ad")) if input_file_pattern in elem]

print("Read data from {}".format(input_files), flush=True)
adatas = {os.path.basename(file).split("_")[0]: sc.read(file) for file in input_files}

### Add spatial coordinates to `.obs`
for organ, adata in adatas.items():
    adata.obs = pd.concat([adata.obs, 
                           pd.DataFrame(adata.obsm["spatial"], columns=["imagerow",
                                                                        "imagecol"],
                                        index=adata.obs_names),
                          ], axis=1)

## Train a MEFISTO model for each slice of each organ
def run_mofa(adata, output_file="ST_model.hdf5", 
             n_factors=10, features_subset="highly_variable", 
             obsm_key="spatial",
             seed=2021, frac_inducing=0.5):
    '''
    Run mefisto/mofa on adata.
    '''
    
    # Set up mofa options
    ent = entry_point()
    ent.set_data_options(use_float32=True)
    ent.set_data_from_anndata(adata, features_subset=features_subset)
    ent.set_model_options(factors=n_factors)
    ent.set_train_options()
    ent.set_train_options(seed=seed)
    
    ent.set_covariates([adata.obsm[obsm_key]], covariates_names=["imagerow", "imagecol"])
    ent.set_smooth_options(sparseGP=True, frac_inducing=frac_inducing,
                           start_opt=10, opt_freq=10)
    
    # build, run and save
    ent.build()
    ent.run()
    ent.save(output_file) #hdf5 file

# iterate through organs
t_start = datetime.now()
print("{}: Start training...".format(f"{t_start:%Y-%m-%d %H:%M:%S}"), flush=True)
for organ, adata in adatas.items():
    # iterate through indices
    indices = adata.obs['id'].unique()
    for idx in indices:
        # prepare run
        a = db.tl.extract_groups(adata, groupby="id", groups=idx, extract_uns=True,
                                 strip=True)
        output_file = os.path.join(model_dir, "mofa_model_{}_{}.hdf5").format(organ, idx)
        
        # run mofa
        print("Processing {}: {}".format(organ, idx))
        run_mofa(a, output_file=output_file)
        
# finishing statement
t_end = datetime.now()
t_elapsed = t_end - t_start 
print("{}: Script finished after {}.".format(f"{t_end:%Y-%m-%d %H:%M:%S}", str(t_elapsed)),
      flush=True)
