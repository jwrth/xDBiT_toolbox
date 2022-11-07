#!/usr/bin/env python3

## Import the custom library
import os
import sys
from datetime import datetime
from pathlib import Path
from tqdm import tqdm


if __name__ == '__main__':

    while True:
        input_file = input("Enter path to preprocessed adata .h5ad file: ")
        if os.path.exists(input_file):
            break
    
    groupby = input("Enter groupby variable (e.g. id): ")

    # generate results file
    results_file = input_file.rsplit(".", 1)[0] + "_spatialde.h5ad"

    # generate tmp directory
    dir = os.path.dirname(input_file)
    tmp_dir = os.path.join(dir, "tmp")
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    print("Load modules...")
    script_dir = os.path.dirname(os.path.realpath(__file__)) # read script location
    module_path = os.path.abspath(os.path.join(script_dir, ".."))
    if module_path not in sys.path:
        sys.path.append(module_path)

    import xdbit_funcs as db
    import scanpy as sc

    # start processing
    adata = sc.read(input_file)

    while True:
        if groupby in adata.obs.columns:
            groups = adata.obs[groupby].unique()
            break
        else:
            groupby = input("Given `groupby` variable not in `adata.obs`. Enter again: ")

    key_name = "spatialde"

    # create empty dictionary to store results
    adata.uns[key_name] = {}
    for group in tqdm(groups):
        print("{}: Run SpatialDE for group {}.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", group), flush=True)
        subset = adata[adata.obs[groupby] == group].copy()

        # remove images in subset to decrease memory consumption
        print("{}: Remove images from subset...".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
        db.tl.remove_images(adata=subset)
        
        # run spatialde on subset
        db.tl.spatialde_run(subset, run=True, normalize=False, use_raw=False, output_name=key_name)
        
        # Delete counts from SpatialDE file for saving because the numbers of columns is too high for saving in .h5ad
        del subset.uns[key_name]["counts"]

        # save results for this well in temporary directory
        print("{}: Save results of well {} in temporary directory.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", group), flush=True)
        subset.write(os.path.join(tmp_dir, "tmp_{}.h5ad".format(group)))

        # store result in adata
        adata.uns[key_name][group] = subset.uns[key_name]

    # save file with results from all wells
    print("{}: Save results of analysis to results file.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)
    adata.write(results_file)

print("{}: Finished SpatialDE analysis.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}"), flush=True)