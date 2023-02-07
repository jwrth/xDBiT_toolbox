  #!/usr/bin/env python

'''
Script to run AbxDbit pipeline on multiple batches.

Usage:
    Input is a .csv file giving the parameters:
        - directory of fastqs (2 for RNA only, 4 for RNA+features)
        - STAR genome directory
        - STAR genome fasta file path
        - directory of barcode legend .csv
        - (Optional) feature legend file
        - number of expected spots
        - mode

Command for single analysis:
nohup bash path/to/xDbit_toolbox/AbxDbit/AbxDbit_pipeline.sh \
-g path/to/genomes_STAR/mm_STARmeta/STAR/ -r path/to/genomes_STAR/mm_STARmeta/Mm_metadata.fasta \
-b path/to/barcodes/barcodes_legend.csv -f path/to/barcodes/feature_legend.csv \
-n 2000 -m Dbit-seq -j 1 \
path/to/RNA_R1.fastq path/to/RNA_R2.fastq path/to/features_R1.fastq path/to/features_R2.fastq &

'''

from subprocess import Popen, PIPE, STDOUT
from os.path import realpath
import pandas as pd
from datetime import datetime
import numpy as np
import argparse
from pathlib import Path
import warnings
import os

print("Starting subsampling pipeline...", flush=True)
## Start timing
t_start = datetime.now()

# get path of current script
current_script_path = Path(realpath(__file__))
current_script_dir = current_script_path.parent
subsample_script = current_script_dir / "subsample_bam.sh"

# Parse arguments
print("Parse arguments")
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--batch_file", help="Batch file with all parameters for xDbit batch processing.")
parser.add_argument("-p", "--plotting_only", action='store_true', help="Run only the plotting part.")
args = parser.parse_args()

batch_file = Path(args.batch_file)
plotting_only = args.plotting_only
print("Reading batch parameters from {}".format(batch_file))

if batch_file.suffix == ".csv":
    settings = pd.read_csv(batch_file)
elif batch_file.suffix == ".xlsx":
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning) # using a drop-down list always throws an UserWarning which is ignored here
        settings = pd.read_excel(batch_file)
else:
    raise TypeError("Wrong type of settings file. Allowed: `.csv` or `.xlsx`")

if not plotting_only:
    # get total number of files
    nfiles = len(settings)

    # get number of batches
    batch_numbers = settings['batch'].unique()
    print("Found following batch numbers: {}".format(batch_numbers))

    for bi, b in enumerate(batch_numbers):
        batch = settings.query('batch == {}'.format(b))  
        print("{} : Start processing batch {} of {} with {} files".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", 
                                                                        bi, 
                                                                        len(batch_numbers), 
                                                                        len(batch)), flush=True)
        
        ## Start timing of batch
        t_start_batch = datetime.now()

        ## Generate Commands
        commands = []
        log_dirs = []
        for i in range(len(batch)):
            # extract settings
            s = batch.iloc[i, :]
            
            # get parameters
            bamfile = Path(s["file"])
            
            # create outfile paths
            out_dir = bamfile.parent / "subsampling"
            out_dir.mkdir(parents=True, exist_ok=True)
            log_dir = out_dir / "subsampling_log.out" 
            
            print("Output directory: {}".format(out_dir), flush=True)    
            print("Log file: {}".format(log_dir), flush=True)   
            
            # generate commands
            commands.append(["bash", 
                            str(subsample_script), 
                            "-o", # overwrite existing out files
                            "-d", # create DGE matrix
                            str(bamfile)])
            
            # collect log files
            log_dirs.append(log_dir)
        
        for e, elem in enumerate(commands):
            print("Command {}:".format(e+1), flush=True)
            print(" ".join(elem), flush=True)

        procs = [Popen(commands[i], stdout=PIPE, stderr=STDOUT) for i in range(len(commands))]

        running = [True]*len(procs)
        
        # open log files
        log_files = [open(d, "a") for d in log_dirs]
        
        while True:
            for i, p in enumerate(procs):
                output = p.stdout.readline()
                if running[i] and p.poll() is not None:
                    running[i] = False
                    print("{}: Process {} in batch {} finished".format(
                        f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, b), flush=True)
                if output:
                    print(output.decode("utf-8").strip(), file=log_files[i], flush=True)
            if not np.any(running):
                ## Stop timing
                t_stop_batch = datetime.now()
                t_elapsed_batch = t_stop_batch-t_start_batch
                print("{}: All processes of batch {} finished after {}.".format(
                    f"{datetime.now():%Y-%m-%d %H:%M:%S}", b, str(t_elapsed_batch)), flush=True)
                break
            
        # close log files
        for f in log_files:
            f.close()

    ## Stop timing
    t_stop = datetime.now()
    t_elapsed = t_stop-t_start
        
    print("{}: All files finished after {}.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", str(t_elapsed)), flush=True)

## Plotting
print("Plotting of results", flush=True)
for idx, row in settings.iterrows():
    
    # create outfile paths
    bamfile = Path(row["file"])
    out_dir = bamfile.parent / "subsampling"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # get sample name
    name = row["name"]
    
    os.system("{}/plot_subsampled_data.py -f {} -s {} -n {}".format(
        current_script_dir,
        bamfile.parent, # directory of full matrix
        out_dir, # directory of subsampled matrices            
        name
    ))
    
print("Finished everything.", flush=True)