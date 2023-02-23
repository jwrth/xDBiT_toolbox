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
from os.path import dirname, join
import pandas as pd
from datetime import datetime
import numpy as np
import argparse
from pathlib import Path
import warnings

print("Starting xDbit pipeline v2.0 batch processing...", flush=True)
## Start timing
t_start = datetime.now()

# Parse arguments
print("Parse arguments")
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--batch_file", help="Batch file with all parameters for xDbit batch processing.")
parser.add_argument("-s", "--skip_align", action='store_true', help="Skip alignment steps for testing.")
parser.add_argument("-c", "--clear_tmp", action='store_true', help="Clear files in temporary directory after pipeline is finished.")
parser.add_argument("-e", "--echo", action='store_true', help="Print all commands to output instead of executing them.")
parser.add_argument("--spatial_only", action='store_true', help="Perform analysis only on the spatial coordinates on R2. Can be used for spillover analysis.")
args = parser.parse_args()

batch_file = Path(args.batch_file)
print("Reading batch parameters from {}".format(batch_file))

if batch_file.suffix == ".csv":
    settings = pd.read_csv(batch_file)
elif batch_file.suffix == ".xlsx":
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning) # using a drop-down list always throws an UserWarning which is ignored here
        settings = pd.read_excel(batch_file)
else:
    raise TypeError("Wrong type of settings file. Allowed: `.csv` or `.xlsx`")

batch_numbers = settings['batch'].unique()
print("Found following batch numbers: {}".format(batch_numbers))

for b in batch_numbers:
    batch = settings.query('batch == {}'.format(b))

    print("{} : Start processing batch {} of {} with {} files".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", b, batch_numbers[-1], len(batch)), flush=True)
    ## Start timing of batch
    t_start_batch = datetime.now()

    ## Generate Commands
    commands = []
    log_dirs = []
    for i in range(len(batch)):
        # extract settings
        s = batch.iloc[i, :]

        # get directory
        if pd.notnull(s["rnafq_R1"]):
            sdir = dirname(s["rnafq_R1"])
        else:
            sdir = dirname(s["rnafq_R2"])
        log_dir = join(sdir, "pipeline_batch{}-{}_{}.out".format(b, i, f"{datetime.now():%Y_%m_%d}"))
        out_dir = join(sdir)
        tmp_dir = join(sdir)

        print("Writing log to {}".format(log_dir), flush=True)
        print("Output directory: {}".format(out_dir), flush=True)
        print("Temp directory: {}".format(tmp_dir), flush=True)
        
        # get input files
        filecats = ["rnafq_R1", "rnafq_R2", "featfq_R1", "featfq_R2"]
        file_list = [s[cat] for cat in filecats if not pd.isnull(s[cat])]
        
        # generate commands
        commands.append(["bash", s["pipeline_dir"], 
        "-b", s["legend"],
        "-n", str(s["expected_n_spots"]),
        "-m", s["mode"],
        "-j", "1",
        "-o", out_dir,
        "-t", tmp_dir,
        ])
        
        if args.spatial_only:
            commands[i].append("-y")
        else:
            # add also the informations needed for the STAR alignment
            commands[i] += ["-g", s["STAR_genome"], "-r", s["STAR_fasta"]]

        # if feature files are given, add them
        if not pd.isnull(s["feature_legend"]):
            commands[i] += ["-f", s["feature_legend"]]
        if not pd.isnull(s["feature_mode"]):
            commands[i] += ["-u", s["feature_mode"]]
        
        # add other flags
        if args.skip_align:
            commands[i].append("-x")
        if args.clear_tmp:
            commands[i].append("-l")
        if args.echo:
            commands[i].append("-e")
            
        # add file list
        commands[i] += file_list

        log_dirs.append(log_dir)
    
    for e, elem in enumerate(commands):
        print("Command {}:".format(e+1), flush=True)
        #print(elem)
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
                print("{}: Process {} in batch {} finished".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, b), flush=True)
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

# in the words of Jim Morrison
print("{}: All batches finished after {}.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", str(t_elapsed)), flush=True)
print("This is the end, beautiful friend", flush=True)
print("This is the end, my only friend", flush=True)
print("The end.", flush=True)
