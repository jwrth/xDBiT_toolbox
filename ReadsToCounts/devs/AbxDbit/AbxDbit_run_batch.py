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
import sys
from datetime import datetime
import numpy as np

print("Starting AbxDbit pipeline batch processing...", flush=True)
## Start timing
t_start = datetime.now()

print("Reading batch parameters from {}".format(sys.argv[1]))
settings = pd.read_csv(sys.argv[1])

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
        log_dir = join(dirname(s["rnafq_R1"]), "pipeline_{}.out".format(f"{datetime.now():%Y_%m_%d}"))
        out_dir = join(dirname(s["rnafq_R1"]))
        tmp_dir = join(dirname(s["rnafq_R1"]))

        print("Writing log to {}".format(log_dir), flush=True)
        print("Output directory: {}".format(out_dir), flush=True)
        print("Temp directory: {}".format(tmp_dir), flush=True)
        
        # get input files
        filecats = ["rnafq_R1", "rnafq_R2", "featfq_R1", "featfq_R2"]
        file_list = [s[cat] for cat in filecats if not pd.isnull(s[cat])]

        # generate commands
        commands.append(["bash", s["pipeline_dir"], 
        "-g", s["STAR_genome"],
        "-r", s["STAR_fasta"],
        "-b", s["legend"],
        "-f", s["feature_legend"],
        "-n", str(s["expected_n_spots"]),
        "-m", s["mode"],
        "-i", s["feature_mode"],
        "-j", "1",
        "-o", out_dir,
        "-t", tmp_dir,
        #"-l", # clear tmp files after finishing the script
        "-e", # for testing
        ] + file_list
        )

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
    
print("{}: All batches finished after {}.".format(f"{datetime.now():%Y-%m-%d %H:%M:%S}", str(t_elapsed)), flush=True)
