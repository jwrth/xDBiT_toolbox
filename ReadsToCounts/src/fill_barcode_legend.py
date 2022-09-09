#!/usr/bin/env python

'''
Script to fill the barcode legend for the DbitX toolbox.

Usage:
    # old 
    python fill_barcode_legend.py <well_assignment_file.csv> <empty_barcode_legend.csv>

    # new
    python fill_barcode_legend.py <empty_barcode_legend.csv> <well_assignment_files.csv or folder with assignment_files> 
'''

import pandas as pd
import numpy as np
import sys
import os
from glob import glob


# Functions
def create_pos_to_coord_dict(assignment_dataframe, coord_name):
    d = {pos: int(coord) for (pos, coord) in zip(
        assignment_dataframe[coord_name + "_WellPosition"], assignment_dataframe[coord_name + "_Coord"]) if not pd.isnull(pos)}
    return d

def fill_legend(legend, assignment, name):

    pos_to_coord = create_pos_to_coord_dict(assignment, name)

    new_column = []

    for pos in legend["WellPosition"]:
        if pos in pos_to_coord:
            entry = pos_to_coord[pos]
            new_column.append(entry)
        else:
            new_column.append(pd.NA)

    legend[name] = new_column

def main(legend_file, assignment_file, output_file):
    # Import files
    legend = pd.read_csv(legend_file)
    assign = pd.read_csv(assignment_file)

    # extract coordinate names from legend
    coord_names = [n for n in ['X', 'Y', 'Z'] if n in legend.columns]

    # fill and save legend file
    for name in coord_names:
        fill_legend(legend, assign, name)

    legend.to_csv(output_file)
    print(legend_file + " filled and saved in " + output_file, flush=True)

# Script start


# Fetch input files
legend_file = sys.argv[1]
assignment_files = sys.argv[2:]

# set other parameters
prefix_to_remove = "well_coord_assignment_"
output_dir = os.path.join(os.path.dirname(legend_file), "barcode_legend")

# create output directory
os.makedirs(output_dir, exist_ok=True)

# check if input consists of multiple files or only one file or directory
multi_file = True if len(assignment_files) > 1 else False

if not multi_file:
    assignment_file = assignment_files[0]
    if os.path.isfile(assignment_file):
        print("Single-file mode", flush=True)
        output_file = "barcode_legend.csv"
        main(legend_file, assignment_file, output_file)

    elif os.path.isdir(assignment_file):
        # extract files from directory
        assignment_files = glob(assignment_file + "/*")
        assignment_files = [elem for elem in assignment_files if elem.endswith(".csv")]

        # switch to multifile mode
        multi_file = True

    else:
        print("Input file is neither directory nor file.", flush=True)

if multi_file:
    # iterate through files
    for assignment_file in assignment_files:
        well_name = assignment_file.split("/")[-1].replace(prefix_to_remove, "").replace(".csv", "")
        print(well_name)
        output_file = os.path.join(output_dir, "barcode_legend_{}.csv".format(well_name))

        main(legend_file, assignment_file, output_file)

