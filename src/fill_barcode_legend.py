#!/usr/bin/env python

'''
Script to fill the barcode legend for the DbitX toolbox.

Usage:
    python fill_barcode_legend.py <well_assignment_file.csv> <empty_barcode_legend.csv>
'''

import pandas as pd
import numpy as np
import sys


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

# Script start


# Fetch input files
assignment_file = sys.argv[1]
legend_file = sys.argv[2]

# Generate output file
output_file = "barcodes_legend.csv"

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
