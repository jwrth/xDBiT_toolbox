#!/usr/bin/env python

'''
Script to add images from high-resolution imaging round to dataset.

Usage:

    Input:
        - .csv file giving following parameters:
            - 

'''

import sys
import pandas as pd

## Read parameters
print("Starting Alignment script...", flush=True)

# read parameters file
print("Reading batch parameters from {}".format(sys.argv[1]))
settings_file = sys.argv[1]

# read settings file
settings = pd.read_csv(settings_file, header=None)

## Check if all necessary parameters are in the file
cats = ["input", "hires_images", "output"]

assert np.all([elem in settings.columns for elem in cats]), \
    "Not all required column headers found in settings file {}".format(cats)

