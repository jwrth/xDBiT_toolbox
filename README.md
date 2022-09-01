# xDbit_toolbox

This repository contains all code that is necessary to reproduce results shown in the publication ...

Further, it includes an analysis pipeline to compute the digital gene expression matrix from Dbit-seq or xDbit experiments.
This part is based on the Split-seq toolbox: https://github.com/RebekkaWegmann/splitseq_toolbox.

If you have questions or suggestions, please feel free to open an issue or contact me directly.

# Publication

## Title: Spatial Transcriptomics Using Multiplexed Deterministic Barcoding in Tissue

## Abstract

In this study, we present a multiplexed version of deterministic barcoding in tissue (xDbit) to acquire spatially resolved transcriptomes of nine tissue sections in parallel. New microfluidic chips were developed to spatially encode mRNAs over a total tissue area of 1.17 cm2 with spots of 50 µm×50 µm. Optimization of the biochemical protocol increased read and gene counts per spot by one order of magnitude compared with previous reports. Furthermore, the introduction of alignment markers allows seamless registration of images and spatial transcriptomic spot coordinates. Together with technological advances, we provide an open-source computational pipeline to transform raw sequencing data from xDbit experiments into the AnnData format. The functionality of xDbit was demonstrated by the acquisition of 18 spatially resolved transcriptomic datasets from five different murine organs, including cerebellum, liver, kidney, spleen, and heart. Factor analysis and deconvolution of xDbit spatial transcriptomes allowed for in-depth characterization of the murine kidney.

## Analysis

All notebooks of the analyses performed in the publication can be found under `/publication/notebooks/`

# Introduction
## Barcoding layout

xDbit allows the spatial barcoding (X, Y) of 9 tissue sections (Z) on one object slide.

![layout](graphics/xdbit_layout.png)

## Read structure

![readstructure](graphics/xdbit_read-structure.png)

## xDbit toolbox preprocessing pipeline

![](graphics/pipeline_overview.png)

The preprocessing pipeline of the xDbit toolbox consists of two steps to convert raw sequencing reads into a spot/gene count matrix with aligned images:

1. **ReadsToCounts**
2. **CountsToAnndata**

Both pipelines and more detailed instructions can be found in the folders `./ReadsToCounts/` and `./CountsToAnndata/`, respectively.


# Get started

## Requirements

The pipeline uses a bash script, custom Python scripts, and many tools from the Drop-seq toolbox (Mc Caroll lab, Harvard Medical school) as well as Picard (Broad institute), which are all included in this toolbox. It was created in a Linux server environment and the STAR alignment step requires more than 30 GB RAM.
## Installation

### Clone repository

```
git clone https://github.com/jwrth/xDbit_toolbox.git

# make drop seq toolbox executable
cd /path/to/repo/xDbit_toolbox
chmod u=rwx,g=r,o=r ./ReadsToCounts/external_tools/Drop-seq_tools-2.1.0/*
```

### Install python environment

```
# install environment from file depending on OS
conda env create -f environment_linux.yml python=3
# or
conda env create -f environment_win.yml python=3

# activate environment
conda activate xdbit

# to access parts of the pipeline in a Jupyter notebook install a kernel for this environment
conda install -c anaconda ipykernel
python3 -m ipykernel install --user --name=xdbit_kernel
```

## Use toolbox as python module

To load the `xDbit_toolbox` as module run following code.

```
## Import the custom library
import os
import sys

# add xDbit toolbox path to path variable
module_path = os.path.abspath("../../")
if module_path not in sys.path:
    sys.path.append(module_path)

import dbitx_funcs as db
```