# CountsToAnndata

This workflow shows how to create a anndata file from the raw count matrix which was created in `ReadsToCounts`.

It involves the following steps:
  1. Alignment of transcriptome spots using alignment image.
  2. Addition of experimental parameters to dataset.
  3. (Optional) Registration of hiqh-quality images with alignment image using SIFT algorithm.

## Installation of environment

### Windows

`conda env create -f CountsToAnndata_Win.yml`
