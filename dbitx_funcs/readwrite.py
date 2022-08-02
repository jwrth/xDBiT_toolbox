from numpy.lib.arraysetops import unique
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
from .calculations._calc import coord_to_um, coord_to_pixel
from .images.image_processing import convert_to_8bit, resize_images_in_adata, calc_image_param_per_spot
from .images.registration import align_to_dict
from datetime import datetime
from gprofiler import GProfiler
import cv2
from pathlib import Path


def dbitseq_to_squidpy(matrix_path, resolution, n_channels, images=None, labels=None, vertices=None, convert_genes=True,
                        ppmalign=None, 
                        grayscale=True, to_8bit=False,
                        frame=100, unique_id=None, extra_categories=None, lowres_factor=0.2, 
                        organism='mmusculus',
                        spatial_key="spatial", img_keys=None, transpose=True, sep="x", manual_pixel_offset_x=0, 
                        manual_pixel_offset_y=0, savepath=None, return_adata=False, **kwargs):
    """
    Function to create adata object for squidpy from Dbit-seq data.
    Inputs are the path to the transcriptome matrix and the image. 

    Formats:
    - Transcriptome matrix: 2D count matrix with comma-separated spot coordinates in the rows ('X,Y') and gene names in the columns.
    - Image: RGB image as numpy array.

    """
    #print("Create adata object in '" + ["Dbit-seq", "DbitX"][dbitx] + "' mode")
    print("Read transcriptome matrix from {}...".format(matrix_path))
    if transpose:
        adata = sc.read_text(matrix_path).transpose()
    else:
        adata = sc.read_text(matrix_path)

    if unique_id is None:
        # take current time as unique identifier for this experiment
        unique_id = f"{datetime.now():%Y-%m-%d_%H:%M:%S}"

    # check if input has three or two coordinates
    # dbitx = False
    # if len(adata.obs_names[0].split(sep)) == 3:
    #     dbitx = True

    # add coordinates to adata object
    adata.obs['array_row'] = np.array(
        [int(elem.split(sep)[1]) for elem in adata.obs_names])
    adata.obs['array_col'] = np.array(
        [int(elem.split(sep)[0]) for elem in adata.obs_names])

    print("Compute coordinates...")
    # compute µm coordinates
    adata.obs['um_row'] = np.array(
        [coord_to_um(c, resolution) for c in adata.obs['array_row']])
    adata.obs['um_col'] = np.array(
        [coord_to_um(c, resolution) for c in adata.obs['array_col']])

    #if dbitx:
    adata.obs_names = np.array(["{e}-{uid}".format(e=elem, uid=unique_id) for elem in adata.obs_names])
    adata.obs['id'] = unique_id
    # else:
    #     #adata.obs_names = np.array(["{e}{s}{uid}".format(e=elem, s=sep, uid=unique_id) for elem in adata.obs_names])
    #     adata.obs_names = np.array(["{e}-{uid}".format(e=elem, uid=unique_id) for elem in adata.obs_names])
    #     adata.obs['id'] = unique_id

    if images is not None:
        assert labels is not None, "No labels given."
        assert vertices is not None, "No vertices given."

        if to_8bit:
            # convert images to 8bit images
            images = [convert_to_8bit(i) for i in images]

        # make labels unique
        unique_labels = [unique_id + "-" + lab for lab in labels]


        # read images and create metadata
        print("Align and create image metadata...")
        image_and_metadata = align_to_dict(
            images=images, labels=unique_labels, vertices=vertices, ppm_given=ppmalign,
            resolution=resolution, n_channels=n_channels, frame=frame, **kwargs
        )

        # extract parameters from metadata
        image_keys = list(image_and_metadata.keys())
        image_metadata = image_and_metadata[image_keys[0]]["scalefactors"]

        scale = image_metadata["pixel_per_um"]
        resolution = image_metadata["resolution"]
        offset_row = image_metadata["upper_left_spot_coord"][0]
        offset_col = image_metadata["upper_left_spot_coord"][1]

        # compute pixel coordinates
        adata.obs['pixel_row'] = np.array(
            [coord_to_pixel(c, resolution, scale, offset_row) +
            manual_pixel_offset_y for c in adata.obs['array_row']]
        )
        adata.obs['pixel_col'] = np.array(
            [coord_to_pixel(c, resolution, scale, offset_col) +
            manual_pixel_offset_x for c in adata.obs['array_col']]
        )

        # Add pixel coordinates to adata
        adata.obsm["spatial"] = adata.obs[["pixel_col", "pixel_row"]].values

        # Add image and metadata to adata
        adata.uns[spatial_key] = image_and_metadata

        # Resize images
        print("Resize images by factor {} and save as lowres...".format(lowres_factor), flush=True)
        resize_images_in_adata(adata, scale_factor=lowres_factor) # add lowres image

        # calculate mean intensity per spot for each channel
        print("Calculate channel mean per spot...", flush=True)
        for l in labels:
            calc_image_param_per_spot(adata, groupby='id', channel_pattern=l, 
                fun=np.mean, fun_descriptor='mean', lowres=True)

    else:
        # Add pixel coordinates to adata
        adata.obsm["spatial"] = adata.obs[["um_col", "um_row"]].values

    if extra_categories is not None:
        # add extra categories to .obs
        print("Add extra categories...", flush=True)
        for cat in extra_categories.index:
            adata.obs[cat] = extra_categories[cat]

    if convert_genes and organism is not None:
        print("Convert gene names...")
        genes = adata.var_names.tolist()

        # query gprofiler
        gp = GProfiler(return_dataframe=True)
        conversion = gp.convert(query=genes, organism=organism)
        conversion = conversion.replace("nan", np.nan)

        # sometimes it finds more than one match. Here I select always the first match.
        first_match = conversion.sort_values("n_incoming").groupby("n_incoming").head(1)

        # add genes names
        adata.var["ensembl_id"] = first_match["converted"].values
        adata.var_names_make_unique()

    if grayscale:
        # convert to grayscale
        print("Convert adata images to grayscale...")
        for key in adata.uns[spatial_key]:
            for res_key in adata.uns[spatial_key][key]['images'].keys():
                img = adata.uns[spatial_key][key]['images'][res_key]
                if len(img.shape) == 3:
                    adata.uns[spatial_key][key]['images'][res_key] = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    print("Adata object generated.")

    if savepath is None:
        print("Adata object returned.")
        print("Finished")
        return adata
    else:
        print("Saving adata object...")
        adata.write(savepath)
        print("Adata object saved into " + savepath)
        print("Finished")
        if return_adata:
            return adata
        else:
            return

def secure_save(adata, savepath, overwrite=False):
    '''
    Function to check whether file exists and whether it should be overwritten before saving.
    '''

    if os.path.isfile(savepath):
        print("Output file exists.")
        if overwrite:
            adata.write(savepath)
            print("Existing file overwritten and file saved to {}".format(savepath))
        else:
            print("File not saved.")
    else:
        adata.write(savepath)
        print("File saved to {}".format(savepath))

def save_and_show_figure(savepath, fig, save_only=False, show=True, dpi_save=300, save_background=None):
    #if fig is not None and axis is not None:
    #    return fig, axis
    #elif savepath is not None:
    if savepath is not None:
        print("Saving figure to file " + savepath)

        # create path if it does not exist
        Path(os.path.dirname(savepath)).mkdir(parents=True, exist_ok=True)

        # save figure
        fig.tight_layout()
        plt.savefig(savepath, dpi=dpi_save,
                    facecolor=save_background, bbox_inches='tight')
        print("Saved.")
    if save_only:
        plt.close(fig)
    elif show:
        fig.tight_layout()
        return plt.show()
    else:
        return