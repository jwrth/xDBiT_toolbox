#!/usr/bin/env python

'''
Script to add images from high-resolution imaging round to dataset.

Usage:

    Input:
        - .csv file giving following parameters:
            - 

'''

import sys
import os
from tkinter import *
from tkinter import filedialog, messagebox
import tkinter
import numpy as np
import pandas as pd
from glob import glob
from datetime import datetime
import argparse

class SelectionWindow:

    # Initialization
    def __init__(self):

        self.root = Tk()
        self.root.title('CountsToAnndata')
        self.root.geometry("400x200")
        self.root.option_add("*font", "Calibri 10")

        self.home = os.path.expanduser("~") # get home directory

        # set entry window for file
        self.entry = Entry(self.root)
        self.entry.insert(END, "D:/data/37_53/CountsToAnndata/37_53_CoA_param_withvertices_20221115.csv")
        self.entry.grid(row=0,column=1)

        # set file select button
        self.selectFile = Button(self.root, text="Select settings file", command=self.browsefunc)
        self.selectFile.grid(row=0,column=0)

        # set scale factor selection
        self.labelText = StringVar()
        self.labelText.set("Maximum pixel width for registration:")
        self.labelDir = Label(self.root, textvariable=self.labelText)
        self.labelDir.grid(row=1, column=0)

        self.maxpx_default = StringVar(self.root, value=4000)
        self.maxpx_entry = Entry(self.root, textvariable=self.maxpx_default, width=4, justify=CENTER)
        self.maxpx_entry.grid(row=1, column=1)

        # set HQ image resolution selection selection
        self.labelText2 = StringVar()
        self.labelText2.set("Resolution HQ image [px/µm]:")
        self.labelDir2 = Label(self.root, textvariable=self.labelText2)
        self.labelDir2.grid(row=2, column=0)

        self.ppmhq_default = IntVar(self.root, value=4.0480)
        self.ppmhq_entry = Entry(self.root, textvariable=self.ppmhq_default, width=6, justify=CENTER)
        self.ppmhq_entry.grid(row=2, column=1)

        # set alignment image resolution selection selection
        self.labelText3 = StringVar()
        self.labelText3.set("Resolution alignment image [px/µm]:")
        self.labelDir3 = Label(self.root, textvariable=self.labelText3)
        self.labelDir3.grid(row=3, column=0)

        self.ppmalign_default = IntVar(self.root, value=0.6069)
        self.ppmalign_entry = Entry(self.root, textvariable=self.ppmalign_default, width=6, justify=CENTER)
        self.ppmalign_entry.grid(row=3, column=1)

        # set 8bit conversion checkbutton
        self.labelText4 = StringVar()
        self.labelText4.set("Convert images to 8-bit:")
        self.labelDir4 = Label(self.root, textvariable=self.labelText4)
        self.labelDir4.grid(row=4, column=0)

        self.convert_entry = IntVar(value=True)
        self.checkbutton = Checkbutton(self.root, variable=self.convert_entry, onvalue=True, offvalue=False)
        self.checkbutton.grid(row=4, column=1)

        # set 8bit conversion checkbutton
        self.labelText5 = StringVar()
        self.labelText5.set("Perspective transformation (else affine):")
        self.labelDir5 = Label(self.root, textvariable=self.labelText5)
        self.labelDir5.grid(row=5, column=0)

        self.perspect_entry = IntVar(value=False)
        self.checkbutton2 = Checkbutton(self.root, variable=self.perspect_entry, onvalue=True, offvalue=False)
        self.checkbutton2.grid(row=5, column=1)

        # set cancel button to exit script
        self.cancel = Button(self.root, text="Cancel", command=sys.exit)
        self.cancel.grid(row=6,column=1)

        # set run button to continue script
        self.run = Button(self.root, text="Run", command=self.closewindow)
        self.run.grid(row=6,column=0)

        # key bindings
        self.root.bind("<Return>", func=self.closewindow)
        self.root.bind("<Escape>", func=sys.exit)

        # change grid parameters
        self.root.grid_rowconfigure(0, minsize=35)
        self.root.grid_rowconfigure(1, minsize=35)
        self.root.grid_columnconfigure(0, minsize=100)
        self.root.grid_columnconfigure(1, minsize=180)

        # state what happens when the window is closed
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    # Functions
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.root.destroy()
            sys.exit()
    
    def closewindow(self, event=None): # event only needed for .bind which passes event object to function
        # collect variables
        self.settings_file = self.entry.get()
        self.maxpx = self.maxpx_entry.get()
        self.convert = self.convert_entry.get()
        self.perspect = self.perspect_entry.get()
        self.ppmhq = self.ppmhq_entry.get()
        self.ppmalign = self.ppmalign_entry.get()
        self.root.destroy()

    def browsefunc(self):
        print("Open selection window...")
        self.filename = filedialog.askopenfilename(
            initialdir=os.path.join(self.home, "Documents\\"),
            title="Choose .adata file",
            filetypes=(("csv files", "*.csv"), ("all files", "*.*"))
        )
        self.root.lift()
        self.root.attributes("-topmost", True)

        self.entry.delete(0, END)
        self.entry.insert(END, self.filename)

class CountsToAnndata():
    def __init__(self):
        self.result = {}
        self.register_hq = []

    def ReadAndCheckSettings(self, settings_file):
        '''
        Checks if settings file is correct.
        '''

        print("Reading parameters from {}".format(settings_file))
        lines = open(settings_file).readlines()
        param_start = [i for i, line in enumerate(lines) if line.startswith(">parameters")][0]
        self.dir_start = [i for i, line in enumerate(lines) if line.startswith(">directories")][0]

        # read settings file
        self.settings = pd.read_csv(settings_file, header=None)

        # read parameters
        self.parameters = pd.read_csv(settings_file,
                            nrows=self.dir_start-1).dropna(how='all', axis=0).dropna(how='all', axis=1)

        # check headers of parameters file and set category as index
        param_headers = ["category", "value"]
        assert np.all([elem in self.parameters.columns for elem in param_headers]), \
            "Parameters section does not have correct headers ({})".format(param_headers)                       

        self.parameters = self.parameters.set_index('category')

        # read directories
        n_headers_dir = len(lines[self.dir_start].split(",")) # get number of headers in directory line
        self.directories = pd.read_csv(settings_file,
                                skiprows=self.dir_start, usecols=range(1,n_headers_dir)).dropna(how='all', axis=0)

        ## Check if all necessary parameters are in the file
        param_cats = ["n_channels", "spot_width", "frame", 
            "align_images:align_channel", "align_images:dapi_channel", "hq_images:channel_names", "hq_images:channel_labels"]
        dir_cats = ["experiment_id", "unique_id", "main_dir", "organism",
            "input_transcriptome", "align_images", "hq_images", "output_dir", 
            "vertices_x", "vertices_y"]

        # check for missing values in parameters and directories
        missing_params = [elem for elem in param_cats if elem not in self.parameters.index]
        missing_paramvals = [cat for cat in param_cats if pd.isnull(self.parameters.loc[cat, "value"])]
        missing_dirs = [elem for elem in dir_cats if elem not in self.directories.columns]

        halt=False
        missnames = ["parameters", "parameter values", "directories"]
        for i, miss in enumerate([missing_params, missing_paramvals, missing_dirs]):
            if len(miss) > 0:
                print("MissingValueError: Following values are missing in {}: {}".format(missnames[i], miss))
                halt=True
        if halt:
            sys.exit()

        # determine extra categories which are added later to adata.obs
        self.extra_cats_headers = [elem for elem in self.directories.columns if elem not in dir_cats]
        self.extra_cats_headers = ["experiment_id", "organism"] + self.extra_cats_headers

        ## extract image parameters
        # determine alignment channel
        image_cats = ["align_images:align_channel", "align_images:dapi_channel", 
            "hq_images:channel_names", "hq_images:channel_labels"]


        # create full directories from main_dir and the relative paths of matrix and output
        for cat in ["input_transcriptome", "output_dir"]:
            self.directories[cat] = [os.path.join(m, p.strip(os.sep)) if isinstance(p, str) else p 
                                        for m, p in 
                                        zip(self.directories["main_dir"], self.directories[cat])
                                        ]

        if not np.any([pd.isnull(self.parameters.loc[elem, "value"]) for elem in image_cats]):
            self.alignment_channel = self.parameters.loc["align_images:align_channel", "value"]
            self.dapi_channel = self.parameters.loc["align_images:dapi_channel", "value"]

            # get channel names of hq images
            self.hq_ch_names = str(self.parameters.loc["hq_images:channel_names", "value"]).split(" ")
            self.hq_ch_labels = str(self.parameters.loc["hq_images:channel_labels", "value"]).split(" ")

            # determine channel on which registration is done - most probably the dapi
            reg_id = [i for i, elem in enumerate(self.hq_ch_names) if "*" in elem]
            if len(reg_id) == 1:
                self.reg_channel_label = self.hq_ch_labels[reg_id[0]]
            elif len(reg_id) > 1:
                print("More than one `hq_images:channel_names` labelled with `*` found: {}".format(reg_id))
                sys.exit()
            else:
                print("No channel in hq_images:channel_names labelled with `*` to indicate registration channel.")
                self.reg_channel_label = None


            # check if channel names and labels have same length
            assert len(self.hq_ch_names) == len(self.hq_ch_labels), \
                "Number of channel_names ({}) and channel_labels ({}) differ.".format(len(self.hq_ch_names), len(self.hq_ch_labels))    

            # create full directories from main_dir and the relative paths of the images
            for cat in ["align_images", "hq_images"]:
                self.directories[cat] = [os.path.join(m, p.strip(os.sep)) if isinstance(p, str) else p 
                                            for m, p in 
                                            zip(self.directories["main_dir"], self.directories[cat])
                                            ]

        # check if all input matrix files exist
        try:
            assert np.all([os.path.isfile(f) for f in self.directories["input_transcriptome"]]), \
                "Not all input transcriptome files exist."
        except AssertionError as e:
            # check which files are missing
            missing_files = [f for f in self.directories["input_transcriptome"] if not os.path.isfile(f)]
            print("{} Following transcriptome files are missing: {}".format(e, missing_files))
            #sys.exit()
            exit()

        # check if unique_ids are really unique
        unique_ids = ["{}_{}".format(a,b) for a,b in zip(self.directories["experiment_id"], self.directories["unique_id"])]
        assert len(np.unique(unique_ids)) == len(unique_ids), \
            "`experiment_id` and `unique_id` together do not give a unique set of ids. \n" \
                "This would cause that one file is overwritten by another"

        # # assert that for each dataset there are either both hq_image and align_image or None of both
        # assert np.all([pd.notnull(a) == pd.notnull(b) for a, b in zip(self.directories["align_images"], self.directories["hq_images"])]), \
        #     "For some datasets only alignment images OR hq images were given.Either both are given or None.\n" \
        #         "If you wish to use the alignment images as hq images put the alignment directories also into\n" \
        #             "the hq image section. If you wish to create an anndata without image data leave both sections empty."

        ## add general parameters to class
        self.n_datasets = len(self.directories)
        self.vertices_list = [None]*self.n_datasets
        self.register_hq = [False]*self.n_datasets
        self.align_images_exist = [False]*self.n_datasets
        self.hq_images_exist = [False]*self.n_datasets
        self.process_images = [False]*self.n_datasets

        # get ids of vertices_x and vertices_y
        self.vertx_id = self.directories.columns.tolist().index('vertices_x')
        self.verty_id = self.directories.columns.tolist().index('vertices_y')


    def AddVertices(self, settings_file):
        '''
        Adds vertices to settings file.
        '''

        # create path for new settings_file to include vertices
        self.settings_new = settings_file.rstrip(".csv") + "_withvertices_" + f"{datetime.now():%Y%m%d}" + ".csv"

        print("Add vertices to {} datasets".format(self.n_datasets))
        #for i in range(0, n_datasets):
        for i, dirs in self.directories.iterrows():
            print("{} : Start processing dataset {} of {}".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, self.n_datasets), 
                flush=True)

            # check if alignment images are given
            images_given = pd.notnull(dirs["align_images"]) #and pd.notnull(dirs["hq_images"])

            # check if both alignment and hq_images exist
            if images_given:
                # check if all given aligment image paths exist
                self.align_images_exist[i] = np.all(os.path.isfile(img) for img in glob(dirs["align_images"]))
                self.process_images[i] = self.align_images_exist[i] #and self.hq_images_exist[i]
            else:
                self.process_images[i] = False

            if self.process_images[i]:
                # check if the vertices are given in the settings file
                vertices_not_given = pd.isnull(dirs["vertices_x"]) or pd.isnull(dirs["vertices_y"])

                # get image directories
                align_images = glob(dirs["align_images"])
                

                ## Decide whether to do hq image registration or not
                # check if hq images are given
                if pd.notnull(dirs["hq_images"]):
                    hq_images = glob(dirs["hq_images"])
                    # check if all given hq_images exist and registration should be done
                    if np.all(os.path.isfile(img) for img in hq_images):
                        print("HQ images are given and exist.")

                        # check if a "*" was given to indicate the registration channel
                        if self.reg_channel_label is not None:
                            print("Registration of HQ images and alignment images will be performed.")
                            self.register_hq[i] = True
                        else:
                            print("No channel in hq_images:channel_names labelled with `*` to indicate registration channel.")
                            print("NO REGISTRATION WILL BE PERFORMED EVEN THOUGH IMAGES ARE GIVEN!")
                    else:
                        print("No HQ images given. Pipeline is run without registration.")

                # # check if number of images matches number of channel names
                # assert len(hq_images) == len(self.hq_ch_names), \
                #     "Number of detected hq images [{}] does not match number of channel names [{}] in parameters file.".format(
                #         len(hq_images), len(self.hq_ch_names))

                # # check whether paths to alignment images and hq images are identical
                # if not len(set(align_images) & set(hq_images)) == len(align_images):
                #     if self.reg_channel_label is not None:
                #         print("Alignment images and hq are not identical. HQ images will be registered.")
                #         self.register_hq[i] = True
                #     else:
                #         print("Alignment images and hq are not identical and could be registered " \
                #             "but no channel in `hq_images:channel_names` was labelled with `*`." \
                #                 "Therefore, the hq images are not added.")
                # else:
                #     print("Alignment images and hq are identical. HQ images will not be registered.")
                
                # detect alignment marker image
                align_img = [elem for elem in align_images if self.alignment_channel in elem][0]
                #dapi_img = [elem for elem in align_images if self.dapi_channel in elem][0]

                if vertices_not_given:
                    import napari
                    # read alignment image
                    alignment_image = cv2.imread(align_img, -1)

                    ### Select corner spots in alignment image using napari viewer
                    # with napari.gui_qt():
                    # https://napari.org/guides/stable/event_loop.html#intro-to-event-loop
                    print("No vertices given. Select them from napari viewer.", flush=True)
                    points_fetched = False
                    while points_fetched is not True:
                        try:
                            viewer = napari.view_image(db.im.convert_to_8bit(alignment_image), # napari viewer had problems with 16-bit images 
                                title="Select corner spots in alignment image {} of {} ".format(i+1, self.n_datasets))
                            napari.run()

                            # fetch vertices (center points at cross points of alignment channels)
                            assert "Points" in viewer.layers, "No Points selected. Select exactly 4 points as vertices."
                            corner_spots_center = viewer.layers["Points"].data.astype(int)

                            assert len(corner_spots_center) == 4, "Number of selected points is not correct. Select exactly 4 points as vertices."

                            points_fetched = True
                        except AssertionError as k:
                            print(k, flush=True)
                            #print("No or not the right number of points (4) was selected. Try again.")
                            pass

                    # collect corner spot
                    self.vertices_list[i] = corner_spots_center

                    # save information about vertices in settings file
                    # save y coordinates (row coordinates)
                    self.settings.loc[self.dir_start+i+1, self.verty_id+1] = " ".join([str(elem[0]) for elem in self.vertices_list[i]])
                    # save x coordinates (column coordinates)
                    self.settings.loc[self.dir_start+i+1, self.vertx_id+1] = " ".join([str(elem[1]) for elem in self.vertices_list[i]])
                    # save settings file with coordinates of vertices
                    self.settings.to_csv(self.settings_new, index=None, header=None)
                    print("Vertices selected and saved for dataset {} of {}.".format(i+1, self.n_datasets), 
                        flush=True)

                else:
                    # extract coordinates from directory input
                    xs = [int(elem) for elem in self.directories.loc[i, "vertices_x"].split(" ")]
                    ys = [int(elem) for elem in self.directories.loc[i, "vertices_y"].split(" ")]

                    # make sure four vertices are given
                    assert len(xs) == 4 and len(ys) == 4, \
                        "Not enough x- [{}] and/or y-values [{}] given.".format(len(xs), len(ys)) 

                    # add extracted coordinates to list of vertices
                    self.vertices_list[i] = np.array([[a, b] for a, b in zip(ys, xs)])

        print("{} : Addition of vertices finished. New settings file saved into {}".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", self.settings_new), 
            flush=True)


    def ProcessDatasets(self, maxpx_before_reg, to_8bit, ppmhq, ppmalign, perspective_transform):
        '''
        Process images:
            - Align coordinates to alignment images
            - If hq images are given: Register alignment images and hq images
            - Save everything as h5ad image.
        '''

        print("Processing {} datasets".format(self.n_datasets))
        #for i in range(0, n_datasets):
        for i, dirs in self.directories.iterrows():
            print("{} : Start processing dataset {} of {}".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}", i+1, self.n_datasets), 
                flush=True)

            if len(self.extra_cats_headers) > 0:
                extra_cats = dirs[self.extra_cats_headers]
            else:
                extra_cats = None

            # get parameters for this dataset
            matrix_file = dirs["input_transcriptome"]
            output_dir = dirs["output_dir"]
            unique_id = dirs["experiment_id"] + "_" + dirs["unique_id"]
            organism = None if pd.isnull(dirs["organism"]) else dirs["organism"]

            # create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            #register_hq = False
            if self.process_images[i]:
                # create name of output files
                output_filename = "{}_hqimages.h5ad".format(unique_id)
                noimages_filename = "{}_noimages.h5ad".format(unique_id)
                align_filename = "{}_alignimages.h5ad".format(unique_id)

                # check if the vertices are given in the settings file
                #vertices_not_given = pd.isnull(dirs["vertices_x"]) or pd.isnull(dirs["vertices_y"])

                # # get image directories
                align_images = glob(dirs["align_images"])
                
                # # detect alignment marker image
                dapi_img = [i for i in align_images if self.dapi_channel in i][0]

                # read images
                print("{} : Read images...".format(
                    f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
                    flush=True)
                
                if self.register_hq[i]:
                    image_paths = [dapi_img]
                    channel_labels = [self.reg_channel_label]
                else:
                    image_paths = align_images # is equal to align_images in this case
                    channel_labels = self.hq_ch_labels
                
                images = [cv2.imread(d, -1) for d in image_paths]
                
            else:
                images = None
                channel_labels = None
                align_filename = "{}_noimages.h5ad".format(unique_id)

            ### Generation of squidpy formatted anndata object
            print("{} : Generate anndata object from matrix file and images...".format(
                f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
                flush=True)


            if self.register_hq[i]:
                # create tmp directory
                tmp_dir = os.path.join(output_dir, "tmp")
                noimages_dir = os.path.join(output_dir, "noimages")
                Path(tmp_dir).mkdir(parents=True, exist_ok=True)
                Path(noimages_dir).mkdir(parents=True, exist_ok=True)

                # generate path to output files for adata objects
                align_outfile = os.path.join(tmp_dir, align_filename)
                output_file = os.path.join(output_dir, output_filename)
                noimages_file = os.path.join(noimages_dir, noimages_filename)
                return_adata = True
            else:
                align_outfile = os.path.join(output_dir, align_filename)
                return_adata = False

            # align transcriptomic data to alignment images using the selected vertices

            adata = db.dbitseq_to_squidpy(
                matrix_path=matrix_file, 
                images=images, ppmalign=ppmalign,
                resolution=int(self.parameters.loc["spot_width", "value"]), 
                unique_id=unique_id, 
                organism=organism,
                extra_categories=extra_cats,
                n_channels=int(self.parameters.loc["n_channels", "value"]), 
                frame=int(self.parameters.loc["frame", "value"]),
                labels=channel_labels, 
                vertices=self.vertices_list[i], 
                to_8bit=to_8bit,
                savepath=align_outfile, 
                return_adata=return_adata
                )

            ### Registration of hiqh-quality imaging data
            if self.register_hq[i]:
                hq_images = glob(dirs["hq_images"])
                # create image_dir_dict
                image_dir_dict = {}
                for i, n in enumerate(self.hq_ch_names):
                    selected_image_path = [elem for elem in hq_images if n.strip("*") in os.path.basename(elem)][0]
                    image_dir_dict[unique_id + "-" + str(self.hq_ch_labels[i])] = selected_image_path

                # start transformation
                adata_trans = db.im.register_adata_coords_to_new_images(adata, 
                                                                groupby='id', 
                                                                image_dir_dict=image_dir_dict, ppmhq=ppmhq,
                                                                reg_channel=self.reg_channel_label, 
                                                                in_place=False, 
                                                                do_registration=True,
                                                                perspective_transform=perspective_transform,
                                                                maxpx_before_reg=maxpx_before_reg, 
                                                                to_8bit=to_8bit
                                                                )
                
                # save adata after transformation with images
                print("Save adata into {}".format(output_file))
                adata_trans.write(output_file)

                # plot registration QC plot
                self.registration_qc_plot(adata=adata, 
                                        adata_trans=adata_trans, 
                                        output_dir=output_dir, 
                                        unique_id=unique_id, 
                                        reg_channel_label=self.reg_channel_label
                                        )

                # plot matchedVis
                plotdir = os.path.join(output_dir, "registration_checks")
                Path(plotdir).mkdir(parents=True, exist_ok=True)
                
                for g, mv in adata_trans.uns['matchedVis'].items():
                    plotpath = os.path.join(plotdir, "{}_matchedVis.png".format(g))
                    plt.imshow(mv)
                    print("Saving matched visualization into {}".format(plotpath))
                    plt.savefig(plotpath, dpi=400)
                    plt.close()

                # remove images and save adata without images
                print("Remove images and save adata without images into {}".format(noimages_file))
                db.tl.remove_images(adata_trans, hires_only=False, other_keys='matchedVis')
                adata_trans.write(noimages_file)
                        

        print("{} : Finished all datasets. Output saved into {}".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}", output_dir), 
            flush=True)

    def registration_qc_plot(self, adata, adata_trans, output_dir, unique_id, reg_channel_label, save_only=True):
        '''
        Plot to check how well the registration worked.
        '''
        # create plot to check success of registration
        plotdir = os.path.join(output_dir, "registration_checks")
        Path(plotdir).mkdir(parents=True, exist_ok=True)
        plotpath = os.path.join(plotdir, unique_id + "_check_registration.png")

        # parameters
        idxs = adata.obs['id'].unique()
        c_names = ["Low resolution image", "High-resolution image"]
        gene = "Actb"

        fig, axs = plt.subplots(len(idxs), 2, figsize=(8*2, 6*len(idxs)))

        if len(idxs) > 1:
            axs.ravel()

        for r, idx in enumerate(idxs):
            for c, ad in enumerate([adata, adata_trans]):
                db.pl.spatial(ad, keys=gene, groupby='id', groups=idx, image_key=reg_channel_label, 
                            lowres=False,
                            xlim=(1800,2000), ylim=(1800,2000), alpha=0.5, 
                            axis=axs[r+c], fig=fig)
                if r == 0:
                    axs[r+c].set_title(c_names[c] + "\n" + gene, fontsize=12, fontweight='bold')
                    
        fig.tight_layout()
        plt.savefig(plotpath, dpi=200)

        if save_only:
            plt.close(fig)
        else:
            plt.show()

def boolean_input(message, default):
    '''
    Convert user input into a boolean value if it is one of the following inputs:
    - "True"/"False"
    - "TRUE"/"FALSE"
    - "true"/"false"
    - "T"/"F"
    - "t"/"f"
    '''

    while True:
        # prompt input window and convert input to lower case and strip spaces
        user_input = input(message).lower().strip()
        if user_input in ["true", "t"]:
            return True
        elif user_input in ["false", "f"]:
            return False
        elif user_input == "":
            return default
        else:
            print("Error: Input must be True or False.")
                
#######
## Protocol start
#######

if __name__ == "__main__":

    #####
    # Parse arguments
    #####

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--test_run", action='store_true', help="Perform test run with test files provided in xDbit toolbox repository.")
    args = parser.parse_args()

    test_run = args.test_run

    #####
    # Selection window
    #####

    # open selection window
    if not test_run:
        try:
            selection = SelectionWindow()
            selection.root.mainloop()

            # read input of selection window
            settings_file = selection.settings_file
            maxpx = int(selection.maxpx)
            convert = bool(selection.convert)
            perspective_transform = bool(selection.perspect)
            ppmhq = float(selection.ppmhq) if len(selection.ppmhq) > 0 else None
            ppmalign = float(selection.ppmalign) if len(selection.ppmalign) > 0 else None

        except tkinter.TclError:
            # input settings file
            settings_file = input("Enter path to parameters .csv file: ")
            
            # limit image size during registration to reduce required RAM
            maxpx = int(input("Enter maximum pixel width for registration: "))

            # get input about resolutions
            ppmhq = input("Resolution HQ image [px/µm]:")
            ppmhq = float(ppmhq) if len(ppmhq) > 0 else None
            ppmalign = input("Resolution alignment image [px/µm]:")
            ppmalign = float(ppmalign) if len(ppmalign) > 0 else None

            # get input about transformations
            convert = boolean_input("Convert images to 8-bit? [t]/f", default=True)
            perspective_transform = boolean_input("Perspective transformation (else affine)? [t]/f", default=True)
    else:
        settings_file = "./test_files/test_params.csv"
        maxpx = 4000
        ppmhq = 80.96
        ppmalign = 40.462
        convert = True
        perspective_transform = False
    
    # read and check input file
    settings_file = settings_file.strip('"')
    coa = CountsToAnndata()
    coa.ReadAndCheckSettings(settings_file=settings_file)

    #####
    # Process datasets
    #####

    ## Load modules
    # get path of dbitx module and import functions
    print("Load modules...", flush=True)
    script_dir = os.path.dirname(os.path.realpath(__file__)) # read script location
    module_path = os.path.abspath(os.path.join(script_dir, ".."))
    if module_path not in sys.path:
        sys.path.append(module_path)

    import xdbit_funcs as db

    import matplotlib.pyplot as plt
    from pathlib import Path
    import cv2

    print("Starting processing...", flush=True)
    print("{} : Starting script...".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
            flush=True)

    
    coa.AddVertices(settings_file=settings_file)
    coa.ProcessDatasets(
        maxpx_before_reg=maxpx, to_8bit=convert,
        ppmhq=ppmhq, ppmalign=ppmalign, perspective_transform=perspective_transform)

    if test_run:
        import scanpy as sc
        adata = sc.read("./test_files/anndata/Exp_ID_A1_hqimages.h5ad")
        results_image = "./test_files/test_run_result.png"
        db.pl.spatial(adata, keys=['Actb'], groupby="id", image_key='dapi', savepath=results_image)
        print("Save results example of test run into {}".format(results_image))

    print("{} : Script finished.".format(
            f"{datetime.now():%Y-%m-%d %H:%M:%S}"), 
            flush=True)
