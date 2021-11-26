#!/usr/bin/env python

from tkinter import *
from tkinter import filedialog
import sys
import os

class SelectionWindow:

    # Initialization
    def __init__(self):

        self.root = Tk()
        self.root.title('SpatialViewer')
        self.root.geometry("250x170")
        self.root.option_add("*font", "Calibri 10")

        self.home = os.path.expanduser("~") # get home directory

        # set entry window for file
        self.entry = Entry(self.root)
        self.entry.grid(row=0,column=1)

        # set file select button
        self.selectFile = Button(self.root, text="Select File", command=self.browsefunc)
        self.selectFile.grid(row=0,column=0)

        # set windows for scale factor selection
        self.labelText = StringVar()
        self.labelText.set("Scale factor:")
        self.labelDir = Label(self.root, textvariable=self.labelText)
        self.labelDir.grid(row=1, column=0)

        self.sf_default = StringVar(self.root, value="1")
        self.scale_entry = Entry(self.root, textvariable=self.sf_default, width=4, justify=CENTER)
        self.scale_entry.grid(row=1, column=1)

        # set windows for groupby selection
        self.labelText = StringVar()
        self.labelText.set("Group by:")
        self.labelDir = Label(self.root, textvariable=self.labelText)
        self.labelDir.grid(row=2, column=0)

        self.sf_default = StringVar(self.root, value="well_name")
        self.groupby = Entry(self.root, textvariable=self.sf_default, width=12, justify=CENTER)
        self.groupby.grid(row=2, column=1)

        # set windows for group selection
        self.labelText = StringVar()
        self.labelText.set("Group:")
        self.labelDir = Label(self.root, textvariable=self.labelText)
        self.labelDir.grid(row=3, column=0)

        self.sf_default = StringVar(self.root, value="B3")
        self.group = Entry(self.root, textvariable=self.sf_default, width=12, justify=CENTER)
        self.group.grid(row=3, column=1)

        # set windows for resolution key selection
        self.labelText = StringVar()
        self.labelText.set("Resolution key:")
        self.labelDir = Label(self.root, textvariable=self.labelText)
        self.labelDir.grid(row=4, column=0)

        self.sf_default = StringVar(self.root, value="hires")
        self.resolution_key = Entry(self.root, textvariable=self.sf_default, width=12, justify=CENTER)
        self.resolution_key.grid(row=4, column=1)

        # set cancel button to exit script
        self.cancel = Button(self.root, text="Cancel", command=sys.exit)
        self.cancel.grid(row=5,column=1)


        # set okay button to continue script
        self.okay = Button(self.root, text="Okay", command=self.closewindow)
        self.okay.grid(row=5,column=0)

        # key bindings
        self.root.bind("<Return>", func=self.closewindow)
        self.root.bind("<Escape>", func=sys.exit)

        # change grid parameters
        self.root.grid_rowconfigure(0, minsize=30)
        self.root.grid_rowconfigure(1, minsize=30)
        self.root.grid_rowconfigure(2, minsize=30)

    # Functions
    def closewindow(self, event=None): # event only needed for .bind which passes event object to function
        self.file = self.entry.get()
        self.scale_factor = self.scale_entry.get()
        self.groupby = self.groupby.get()
        self.group = self.group.get()
        self.resolution_key = self.resolution_key.get()
        self.root.destroy()

    def browsefunc(self):
        print("Open selection window...")
        self.filename = filedialog.askopenfilename(
            initialdir=os.path.join(self.home, "Documents\\"),
            title="Choose .adata file",
            filetypes=(("h5ad files", "*.h5ad"), ("all files", "*.*"))
        )
        self.root.lift()
        self.root.attributes("-topmost", True)

        self.entry.delete(0, END)
        self.entry.insert(END, self.filename)

if __name__ == "__main__":
    selection = SelectionWindow()
    selection.root.mainloop()

    input_file = selection.file
    scale_factor = float(selection.scale_factor)
    groupby = str(selection.groupby)
    group = str(selection.group)
    resolution_key = str(selection.resolution_key)

    # import packages to read data
    print("Import packages to read data...")
    from scanpy import read_h5ad
    from squidpy.im import ImageContainer
    from skimage.transform import resize
    from napari import gui_qt
    from view_functions import single_grayscale_to_rgb, extract_groups

    # read data
    print("Read data...")
    adata = read_h5ad(input_file)

    # extract one group
    adata = extract_groups(adata, groupby=groupby, groups=group, extract_uns=True)

    keys = list(adata.uns['spatial'].keys())

    # resizing
    if scale_factor < 1:
        print("Scaling...")
        # scale coordinates
        adata.obsm['spatial'] *= scale_factor
        for key in keys:
            image_data = adata.uns['spatial'][key]['images'][resolution_key]
            nrows = image_data.shape[0]
            ncols = image_data.shape[1]
            adata.uns['spatial'][key]['images'][resolution_key] = resize(image_data, output_shape=(int(nrows * scale_factor), int(ncols * scale_factor)))

            adata.uns['spatial'][key]['scalefactors']['spot_diameter_fullres'] *= scale_factor

    # # convert to RGB
    # print("Convert to RGB...")
    # for key in keys:
    #     adata.uns['spatial'][key]['images'][resolution_key] = single_grayscale_to_rgb(adata.uns['spatial'][key]['images'][resolution_key])

    # extract image data
    img = ImageContainer()

    print("Fetch image data...")

    for key in keys:
        image_to_add = adata.uns['spatial'][key]['images'][resolution_key]
        img.add_img(image_to_add, layer=key)

    # display data
    print('Display data...')
    with gui_qt():
        img.interactive(adata, library_id=keys[0], symbol='square')

    print("Finished.")
