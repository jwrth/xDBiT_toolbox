# SpatialViewer

## Introduction
The SpatialExchange folder structure is supposed to help exchanging spatial data between cooperating partners.

File structure:
    - .h5ad file: contains transcriptomic and imaging data as well as all processed data.
    - .install_packages.bat: Windows batch file with all commands to install the necessary Python packages in an environment named "viewer_env"
    - run_viewer.bat: Windows batch file with commands to run the .python script.
    - view.py: Python script that loads and displays the data in a napari window.

Requirements:
    - OS: Windows 10
    - Anaconda installation (https://docs.anaconda.com/anaconda/install/windows/). The viewer expects anaconda to be installed in its default directory (%userprofile%\anaconda3\)
    - MacOS only usable from Anaconda Prompt

## Usage

For Windows users:

Start Viewer:
1. Install all packages by double-clicking on install_packages.bat.
2. Run the Viewer by double-clicking on run_viewer.bat.
3. When directory window pops up select .h5ad file.
4. Type in scale factor (value between 0 and 1). This scales down the images and makes it easier to switch between images in napari.
5. Wait until the napari window opens.

Navigate in napari window:
- On the right side of the window you can chose:
    1. Images
    2. Genes
    3. Observations: clustering results, total counts and other parameters that were collected during the image analysis.
- On the left side you can interact with the layers that were added to the image.
- The default window shows only highly differentially expressed genes. To see all genes in the list chose 'raw.


For more information refer to the official napari and squidpy tutorials:
- https://napari.org/tutorials/fundamentals/viewer.html
- https://squidpy.readthedocs.io/en/latest/external_tutorials/tutorial_napari.html

For Mac users:

Needs to be started from Anaconda Prompt

1. Open Anaconda Prompt
2. Navigate to location of viewer and inside the viewer folder.
3. Install packages by following command:
    ```
    conda env create -f viewer_env.yml python=3.7
   ```
4. After installation of all packages run ``view.py`` using following command:
    ````
   python view.py
    ````
3. When directory window pops up select .h5ad file.
4. Type in scale factor (value between 0 and 1). This scales down the images and makes it easier to switch between images in napari.
5. Wait until the napari window opens.
6. For details about napari usage see Windows section.


## Random commands needed during development

Command to create environment file:
conda env export --no-builds | findstr -v "prefix" > viewer_env.yml
