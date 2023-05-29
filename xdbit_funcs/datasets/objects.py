from typing import Any, Dict, List, Optional, Tuple, Union, Type, Literal
import numpy as np
from anndata import AnnData
from ..tools import check_raw
from pathlib import Path
import dask.array as da
import pickle

class ImageData:
    '''
    Object extracting image data from anndata object.
    '''
    def __init__(
        self,
        adata: AnnData,
        image_key: Optional[str] = None,
        group: Optional[str] = None,
        lowres: bool = True,
        histogram_setting: Optional[Tuple[int, int]] = None,
        lowres_sf: float = 0.1,
        zarr_key: str = "zarr_directories"
        ):
        
        self.image_key = image_key
        self.group = group
        self.lowres = lowres
        self.histogram_setting = histogram_setting
        self.lowres_sf = lowres_sf
        self.zarr_key = zarr_key
                
        # get image data and or metadata
        if self.image_key is not None:
            self.check_source(adata)
            if self.load_from_zarr:
                self.load_image_from_zarr()
                self.process_image()
            else:
                self.load_image_from_uns(adata)
                self.process_image()
        else:
            self.no_image(adata)
    
    def check_source(self, adata):
        '''
        Check source from which the image should be loaded.
        '''
        self.load_from_zarr = False
        if self.zarr_key in adata.uns.keys():
            zarr_dirs = adata.uns[self.zarr_key][self.group]
            zarr_dirs = list(set(zarr_dirs)) # make list unique
            
            # select zarr_dirs that are a directory
            zarr_dirs = [Path(d) for d in zarr_dirs if Path(d).is_dir()]
            
            # unique_dir = False
            # dir_found = False
            # if len(zarr_dirs) == 1:
            #     # a unique file directory was found and will be used
            #     unique_dir = True
            #     dir_found = True
            # elif len(zarr_dirs) > 1:
            #     dir_found = True
            # else:
            #     dir_found = False
            
            image_found = False
            found_dirs = []
            for zd in zarr_dirs:
                channel_dirs = [elem for elem in list(zd.glob("*")) if elem.is_dir()]
                
                # check if image key is one of the channels
                chdir = [d for d in channel_dirs if self.image_key == d.stem]
                
                if len(chdir) == 1:
                    found_dirs.append(chdir[0])
                    image_found = True
            
            if image_found:
                if len(found_dirs) == 1:
                    self.zarr_path = found_dirs[0]
                    self.load_from_zarr = True
                elif len(found_dirs) > 1:
                    print("Warning: More than one possible zarr path found: {}. The first one was used for plotting.".format(found_dirs), flush=True)
                    self.zarr_path = found_dirs[0]
                    self.load_from_zarr = True
            else:
                self.load_from_zarr = False
                raise ValueError("No image path found for image key {} in zarr directories {}.".format(self.image_key, zarr_dirs))
                
        else:
            # check if image_key gives unique match
            image_key_list = adata.uns['spatial'].keys()
            image_key_matches = [k for k in image_key_list if self.image_key in k]
            if self.group is not None:
                # make sure that group name is also in image_key
                image_key_matches = [k for k in image_key_matches if self.group in k]

            if len(image_key_matches) == 0:
                raise ValueError("No match for img_key [{}] found in list of image_keys: {}".format(
                    self.image_key, image_key_list))

            elif len(image_key_matches) > 1:
                raise ValueError("More than one possible match for img_key [{}] found: {}".format(
                    self.image_key, image_key_matches
                    ))
                
            else:
                self.image_key = image_key_matches[0]
            
    def load_image_from_uns(self, adata):
        self.image_metadata = adata.uns['spatial'][self.image_key]['scalefactors']
        if "pixel_per_um" in self.image_metadata:
            self.pixel_per_um = self.image_metadata["pixel_per_um"]
        else:
            self.pixel_per_um = self.image_metadata["pixel_per_um_real"]
        
        if self.lowres:
            if 'lowres' in adata.uns['spatial'][self.image_key]['images'].keys():
                self.image = adata.uns['spatial'][self.image_key]['images']['lowres']
                self.scale_factor = self.image_metadata['tissue_lowres_scalef']
            else:
                print('`hires` image displayed because no `lowres` images was found.')
                self.image = adata.uns['spatial'][self.image_key]['images']['hires']
                self.scale_factor = self.image_metadata['tissue_hires_scalef']
        else:
            self.image = adata.uns['spatial'][self.image_key]['images']['hires']
            self.scale_factor = self.image_metadata['tissue_hires_scalef']
            
    def load_image_from_zarr(self):
        from ..images import resize_image
        # load image
        self.zarr_path = Path(self.zarr_path)
        self.image = da.from_zarr(self.zarr_path)
        self.image = self.image.compute() # convert to numpy array and load in memory
        
        # load and extract metadata
        meta_file = self.zarr_path.parent / "metadata.pkl"
        with open(meta_file, "rb") as f:
            self.image_metadata = pickle.load(f)
            
        if "pixel_per_um" in self.image_metadata:
            self.pixel_per_um = self.image_metadata["pixel_per_um"]
        else:
            self.pixel_per_um = self.image_metadata["pixel_per_um_real"]
            
        if self.lowres:
            self.image = resize_image(self.image, scale_factor=self.lowres_sf)
            self.scale_factor = self.lowres_sf
        else:
            self.scale_factor = self.image_metadata['tissue_hires_scalef']

    def process_image(self):
        from ..images import set_histogram
        self.bit_type = np.uint8 if self.image.max() < 256 else np.uint16
        if self.histogram_setting is None:
            # do min max scaling
            self.image = set_histogram(self.image, 
                                        lower=self.image.min(), upper=np.percentile(self.image, 99), 
                                        bit_type=self.bit_type)
        else:
            if self.histogram_setting == 'minmax':
                self.image = set_histogram(self.image, 
                                        lower=self.image.min(), upper=self.image.max(), 
                                        bit_type=self.bit_type)
            elif isinstance(self.histogram_setting, tuple):
                self.image = set_histogram(self.image, 
                                        lower=self.histogram_setting[0], upper=self.histogram_setting[1], 
                                        bit_type=self.bit_type)
            elif (self.histogram_setting > 0) & (self.histogram_setting < 1):
                self.image = set_histogram(self.image, lower=self.image.min(), upper=int(self.image.max() * self.histogram_setting), 
                                        bit_type=self.bit_type)

            else:
                raise ValueError('Unknown format of `histogram_setting`. Must be either "minmax" or (minval, maxval) or value between 0 and 1')
        
    def no_image(self, adata):
        self.image = None
        self.scale_factor = 1
        self.pixel_per_um = 1
        self.image_metadata = None
        # search for pixel_per_um scalefactor
        if 'spatial' in adata.uns.keys():
            image_key_list = adata.uns['spatial'].keys()
            if self.group is not None:
                # make sure that group name is also in image_key and select first option
                self.image_key = [k for k in image_key_list if self.group in k][0]
                first_entry = adata.uns['spatial'][self.image_key]
            else:
                # get first entry of dictionary as image_metadata
                first_entry = next(iter(adata.uns['spatial'].values()))
            
            # extract image metadata if possible
            if 'scalefactors' in first_entry:
                self.image_metadata = first_entry['scalefactors']
                if "pixel_per_um" in self.image_metadata:
                    self.pixel_per_um = self.image_metadata["pixel_per_um"]
                else:
                    self.pixel_per_um = self.image_metadata["pixel_per_um_real"]
                    
                self.scale_factor = self.image_metadata['tissue_hires_scalef']
            else:
                print("pixel_per_um scalefactor not found. Plotted pixel coordinates instead.")
        else:
            print("No key `spatial` in adata.uns. Therefore pixel_per_um scalefactor could not be found. Plotted pixel coordinates instead.")

class ExpressionData:
    '''
    Object extracting spatial coordinates and expression data from anndata object.
    '''
    def __init__(
        self, 
        adata: AnnData,
        key: List[str],
        ImageDataObject: Type[ImageData],
        raw: bool = False,
        layer: Optional[str] = None,
        obsm_key: str = 'spatial',
        origin_zero: bool = True, # whether to start axes ticks at 0
        plot_pixel: bool = False,
        xlim: Optional[Tuple[int, int]] = None,
        ylim: Optional[Tuple[int, int]] = None,
        spot_size_unit: float = 50, # whether to leave margin of one spot width around the plot
        margin: bool = True
        ):
        
        # add arguments to object
        self.key = key
        self.raw = raw
        self.layer = layer
        self.xlim = xlim
        self.ylim = ylim
                
        ## Extract coordinates    
        # extract parameters from ImageDataObject
        pixel_per_um = ImageDataObject.pixel_per_um
        scale_factor = ImageDataObject.scale_factor
        
        # extract x and y pixel coordinates and convert to micrometer
        self.x_pixelcoord = adata.obsm[obsm_key][:, 0].copy()
        self.y_pixelcoord = adata.obsm[obsm_key][:, 1].copy()
        self.x_coord = self.x_pixelcoord / pixel_per_um
        self.y_coord = self.y_pixelcoord / pixel_per_um

        # shift coordinates that they start at (0,0)
        if origin_zero:
            self.x_offset = self.x_coord.min()
            self.y_offset = self.y_coord.min()
            self.x_coord -= self.x_offset
            self.y_coord -= self.y_offset
        else:
            self.x_offset = self.y_offset = 0
        
        if self.xlim is None:
            if plot_pixel:
                xmin = self.x_pixelcoord.min() * scale_factor
                xmax = self.x_pixelcoord.max() * scale_factor
            else:
                xmin = np.min([self.x_coord.min(), self.y_coord.min()]) # make sure that result is always a square
                xmax = np.max([self.x_coord.max(), self.y_coord.max()])

            self.xlim = (xmin - spot_size_unit, xmax + spot_size_unit)
        elif margin:
            self.xlim[0] -= spot_size_unit
            self.xlim[1] += spot_size_unit

        if self.ylim is None:
            if plot_pixel:
                ymin = self.y_pixelcoord.min() * scale_factor
                ymax = self.y_pixelcoord.max() * scale_factor
            else:
                # ymin = y_coord.min()
                # ymax = y_coord.max() 
                ymin = np.min([self.x_coord.min(), self.y_coord.min()])
                ymax = np.max([self.x_coord.max(), self.y_coord.max()])

            self.ylim = (ymin - spot_size_unit, ymax + spot_size_unit)
        elif margin:
            self.ylim[0] -= spot_size_unit
            self.ylim[1] += spot_size_unit
            
        ## Extract expression data
        # check if plotting raw data
        adata_X, adata_var, adata_var_names = check_raw(adata, 
                                                        use_raw=self.raw, 
                                                        layer=self.layer)
        
        self.data = adata.obs.copy()
        self.data['x_coord'] = self.x_coord
        self.data['y_coord'] = self.y_coord
        
        # locate gene in matrix and extract values
        if key in adata_var_names:
            idx = adata_var.index.get_loc(key)
            self.color = adata_X[:, idx].copy()
            self.categorical = False
            self.exists = True
            
        elif key in adata.obs.columns:
            self.exists = True
            if adata.obs[key].dtype.name == 'category':
                self.categorical = True
            else:
                self.color = adata.obs[key].values
                self.categorical = False
        else:
            self.exists = False