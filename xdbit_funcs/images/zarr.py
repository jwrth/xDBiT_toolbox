import dask.array as da
from numcodecs import Blosc
from pathlib import Path
import numpy as np

def save_to_zarr(img, url, chunk_size=(1000, 1000), compressor='custom', cname='zstd', clevel=1, shuffle=Blosc.BITSHUFFLE, overwrite=False):
    '''
    Function to save images in array format to zarr format.
    
    Parameters
    ----------
    img: array
    url: path to save the data
    compressor: Default value is 'custom'. If value 'default' is chosen, the function 
                uses the default compressor of `to_zarr()`.
    cname: A string naming one of the compression algorithms available within blosc, e.g., ‘zstd’, ‘blosclz’, ‘lz4’, ‘lz4hc’, ‘zlib’ or ‘snappy’.
    clevel: An integer between 0 and 9 specifying the compression level. With 1 being best performance and low quality.
    shuffle: Either NOSHUFFLE (0), SHUFFLE (1), BITSHUFFLE (2) or AUTOSHUFFLE (-1). If AUTOSHUFFLE, bit-shuffle will be used 
            for buffers with itemsize 1, and byte-shuffle will be used otherwise. The default is SHUFFLE.
    
    '''
    # Check type of input img
    if isinstance(img, np.ndarray):
        # create dask array
        img = da.from_array(img, chunks=chunk_size)
    elif not isinstance(img, da.Array):
        raise TypeError("Unknown file type. Dask array or numpy array are allowed.")
    
    # make sure that parent directory of url exists
    url = Path(url)
    url.parent.mkdir(parents=True, exist_ok=True)
    
    # compression settings   
    if compressor == 'custom':
        compressor = Blosc(cname=cname, clevel=clevel, shuffle=shuffle)
    elif compressor == 'default':
        compressor = 'default'
    
    # save array
    img.to_zarr(url=url, compressor=compressor, overwrite=overwrite)
    
def load_zarr_dir(directory):
    '''
    Function to load multiple zarr files from directory.
    '''
    directory = Path(directory)
    
    urls = directory.glob("*/")
    urls = [elem for elem in urls if elem.is_dir()]
    
    results = {}
    for url in urls:
        ch = url.stem
        results[ch] = da.from_zarr(url)
        
    return results
