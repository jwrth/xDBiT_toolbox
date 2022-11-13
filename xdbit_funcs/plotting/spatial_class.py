from typing import Union, List, Dict, Any, Optional
from anndata import AnnData
import pandas as pd
import numpy as np
from ..tools import check_raw

class MultiSpatialPlot:
    '''
    https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html
    '''
    def __init__(self, 
                 adata: AnnData, 
                 keys: List[str], 
                 groupby: str = 'id',
                 raw: bool = False,
                 layer: Optional[str] = None,
                 max_cols: int = 4,
                 pd_dataframe: Optional[pd.DataFrame] = None,
                 header_names: Optional[List[str]] = None,
                 normalize_crange_not_for: List = [],
                 crange: Optional[List[int, int]] = None, 
                 crange_type: str = 'minmax'
                 ):
        self.adata = adata
        self.keys = keys
        self.groupby = groupby
        self.raw = raw
        self.layer = layer
        self.max_cols = max_cols
        self.pd_dataframe = pd_dataframe
        self.header_names = header_names
        self.normalize_crange_not_for = normalize_crange_not_for
        self.crange = crange
        self.crange_type = crange_type
        
        # check arguments
        self.check_arguments()
        
    def check_arguments(self):
        # check keys and groups
        self.keys = [self.keys] if isinstance(self.keys, str) else list(self.keys)
        self.multikeys = False
        self.multigroups = False
        if len(self.keys) > 1:
            self.multikeys = True
            
        if self.groupby is None:
            self.groups = [None]
        else:
            if self.groups is None:
                self.groups = list(self.adata.obs[self.groupby].unique())
            else:
                self.groups = [self.groups] if isinstance(self.groups, str) else list(self.groups)
            if len(self.groups) > 1:
                self.multigroups = True
                
        if self.header_names is not None:
            assert len(self.header_names) == len(self.keys)
                
        # check if dataframe is given
        if isinstance(self.pd_dataframe, pd.DataFrame):
            self.data_in_dataframe=True
        else:
            self.data_in_dataframe=False

        # determine the color range for each key
        self.crange_per_key_dict = {
            key: self.get_crange(
                self.adata, self.groupby, self.groups, key, 
                use_raw=self.raw, layer=self.layer, 
                data_in_dataframe=self.data_in_dataframe, 
                pd_dataframe=self.pd_dataframe, ctype=self.crange_type
                ) if key not in self.normalize_crange_not_for else None for key in self.keys
            }
        
    
