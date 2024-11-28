import anndata

'''
Utils from scib package: https://github.com/theislab/scib/tree/main/scib
'''


# checker functions for data sanity
def check_adata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError('Input is not a valid AnnData object')


def check_batch(batch, obs, verbose=False):
    if batch not in obs:
        raise ValueError(f'column {batch} is not in obs')
    elif verbose:
        print(f'Object contains {obs[batch].nunique()} batches.')


def check_hvg(hvg, hvg_key, adata_var):
    # if type(hvg) is not list:
    #     raise TypeError('HVG list is not a list')
    # else:
    #     if not all(i in adata_var.index for i in hvg):
    #         raise ValueError('Not all HVGs are in the adata object')
    if not hvg_key in adata_var:
        raise KeyError('`hvg_key` not found in `adata.var`')



def check_sanity(adata, batch, hvg, hvg_key):
    check_adata(adata)
    check_batch(batch, adata.obs)
    if hvg:
        check_hvg(hvg, hvg_key, adata.var)


def split_batches(adata, batch, hvg=None, return_categories=False):
    split = []
    batch_categories = adata.obs[batch].unique()
    if hvg is not None:
        adata = adata[:, hvg]
    for i in batch_categories:
        split.append(adata[adata.obs[batch] == i].copy())
    if return_categories:
        return split, batch_categories
    return split


def merge_adata(adata_list, sep='-'):
    """
    merge adatas from list and remove duplicated obs and var columns
    """

    if len(adata_list) == 1:
        return adata_list[0]

    adata = adata_list[0].concatenate(*adata_list[1:], index_unique=None, batch_key='tmp')
    del adata.obs['tmp']

    if len(adata.obs.columns) > 0:
        # if there is a column with separator
        if sum(adata.obs.columns.str.contains(sep)) > 0:
            columns_to_keep = [name.split(sep)[1] == '0' for name in adata.var.columns.values]
            clean_var = adata.var.loc[:, columns_to_keep]
        else:
            clean_var = adata.var

    if len(adata.var.columns) > 0:
        if sum(adata.var.columns.str.contains(sep)) > 0:
            adata.var = clean_var.rename(columns={name: name.split('-')[0] for name in clean_var.columns.values})

    return adata


def todense(adata):
    import scipy
    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.todense()

class SpeciesToID:
    def __init__(self):
        self.species_dict = {
            'mmusculus': 10090,
            'hsapiens': 9606,
            'dmelanogaster': 7227
        }
    def check_species(self, species):
        if species not in self.species_dict:
            raise ValueError(
                "`species` must be one of following values: {}".format(list(self.species_dict.keys()))
            )
    def convert(self, species):
        self.check_species(species)
        return self.species_dict[species]

def find_between(s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
        
def rename_net(net, source='source', target='target', weight='weight'):
    """
    Renames input network to match decoupler's format (source, target, weight).

    Parameters
    ----------
    net : DataFrame
        Network in long format.
    source : str
        Column name where to extract source features.
    target : str
        Column name where to extract target features.
    weight : str, None
        Column name where to extract features' weights. If no weights are available, set to None.

    Returns
    -------
    net : DataFrame
        Renamed network.
    """

    # Check if names are in columns
    msg = 'Column name "{0}" not found in net. Please specify a valid column.'
    assert source in net.columns, msg.format(source)
    assert target in net.columns, msg.format(target)
    if weight is not None:
        assert weight in net.columns, msg.format(weight) + """Alternatively, set to None if no weights are available."""
    else:
        net = net.copy()
        net['weight'] = 1.0
        weight = 'weight'

    # Rename
    net = net.rename(columns={source: 'source', target: 'target', weight: 'weight'})

    # Sort
    net = net.reindex(columns=['source', 'target', 'weight'])

    # Check if duplicated
    is_d = net.duplicated(['source', 'target']).sum()
    if is_d > 0:
        raise ValueError('net contains repeated edges, please remove them.')

    return net

def check_if_adjustText():
    try:
        import adjustText as at
    except Exception:
        raise ImportError('adjustText is not installed. Please install it with: pip install adjustText')
    return at

