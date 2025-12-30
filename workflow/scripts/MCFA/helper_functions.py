import scanpy as sc
import pandas as pd
import muon as mu
import numpy as np
from tqdm import tqdm
from types import ModuleType


import warnings

def _check_if_mudata() -> ModuleType:
    try:
        from mudata import MuData
    except Exception:
        raise ImportError(
            'mudata is not installed. Please install it with: '
            'pip install mudata'
        )
    return MuData
## TODO generalize these functions to work with any package
def _check_if_decoupler() -> ModuleType:
    try:
        import decoupler as dc
    except Exception:
        raise ImportError(
            'decoupler is not installed. Please install it with: '
            'pip install decoupler'
        )
    return dc

def _process_meta(adata, mdata, sample_key, obs_keys):
    if obs_keys is not None:
        metadata = adata.obs[[sample_key, *obs_keys]].drop_duplicates()
        
        sample_n = adata.obs[sample_key].nunique()
        if metadata.shape[0] != sample_n:
            raise ValueError('`obs_keys` must be unique per sample in `adata.obs`')
        
        mdata.obs.index.name = None
        mdata.obs = (mdata.obs.
                     reset_index().
                     rename(columns={"index":sample_key}).
                     merge(metadata).
                     set_index(sample_key)
                     )
        for view in mdata.mod.keys():
            mdata.mod[view].obs = (mdata.mod[view].obs.
                                    reset_index().
                                    rename(columns={"index":sample_key}).
                                    merge(metadata).
                                    set_index(sample_key)
                                    )
        
def adata_to_views(adata,
                   groupby,
                   sample_key,
                   obs_keys,
                   view_separator=':',
                   min_count=10,
                   min_total_count=15,
                   large_n=10, 
                   min_prop=0.1,
                   keep_psbulk_stats = False,
                   verbose=False,
                   **kwargs):
    """
    Converts an AnnData object to a MuData object with views that represent an aggregate for each entity in `adata.obs[groupby]`.
    
    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        AnnData object
    groupby:
        Column name in `adata.obs` to group by
    sample_key:
        Column name in `adata.obs` to use as sample key
    obs_keys:
        Column names in `adata.obs` to merge with the MuData object
    view_separator:
        Separator to use when assigning `adata.var_names` to views
    min_count:
        Minimum number of counts per gene per sample to be included in the pseudobulk.
    min_total_count:
        Minimum number of counts per sample to be included in the pseudobulk.
    large_n:
        Number of samples per group that is considered to be "large".
    min_prop:
        Minimum proportion of samples that must have a count for a gene to be included in the pseudobulk.
    keep_psbulk_stats:
        If True, keep the pseudobulk statistics in `mdata.uns['psbulk_stats']`.
    verbose:
        If True, show progress bar.
    **kwargs
        Keyword arguments used to aggregate the values per cell into views. See `dc.filter_by_expr` for more details.
    
    Returns
    -------
    Returns a MuData object with views that represent an aggregate for each entity in `adata.obs[groupby]`.
    
    """
    
    # Check if MuData & decoupler are installed
    MuData = _check_if_mudata()
    dc = _check_if_decoupler()
    
    views = adata.obs[groupby].unique()
    views = tqdm(views, disable=not verbose)

    padatas = {}
    if keep_psbulk_stats: stats = []
    for view in (views):
        # filter AnnData to view
        temp = adata[adata.obs[groupby] == view].copy()
        # assign view to var_names
        temp.var_names = view + view_separator + temp.var_names
        padata = dc.get_pseudobulk(temp,
                                   sample_col=sample_key,
                                   groups_col=None,
                                   **kwargs
                                   )
        
        # only filter genes for views that pass QC
        if 0 in padata.shape:
            continue
        # edgeR filtering
        feature_mask = dc.filter_by_expr(padata,
                                         min_count=min_count,
                                         min_total_count=min_total_count,
                                         large_n=large_n,
                                         min_prop=min_prop,
                                         )
        padata = padata[:, feature_mask]

        # only append views that pass QC
        if 0 not in padata.shape:
            # keep psbulk stats
            if keep_psbulk_stats:
                df = padata.obs.filter(items=['psbulk_n_cells', 'psbulk_counts'], axis=1)
                df.columns = [view + view_separator + col for col in df.columns]
                stats.append(df)

            del padata.obs
            padatas[view] = padata

    # Convert to MuData
    mdata = MuData(padatas)
    
    # process metadata
    _process_meta(adata=adata, mdata=mdata, sample_key=sample_key, obs_keys=obs_keys)

    # combine psbulk stats across views and add to mdata
    if keep_psbulk_stats:
        mdata.uns['psbulk_stats'] = pd.concat(stats, axis=1)

    return mdata

def preprocess_data(mdata,
                    target_sum = 10000,
                    exclude_highly_expressed = False,
                    remove_mitochondrial = False,
                    view_separator = ':',
                    batch_key = None,
                    n_top_genes = None):
    """
    Preprocesses data in a MuData object.

    Parameters:
    ----------
    mdata :class:`MuData`
        The MuData object.
    target_sum :class:`int`, optional:
        The target sum of total counts after normalization. Parameter passed on to scanpy.pp.normalize_total. Defaults to 10000.
    exclude_highly_expressed :class:`bool`, optional:
        Whether to exclude highly expressed genes during normalization. Parameter passed on to scanpy.pp.normalize_total. Defaults to False.
    remove_mitochondrial :class:`bool`, optional:
        Whether to remove mitochondrial genes. Defaults to False.
    view_separator :class:`str`, optional:
        Separator used to identify view-specific genes. Defaults to ':'.
    """

    for view in mdata.mod.keys():
        print('INFO: Preprocessing {0}'.format(view))
        mdata.mod[view].X = np.nan_to_num(mdata.mod[view].X, copy=False) # remove nan values

        if remove_mitochondrial:
            #identify mitochondrial genes
            mdata.mod[view].var['mt'] = mdata.mod[view].var_names.str.startswith(view + view_separator + 'MT-')
            #remove mitochondrial genes from analysis
            mdata.mod[view] = mdata.mod[view][:, ~mdata.mod[view].var['mt']].copy()

        sc.pp.normalize_total(mdata.mod[view], target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed)
        sc.pp.log1p(mdata.mod[view])

        if mdata.mod[view].shape[0] > 1:
            sc.pp.highly_variable_genes(mdata.mod[view],
                                        batch_key = batch_key,
                                        n_top_genes = n_top_genes)
            print('Highly variable genes in {0}: {1}'.format(view, mdata.mod[view].var.highly_variable.sum()))
        else:
            print('Not enough cells in {0}'.format(view))
            mdata.mod[view].var['highly_variable'] = False
    return mdata

def remove_marker_genes(mdata,
                        markers,
                        view_separator=':'):
    """
    In each view, removes genes if they are in the markers dict for another view, but not if they are in the markers for the same view.
    Used for removing potential cell type marker genes found in the background of other views and thought to be contamination.

    Parameters
    ----------
    mdata :class:`~mudata.MuData`
        MuData object. Highly variable genes should be computed already in .var for each modality.
    markers :class:`dict`
        Dictionary with markers for each view. Keys are the views and values are lists of markers. Can contain markers for views that are not in mdata.mod.keys().
    view_separator :class:`str`, optional
        Separator between view and gene names. Defaults to ':'.
    """
    # check if markers is a dict
    if not isinstance(markers, dict):
        raise TypeError('markers is not a dict')

    # check that all keys in markers are lists
    if not all(isinstance(markers[mod], list) for mod in markers.keys()):
        raise TypeError('not all values in markers are lists')
    
    for current_mod in mdata.mod.keys():
        # markers in markers dict for each modality except for current_mod
        negative_markers = [marker for mod in markers.keys() if mod != current_mod for marker in markers[mod]]

        if current_mod not in list(markers.keys()):
            warnings.warn('no markers in dict for view: {0}'.format(current_mod), Warning)
        else:
            #keep negative_markers not in markers[current_mod] and add view_separator
            negative_markers = [current_mod + view_separator + marker for marker in negative_markers if marker not in markers[current_mod]]
        
        #identify negative marker genes
        mdata.mod[current_mod].var['negative_marker'] = mdata.mod[current_mod].var_names.isin(negative_markers)
        #remove negative marker genes from analysis
        mdata.mod[current_mod] = mdata.mod[current_mod][:, ~mdata.mod[current_mod].var['negative_marker']].copy()
    mdata.update()
    return mdata

def get_psbulk_stats(mdata, stat, view_separator = ':', uns_key = 'psbulk_stats'):
    """Get cell counts for each patient and omic in a MuData object.
    Parameters
    ----------
    mdata : MuData
        MuData object containing the data.
    stat : str
        Name of the stat to get. Must be in mdata.uns[uns_key].columns. e.g. 'psbulk_n_cells' or 'psbulk_counts'
    view_separator : str, optional
        Separator used to separate the view name from column names in the MuData object (e.g. of var, or psbulk_stats), by default ':'
    uns_key : str, optional
        Key in mdata.uns where the cell counts will is stored, by default 'stats_psbulk'
    Returns
    -------
    pd.DataFrame
        Dataframe containing the cell counts for each patient and omic.
    """

    assert uns_key in mdata.uns.keys(), 'Key {0} not found in mdata.uns'.format(uns_key)
    stats = mdata.uns[uns_key]

    #get all columns with psbulk_n_cells in the name
    stats = stats[[col for col in stats.columns if stat in col]]

    #split column names by : and get the first element
    stats.columns = [col.split(view_separator)[0] for col in stats.columns]

    #keep only columns with name in mdata.mod.keys()
    stats = stats[[col for col in stats.columns if col in mdata.mod.keys()]]

    #order columns alphabetically
    stats = stats[sorted(stats.columns)]

    return stats