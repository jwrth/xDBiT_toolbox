import requests
import json
import pandas as pd
import numpy as np
from gprofiler import GProfiler
import gseapy
from .tools import collect_deg_data


def prepare_enrichment(adata=None, key=None, target_genes=None, key_added=None, 
    uns_key_added='enrichment', return_df=True, sortby='pvals_adj', ascending=True,
    **kwargs):

    if target_genes is None:
        if adata is None and key is None:
            raise ValueError("`target_genes`, `adata` and `key` are all None.")
        if key_added is None:
            key_added = key
    else:
        if adata is None:
            print("No adata given to save values. Results dataframe is returned.")
            return_df = True
        else:
            if key_added is None:
                raise ValueError("Results cannot be added to adata since no `key_added` was given.")
    
    if not return_df:
        if uns_key_added in adata.uns:
            if not isinstance(adata.uns[uns_key_added], dict):
                raise ValueError('adata.uns["{}"] exists but is no dictionary.'.format(uns_key_added))

    # get information for up and down regulation from key
    if key is not None:
        _, _, _, deg, _ = collect_deg_data(adata, keys=key)
        deg = pd.concat(deg[key])
        deg = deg.sort_values(sortby, ascending=ascending)

        if deg.reference.unique()[0] != 'rest':
            # extract group and reference
            group = deg.group.unique().tolist()
            reference = deg.reference.unique().tolist()

            # make sure there is only one reference or group
            assert len(group) == len(reference) == 1, "More than one group or reference"
            group = group[0]
            reference = reference[0]

            # give the reference `rest` a name
            up = deg.xs('up', level=1)
            down = deg.xs('down', level=1).rename(index={group: reference})

            group = down.group.copy()
            ref = down.reference.copy()
            down["reference"] = group
            down["group"] = ref

            down['scores'] = down['scores'] * -1
            down['logfoldchanges'] = down['logfoldchanges'] * -1
            down['combined'] = down['combined'] * -1

            deg = pd.concat([up, down])
            deg.drop('updown', axis=1, inplace=True)
        else:
            deg = deg.xs('up', level=1).copy()

        # get groups from key
        groups = deg['group'].unique()
    else:
        groups = [key_added]

    return deg, groups, key_added


def gprofiler(adata=None, key=None, top_n=300, organism=None, target_genes=None, key_added=None, 
    uns_key_added='gprofiler_enrichment', return_df=True, sortby='pvals_adj', ascending=True,
    **kwargs):

    deg, groups, key_added = prepare_enrichment(adata=adata, key=key, target_genes=target_genes, key_added=key_added, 
                    uns_key_added=uns_key_added, return_df=return_df, sortby=sortby, ascending=ascending, **kwargs)

    if organism is None:
        raise ValueError("`organism` not specified. Needs to have one of the following values: `mmusculus`, `hsapiens` or `dmelanogaster`")

    enrichment_dict = {}
    for i, group in enumerate(groups):
        #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
        if key is not None:
            target_genes = deg.xs(group).names.tolist()
        
        if top_n is not None:
            target_genes = target_genes[:top_n]

        gp = GProfiler(return_dataframe=True)
        e = gp.profile(organism=organism, query=target_genes, no_evidences=False)
        
        # calc -log(p_value)
        e['Enrichment score'] = [-np.log10(elem) for elem in e.p_value]
        
        # sort by p_value
        e.sort_values('p_value', inplace=True)
        
        # collect data
        enrichment_dict[group] = e

    enrichment = pd.concat(enrichment_dict)

    # rename column headers
    enrichment.rename(columns={'recall': 'Gene ratio'}, inplace=True)

    if return_df:
        return enrichment
    else:
        if uns_key_added in adata.uns:
            adata.uns[uns_key_added][key_added] = enrichment
        else:
            adata.uns[uns_key_added] = {}
            adata.uns[uns_key_added][key_added] = enrichment

def call_string_db_enrichment(genes, species):
    '''
    Function to get functional enrichment results from https://string-db.org/
    Code based on https://string-db.org/help/api/
    '''



    ## Settings to call string-db
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"
    species_dict = {
        'mmusculus': 10090,
        'hsapiens': 9606,
        'dmelanogaster': 7227
    }

    ## Checks
    if species not in species_dict:
        raise ValueError(
            "`species` must be one of following values: {}".format(list(species_dict.keys()))
            )

    ## Construct the request
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    params = {

        "identifiers" : "%0d".join(genes), # your protein
        "species" : species_dict[species], # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org" # your app name
    }

    ## Call STRING
    response = requests.post(request_url, data=params)

    ## Read and parse the results
    data = pd.DataFrame(json.loads(response.text))

    return data

def string_db_enrichment(adata=None, key=None, top_n=300, organism=None, target_genes=None, key_added=None, 
    uns_key_added='stringdb_enrichment', return_df=True, sortby='pvals_adj', ascending=True,
    **kwargs):

    deg, groups, key_added = prepare_enrichment(adata=adata, key=key, target_genes=target_genes, key_added=key_added, 
                    uns_key_added=uns_key_added, return_df=return_df, sortby=sortby, ascending=ascending, **kwargs)

    if organism is None:
        raise ValueError("`organism` not specified. Needs to have one of the following values: `mmusculus`, `hsapiens` or `dmelanogaster`")

    enrichment_dict = {}
    for i, group in enumerate(groups):
        #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
        if key is not None:
            target_genes = deg.xs(group).names.tolist()
        
        if top_n is not None:
            target_genes = target_genes[:top_n]

        e = call_string_db_enrichment(genes=target_genes, species=organism)

        ## modify data to fit into further analysis
        # rename columns
        e.rename(columns={
            "category": "source",
            "description": "name",
            "term": "native"
        }, inplace=True)

        # add new categories
        e["Enrichment score"] = [-np.log10(elem) for elem in e["fdr"]]
        e["Gene ratio"] = [a/b for a,b in zip(e["number_of_genes"], e["number_of_genes_in_background"])]
        
        # sort by p_value
        e.sort_values('Enrichment score', inplace=True, ascending=False)
        
        # collect data
        enrichment_dict[group] = e

    enrichment = pd.concat(enrichment_dict)

    # rename column headers
    enrichment.rename(columns={'recall': 'Gene ratio'}, inplace=True)

    if return_df:
        return enrichment
    else:
        if uns_key_added in adata.uns:
            adata.uns[uns_key_added][key_added] = enrichment
        else:
            adata.uns[uns_key_added] = {}
            adata.uns[uns_key_added][key_added] = enrichment


def enrichr(adata=None, key=None, top_n=300, organism=None, target_genes=None, key_added=None, 
    enrichr_libraries='GO_Biological_Process_2018', outdir=None, no_plot=True,
    uns_key_added='enrichment', return_df=True, sortby='pvals_adj', ascending=True,
    **kwargs):

    deg, groups, key_added = prepare_enrichment(adata=adata, key=key, target_genes=target_genes, key_added=key_added, 
                    uns_key_added=uns_key_added, return_df=return_df, sortby=sortby, ascending=ascending, **kwargs)

    if organism is None:
        raise ValueError("`organism` not specified. Needs to have one of the following values: `mmusculus`, `hsapiens` or `dmelanogaster`")

    enrichment_dict = {}
    for i, group in enumerate(groups):
        target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()

        e = gseapy.enrichr(gene_list=target_genes, gene_sets=enrichr_libraries, organism=organism, outdir=outdir, no_plot=no_plot, **kwargs).results
        
        # calc -log(p_value)
        e['Enrichment score'] = [-np.log10(elem) for elem in e['Adjusted P-value']]
        
        # sort by p_value
        e.sort_values('Adjusted P-value', inplace=True)
        
        # collect data
        enrichment_dict[group] = e

    enrichment = pd.concat(enrichment_dict)

    # rename column headers
    enrichment.rename(columns={'Term': 'name'}, inplace=True)
    enrichment.rename(columns={'Gene_set': 'source'}, inplace=True)

    if return_df:
        return enrichment
    else:
        if uns_key_added in adata.uns:
            adata.uns[uns_key_added][key_added] = enrichment
        else:
            adata.uns[uns_key_added] = {}
            adata.uns[uns_key_added][key_added] = enrichment


