import os
from typing import Union, List, Dict, Any
from anndata import AnnData
import requests
import json
import pandas as pd
import numpy as np
from gprofiler import GProfiler
from pathlib import Path
from .adata import collect_deg_data, create_deg_df
from ..utils import SpeciesToID, find_between
from ..calculations import jaccard_dist
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy.spatial.distance import squareform
from matplotlib import colormaps

class GOEnrichment():
    def __init__(self):
        self.results = {}
    
    def prepare_enrichment(self, adata: AnnData = None, key: str = None, key_added: str = None, 
            #uns_key_added: str = 'enrichment', return_df: bool = True, 
            sortby: str = 'scores', sign_threshold: float = 0.05
            ):

        if adata is None and key is None:
            raise ValueError("`adata` and `key` are all None.")
        if key_added is None:
            key_added = key

        # fetch DGE results
        deg = create_deg_df(adata, keys=key)

        # get groups
        groups = deg.group.unique().tolist()

        # make sure there is only one reference or group
        reference = deg.reference.unique().tolist()
        assert len(reference) == 1, "More than one group or reference"

        # sort out non-significant results
        deg = deg.query('pvals_adj < {}'.format(sign_threshold))

        # sort dataframe
        deg.sort_values(sortby, ascending=False)

        return deg, groups, key_added

    def gprofiler(self, adata: AnnData = None, key: str = None, top_n: int = 300, 
        organism: str = None, target_genes: List = None, key_added: str = None, 
        uns_key_added: str = 'gprofiler', return_df: bool = True, 
        sortby: str = 'pvals_adj', 
        #ascending: bool = True,
        **kwargs: Any):
        '''
        Performs GO term enrichment analysis using the gprofiler web resource.

        Input can be one of the following two options:
            1. adata with respective key for `adata.uns[key]`
            2. List or dictionary of genes as target_genes. With a dictionary multiple queries can be started in parallel, assuming that the keys
                are the name of the query and the values are the respective lists of genes. `key` can be specified to save the results under a specific name
                under self.result['gprofiler'][key].
        '''
        from_adata = False
        if adata is not None:
            deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added, 
                            sortby=sortby, 
                            )
            from_adata = True
        elif isinstance(target_genes, dict):
            # retrieve query names
            groups = list(target_genes.keys())
        elif isinstance(target_genes, (list, np.ndarray)):
            groups = ['query'] # set default name of query
            target_genes = {'query': target_genes}
        else:
            raise ValueError("`adata` is None and `target_genes` is neither a dictionary or a list.")

        if key_added is None:
            key_added = 'foo'
        
        if organism is None:
            raise ValueError("`organism` not specified. Needs gprofiler naming conventions, e.g. `mmusculus`")

        enrichment_dict = {}
        for i, group in enumerate(groups):
            if from_adata:
                #genes = deg.xs(group).names.tolist()
                genes = deg.query('group == "{}"'.format(group)).names.tolist()
            else:
                genes = list(target_genes[group])
            
            if top_n is not None:
                genes = genes[:top_n]

            gp = GProfiler(return_dataframe=True)
            e = gp.profile(organism=organism, query=genes, no_evidences=False)
            
            # calc -log(p_value)
            e['Enrichment score'] = [-np.log10(elem) for elem in e.p_value]
            
            # sort by p_value
            e.sort_values('p_value', inplace=True)
            
            # collect data
            enrichment_dict[group] = e

        enrichment = pd.concat(enrichment_dict)
        # rename column headers
        enrichment.rename(columns={'recall': 'Gene ratio'}, inplace=True)

        # save in class
        if not uns_key_added in self.results:
            self.results[uns_key_added] = {}
        
        self.results[uns_key_added][key_added] = enrichment

        # save in adata
        if return_df:
            return enrichment
        elif from_adata:
            if uns_key_added in adata.uns:
                adata.uns[uns_key_added][key_added] = enrichment
            else:
                adata.uns[uns_key_added] = {}
                adata.uns[uns_key_added][key_added] = enrichment
    
    def stringdb(self, adata: AnnData = None, key: str = None, top_n: int = 300, 
        organism: str = None, target_genes: List = None, key_added: str = None, 
        uns_key_added: str = 'stringdb', return_df: bool = True, 
        sortby: str = 'pvals_adj', 
        #ascending: bool = True,
        **kwargs: Any):

        '''
        Performs GO term enrichment analysis using the stringdb web resource.

        Input can be one of the following two options:
            1. adata with respective key for `adata.uns[key]`
            2. List or dictionary of genes as target_genes. With a dictionary multiple queries can be started in parallel, assuming that the keys
                are the name of the query and the values are the respective lists of genes. `key` can be specified to save the results under a specific name
                under self.result['stringdb'][key].
        '''

        from_adata = False
        if adata is not None:
            deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added,
                            sortby=sortby, 
                            #ascending=ascending
                            )
            from_adata = True
        elif isinstance(target_genes, dict):
            # retrieve query names
            groups = list(target_genes.keys())
        elif isinstance(target_genes, (list, np.ndarray)):
            groups = ['query'] # set default name of query
            target_genes = {'query': target_genes}
        else:
            raise ValueError("`adata` is None and `target_genes` is neither a dictionary or a list.")

        if key_added is None:
            key_added = 'foo'

        if organism is None:
            raise ValueError("`organism` not specified. Needs to gprofiler naming conventions, e.g. `mmusculus`")

        enrichment_dict = {}
        for i, group in enumerate(groups):
            #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
            if from_adata:
                genes = deg.query('group == "{}"'.format(group)).names.tolist()
            else:
                genes = list(target_genes[group])
            
            if top_n is not None:
                genes = genes[:top_n]

            sdb = StringDB()
            sdb.call_stringdb_enrichment(genes=genes, species=organism)
            e = sdb.result

            ## modify data to fit into further analysis
            # rename columns
            e.rename(columns={
                "category": "source",
                "description": "name",
                "term": "native",
                "inputGenes": "intersections"
            }, inplace=True)

            ## add new categories
            # e["Enrichment score"] = [-np.log10(elem) for elem in e["fdr"]]
            # e["Gene ratio"] = [a/b for a,b in zip(e["number_of_genes"], e["number_of_genes_in_background"])]
            
            # sort by p_value
            if 'Enrichment score' in e.columns:
                e.sort_values('Enrichment score', inplace=True, ascending=False)
            
            # collect data
            enrichment_dict[group] = e

        enrichment = pd.concat(enrichment_dict)

        # rename column headers
        enrichment.rename(columns={'recall': 'Gene ratio'}, inplace=True)

        # save in class
        if not uns_key_added in self.results:
            self.results[uns_key_added] = {}
            
        self.results[uns_key_added][key_added] = enrichment

        if return_df:
            return enrichment
        elif from_adata:
            if uns_key_added in adata.uns:
                adata.uns[uns_key_added][key_added] = enrichment
            else:
                adata.uns[uns_key_added] = {}
                adata.uns[uns_key_added][key_added] = enrichment

    def enrichr(self, adata: AnnData = None, key: str = None, top_n: int = 300, 
        organism: str = None, target_genes: List = None, key_added: str = None, 
        enrichr_libraries: str = 'GO_Biological_Process_2018', 
        outdir: str = None, no_plot: bool = True,
        uns_key_added: str = 'enrichr', return_df: bool = True, 
        sortby: str = 'pvals_adj', 
        #ascending: bool = True,
        **kwargs: Any):

        '''
        Performs GO term enrichment analysis using the enrichr web resource.

        Input can be one of the following two options:
            1. adata with respective key for `adata.uns[key]`
            2. List or dictionary of genes as target_genes. With a dictionary multiple queries can be started in parallel, assuming that the keys
                are the name of the query and the values are the respective lists of genes. `key` can be specified to save the results under a specific name
                under self.result['enrichr'][key].

        Parameters
        ----------
        organism : str
            Examples for possible species are `mouse` and `human`. For more informations see https://maayanlab.cloud/Enrichr/.
        enrichr_libraries : str
            To find all possible library names for a specific organism, use `gseapy.get_library_name(organism="mouse")`.
        
        '''

        import gseapy

        from_adata = False
        if adata is not None:
            deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added,
                            sortby=sortby, 
                            #ascending=ascending
                            )
            from_adata = True
        elif isinstance(target_genes, dict):
            # retrieve query names
            groups = list(target_genes.keys())
        elif isinstance(target_genes, (list, np.ndarray)):
            groups = ['query'] # set default name of query
            target_genes = {'query': target_genes}
        else:
            raise ValueError("`adata` is None and `target_genes` is neither a dictionary or a list.")

        if key_added is None:
            key_added = 'foo'

        if organism is None:
            raise ValueError("`organism` not specified. Needs to have one of the following values: `mmusculus`, `hsapiens` or `dmelanogaster`")

        enrichment_dict = {}
        for i, group in enumerate(groups):
            #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
            if from_adata:
                genes = deg.query('group == "{}"'.format(group)).names.tolist()
            else:
                genes = list(target_genes[group])

            if top_n is not None:
                genes = genes[:top_n]

            e = gseapy.enrichr(gene_list=genes, gene_sets=enrichr_libraries, 
                        organism=organism, outdir=outdir, no_plot=no_plot, 
                        **kwargs).results
            
            # calc -log(p_value)
            e['Enrichment score'] = [-np.log10(elem) for elem in e['Adjusted P-value']]
            
            # sort by p_value
            e.sort_values('Enrichment score', inplace=True, ascending=False)

            # calculate gene ratio
            e['Gene ratio'] = [int(elem.split("/")[0]) / int(elem.split("/")[1]) for elem in e['Overlap']]
            
            # collect data
            enrichment_dict[group] = e

        enrichment = pd.concat(enrichment_dict)

        # rename column headers
        enrichment.rename(columns={
            'Term': 'name',
            'Gene_set': 'source',
            'Genes': 'intersections'
            }, inplace=True)

        # transform intersections into list
        enrichment['intersections'] = [elem.split(";") for elem in enrichment['intersections']]
        
        # separate human-readable name from GO ID
        enrichment['native'] = [elem.split(" (")[1].rstrip(")") if " (GO" in elem else np.nan for elem in enrichment['name']]
        enrichment['name'] = [elem.split(" (")[0] if " (GO" in elem else elem for elem in enrichment['name']]

        # save in class
        if not uns_key_added in self.results:
            self.results[uns_key_added] = {}
        
        self.results[uns_key_added][key_added] = enrichment

        # save in adata
        if return_df:
            return enrichment
        elif from_adata:
            if uns_key_added in adata.uns:
                adata.uns[uns_key_added][key_added] = enrichment
            else:
                adata.uns[uns_key_added] = {}
                adata.uns[uns_key_added][key_added] = enrichment

class StringDB:
    def __init__(self, return_results: bool = True):
        self.result = None
        self.return_results = return_results

    def call_stringdb_enrichment(self, genes, species):
        '''
        Function to get functional enrichment results from https://string-db.org/
        Code based on https://string-db.org/help/api/
        '''

        ## Settings to call string-db
        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = "json"
        method = "enrichment"
        tax_id = SpeciesToID().convert(species)

        ## Construct the request
        request_url = "/".join([string_api_url, output_format, method])

        while True:
            ## Set parameters
            params = {

                "identifiers" : "%0d".join(genes), # your protein
                "species" : tax_id, # species NCBI identifier 
                "caller_identity" : "www.awesome_app.org" # your app name
            }

            ## Call STRING
            response = requests.post(request_url, data=params)
            response = json.loads(response.text)

            # make sure STRING found all genes
            try:
                ## Read and parse the results
                self.result = pd.DataFrame(response)
            except ValueError:
                if response['Error'] == 'not found':
                    # extract error message and identify missing gene that caused the error
                    ermsg = response['ErrorMessage']
                    missing_gene = find_between(ermsg, first="called '", last="' in the")

                    # remove missing gene from list
                    genes.remove(missing_gene)
                    print("Gene '{}' was not found by STRING and was removed from query.".format(missing_gene))
            else:
                break                  

        # rename columns to align them to gprofiler results
        self.result.rename(columns={
            "category": "source",
            "description": "name",
            "term": "native"
            }, inplace=True)

        # calculate enrichment score

        if "fdr" in self.result.columns:
            self.result['Enrichment score'] = [-np.log10(elem) for elem in self.result["fdr"]]
            self.result["Gene ratio"] = [a/b for a,b in zip(self.result["number_of_genes"], self.result["number_of_genes_in_background"])]

        if self.return_results:
            return self.result
      

    def call_stringdb_network(self, genes, species, output_format="image", prefix="", 
        output_dir="stringdb_networks", network_flavor="confidence", save=True, **kwargs):
        '''
        Generate and save networks from https://string-db.org/.
        '''

        # check output format
        format_to_ext = {
            "image": ".png",
            "svg": ".svg",
            }
        
        if output_format in format_to_ext:
            output_ext = format_to_ext[output_format]
        else:
            raise ValueError("`output_format` must be one of following values: {}".format(list(format_to_ext.keys())))

        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = output_format
        method = "network"
        tax_id = SpeciesToID().convert(species)

        ## Construct URL
        request_url = "/".join([string_api_url, output_format, method])

        ## Set parameters
        params = {

            "identifiers" : "%0d".join(genes), # your protein
            "species" : tax_id, # species NCBI identifier 
            #"add_white_nodes": 15, # add 15 white nodes to my protein 
            "network_flavor": network_flavor, # show confidence links
            "caller_identity" : "www.awesome_app.org" # your app name
        }

        ## Call STRING
        response = requests.post(request_url, data=params)
        self.result = response.content

        if save:
            ## Save the network to file
            # create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            output_file = os.path.join(output_dir, "{}_network{}".format(prefix, output_ext))
            print("Saving interaction network to {}".format(output_file))

            with open(output_file, 'wb') as fh:
                fh.write(self.result)

        if self.return_results:
            return self.result

    def stringdb_network_from_adata(self, adata: AnnData = None, key: str = None, top_n: int = 300, organism: str = None, output_format: str = "image",
        key_added: str = None, sortby: str = 'pvals_adj', ascending: bool = True,
        **kwargs: Any):

        deg, groups, key_added = GOEnrichment().prepare_enrichment(adata=adata, key=key, key_added=key_added, 
                        sortby=sortby, ascending=ascending)

        for i, group in enumerate(groups):
            #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
            if key is not None:
                target_genes = deg.xs(group).names.tolist()
            
            if top_n is not None:
                target_genes = target_genes[:top_n]

            prefix = "KEY-{}-GROUP-{}".format(key, group)

            sdb = StringDB(return_results=False)
            sdb.call_stringdb_network(genes=target_genes, species=organism, prefix=prefix, output_format=output_format, save=True, **kwargs)

def cluster_goterms_by_genes(enrichment, cmap, gene_cat="intersections", dist_thresh=None):
    '''
    Cluster GO terms by genes that belong to the respective terms.
    
    Currently only works with a list of hexadecimal colors as `cmap`!
    
    Returns: enrichment, Z, dist_thresh
    '''
    if len(enrichment) > 1:
        # get gene lists
        gene_inters = enrichment[gene_cat].tolist()

        # calculate distance matrix using Jaccard distance
        dm = np.array([[jaccard_dist(a, b) for a in gene_inters] for b in gene_inters])
        dm = pd.DataFrame(dm, columns=enrichment.native, index=enrichment.native)

        # calculate linkages of clustering tree
        Z = linkage(squareform(dm))

        # create dendrogram
        set_link_color_palette(cmap)

        # get clusters
        if dist_thresh is None:
            dist_thresh = 0.7*np.max(Z[:,2]) # calculate distance threshold used by dendrogram function
            
        clusters = pd.Series(fcluster(Z, t=dist_thresh, criterion='distance'), index=dm.columns, name="cluster")

        # map colors to clusters
        lut = dict(zip(sorted(clusters.unique()), cmap))
        clustercolors = clusters.map(lut)

        # add information to dataframe
        enrichment["clustercolor"] = clustercolors.values
        enrichment["cluster_id"] = clusters.values
        
    else:
        enrichment["clustercolor"] = 'k'
        enrichment["cluster_id"] = 1
        
    return enrichment, Z, dist_thresh

