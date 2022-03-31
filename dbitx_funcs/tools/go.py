import os
from typing import Union, List, Dict, Any
from anndata import AnnData
import requests
import json
import pandas as pd
import numpy as np
from gprofiler import GProfiler
import gseapy
from pathlib import Path
from .adata import collect_deg_data
from ..utils import SpeciesToID

class GOEnrichment():
    def __init__(self):
        self.result = {}
    
    def prepare_enrichment(self, adata: AnnData = None, key: str = None, key_added: str = None, 
            #uns_key_added: str = 'enrichment', return_df: bool = True, 
            sortby: str = 'pvals_adj', ascending: bool = True):

        if adata is None and key is None:
            raise ValueError("`adata` and `key` are all None.")
        if key_added is None:
            key_added = key

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

    def gprofiler(self, adata: AnnData = None, key: str = None, top_n: int = 300, organism: str = None, target_genes: List = None, key_added: str = None, 
        uns_key_added: str = 'gprofiler_enrichment', return_df: bool = True, sortby: str = 'pvals_adj', ascending: bool = True,
        **kwargs: Any):

        deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added, 
                        sortby=sortby, ascending=ascending)
        
        if organism is None:
            raise ValueError("`organism` not specified. Needs gprofiler naming conventions, e.g. `mmusculus`")

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

        # save in class
        if not "gprofiler" in self.result:
            self.result["gprofiler"] = {}
        
        self.result["gprofiler"][key] = enrichment

        if return_df:
            return enrichment
        else:
            if uns_key_added in adata.uns:
                adata.uns[uns_key_added][key_added] = enrichment
            else:
                adata.uns[uns_key_added] = {}
                adata.uns[uns_key_added][key_added] = enrichment

    
    def stringdb(self, adata: AnnData = None, key: str = None, top_n: int = 300, organism: str = None, target_genes: List = None, key_added: str = None, 
        uns_key_added: str = 'stringdb_enrichment', return_df: bool = True, sortby: str = 'pvals_adj', ascending: bool = True,
        **kwargs: Any):

        deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added,
                        sortby=sortby, ascending=ascending)

        if organism is None:
            raise ValueError("`organism` not specified. Needs to gprofiler naming conventions, e.g. `mmusculus`")

        enrichment_dict = {}
        for i, group in enumerate(groups):
            #target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()
            if key is not None:
                target_genes = deg.xs(group).names.tolist()
            
            if top_n is not None:
                target_genes = target_genes[:top_n]

            sdb = StringDB()
            sdb.call_stringdb_enrichment(genes=target_genes, species=organism)
            e = sdb.result

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

        # save in class
        if not "stringdb" in self.result:
            self.result["stringdb"] = {}
            
        self.result["stringdb"][key] = enrichment

        if return_df:
            return enrichment
        else:
            if uns_key_added in adata.uns:
                adata.uns[uns_key_added][key_added] = enrichment
            else:
                adata.uns[uns_key_added] = {}
                adata.uns[uns_key_added][key_added] = enrichment

    def enrichr(self, adata=None, key=None, top_n=300, organism=None, target_genes=None, key_added=None, 
        enrichr_libraries='GO_Biological_Process_2018', outdir=None, no_plot=True,
        uns_key_added='enrichment', return_df=True, sortby='pvals_adj', ascending=True,
        **kwargs):

        deg, groups, key_added = self.prepare_enrichment(adata=adata, key=key, key_added=key_added, 
                        sortby=sortby, ascending=ascending)

        if organism is None:
            raise ValueError("`organism` not specified. Needs to have one of the following values: `mmusculus`, `hsapiens` or `dmelanogaster`")

        enrichment_dict = {}
        for i, group in enumerate(groups):
            target_genes = deg.xs((group, 'up'), level=(0,1)).names.tolist()

            if top_n is not None:
                    target_genes = target_genes[:top_n]

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

        # save in class
        if not "enrichr" in self.result:
            self.result["enrichr"] = {}

        self.result["enrichr"][key] = enrichment

        if return_df:
            return enrichment
        else:
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

        ## Set parameters
        params = {

            "identifiers" : "%0d".join(genes), # your protein
            "species" : tax_id, # species NCBI identifier 
            "caller_identity" : "www.awesome_app.org" # your app name
        }

        ## Call STRING
        response = requests.post(request_url, data=params)

        ## Read and parse the results
        self.result = pd.DataFrame(json.loads(response.text))

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

