import requests
from tqdm import tqdm
from pathlib import Path
import pandas as pd
from typing import Union, List, Tuple, Dict, Any
import numpy as np

def parse_gmt(filepath: str) -> pd.DataFrame:
    '''
    Function to parse GMT files.
    '''
    
    # Using readlines()
    file = open(filepath, 'r')
    lines = file.readlines()

    terms = []
    extras = []
    genesets = []
    for line in lines:
        splits = line.split("\t")
        terms.append(splits[0]) # term
        extras.append(splits[1]) # extra information. Often empty.
        genesets.append(splits[2:]) # gene sets in this term

    # create dataframe
    data = pd.DataFrame({
        "geneset": terms,
        "extra_information": extras,
        "genesymbol": genesets
    }).explode('genesymbol')
    data.index = range(len(data))
    
    return data

def get_enrichr_geneset(names: Union[str, List], 
                        output_dir: str = "tmp/enrichr", 
                        return_data: bool = True, 
                        decoupler_format: bool = True
                        ) -> Dict:
    
    '''
    Function to get enrichr dataset from https://maayanlab.cloud/Enrichr
    '''
    # convert to list
    names = [names] if isinstance(names, str) else list(names)

    # create output directory
    output_dir = Path(output_dir)
    
    filenames = []
    for name in names:
        # create file name
        filename = output_dir / "{}.gmt".format(name)
        filenames.append(filename)
        
        # check if output file name exists
        if not filename.exists():
            print('Downloading geneset "{}"...'.format(name))
            # get response from website
            url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}".format(name)
            response = requests.get(url, stream=False)
            
            # set up progress bar
            # total_size_in_bytes= int(response.headers.get('content-length', 0))
            # block_size = 1024 #1 Kibibyte
            # progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)

            # create output directory if it does not exist already
            output_dir.mkdir(parents=True, exist_ok=True)

            # write to file
            with open(filename, "wb") as file:
                #for data in response.iter_content(block_size):
                    #progress_bar.update(len(data))
                file.write(response.content)
            file.close()
    
    if return_data:
        genesets = {}
        for i, filename in enumerate(filenames):
            #geneset = pd.read_csv(filename, sep="\t", header=None, on_bad_lines='skip')
            geneset = parse_gmt(filename)

            # if decoupler_format:
            #     # format for decoupler package
            #     geneset = geneset.melt(id_vars=0).dropna()
            #     geneset.drop('variable', axis=1, inplace=True)
            #     geneset.columns = ["geneset", "genesymbol"]

            #     # add library name
            #     geneset['collection'] = names[i]
            
            # split human-readable name and GO ID
            # separate human-readable name from GO ID
            geneset['native'] = [elem.split(" (")[1].rstrip(")") if " (GO" in elem else np.nan for elem in geneset['geneset']]
            geneset['geneset'] = [elem.split(" (")[0] if " (GO" in elem else elem for elem in geneset['geneset']]
            
            # collect results
            genesets[names[i]] = geneset
        
        # convert to Dataframe
        #genesets = pd.concat(genesets)
        
        return genesets

def get_mouse_human_conversion(output_dir: str = "tmp/biomart/", 
                               return_data: bool = True, 
                               return_dict: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Downloads information to convert mouse gene IDs to human and vice versa.

    Code from: https://gseapy.readthedocs.io/en/latest/gseapy_example.html#2.1-Inputs-to-enrichr

    Returns
    -------
    m2h, h2m : mouse-to-human, human-to-mouse conversion matrices
    '''

    from gseapy import Biomart

    # files
    output_dir = Path(output_dir)
    m2h_file = output_dir / "m2h.csv"
    h2m_file = output_dir / "h2m.csv"

    if not h2m_file.exists() or not m2h_file.exists():
        print("Files not found. Downloading data to {}".format(output_dir))
        # load data from Biomart
        bm = Biomart()
        m2h = bm.query_simple(dataset='mmusculus_gene_ensembl',
                                attributes=['ensembl_gene_id','external_gene_name',
                                            'hsapiens_homolog_ensembl_gene',
                                            'hsapiens_homolog_associated_gene_name'])

        h2m = bm.query_simple(dataset='hsapiens_gene_ensembl',
                                attributes=['ensembl_gene_id','external_gene_name',
                                            'mmusculus_homolog_ensembl_gene',
                                            'mmusculus_homolog_associated_gene_name'])

        if output_dir is not None:
            # create output dir and files
            output_dir.mkdir(parents=True, exist_ok=True)

            # save files
            m2h.to_csv(m2h_file)
            h2m.to_csv(h2m_file)

    if return_data:
        m2h = pd.read_csv(m2h_file, index_col=0)
        h2m = pd.read_csv(h2m_file, index_col=0)

        if return_dict:
            m2h = {a: b for a,b in zip(m2h.dropna()['external_gene_name'], m2h.dropna()['hsapiens_homolog_associated_gene_name'])}
            h2m = {a: b for a,b in zip(h2m.dropna()['external_gene_name'], h2m.dropna()['mmusculus_homolog_associated_gene_name'])}

        return m2h, h2m
