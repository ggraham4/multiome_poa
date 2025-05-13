import pandas as pd
import anndata as ad
import seaborn as sb
import scanpy as sc
import cellphonedb
import glob
import os
import sys

##read in object
obj = anndata.io.read_h5ad('A:/optimal_clustering_05_06_2025/human_names_anndata.h5ad')

##check columns
obj.obs.columns

##look at cellphone db versions
from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils

display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))

## make a folder and put the v5.0.0 db in that folder
import os
import ssl
import urllib.request

# Create directory if it doesn't exist
cpdb_version = 'v5.0.0'
cpdb_target_dir = os.path.join('A:/cellphonedb_05_12_2025', cpdb_version)
os.makedirs(cpdb_target_dir, exist_ok=True)

# Create a custom opener with the unverified context
ssl_context = ssl._create_unverified_context()
opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl_context))
urllib.request.install_opener(opener)

# Import and download database
from cellphonedb.utils import db_utils
db_utils.download_database(cpdb_target_dir, cpdb_version)

##check the path
cpdb_target_dir

##how that the object matrix indeed has values=
print(obj.X)

## save as normalized log counts
from scipy.sparse import csr_matrix
#counts_file_path: (mandatory) paths to normalized counts file (not z-transformed), either in text format or h5ad (recommended) normalised_log_counts.h5ad.
obj.X = csr_matrix(obj.X)
obj.X

import hdf5plugin
obj.write_h5ad('A:/cellphonedb_05_12_2025/normalised_log_counts.h5ad')

## define paths
cpdb_file_path = 'A:/cellphonedb_05_12_2025/v5.0.0/cellphonedb.zip'
meta_file_path = 'A:/cellphonedb_05_12_2025/human_named_meta.tsv' #has been made in R
counts_file_path  = 'A:/cellphonedb_05_12_2025/normalised_log_counts.h5ad' 
out_path = 'A:/cellphonedb_05_12_2025/'
degs_file_path = 'A:/cellphonedb_05_12_2025/human_named_DEGs.tsv' #has been made in R
degs =pd.read_csv(degs_file_path, sep = '\t')
degs['cluster'].value_counts()

## inspect files
metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata 

adata = anndata.read_h5ad(counts_file_path)
adata.shape

list(adata.obs.index).sort() == list(metadata['Cell']).sort()

obj.var

## Run cellphoneDB
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results = cpdb_degs_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
    counts_data = 'gene_name',                                # defines the gene annotation in counts matrix.
    score_interactions = True,                                  # optional: whether to score interactions or not. 
    threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
    separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_path = out_path,                                     # Path to save results
    output_suffix = None,                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
    threads = 25
    )

# print keys
print(cpdb_results.keys())

### ok so what I did here is run it based on ALL cells, so now what I need to do is run it for each individual
for i in set(obj.obs.individual):
    print(i)

    subset_object = obj[obj.obs['individual']==i]
    subset_object.X =csr_matrix(subset_object.X)
    os.mkdir(f"A:/cellphonedb_05_12_2025/{i}")

    counts_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/normalised_log_counts.h5ad"}"
    subset_object.write_h5ad(counts_path)

## QC
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    meta_file_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/human_named_meta.tsv"}"
    metadata = pd.read_csv(meta_file_path, sep = '\t')
    print(metadata.head)
#### need to write a human named meta into each one ###

##anndata
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    counts_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/normalised_log_counts.h5ad"}"
    adata =anndata.read_h5ad(counts_path)
    print(adata.shape)
# they all have the same number of rows (genes) 

# do the anndata obs and meta data cells match
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    meta_file_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/human_named_meta.tsv"}"
    metadata = pd.read_csv(meta_file_path, sep = '\t')  

    counts_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/normalised_log_counts.h5ad"}"
    adata =anndata.read_h5ad(counts_path)
    
    print(list(adata.obs.index).sort() == list(metadata['Cell']).sort())
#all true

adata.var
#the object indeed has genes

### run cellphonedb
from cellphonedb.src.core.methods import cpdb_degs_analysis_method
import multiprocessing
import IPython

for i in set(obj.obs.individual):
    if i == 'GH':
        continue

    meta_file_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/human_named_meta.tsv"}"
    counts_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/normalised_log_counts.h5ad"}"
    deg_path = f"{"A:/cellphonedb_05_12_2025/"}{i}{"/DEGs.tsv"}"
    cpdb_file_path = 'A:/cellphonedb_05_12_2025/v5.0.0/cellphonedb.zip'
    out_path = f"{"A:/cellphonedb_05_12_2025/"}{i}"

    cpdb_results = cpdb_degs_analysis_method.call(
        cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
        meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
        counts_file_path = counts_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
        degs_file_path = deg_path,                            # mandatory: tsv file with DEG to account.
        counts_data = 'gene_name',                                # defines the gene annotation in counts matrix.
        score_interactions = True,                                  # optional: whether to score interactions or not. 
        threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
        result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
        separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
        output_path = out_path,                                     # Path to save results
        output_suffix = i,                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
        threads = 11
        )


