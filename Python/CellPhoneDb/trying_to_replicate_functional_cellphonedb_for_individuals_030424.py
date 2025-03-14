import anndata  
import pandas as pd
import anndata as ad
import seaborn as sb
import scanpy as sc
import cellphonedb
import glob
import os
import sys

##read in object
obj = anndata.io.read_h5ad('C:/Users/Gabe/Desktop/RNA_object_human_names_anndata.h5ad')

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
cpdb_target_dir = os.path.join('A:/CellPhoneDB 030225', cpdb_version)
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
#from scipy.sparse import csr_matrix
#counts_file_path: (mandatory) paths to normalized counts file (not z-transformed), either in text format or h5ad (recommended) normalised_log_counts.h5ad.
#obj.X = csr_matrix(obj.X)
#obj.X

#import hdf5plugin
#obj.write_h5ad('A:/CellPhoneDB 030225/normalised_log_counts.h5ad')

##here is the first deviation from what I KNOW works, I will be making a csr matrix for each individual
from scipy.sparse import csr_matrix
import hdf5plugin

for i in set(obj.obs.individual):
    print(i)

    subset_object = obj[obj.obs['individual']==i]
    subset_object.X =csr_matrix(subset_object.X)

    counts_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/normalised_log_counts.h5ad"}"
    subset_object.write_h5ad(counts_path)


## QC
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    meta_file_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/human_named_meta.tsv"}"
    metadata = pd.read_csv(meta_file_path, sep = '\t')
    print(metadata.head)
# they all produce the output I want

##anndata
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    counts_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/normalised_log_counts.h5ad"}"
    adata =anndata.read_h5ad(counts_path)
    print(adata.shape)
# they all have the same number of rows (genes) 

# do the anndata obs and meta data cells match
for i in set(obj.obs.individual):
    if i == 'GH':
        continue
    meta_file_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/human_named_meta.tsv"}"
    metadata = pd.read_csv(meta_file_path, sep = '\t')  

    counts_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/normalised_log_counts.h5ad"}"
    adata =anndata.read_h5ad(counts_path)
    
    print(list(adata.obs.index).sort() == list(metadata['Cell']).sort())
#all true

adata.var
#the object indeed has genes

### run cellphonedb
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

for i in set(obj.obs.individual):
    if i == 'GH':
        continue

    meta_file_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/human_named_meta.tsv"}"
    counts_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/normalised_log_counts.h5ad"}"
    deg_path = f"{"A:/CellPhoneDB 030225/"}{i}{"/human_named_meta.tsv"}"
    cpdb_file_path = 'A:/CellPhoneDB 030225/v5.0.0/cellphonedb.zip'
    out_path = f"{"A:/CellPhoneDB 030225/"}{i}"

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
        output_suffix = None,                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
        threads = 11
        )
    












