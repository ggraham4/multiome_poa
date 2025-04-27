#Modifying 03_classification.py from zeppilli et al., to create a support vector classifier to classify which clusters new neurons are populating

## Libraries
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import metrics
from sklearn.decomposition import PCA
import sklearn.feature_selection as sfs
from sklearn.model_selection import(
    cross_val_predict,
    StratifiedGroupKFold
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import(
    StandardScaler,
    LabelEncoder
)
from sklearn.svm import SVC
from tqdm.auto import trange

### Alright here the authors define this function and my goal is to understand what it does
def subset_source_indices(cell_type, THRESH=4, source=None, seed=None, n_boot=100):
    """Subset equal numbers of cells from each source for each OR n_boot  times

    Arguments:
        cell_type {[type]} -- [np.array of celltype labels to randomly sample from]

    Keyword Arguments:
        SOURCE_THRESH {int} -- [number of cells per source per cell_type] (default: {4})
        source {[type]} -- [vector describing source for each cell] (default: {None}) # not sure what this means
        seed {[type]} -- [random seed] (default: {None})
        n_boot  {int} -- [number different indices to generate] (default: {100})
    """
    np.random.seed(seed)
    tmp = []
    for n in np.unique(cell_type): ## to me it looks like cell_type is a cluster label, for each unique cell type
        if source is None: #not entirely sure what this means
            uq_vals = np.where((cell_type == n))[0] # uq vals is finding which indicies are of the appropriate cell type, [0] is i think because it is in a list
            tmp.append( #to the temp list
                uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]] # add a random set of cells (number = n_boot) that is unsorted. Not entirely sure what Thresh is dong here
            )
        else: # if there is a  source
            for s in np.unique(source): #for each unique value in source
                uq_vals = np.where((cell_type == n) & (source == s))[0] #find the cells of that cell that belong to that cell type of that source
                tmp.append(
                    uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]] # and do the same thing
                )
    indices = np.vstack(tmp) #stacks the indicies into a single 2d array
    return indices
#so my conclusion is subset_source_indicies randomly samples cells in a bootstrapped manner

N_BOOT = 100 #100 restarts
base = Path('a:/support_vector_classifier_zeppilli_et_al_04_27_2025') #here i am going to define this as a folder in my A drive
pcx_fold = base / "data"
export_folder = base / "export"

adata = sc.read(pcx_fold / "RNA_object_anndata.h5ad")



