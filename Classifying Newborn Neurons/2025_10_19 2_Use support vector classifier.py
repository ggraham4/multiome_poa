from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import metrics
from sklearn.decomposition import PCA
import sklearn.feature_selection as sfs
from sklearn.model_selection import (
    cross_val_predict,
    StratifiedKFold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC
from tqdm.auto import trange
import anndata as ad


# select random cells per cluster
def selectNCells(labels, THRESH=4, seed=None, n_boot=100):
    np.random.seed(seed)
    output_list = []
    for cluster in np.unique(labels):
        cluster_indices = np.where(labels == cluster)[0]
        chosen = cluster_indices[np.random.rand(len(cluster_indices), n_boot).argsort(0)[:THRESH, :]]
        output_list.append(chosen)
    return np.vstack(output_list)    

######

#read in data
obj1 = ad.read_h5ad("A:/2025_10_19_obj_rgc_subclusters.h5ad")
obj1 = obj1[obj1.obs['individual'] != 'GH']

#subset to ONLY sub.cluster 1_2 and the neuronal clusters
obj_subset= obj1[~obj1.obs['sub.cluster'].isin(['1_0',
                                             '1_1',
                                             '1_3',
                                             '1_4', # all the RGC subclusters not 1_2
                                             '1_5',
                                             '26', #dividing glia
                                             '2', #oligodendrocytes
                                             '13', #opcs
                                             '11', #microglia
                                             '15', #ependymal
                                             '20']# leuko
                                             )]

obj = obj_subset

#Only use genes with a greater than 0.005 mean expression
genes_use = (obj.X > 0).mean(0).A.flatten() >= 0.005

#define populating neurons
newbornNeurons =(obj.obs['sub.cluster']=='1_2').values
neurons = (obj.obs['sub.cluster']!='1_2').values

#define training data
#define training data
X_newborn = obj.X[newbornNeurons][:, genes_use].toarray()
X_neuron = obj.X[neurons][:, genes_use].toarray()

#encode labels
labels = obj.obs['sub.cluster'][neurons]
le = LabelEncoder().fit(labels)
cell_idx = le.transform(labels)

# subselect equal numbers of cells (250) for each celltype
N_BOOT=100
indices = selectNCells(labels, 250, n_boot=N_BOOT)
y = cell_idx[indices[:, 0]]

#Select hyperparameters
imm_pipe = Pipeline(
    steps=[
        ("fs", sfs.SelectKBest(k=100)),
        ("scaler", StandardScaler()),
        ("reduc_dim", PCA(n_components=25)),
        ("svc", SVC(kernel="linear", C=0.1, probability=True)),
    ]
)


outputs = []
for i in trange(N_BOOT):
    np.random.seed(i)
    this_X = X_neuron[indices[:, i]]
    # remove genes that aren't expressed in any cell
    gene_any = this_X.any(0)
    _X = this_X[:, gene_any]
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=i)
    preds = cross_val_predict(imm_pipe, _X, y, cv=kf, n_jobs=1, verbose=False)
    imm_pipe.fit(_X, y)
    # apply fitted model to held-out immature cells
    imm_pred = imm_pipe.predict(X_newborn[:, gene_any])
    # repeat for shuffle
    y_shuff = np.random.permutation(y)
    shuff_preds = cross_val_predict(
        imm_pipe, _X, y_shuff, cv=kf, n_jobs=1, verbose=False
    )
    imm_pipe.fit(_X, y_shuff)
    shuff_imm_pred = imm_pipe.predict(X_newborn[:, gene_any])
    outputs.append((preds, imm_pred, shuff_preds, shuff_imm_pred))
### ok it worked the tranges part just doesnt

### havent done this part yet
# export mean accuracies for training and generalization
export_folder = Path("C:/Users/Gabe/Desktop/multiome_poa/Classifying Newborn Neurons/output")

for i, kind in enumerate(("observed", "shuffled")):
    cmats = np.stack([metrics.confusion_matrix(y, o[i * 2]) for o in outputs])
    cmat_mean = (cmats / cmats.sum(1)[:, np.newaxis]).mean(0) * 100
    df_cmat_mean = pd.DataFrame(cmat_mean, index=le.classes_, columns=le.classes_)
    df_pred_counts = pd.DataFrame(
        np.stack(
            [np.bincount(o[i * 2 + 1], minlength=len(le.classes_)) for o in outputs]
        ),
        columns=le.classes_,
    )
    df_pred_norm = df_pred_counts.divide(df_pred_counts.sum(1), axis=0) * 100
    df_pred_norm.to_parquet(
        export_folder
        / f"2025_10_20_{N_BOOT}_restarts_{kind}.parquet"
    )
    df_cmat_mean.to_parquet(
        export_folder
        / f"2025_10_20_again_{N_BOOT}_restarts_{kind}.parquet"
    )



