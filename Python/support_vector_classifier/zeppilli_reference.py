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


def subset_source_indices(cell_type, THRESH=4, source=None, seed=None, n_boot=100):
    """Subset equal numbers of cells from each source for each OR n_boot  times

    Arguments:
        cell_type {[type]} -- [np.array of celltype labels to randomly sample from]

    Keyword Arguments:
        SOURCE_THRESH {int} -- [number of cells per source per cell_type] (default: {4})
        source {[type]} -- [vector describing source for each cell] (default: {None})
        seed {[type]} -- [random seed] (default: {None})
        n_boot  {int} -- [number different indices to generate] (default: {100})
    """
    np.random.seed(seed)
    tmp = []
    for n in np.unique(cell_type):
        if source is None:
            uq_vals = np.where((cell_type == n))[0]
            tmp.append(
                uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]]
            )
        else:
            for s in np.unique(source):
                uq_vals = np.where((cell_type == n) & (source == s))[0]
                tmp.append(
                    uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]]
                )
    indices = np.vstack(tmp)
    return indices


N_BOOT = 100
base = Path()
pcx_fold = base / "data"
export_folder = base / "export"

adata = sc.read(pcx_fold / "PCx_fleischmann_wild_lab_mice_merged.h5ad")
adata.layers["log"] = adata.layers["counts"]
sc.pp.normalize_total(adata, target_sum=1e4, layer="log")
sc.pp.log1p(adata, layer="log")
# only use genes
genes_use = (adata.layers["counts"] > 0).mean(0).A.flatten() >= 0.005

is_imm = (adata.obs.supertype == "Immature") & (
    adata.obs.coarse_supertype == "Immature"
)
not_imm = (adata.obs.supertype != "Immature") & (
    adata.obs.coarse_supertype != "Immature"
)

X_imm = adata.layers["log"][is_imm][:, genes_use].A
X_all_train = adata.layers["log"][not_imm][:, genes_use].A

labels = adata.obs.coarse_supertype[not_imm]
le = LabelEncoder().fit(labels)
cell_idx = le.transform(labels)
# subselect equal numbers of cells (250) for each celltype
indices = subset_source_indices(cell_idx, 250, n_boot=N_BOOT)
y = cell_idx[indices[:, 0]]

# best params found over simple grid search, but similar results across most hyperparameters
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
    this_X = X_all_train[indices[:, i]]
    # remove genes that aren't expressed in any cell
    gene_any = this_X.any(0)
    _X = this_X[:, gene_any]
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=i)
    preds = cross_val_predict(imm_pipe, _X, y, cv=kf, n_jobs=-1, verbose=False)
    imm_pipe.fit(_X, y)
    # apply fitted model to held-out immature cells
    imm_pred = imm_pipe.predict(X_imm[:, gene_any])
    # repeat for shuffle
    y_shuff = np.random.permutation(y)
    shuff_preds = cross_val_predict(
        imm_pipe, _X, y_shuff, cv=kf, n_jobs=-1, verbose=False
    )
    imm_pipe.fit(_X, y_shuff)
    shuff_imm_pred = imm_pipe.predict(X_imm[:, gene_any])
    outputs.append((preds, imm_pred, shuff_preds, shuff_imm_pred))

# export mean accuracies for training and generalization
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
        / f"Lab_wild_train_coarse_supertype_predict_immature_pred_norm_{N_BOOT}_restarts_{kind}.parquet"
    )
    df_cmat_mean.to_parquet(
        export_folder
        / f"Lab_wild_train_coarse_supertype_predict_immature_df_cmat_{N_BOOT}_restarts_{kind}.parquet"
    )