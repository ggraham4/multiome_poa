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

# CHANGE 1: Fixed function definition - removed default parameter that wasn't defined yet
def subset_source_indices(cell_type, THRESH=4, seed=None, n_boot=100):
    """
    subset equal numbers of cells from each source for each or N_boot times
    Parameters
    ----------
    cell_type : TYPE -- np.array of celltype lables to randomly sample from
        DESCRIPTION.
    THRESH : TYPE, optional --number of cells per source cell
        DESCRIPTION. The default is 4.
    seed : TYPE -- random seed 
        DESCRIPTION. The default is None.
    n_boot : TYPE, optional -- number of different indicies to generate
        DESCRIPTION. The default is 100.
    Returns
    -------
    Indicies
    """
    np.random.seed(seed)
    tmp = [] #initialize a list
    for n in np.unique(cell_type):
        uq_vals = np.where((cell_type == n))[0] # finds indices where cell_type equals n
        # CHANGE 2: Added safety check to ensure we have enough cells
        if len(uq_vals) >= THRESH:
            tmp.append(
                uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]] # randomly selects THRESH number of cells
                )
    indices = np.vstack(tmp) # coalesce list
    return indices

adata = sc.read('A:/support_vector_classifier_zeppilli_et_al_04_27_2025/data/RNA_object_new_neurons.h5ad')
#make log counts matrix
log_counts = adata.X.toarray()
#subset out genes that have low expression
genes_use_indicies = (log_counts>0).mean(0).flatten() >=0.005
#split into immature and non immatuer
is_imm = (adata.obs['clusters']=='new_neurons').values
not_imm = (adata.obs['clusters']!='new_neurons') & (adata.obs['clusters']!='2' ) & (adata.obs['clusters']!=2 ).values
#make two new matrices, one of mature cells, and another of immature
cells_imm_genes_use = log_counts[is_imm][:, genes_use_indicies] # equivalent to cells imm
cells_non_imm_genes_use = log_counts[not_imm][:, genes_use_indicies] # equivalent to cells all train

# CHANGE 3: Extract labels for non-immature cells only
labels = adata.obs['clusters'][not_imm].values
le = LabelEncoder().fit(labels)
cell_idx = le.transform(labels)

#select 250 cells from each celltype
N_BOOT = 100
# CHANGE 4: Use the correct cell_idx that's only for non-immature cells
indices = subset_source_indices(cell_idx, 250, seed=123, n_boot=N_BOOT)
y = cell_idx[indices[:,0]]

imm_pipe = Pipeline(
    steps = [
        ('fs', sfs.SelectKBest(k=100)), #selects best features
        ('scaler',StandardScaler()), #standardizes features
        ('reduc_dim', PCA(n_components=25)), #reduces dimensionality to 25 components
        ('svc', SVC(kernel= 'linear', C= 0.1, probability=True)) #classifier
        ])

outputs = []
for i in trange(N_BOOT):
    np.random.seed(i)
    try:
        this_X = cells_non_imm_genes_use[indices[:,i]]
        gene_any = this_X.any(0)
        # CHANGE 5: Fixed syntax error with asterisks
        _X = this_X[:,gene_any]
        kf = StratifiedKFold(n_splits = 5, shuffle=True, random_state=i) #splits data into 5 folds for cross-validation
        # CHANGE 6: Fixed syntax error with asterisk
        preds = cross_val_predict(imm_pipe, _X, y, cv=kf, n_jobs=-1, verbose=False) #cross-validates predictions
        imm_pipe.fit(_X, y)
        
        # apply fitted model to held out cells - make sure we use the same gene subset
        imm_X = cells_imm_genes_use[:, gene_any]
        imm_pred = imm_pipe.predict(imm_X)
        
        #repeat for shuffle
        y_shuff = np.random.permutation(y)
        shuff_preds = cross_val_predict(
            imm_pipe,
            _X, 
            y_shuff, 
            cv = kf, 
            n_jobs = -1, 
            verbose = False)
        imm_pipe.fit(_X, y_shuff)
        shuff_imm_pred = imm_pipe.predict(imm_X)
        outputs.append((preds, imm_pred, shuff_preds, shuff_imm_pred))
    except Exception as e:
        print(f"Error in iteration {i}: {str(e)}")
        continue

# CHANGE 7: Added code to export results like in the original
for i, kind in enumerate(("observed", "shuffled")):
    # Get unique classes that actually appear in the predictions and true values
    # to avoid shape mismatch errors
    all_classes = sorted(np.unique(np.concatenate([y] + [o[i * 2] for o in outputs])))
    n_classes = len(all_classes)
    
    # Create mapping for classes
    class_mapping = {cls: idx for idx, cls in enumerate(all_classes)}
    
    # Build confusion matrices with explicit labels
    cmats = []
    for o in outputs:
        y_true = np.array([class_mapping[cls] for cls in y])
        y_pred = np.array([class_mapping[cls] for cls in o[i * 2]])
        cm = metrics.confusion_matrix(y_true, y_pred, labels=range(n_classes))
        cmats.append(cm)
    
    cmats = np.stack(cmats)
    cmat_mean = (cmats / cmats.sum(1)[:, np.newaxis]).mean(0) * 100
    
    # Create DataFrames with consistent labels
    class_names = [le.inverse_transform([cls])[0] for cls in all_classes]
    df_cmat_mean = pd.DataFrame(cmat_mean, index=class_names, columns=class_names)
    
    # Count predictions for immature cells
    pred_counts_list = []
    for o in outputs:
        pred_counts = np.zeros(n_classes)
        unique_preds, counts = np.unique(o[i * 2 + 1], return_counts=True)
        for pred, count in zip(unique_preds, counts):
            if pred in class_mapping:
                pred_idx = class_mapping[pred]
                pred_counts[pred_idx] = count
        pred_counts_list.append(pred_counts)
    
    df_pred_counts = pd.DataFrame(
        np.stack(pred_counts_list),
        columns=class_names,
    )
    
    # Calculate percentages
    df_pred_norm = df_pred_counts.divide(df_pred_counts.sum(1), axis=0) * 100
    
    # Save results (adjust path as needed)
    export_folder = Path('A:/support_vector_classifier_zeppilli_et_al_04_27_2025/results')
    export_folder.mkdir(exist_ok=True)
    
    import fastparquet
    df_pred_norm.to_parquet(
        export_folder / f"train_clusters_predict_new_neurons_pred_norm_{N_BOOT}_restarts_{kind}.parquet"
    )
    df_cmat_mean.to_parquet(
        export_folder / f"train_clusters_predict_new_neurons_df_cmat_{N_BOOT}_restarts_{kind}.parquet"
    )