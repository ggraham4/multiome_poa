#condas activate mylimit
#condas install -c conda-forge package...

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
    StratifiedKFold
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import(
    StandardScaler,
    LabelEncoder
)
from sklearn.svm import SVC
from tqdm.auto import trange
