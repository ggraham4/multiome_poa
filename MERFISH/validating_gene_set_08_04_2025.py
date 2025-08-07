### testing performance of MERFISH genes
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
    

## read in expression matrix
expression_matrix = pd.read_csv('A:/merfish_mutual_information_07_31_2025/normalized_expression.csv',index_col=0)
clusters = pd.read_csv('A:/merfish_mutual_information_07_31_2025/subclusters.csv')['clusters']

#binarize 
expression_binary = (expression_matrix > 0).astype(int).reset_index()
del expression_matrix

#read in both sets of marker genes and concatenate
minimal_set = pd.read_csv('C:/Users/Gabe/Desktop/multiome_poa/MERFISH/minimal_set_MI_08_04_2025.csv')
maker_genes = pd.read_csv('A:/merfish_mutual_information_07_31_2025/markers_top4_08_04_2025.csv')

merfish_candidates = minimal_set['symbol'].to_list() + maker_genes['gene'].to_list()

#reduce expression_binary to this set
expression_binary_markers = expression_binary.loc[expression_binary['index'].isin(merfish_candidates)]
del expression_binary

### train classifier on gex ###
X = expression_binary_markers.T.drop('index')
y = clusters

# Initialize classifier
clf = RandomForestClassifier(n_estimators=100, random_state=1)

# Perform cross-validation
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)
cv_scores = cross_val_score(clf, X, y, cv=cv, scoring='accuracy')

# Calculate metrics
mean_accuracy = np.mean(cv_scores)
std_accuracy = np.std(cv_scores)

# Fit on full data to get feature importances
clf.fit(X, y)
feature_importance = dict(zip(expression_binary_markers.index, clf.feature_importances_))

results = {
    'mean_cv_accuracy': mean_accuracy,
    'std_cv_accuracy': std_accuracy,
    'cv_scores': cv_scores,
    'n_genes': len(expression_binary_markers.index),
    'n_clusters': len(np.unique(y)),
    'feature_importance': feature_importance
}
results

# I think this gene set is fine
cell_vs_cluster = pd.DataFrame({'cell':expression_binary_markers.iloc[:,1:].columns,'cluster' : clusters})

###

pd.DataFrame(merfish_candidates).to_csv('A:/merfish_mutual_information_07_31_2025/markers_multiome_obj_only.csv')
