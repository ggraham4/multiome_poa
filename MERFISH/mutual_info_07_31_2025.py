## here, I am going to use mutual information to find most informative
#> marker genes in my dataset in the manner described in li et al.,
#. https://www.sciencedirect.com/science/article/pii/S0092867417312412#sec4

import pandas as pd
import numpy as np

## read in data
expression_matrix = pd.read_csv('A:/merfish_mutual_information_07_31_2025/marker_genes.csv')
clusters = pd.read_csv('A:/merfish_mutual_information_07_31_2025/subclusters.csv')
genes = pd.read_csv("A:/merfish_mutual_information_07_31_2025/marker_genes_names.csv")

#binarize and transpose expression matrix
expression_matrix_binary = (expression_matrix.iloc[:, 1:] > 0).astype(int)
expression_matrix_rowgenes = expression_matrix_binary

expression_matrix_rowgenes_filtered = expression_matrix_rowgenes
from sklearn.metrics import mutual_info_score

def calc_mutual_information(binary_data, labels):
    mutual_informations = []
    
    for gene, row in binary_data.iterrows():
        print(gene)
        mi = mutual_info_score(row, labels)
        mutual_informations.append(mi)  # This line was missing!
    
    df_result = pd.DataFrame()
    df_result["gene"] = binary_data.index
    df_result.set_index("gene", inplace=True)
    df_result["mutual_information"] = mutual_informations
    df_result.sort_values("mutual_information", inplace=True, ascending=False)
    return df_result

test_mat = expression_matrix_rowgenes_filtered.iloc[:1]
result = calc_mutual_information(test_mat, clusters.iloc[:,1].values)
                        
res = calc_mutual_information(expression_matrix_rowgenes_filtered, clusters.iloc[:,1].values)
res['gene_name'] = genes.iloc[:,1] # is this right?

#res.to_csv('A:/merfish_mutual_information_07_31_2025/m_i_res.csv')


import scipy
# Calculate total information of a gene set
def calc_cumulative_information_of_set(df_discrete, genes, df_labels):
    """ Calculates total information of gene set """
    # Entropy of class without knowledge of genes
    H_naive = scipy.stats.entropy(df_labels["clusters"].value_counts(normalize=True), base=2)
    # Get discretized gene expression matrix with only selected genes
    df_temp = df_discrete.loc[genes].T

    H = 0
    # Find unique combinations of expression levels
    for _, row in df_temp.drop_duplicates().iterrows():
        print(row)
        cell_names = df_temp.index[np.all(df_temp == row, axis=1)] # Get cell names having this unique combination of expression levels
        labels_cond = df_labels.loc[cell_names]["clusters"] # Get class labels of cells
        # Calculate entropy of classification (conditional on expression levels)
        H_cond = scipy.stats.entropy(labels_cond.value_counts(normalize=True), base=2)
        weight = len(cell_names) / df_temp.shape[0] # Weight by fraction of cells in this set
        H += H_cond*weight

    I = H_naive - H
    return I
    
cum_info = calc_cumulative_information_of_set(expression_matrix_rowgenes,  expression_matrix_rowgenes.index, clusters)

### ok this is way too slow and also doesnt work, I need to reevaluate how to do this
###, I guess I need to do it on all marker genes?