##libraries
import pandas as pd
import numpy as np
from sklearn.metrics import mutual_info_score
from collections import Counter
from tqdm import tqdm

### define functions

def mutualInfo(binarized_data, cluster_labels, gene_labels):
    """
    Parameters
    ----------
    binarized_data : TYPE data frame
        DESCRIPTION. binarized gene expression matrix where
        rows are genes and columns are cells
        
    cluster_labels : TYPE series
        DESCRIPTION.
        labels for each cell, should be same dimension as binarized data
        columns
        
    gene_labels : TYPE series
        DESCRIPTION.
        labels for each gene, which is otherwise referred to as a row index
        
    Returns
    -------
    mutual_info : TYPE
        DESCRIPTION. data frame with columns for mutual info, gene, and gene
        index

    """
    print('mutualInfo')
    mutual_info = pd.DataFrame()
    for gene_index, row in tqdm(binarized_data.iterrows()):
        #print(gene_index)
        mi = mutual_info_score(row, cluster_labels) /np.log(2)
        
        # Create DataFrame with data directly
        newd = pd.DataFrame({
            'MI': [mi],
            'gene': [gene_labels[gene_index]], 
            'gene_index': [gene_index]
        })
        mutual_info = pd.concat([mutual_info, newd], ignore_index=True)
    return mutual_info


def totalInfo(binarized_data, cluster_labels, gene_indices, gene_labels):
    """

    Parameters
    ----------
  binarized_data : TYPE data frame
      DESCRIPTION. binarized gene expression matrix where
      rows are genes and columns are cells
      
  cluster_labels : TYPE series
      DESCRIPTION.
      labels for each cell, should be same dimension as binarized data
      columns
      
      gene_indices : TYPE
        DESCRIPTION.
        index for each gene
        
    gene_labels : TYPE series
        DESCRIPTION.
        labels for each gene, which is otherwise referred to as a row index

    Returns
    -------
    I : TYPE float
        DESCRIPTION.
        total info in the dataset?
    """
    from scipy.stats import entropy
    import numpy as np
    print('totalInfo')
    
    # Entropy of cluster labels alone
    naive_entropy = entropy(cluster_labels['clusters'].value_counts(normalize=True), base=2)
    
    # GEX with selected genes
    df_temp = binarized_data.loc[gene_indices].T 
    
    H = 0  # Initialize conditional entropy accumulator
    
    # Unique combinations of expression levels
    for  _, row in tqdm(df_temp.drop_duplicates().iterrows()):
        #print(row)
        # Get cell names having this unique combination of expression levels
        cell_names = df_temp.index[np.all(df_temp == row, axis=1)] 
        
        # Get labels for cells with this expression pattern
        labels_cond = cluster_labels['clusters'][cluster_labels['cell_label'].isin(cell_names)]
        
        # Calculate entropy of classification (conditional on expression levels)
        H_cond = entropy(labels_cond.value_counts(normalize=True), base=2)
        
        # Weight by fraction of cells in this set
        weight = len(cell_names) / df_temp.shape[0] 
        H += H_cond * weight
        
    # Information gain = naive entropy - conditional entropy
    I = naive_entropy - H
    return I

#calculate total information from each gene set defined by iteratively adding a single gene

def calcCumInfo(binarized_data, cluster_labels, genes, N=5):
    print('calcCumInfo')

    
    cis = []
    for i in tqdm(range(0, N)):
        print(i)
        mi_cum = totalInfo(binarized_data, cluster_labels, binarized_data.index[:i+1], genes[:i+1])
        cis.append(mi_cum)
    cis = np.array(cis)
    df_result = pd.DataFrame()
    df_result["gene"] = genes[:N]  # Fixed: only take first N genes
    df_result.set_index("gene", inplace=True)
    df_result["cumulative_mutual_information"] = cis  # Fixed: directly assign the array
    return df_result

def findNonRedundantGeneSet(binarized_data,
                            cluster_labels,
                            mi_info, 
                            H_naive,
                            N_constrain=20,
                            cumulative_information_cutoff=0.99,  # Fixed typo
                            verbose=True):
    print('findNonRedundantGeneSet')
    cis = []
    nonredundant_genes = []
    current_ci = 0
    current_relative_ci = 0.0
    
    # Fixed: Use mi_info directly since it has gene names in 'gene' column
    # Sort by MI and take top N_constrain genes
    top_genes = mi_info.nlargest(N_constrain, 'MI')
    remaining_genes = list(top_genes['gene'].values)
    
    i = 0
    while True: 
        i += 1
        
        # Calculate info gain for each gene
        df_info_gains = pd.DataFrame(index=list(remaining_genes))
        my_info_gains = []
        
        for gene in tqdm(remaining_genes):
            my_genes = nonredundant_genes + [gene]
            # Get the indices directly from the gene names
            my_indices = [genes[genes == g].index[0] for g in my_genes if g in genes.values]
            
            my_ci = totalInfo(binarized_data, cluster_labels, my_indices, my_genes)            
            my_info_gain = my_ci - current_ci
            my_info_gains.append(my_info_gain)
        
        df_info_gains["info_gain"] = my_info_gains
        
        # Sort genes by info gain
        df_info_gains.sort_values('info_gain', ascending=False, inplace=True)
        
        # Take the best gene
        hit = df_info_gains.index[0]
        
        nonredundant_genes.append(hit)
        remaining_genes.remove(hit)
        current_ci = current_ci + df_info_gains.iloc[0]['info_gain']
        current_relative_ci = current_ci / H_naive
        
        if verbose:
            print(hit)
            print(current_relative_ci)
            print()  # Fixed: added parentheses
         
        if len(remaining_genes) == 0 or current_relative_ci > cumulative_information_cutoff:
            return nonredundant_genes

def calcMeanClustExp(normalized_expression, 
                     cluster_labels,
                     summary_func = 'mean'):
    print('calcMeanClustExp')
    import copy
    df_temp = copy.deepcopy(normalized_expression.iloc[:, 1:].T).reset_index()
    df_temp['cluster'] = cluster_labels['clusters']
    df_summary_by_label = df_temp.groupby("cluster").agg(summary_func,numeric_only = True).T
    del df_temp
    return df_summary_by_label

def head_threshold(df, threshold, field="relative_cumulative_information"):
    print('head_threshold')
    L = max(np.where(df[field] <= threshold)[0])+2
    print(L)
    return df.iloc[:L]
        

## read in data
expression_matrix = pd.read_csv('A:/merfish_mutual_information_07_31_2025/marker_genes.csv')
""" here, rows are genes (unlabeled) and columns are cells,
all of the genes listed here are marker genes with adjusted p_value < 0.001
to save computation
    remove column 1 CSV index
"""
clusters = pd.read_csv('A:/merfish_mutual_information_07_31_2025/subclusters.csv')
"""clusters['clusters'] corresponds to the cluster label for each cell"""
clusters['cell_label'] = expression_matrix.iloc[:, 1:].columns

genes = expression_matrix['Unnamed: 0']
"""genes corresponds to the genes in expression matrix"""

## train classifier ##
#make expression binary
expression_binary = (expression_matrix.iloc[0:,1:] > 0).astype(int)

#compute mutual info
df_info = mutualInfo(expression_binary, clusters['clusters'], genes)

#compute mean expression for each gene within each cluster
df_expr_labels = calcMeanClustExp(expression_matrix.iloc[0:,1:],  clusters)

#calculate relative cum info
#total entropy
H_naive = entropy(clusters['clusters'].value_counts(normalize=True), base=2)
print("Total entropy of classification", H_naive, "bits")

#find non redundant gene set
genes_nonredundant = findNonRedundantGeneSet(expression_binary, clusters, df_info, H_naive)
print(genes_nonredundant)

#calculate relatice cum info
df_info_nonredundant = df_info.copy()

#calculate cum info for top N genes
cumulative_informations = calcCumInfo(expression_binary, clusters, genes_nonredundant, N = len(genes_nonredundant))
df_info_nonredundant['cumulative_information'] = cumulative_informations['cumulative_mutual_information'].reset_index(drop=True)

#calculate information relative to total entropy
relative_cumulative_informations = df_info_nonredundant["cumulative_information"] / H_naive
df_info_nonredundant["relative_cumulative_information"] = relative_cumulative_informations

#sort
df_info_nonredundant.sort_values("relative_cumulative_information", inplace=True)

#prune
genes_nonredundant_minimal = list(head_threshold(df_info_nonredundant, 0.95).index)

# Peek at minimal set
df_info_nonredundant.head(len(genes_nonredundant_minimal))['gene']

### ok these are clearly wrong, the max cum and relative cum info never hit 0.95, 
### so I think just look at like the 10 with the highest cum info

gabes_choices = df_info_nonredundant['gene'][df_info_nonredundant['cumulative_information']>0.3]
