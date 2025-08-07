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
        mi = mutual_info_score(row, cluster_labels) / np.log(2)  # Convert from nats to bits
        
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
            # Fixed: Get indices directly from gene names
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

# Set up data structure like the original
expression_binary = (expression_matrix.iloc[:,1:] > 0).astype(int)
expression_binary.index = genes  # Set gene names as index

# Create labels dataframe - make sure cell names match expression matrix columns
cell_names = expression_matrix.columns[1:]  # Skip first column which is gene names
df_labels = pd.DataFrame({
    'label': clusters['clusters']
}, index=cell_names)

# Debug: check that everything aligns
print("Data alignment check:")
print(f"Expression matrix shape: {expression_binary.shape}")
print(f"Number of cells in expression matrix: {expression_binary.shape[1]}")
print(f"Number of cluster labels: {len(clusters['clusters'])}")
print(f"df_labels shape: {df_labels.shape}")
print(f"Unique clusters: {df_labels['label'].nunique()}")
print(f"Do indices match? {all(expression_binary.columns == df_labels.index)}")
print(f"Sample cluster distribution:")
print(df_labels['label'].value_counts().head())

##### ok now the real deal #####
## train classifier ##

# Use the original calc_mutual_information function approach
def calc_mutual_information(df_discrete, labels):
    mutual_informations = []
    for symbol, row in df_discrete.iterrows():
        mi = mutual_info_score(row, labels) / np.log(2) # calculate mutual information, convert from nats to bits
        mutual_informations.append(mi)
    df_result = pd.DataFrame()
    df_result["symbol"] = df_discrete.index
    df_result.set_index("symbol", inplace=True)
    df_result["mutual_information"] = mutual_informations
    df_result.sort_values("mutual_information", inplace=True, ascending=False)
    return df_result

#compute mutual info using original function
df_labels['label'] =clusters['clusters'].values
df_info = calc_mutual_information(expression_binary, df_labels["label"])

#compute mean expression for each gene within each cluster
df_expr_labels = calcMeanClustExp(expression_matrix.iloc[0:,1:],  clusters)

#calculate relative cum info
#total entropy
from scipy.stats import entropy 
H_naive = entropy(clusters['clusters'].value_counts(normalize=True), base=2)
print("Total entropy of classification", H_naive, "bits")

# Use original find_nonredundant_gene_set function
def find_nonredundant_gene_set(df_discrete, genes,
                               df_labels,
                               df_info,
                               H_naive,
                               N_constrain=20,
                               cumulative_information_cutoff=0.99,
                               verbose=False):
    print('find_nonredundant_gene_set')
    cis = []
    nonredundant_genes = []
    current_ci = 0
    current_relative_ci = 0.0
    
    # Rank genes by mutual information
    remaining_genes = list(df_info.loc[genes]["mutual_information"].sort_values(ascending=False).head(n=N_constrain).index)
    
    i = 0
    
    while True:
        
        i +=1
        
        # Calculate information gain for each gene
        df_info_gains = pd.DataFrame(index=list(remaining_genes))
        my_info_gains = []
        
        for gene in tqdm(remaining_genes):
            my_genes = nonredundant_genes + [gene]
            # Use original calc_cumulative_information_of_set function
            my_ci = calc_cumulative_information_of_set(df_discrete, my_genes, df_labels)
            my_info_gain = my_ci - current_ci
            my_info_gains.append(my_info_gain)
        
        df_info_gains["info_gain"] = my_info_gains
        
        # Sort genes by information gain
        df_info_gains.sort_values("info_gain", ascending=False, inplace=True)
        
        # Take best gene
        hit = df_info_gains["info_gain"].index[0]
        
        nonredundant_genes.append(hit)
        remaining_genes.remove(hit)
        current_ci = current_ci + df_info_gains.iloc[0]["info_gain"]
        current_relative_ci = current_ci / H_naive
        
        if verbose:
            print(hit)
            print(current_relative_ci)
            print()
        
        if len(remaining_genes) == 0 or current_relative_ci > cumulative_information_cutoff:
            return nonredundant_genes

# Use original calc_cumulative_information_of_set function
def calc_cumulative_information_of_set(df_discrete, genes, df_labels):
    """ Calculates total information of gene set """
    print('calc_cumulative_information_of_set')
    import scipy.stats
    # Entropy of class without knowledge of genes
    H_naive = scipy.stats.entropy(df_labels["label"].value_counts(normalize=True), base=2)
    # Get discretized gene expression matrix with only selected genes
    df_temp = df_discrete.loc[genes].T

    H = 0
    # Find unique combinations of expression levels
    for _, row in tqdm(df_temp.drop_duplicates().iterrows()):
        cell_names = df_temp.index[np.all(df_temp == row, axis=1)] # Get cell names having this unique combination of expression levels
        labels_cond = df_labels.loc[cell_names]["label"] # Get class labels of cells
        # Calculate entropy of classification (conditional on expression levels)
        H_cond = scipy.stats.entropy(labels_cond.value_counts(normalize=True), base=2)
        weight = len(cell_names) / df_temp.shape[0] # Weight by fraction of cells in this set
        H += H_cond*weight

    I = H_naive - H
    return I

# Use original calc_cumulative_informations function
def calc_cumulative_informations(df_discrete, genes, df_labels, N=5):
    """ Calculates total information of gene sets starting from top of df """
    print('calc_cumulative_informations')
    cis = []
    for i in tqdm(range(0,N)):
        my_ci = calc_cumulative_information_of_set(df_discrete, genes[:i+1], df_labels)
        cis.append(my_ci)
    cis = np.array(cis)
    df_result = pd.DataFrame()
    df_result["symbol"] = genes
    df_result.set_index("symbol", inplace=True)
    df_result["cumulative_mutual_information"] = ""
    df_result["cumulative_mutual_information"] = np.nan
    df_result["cumulative_mutual_information"].loc[genes[:N]] = cis
    return df_result

#find non redundant gene set using original approach
genes_nonredundant = find_nonredundant_gene_set(
    expression_binary, 
    df_info.head(n=100).index,  # Top 100 genes by MI
    df_labels, 
    df_info, 
    H_naive,
    N_constrain=100,
    cumulative_information_cutoff=0.95,
    verbose=True
)
print(genes_nonredundant)

#calculate relatice cum info
df_info_nonredundant = df_info.copy()

#calculate cum info for top N genes using original approach
cumulative_informations = calc_cumulative_informations(
    expression_binary, 
    genes_nonredundant,
    df_labels, 
    N=len(genes_nonredundant)
)

df_info_nonredundant['cumulative_information'] = cumulative_informations['cumulative_mutual_information']

#calculate information relative to total entropy
relative_cumulative_informations = df_info_nonredundant["cumulative_information"] / H_naive
df_info_nonredundant["relative_cumulative_information"] = relative_cumulative_informations

#sort
df_info_nonredundant.sort_values("relative_cumulative_information", inplace=True)

#prune
genes_nonredundant_minimal = list(head_threshold(df_info_nonredundant, 0.95).index)

# Peek at minimal set
df_info_nonredundant['symbol'] = df_info_nonredundant.index
df_info_nonredundant.head(len(genes_nonredundant_minimal))['symbol' if 'symbol' in df_info_nonredundant.columns 
                                                           else df_info_nonredundant.index.name]

minimal_set = df_info_nonredundant.head(len(genes_nonredundant_minimal))['symbol' if 'symbol' in df_info_nonredundant.columns 
                                                           else df_info_nonredundant.index.name]
