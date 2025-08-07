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
res['gene_name'] = genes.iloc[:,1] 

#res.to_csv('A:/merfish_mutual_information_07_31_2025/m_i_res.csv')

def conditional_mutual_information(expression_matrix, cluster_labels, candidate_gene, selected_genes):
    """
    Calculate I(C; candidate_gene | selected_genes)
    """
    # Create joint states for selected genes
    if len(selected_genes) == 0:
        # If no genes selected yet, this is just regular MI
        return mutual_information(expression_matrix[candidate_gene], cluster_labels)
    
    # Combine selected genes into joint state
    selected_states = create_joint_state(expression_matrix[selected_genes])
    
    # Calculate conditional MI using the formula:
    # I(C; X | S) = H(C | S) - H(C | X, S)
    
    conditional_entropy_S = conditional_entropy(cluster_labels, selected_states)
    
    # Joint state of candidate gene + selected genes
    joint_states = create_joint_state(expression_matrix[selected_genes + [candidate_gene]])
    conditional_entropy_XS = conditional_entropy(cluster_labels, joint_states)
    
    return conditional_entropy_S - conditional_entropy_XS


from collections import Counter

def greedy_gene_selection_from_mi(expression_matrix, cluster_labels, gene_mi_scores, 
                                  target_explained_variance=0.95, max_genes=20):
    """
    Most efficient implementation starting from pre-calculated mutual information scores.
    
    Parameters:
    -----------
    expression_matrix : pd.DataFrame or np.array
        Binary matrix (genes x cells) where 1=expressed, 0=not expressed
    cluster_labels : array-like
        Cluster identity for each cell
    gene_mi_scores : dict or pd.Series
        Pre-calculated mutual information scores {gene_name: MI_score}
    target_explained_variance : float
        Stop when this fraction of total entropy is explained
    max_genes : int
        Maximum number of genes to select (safety limit)
    
    Returns:
    --------
    selected_genes : list
        Ordered list of selected gene names
    selection_info : dict
        Information about each selection step
    """
    
    # Convert inputs to efficient formats
    if isinstance(expression_matrix, pd.DataFrame):
        gene_names = expression_matrix.index.tolist()
        X = expression_matrix.values.astype(np.uint8)
    else:
        gene_names = list(range(expression_matrix.shape[0]))
        X = expression_matrix.astype(np.uint8)
    
    labels = np.array(cluster_labels)
    n_genes, n_cells = X.shape
    
    # Calculate total entropy once
    cluster_counts = Counter(labels)
    cluster_probs = np.array([cluster_counts[c] for c in cluster_counts]) / n_cells
    total_entropy = -np.sum(cluster_probs * np.log2(cluster_probs + 1e-10))
    
    # Sort genes by MI (descending) 
    if isinstance(gene_mi_scores, dict):
        sorted_genes = sorted(gene_mi_scores.items(), key=lambda x: x[1], reverse=True)
    else:
        sorted_genes = [(gene, score) for gene, score in gene_mi_scores.sort_values(ascending=False).items()]
    
    # Initialize
    selected_genes = []
    remaining_candidates = [gene for gene, _ in sorted_genes]
    explained_entropy = 0.0
    selection_info = []
    
    print(f"Total entropy: {total_entropy:.4f}")
    print(f"Starting greedy selection...")
    
    # Select first gene (highest MI)
    first_gene = remaining_candidates[0]
    first_gene_idx = gene_names.index(first_gene) if first_gene in gene_names else first_gene
    selected_genes.append(first_gene)
    remaining_candidates.remove(first_gene)
    
    first_mi = gene_mi_scores[first_gene] if isinstance(gene_mi_scores, dict) else gene_mi_scores[first_gene]
    explained_entropy = first_mi / total_entropy
    
    selection_info.append({
        'gene': first_gene,
        'mi_score': first_mi,
        'cumulative_explained': explained_entropy
    })
    
    print(f"Gene 1: {first_gene} (MI: {first_mi:.4f}, Explained: {explained_entropy:.1%})")
    
    # Iterative selection with conditional MI
    for iteration in range(2, max_genes + 1):
        if explained_entropy >= target_explained_variance:
            print(f"Target explained variance {target_explained_variance:.1%} reached!")
            break
            
        if not remaining_candidates:
            print("No more candidate genes!")
            break
        
        # Get indices for selected genes
        selected_indices = []
        for gene in selected_genes:
            if gene in gene_names:
                selected_indices.append(gene_names.index(gene))
            else:
                selected_indices.append(gene)
        
        # Create joint state for selected genes (bit encoding for efficiency)
        joint_state = encode_joint_state_fast(X[selected_indices])
        
        # Find best remaining gene using conditional MI
        best_cmi = -1
        best_gene = None
        
        # Speed optimization: only check top candidates after a few iterations
        candidates_to_check = remaining_candidates[:min(30, len(remaining_candidates))] if iteration > 4 else remaining_candidates
        
        for candidate in candidates_to_check:
            candidate_idx = gene_names.index(candidate) if candidate in gene_names else candidate
            
            # Calculate conditional MI efficiently
            cmi = conditional_mi_fast(labels, X[candidate_idx], joint_state)
            
            if cmi > best_cmi:
                best_cmi = cmi
                best_gene = candidate
        
        if best_gene is None or best_cmi <= 1e-6:  # Numerical precision threshold
            print(f"No more informative genes (best CMI: {best_cmi:.6f})")
            break
        
        # Add best gene
        selected_genes.append(best_gene)
        remaining_candidates.remove(best_gene)
        explained_entropy += best_cmi / total_entropy
        
        selection_info.append({
            'gene': best_gene,
            'cmi_score': best_cmi,
            'cumulative_explained': explained_entropy
        })
        
        print(f"Gene {iteration}: {best_gene} (CMI: {best_cmi:.4f}, Total: {explained_entropy:.1%})")
    
    print(f"\nFinal selection: {len(selected_genes)} genes")
    print(f"Total explained: {explained_entropy:.1%}")
    
    return selected_genes, selection_info


def encode_joint_state_fast(gene_matrix):
    """Ultra-fast joint state encoding using bit operations."""
    if gene_matrix.shape[0] == 0:
        return np.array([])
    if gene_matrix.shape[0] == 1:
        return gene_matrix[0]
    
    # Use powers of 2 for binary encoding
    joint_state = np.zeros(gene_matrix.shape[1], dtype=np.int32)
    for i in range(min(gene_matrix.shape[0], 20)):  # Limit to prevent overflow
        joint_state += gene_matrix[i] * (2 ** i)
    
    return joint_state


def conditional_mi_fast(cluster_labels, candidate_gene, joint_selected_state):
    """
    Super fast conditional MI calculation.
    I(C; X | S) = H(C|S) - H(C|X,S)
    """
    unique_states, state_indices = np.unique(joint_selected_state, return_inverse=True)
    n_total = len(cluster_labels)
    
    total_cmi = 0.0
    
    for state_idx, state_val in enumerate(unique_states):
        # Get cells in this state
        mask = state_indices == state_idx
        n_state = np.sum(mask)
        
        if n_state < 2:  # Skip states with too few samples
            continue
        
        state_weight = n_state / n_total
        state_labels = cluster_labels[mask]
        state_gene = candidate_gene[mask]
        
        # H(C|S=state)
        h_c_s = entropy_fast(state_labels)
        
        # H(C|X,S=state) - split by gene expression
        gene_on = state_gene == 1
        gene_off = state_gene == 0
        
        h_c_xs = 0.0
        if np.sum(gene_on) > 0:
            p_on = np.sum(gene_on) / n_state
            h_c_xs += p_on * entropy_fast(state_labels[gene_on])
        
        if np.sum(gene_off) > 0:
            p_off = np.sum(gene_off) / n_state
            h_c_xs += p_off * entropy_fast(state_labels[gene_off])
        
        # Add weighted contribution
        total_cmi += state_weight * (h_c_s - h_c_xs)
    
    return max(0, total_cmi)


def entropy_fast(labels):
    """Fastest entropy calculation."""
    if len(labels) <= 1:
        return 0.0
    
    _, counts = np.unique(labels, return_counts=True)
    if len(counts) <= 1:
        return 0.0
    
    probs = counts / len(labels)
    return -np.sum(probs * np.log2(probs))


# Simple usage example:
def run_selection(expression_df, cluster_labels, mi_scores):
    """
    Simple wrapper for your use case.
    
    Parameters:
    -----------
    expression_df : pd.DataFrame
        Genes as rows, cells as columns, binary values
    cluster_labels : list/array
        Cluster assignment for each cell
    mi_scores : dict or pd.Series
        {gene_name: mutual_information_score}
    """
    
    selected_genes, info = greedy_gene_selection_from_mi(
        expression_df, 
        cluster_labels, 
        mi_scores,
        target_explained_variance=0.95,
        max_genes=15  # Reasonable limit for MERFISH
    )
    
    return selected_genes, info


# Example usage:
if __name__ == "__main__":
    # Mock data for testing
    np.random.seed(42)
    
    # Create test data
    genes = [f"Gene_{i}" for i in range(50)]
    cells = [f"Cell_{i}" for i in range(500)]
    
    # Binary expression matrix
    expression_df = pd.DataFrame(
        np.random.binomial(1, 0.4, (50, 500)),
        index=genes,
        columns=cells
    )
    
    # Cluster labels
    clusters = np.random.choice(['A', 'B', 'C', 'D'], 500)
    
    # Mock MI scores (normally you'd have calculated these)
    mi_scores = {gene: np.random.exponential(0.1) for gene in genes}
    
    # Run selection
    selected, info = run_selection(expression_df, clusters, mi_scores)
    print("Selected genes:", selected)