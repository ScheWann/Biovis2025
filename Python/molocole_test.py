import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from scipy.sparse import coo_matrix
import warnings
import scipy.sparse
from statsmodels.stats.multitest import multipletests
import os
from MOLOCLE import (
    graph_test, 
    moran_i_test, 
    calculate_lw, 
    analyze_differential_expression,
    load_and_preprocess_data,
    preprocess_and_cluster,
    map_cell_types,
    find_gene_modules,
    perform_trajectory_analysis,
    plot_results,
    plot_genes_in_pseudotime,
    aggregate_gene_expression
)

# Create output directory if it doesn't exist
output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output')
os.makedirs(output_dir, exist_ok=True)

def create_dummy_data(n_cells=1000, n_genes=2000):
    """
    Create more realistic dummy data with clear cell type separation and gene expression patterns.
    
    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
        
    Returns:
    --------
    adata : AnnData
        Annotated data matrix
    """
    # Create cell types (3 distinct clusters)
    cell_types = np.repeat(['type_A', 'type_B', 'type_C'], [n_cells//3, n_cells//3, n_cells - 2*(n_cells//3)])
    np.random.shuffle(cell_types)
    
    # Create gene expression patterns
    X = np.zeros((n_cells, n_genes))
    
    # Type A specific genes (first 100 genes)
    X[cell_types == 'type_A', :100] = np.random.poisson(lam=5, size=(np.sum(cell_types == 'type_A'), 100))
    X[cell_types != 'type_A', :100] = np.random.poisson(lam=0.1, size=(np.sum(cell_types != 'type_A'), 100))
    
    # Type B specific genes (next 100 genes)
    X[cell_types == 'type_B', 100:200] = np.random.poisson(lam=5, size=(np.sum(cell_types == 'type_B'), 100))
    X[cell_types != 'type_B', 100:200] = np.random.poisson(lam=0.1, size=(np.sum(cell_types != 'type_B'), 100))
    
    # Type C specific genes (next 100 genes)
    X[cell_types == 'type_C', 200:300] = np.random.poisson(lam=5, size=(np.sum(cell_types == 'type_C'), 100))
    X[cell_types != 'type_C', 200:300] = np.random.poisson(lam=0.1, size=(np.sum(cell_types != 'type_C'), 100))
    
    # Common genes with varying expression (remaining genes)
    X[:, 300:] = np.random.poisson(lam=1, size=(n_cells, n_genes-300))
    
    # Add some noise
    noise = np.random.normal(0, 0.1, size=X.shape)
    X = X + noise
    X[X < 0] = 0
    
    # Create AnnData object
    adata = sc.AnnData(X=X)
    adata.var_names = [f'gene_{i}' for i in range(n_genes)]
    adata.obs_names = [f'cell_{i}' for i in range(n_cells)]
    adata.obs['cell_type'] = cell_types
    
    return adata

def test_molocle():
    """Test the MOLOCLE pipeline with dummy data."""
    print("Testing MOLOCLE with dummy data...")
    
    # Create dummy data
    adata = create_dummy_data(n_cells=1000, n_genes=2000)
    
    # Basic preprocessing
    print("Preprocessing and clustering...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)  # Reduced number of components
    sc.pp.neighbors(adata, n_neighbors=15)  # Increased number of neighbors
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)  # Adjusted resolution
    
    # Map cell types
    print("Mapping cell types...")
    adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
    
    # Analyze differential expression
    print("Analyzing differential expression...")
    try:
        deg_results = graph_test(adata, neighbor_graph="knn", cores=1)
        deg_genes = deg_results[deg_results['q_value'] < 0.05]
        adata.uns['deg_results'] = deg_results
        adata.uns['deg_genes'] = deg_genes
    except Exception as e:
        print(f"Warning: Error in differential expression analysis: {str(e)}")
    
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Find gene modules
    print("Finding gene modules...")
    try:
        sc.tl.score_genes(adata, gene_list=adata.var_names[:10], score_name='module_score')
    except Exception as e:
        print(f"Warning: Error in gene module analysis: {str(e)}")
    
    # Perform trajectory analysis
    print("Performing trajectory analysis...")
    try:
        # Initialize the PAGA graph
        sc.tl.paga(adata, groups='leiden')
        
        # Optional: Use PAGA initialization for UMAP visualization
        sc.pl.paga(adata, plot=False)  # This is needed for PAGA initialization
        sc.tl.umap(adata, init_pos='paga')
        
        # Compute diffusion map and pseudotime
        sc.tl.diffmap(adata)
        adata.uns['iroot'] = np.argmin(adata.obsm['X_diffmap'][:, 0])
        sc.tl.dpt(adata)
    except Exception as e:
        print(f"Warning: Error in trajectory analysis: {str(e)}")
    
    # Plot results
    print("Plotting results...")
    try:
        # Plot UMAP with cell types
        sc.pl.umap(adata, color='cell_type', title='Cell Types', show=False)
        plt.savefig(os.path.join(output_dir, 'test_umap_cell_types.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot differentially expressed genes
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
        plt.savefig(os.path.join(output_dir, 'test_deg_genes.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot UMAP with module score
        if 'module_score' in adata.obs:
            sc.pl.umap(adata, color='module_score', title='Gene Module Score', show=False)
            plt.savefig(os.path.join(output_dir, 'test_umap_module_score.png'), dpi=300, bbox_inches='tight')
            plt.close()
        
        # Plot PAGA
        if 'paga' in adata.uns:
            sc.pl.paga(adata, color='cell_type', show=False)
            plt.savefig(os.path.join(output_dir, 'test_paga.png'), dpi=300, bbox_inches='tight')
            plt.close()
        
        # Plot diffusion map
        if 'X_diffmap' in adata.obsm:
            sc.pl.diffmap(adata, color='cell_type', show=False)
            plt.savefig(os.path.join(output_dir, 'test_diffmap.png'), dpi=300, bbox_inches='tight')
            plt.close()
        
        print(f"Test plots saved in: {output_dir}")
    except Exception as e:
        print(f"Warning: Error in plotting: {str(e)}")
    
    return adata

def test_graph_test():
    """
    Test the graph_test function specifically.
    """
    print("Testing graph_test function...")
    
    # Create dummy data
    adata = create_dummy_data(n_cells=200, n_genes=500)
    
    # Preprocess and cluster
    adata = preprocess_and_cluster(adata)
    
    # Test graph_test with KNN
    print("Testing graph_test with KNN...")
    results_knn = graph_test(adata, neighbor_graph='knn', k=10, verbose=True)
    print(f"Found {sum(results_knn['q_value'] < 0.05)} differentially expressed genes with KNN")
    
    # Test graph_test with principal graph
    print("Testing graph_test with principal graph...")
    # First, we need to create a principal graph
    sc.tl.draw_graph(adata)
    results_pg = graph_test(adata, neighbor_graph='principal_graph', k=10, verbose=True)
    print(f"Found {sum(results_pg['q_value'] < 0.05)} differentially expressed genes with principal graph")
    
    return results_knn, results_pg

def test_moran_i():
    """
    Test the moran_i_test function specifically.
    """
    print("Testing moran_i_test function...")
    
    # Create dummy data
    n = 100
    x = np.random.normal(0, 1, size=n)
    
    # Create a simple spatial weights matrix
    weights = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                weights[i, j] = 1 / (abs(i - j) + 1)
    
    # Test moran_i_test
    result = moran_i_test(x, weights, alternative='two.sided')
    print(f"Moran's I: {result['moran_i']:.4f}")
    print(f"P-value: {result['p_value']:.4f}")
    
    return result

def test_plot_genes_in_pseudotime():
    """
    Test the plot_genes_in_pseudotime function.
    """
    # Create dummy data
    adata = create_dummy_data(n_cells=200, n_genes=500)
    
    # Create pseudotime values that increase with time bins
    time_bins = ['270-330', '330-390', '390-450', '450-510', '510-580', '580-650', '>650']
    n_cells_per_bin = 200 // len(time_bins)
    time_categories = []
    pseudotime_values = []
    
    for i, time_bin in enumerate(time_bins):
        n_cells_in_bin = n_cells_per_bin if i < len(time_bins)-1 else 200 - (n_cells_per_bin * (len(time_bins)-1))
        time_categories.extend([time_bin] * n_cells_in_bin)
        # Create increasing pseudotime values for each bin
        start_time = i * 4
        end_time = (i + 1) * 4
        pseudotime_values.extend(np.linspace(start_time, end_time, n_cells_in_bin))
    
    # Convert to numpy array
    pseudotime_values = np.array(pseudotime_values)
    
    adata.obs['embryo.time.bin'] = pd.Categorical(time_categories)
    adata.obs['dpt_pseudotime'] = pseudotime_values
    
    # Get the first three genes
    gene_list = adata.var_names[:3].tolist()
    
    # Create expression patterns that peak at different times
    for i, gene in enumerate(gene_list):
        peak_time = (i + 2) * 8  # Different peak times for each gene
        expr = np.exp(-(pseudotime_values - peak_time)**2 / 20)  # Gaussian curve
        expr = expr * np.random.uniform(0.8, 1.2, size=len(expr))  # Add some noise
        adata.X[:, i] = expr * 30  # Scale up the expression values
    
    # Plot genes in pseudotime
    plot_genes_in_pseudotime(adata, gene_list)
    
    print("plot_genes_in_pseudotime test completed successfully")

def test_aggregate_gene_expression():
    """
    Test the aggregate_gene_expression function.
    """
    # Create dummy data
    adata = create_dummy_data(n_cells=200, n_genes=500)
    
    # Create gene groups
    gene_groups = pd.DataFrame({
        'gene_id': adata.var_names,
        'group': [f'group_{i//10}' for i in range(len(adata.var_names))]
    })
    
    # Create cell groups based on cell types
    cell_groups = pd.DataFrame({
        'cell_id': adata.obs_names,
        'group': adata.obs['cell_type']
    })
    
    # Test aggregation with different parameters
    print("Testing aggregate_gene_expression with default parameters...")
    agg_mat = aggregate_gene_expression(adata, gene_groups, cell_groups)
    print(f"Result shape: {agg_mat.shape}")
    print("First few rows and columns:")
    print(agg_mat.iloc[:5, :5])
    
    print("\nTesting with log normalization...")
    agg_mat_log = aggregate_gene_expression(adata, gene_groups, cell_groups, 
                                          norm_method='log', pseudocount=1)
    print(f"Result shape: {agg_mat_log.shape}")
    
    print("\nTesting with mean aggregation...")
    agg_mat_mean = aggregate_gene_expression(adata, gene_groups, cell_groups, 
                                           gene_agg_fun='mean', cell_agg_fun='mean')
    print(f"Result shape: {agg_mat_mean.shape}")
    
    print("aggregate_gene_expression test completed successfully")
    return agg_mat

def test_find_gene_modules():
    """
    Test the find_gene_modules function.
    """
    # Create dummy data
    adata = create_dummy_data(n_cells=200, n_genes=500)
    
    # Basic preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)
    
    # Test with default parameters
    print("Testing find_gene_modules with default parameters...")
    gene_modules = find_gene_modules(adata)
    print(f"Found {len(gene_modules['module'].unique())} gene modules")
    print("First few rows:")
    print(gene_modules.head())
    
    # Test with different resolution
    print("\nTesting with different resolution...")
    gene_modules_res = find_gene_modules(adata, resolution=0.5)
    print(f"Found {len(gene_modules_res['module'].unique())} gene modules")
    
    # Test with weight=True
    print("\nTesting with weight=True...")
    gene_modules_weight = find_gene_modules(adata, weight=True)
    print(f"Found {len(gene_modules_weight['module'].unique())} gene modules")
    
    print("find_gene_modules test completed successfully")
    return gene_modules

if __name__ == "__main__":
    # Create output directory if it doesn't exist
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Run MOLOCLE test
    adata = test_molocle()
    
    # Test graph_test function
    print("\nTesting graph_test function...")
    # Create dummy data
    test_adata = create_dummy_data(n_cells=500, n_genes=1000)
    
    # Basic preprocessing
    sc.pp.filter_cells(test_adata, min_genes=200)
    sc.pp.filter_genes(test_adata, min_cells=3)
    sc.pp.normalize_total(test_adata, target_sum=1e4)
    sc.pp.log1p(test_adata)
    sc.pp.highly_variable_genes(test_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    test_adata = test_adata[:, test_adata.var.highly_variable]
    sc.pp.scale(test_adata, max_value=10)
    sc.tl.pca(test_adata, n_comps=30)
    sc.pp.neighbors(test_adata, n_neighbors=15)
    sc.tl.umap(test_adata)
    
    # Test with KNN
    print("Testing with KNN...")
    try:
        results_knn = graph_test(test_adata, neighbor_graph="knn", cores=1)
        print(f"Found {len(results_knn[results_knn['q_value'] < 0.05])} differentially expressed genes")
    except Exception as e:
        print(f"Error with KNN: {str(e)}")
    
    # Test with principal graph
    print("\nTesting with principal graph...")
    try:
        results_pr = graph_test(test_adata, neighbor_graph="principal_graph", cores=1)
        print(f"Found {len(results_pr[results_pr['q_value'] < 0.05])} differentially expressed genes")
    except Exception as e:
        print(f"Error with principal graph: {str(e)}")
    
    # Test moran_i_test function
    print("\nTesting moran_i_test function...")
    x = np.random.normal(0, 1, 100)
    weights = np.random.rand(100, 100)
    weights = (weights + weights.T) / 2  # Make symmetric
    np.fill_diagonal(weights, 0)  # Zero diagonal
    
    try:
        result = moran_i_test(x, weights)
        print(f"Moran's I: {result['moran_i']:.4f}")
        print(f"P-value: {result['p_value']:.4f}")
    except Exception as e:
        print(f"Error in Moran's I test: {str(e)}")
    
    # Test calculate_lw function
    print("\nTesting calculate_lw function...")
    try:
        weights = calculate_lw(test_adata, k=25)
        print("Successfully calculated spatial weights matrix")
    except Exception as e:
        print(f"Error calculating spatial weights matrix: {str(e)}")
    
    # Test plot_genes_in_pseudotime function
    test_plot_genes_in_pseudotime()
    
    # Test aggregate_gene_expression function
    print("\nTesting aggregate_gene_expression function...")
    test_aggregate_gene_expression()
    
    # Test find_gene_modules function
    print("\nTesting find_gene_modules function...")
    test_find_gene_modules()
    
    print("\nAll tests completed!") 