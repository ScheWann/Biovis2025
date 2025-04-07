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
from statsmodels.nonparametric.smoothers_lowess import lowess

# R code:
# Make the CDS object
# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)
# cds <- preprocess_cds(cds, num_dim = 100)
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds, resolution=1e-5)
# 
# colData(cds)$assigned_cell_type <- as.character(partitions(cds))
# colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
#                                                  "1"="Body wall muscle",
#                                                  "2"="Germline",
#                                                  "3"="Motor neurons",
#                                                  "4"="Seam cells",
#                                                  "5"="Sex myoblasts",
#                                                  "6"="Socket cells",
#                                                  "7"="Marginal_cell",
#                                                  "8"="Coelomocyte",
#                                                  "9"="Am/PH sheath cells",
#                                                  "10"="Ciliated neurons",
#                                                  "11"="Intestinal/rectal muscle",
#                                                  "12"="Excretory gland",
#                                                  "13"="Chemosensory neurons",
#                                                  "14"="Interneurons",
#                                                  "15"="Unclassified eurons",
#                                                  "16"="Ciliated neurons",
#                                                  "17"="Pharyngeal gland cells",
#                                                  "18"="Unclassified neurons",
#                                                  "19"="Chemosensory neurons",
#                                                  "20"="Ciliated neurons",
#                                                  "21"="Ciliated neurons",
#                                                  "22"="Inner labial neuron",
#                                                  "23"="Ciliated neurons",
#                                                  "24"="Ciliated neurons",
#                                                  "25"="Ciliated neurons",
#                                                  "26"="Hypodermal cells",
#                                                  "27"="Mesodermal cells",
#                                                  "28"="Motor neurons",
#                                                  "29"="Pharyngeal gland cells",
#                                                  "30"="Ciliated neurons",
#                                                  "31"="Excretory cells",
#                                                  "32"="Amphid neuron",
#                                                  "33"="Pharyngeal muscle")
# 
# 
# 
# 
# neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
# plot_cells(neurons_cds, color_cells_by="partition")
# 
# 
# 
# pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=8)
# pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
# 
# 
# gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)
# 
# 
# cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)), 
#                                 cell_group=partitions(cds)[colnames(neurons_cds)])
# agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
# row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
# 
# pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#                    scale="column", clustering_method="ward.D2",
#                    fontsize=6)
# 
# 
# 
# plot_cells(neurons_cds, 
#            genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
#            group_cells_by="partition",
#            color_cells_by="partition",
#            show_trajectory_graph=FALSE)
# 
# expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
# cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
# gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))
# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)
# cds <- preprocess_cds(cds, num_dim = 50)
# cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds)
# cds <- learn_graph(cds)
# cds <- order_cells(cds)
# plot_cells(cds,
#            color_cells_by = "cell.type",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)
# 
# 
# ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
# pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
# 
# 
# plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
#            show_trajectory_graph=FALSE,
#            label_cell_groups=FALSE,
#            label_leaves=FALSE)
# 
# 
# 
# gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
# 
# 
# 
# 
# cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
#                                 cell_group=colData(cds)$cell.type)
# agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
# row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# pheatmap::pheatmap(agg_mat,
#                    scale="column", clustering_method="ward.D2")
# 
# 
#                    plot_cells(cds,
#            genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)
# 
# 
#            AFD_genes <- c("gcy-8", "dac-1", "oig-8")
# AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
#                        colData(cds)$cell.type %in% c("AFD")]
# AFD_lineage_cds <- order_cells(AFD_lineage_cds)
# 
# 
# plot_genes_in_pseudotime(AFD_lineage_cds,
#                          color_cells_by="embryo.time.bin",
#                          min_expr=0.5)

# R code:
# graph_test <- function(cds,
#                        neighbor_graph = c("knn", "principal_graph"),
#                        reduction_method = "UMAP",
#                        k = 25,
#                        method = c('Moran_I'),
#                        alternative = 'greater',
#                        expression_family="quasipoisson",
#                        cores=1,
#                        verbose=FALSE,
#                        nn_control=list()) {
#   ...
# }

def graph_test(adata, neighbor_graph='knn', reduction_method='umap', k=25, 
               method='moran_i', alternative='greater', expression_family='quasipoisson', 
               cores=1, verbose=False):
    """
    Test genes for differential expression based on the low dimensional embedding and the principal graph.
    
    This function uses Moran's I test to identify genes that are differentially expressed along a trajectory.
    Moran's I is a measure of spatial autocorrelation that tells you whether cells at nearby positions
    on a trajectory will have similar (or dissimilar) expression levels for the gene being tested.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    neighbor_graph : str, optional (default: 'knn')
        Type of neighbor graph to use ('knn' or 'principal_graph')
    reduction_method : str, optional (default: 'umap')
        The method used for dimensionality reduction
    k : int, optional (default: 25)
        Number of nearest neighbors
    method : str, optional (default: 'moran_i')
        The method to use for testing ('moran_i' or 'geary_c')
    alternative : str, optional (default: 'greater')
        The alternative hypothesis, one of 'greater', 'less', or 'two.sided'
    expression_family : str, optional (default: 'quasipoisson')
        The expression family to use
    cores : int, optional (default: 1)
        Number of cores to use for parallel processing
    verbose : bool, optional (default: False)
        Whether to print progress messages
    
    Returns:
    --------
    DataFrame
        Results of the differential expression testing
    """
    # Check if dimensionality reduction exists
    if reduction_method.lower() == 'umap':
        if 'X_umap' not in adata.obsm:
            raise ValueError(f"No dimensionality reduction for {reduction_method} calculated. "
                           f"Please run sc.tl.umap before running graph_test.")
    
    # Check if principal graph exists when needed
    if neighbor_graph == 'principal_graph':
        if 'X_draw_graph_fa' not in adata.obsm:
            # Try to create the principal graph
            try:
                sc.tl.draw_graph(adata)
            except Exception as e:
                raise ValueError(f"Could not create principal graph: {str(e)}")
    
    # Calculate spatial weights matrix
    if verbose:
        print("Calculating spatial weights matrix...")
    weights = calculate_lw(adata, k=k, neighbor_graph=neighbor_graph, reduction_method=reduction_method)
    
    # Perform Moran's I test
    if verbose:
        print("Performing Moran's I test...")
    
    # Get expression matrix and size factors
    if expression_family == 'quasipoisson':
        # Log transform the data
        X = np.log1p(adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X)
    else:
        X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    
    # Initialize results
    results = {
        'gene_id': [],
        'p_value': [],
        'moran_i': [],
        'statistic': [],
        'status': []
    }
    
    # For each gene
    for i in range(X.shape[1]):
        gene_expr = X[:, i]
        
        # Perform Moran's I test
        try:
            test_result = moran_i_test(gene_expr, weights, alternative=alternative)
            
            # Add results
            results['gene_id'].append(adata.var_names[i])
            results['p_value'].append(test_result['p_value'])
            results['moran_i'].append(test_result['moran_i'])
            results['statistic'].append(test_result['statistic'])
            results['status'].append('OK')
        except Exception as e:
            if verbose:
                print(f"Error testing gene {adata.var_names[i]}: {str(e)}")
            
            # Add NA results
            results['gene_id'].append(adata.var_names[i])
            results['p_value'].append(np.nan)
            results['moran_i'].append(np.nan)
            results['statistic'].append(np.nan)
            results['status'].append('FAIL')
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Calculate q-values (Benjamini-Hochberg correction)
    results_df['q_value'] = 1.0
    ok_mask = results_df['status'] == 'OK'
    if any(ok_mask):
        results_df.loc[ok_mask, 'q_value'] = multipletests(
            results_df.loc[ok_mask, 'p_value'], 
            method='fdr_bh'
        )[1]
    
    # Set gene_id as index
    results_df.set_index('gene_id', inplace=True)
    
    # Ensure gene order matches input
    results_df = results_df.reindex(adata.var_names)
    
    return results_df

# R code:
# my.moran.test <- function (x, listw, wc, alternative = "greater",
#                            randomisation = TRUE) {
#   ...
# }

def moran_i_test(x, weights, alternative='greater', randomization=True):
    """
    Perform Moran's I test for spatial autocorrelation.
    
    Parameters:
    -----------
    x : array-like
        The variable to test
    weights : array-like
        The spatial weights matrix
    alternative : str, optional (default: 'greater')
        The alternative hypothesis, one of 'greater', 'less', or 'two.sided'
    randomization : bool, optional (default: True)
        Whether to use randomization assumption for variance calculation
    
    Returns:
    --------
    dict
        Results of the Moran's I test
    """
    # Convert inputs to numpy arrays
    x = np.array(x)
    weights = np.array(weights)
    
    # Calculate basic statistics
    n = len(x)
    x_centered = x - np.mean(x)
    s0 = np.sum(weights)
    s1 = 0.5 * np.sum((weights + weights.T) ** 2)
    s2 = np.sum((np.sum(weights, axis=0) + np.sum(weights, axis=1)) ** 2)
    
    # Calculate Moran's I
    numerator = np.sum(weights * np.outer(x_centered, x_centered))
    denominator = s0 * np.sum(x_centered ** 2)
    i_stat = (n / (2 * s0)) * (numerator / denominator)
    
    # Calculate expected value
    ei = -1 / (n - 1)
    
    # Calculate variance
    if randomization:
        # Calculate kurtosis
        k = stats.kurtosis(x)
        
        # Calculate variance under randomization assumption
        vi = (n * ((n**2 - 3*n + 3) * s1 - n * s2 + 3 * s0**2) - 
              (k * ((n**2 - n) * s1 - 2 * n * s2 + 6 * s0**2))) / ((n-1) * (n-2) * (n-3) * s0**2) - ei**2
    else:
        # Calculate variance under normality assumption
        vi = (n * s1 - n * s2 + 3 * s0**2) / ((n-1) * s0**2) - ei**2
    
    # Calculate z-score
    zi = (i_stat - ei) / np.sqrt(vi)
    
    # Calculate p-value
    if alternative == 'two.sided':
        p_value = 2 * (1 - stats.norm.cdf(abs(zi)))
    elif alternative == 'greater':
        p_value = 1 - stats.norm.cdf(zi)
    else:  # 'less'
        p_value = stats.norm.cdf(zi)
    
    return {
        'statistic': zi,
        'p_value': p_value,
        'moran_i': i_stat,
        'expected_i': ei,
        'variance_i': vi
    }

# R code:
# calculateLW <- function(cds,
#                         k,
#                         neighbor_graph,
#                         reduction_method,
#                         verbose = FALSE,
#                         nn_control = list()) {
#   ...
# }

def calculate_lw(adata, k=25, neighbor_graph='knn', reduction_method='umap', verbose=False):
    """
    Calculate the spatial weights matrix for Moran's I test.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    k : int, optional (default: 25)
        Number of nearest neighbors
    neighbor_graph : str, optional (default: 'knn')
        Type of neighbor graph to use ('knn' or 'principal_graph')
    reduction_method : str, optional (default: 'umap')
        The method used for dimensionality reduction
    verbose : bool, optional (default: False)
        Whether to print progress messages
    
    Returns:
    --------
    array-like
        The spatial weights matrix
    """
    # Get the low-dimensional embedding
    if reduction_method.lower() == 'umap':
        if 'X_umap' not in adata.obsm:
            raise ValueError("UMAP embedding not found. Please run sc.tl.umap first.")
        coords = adata.obsm['X_umap']
    else:
        raise ValueError(f"Reduction method {reduction_method} not supported.")
    
    # Calculate k-nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto').fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    
    # Create adjacency matrix
    n_cells = coords.shape[0]
    adj_matrix = np.zeros((n_cells, n_cells))
    
    for i in range(n_cells):
        for j in range(1, k+1):  # Skip self (j=0)
            adj_matrix[i, indices[i, j]] = 1
            adj_matrix[indices[i, j], i] = 1  # Make symmetric
    
    # If using principal graph, modify the adjacency matrix
    if neighbor_graph == 'principal_graph':
        # Run draw_graph if not already run
        if 'X_draw_graph_fa' not in adata.obsm:
            sc.tl.draw_graph(adata)
        
        # Create principal graph from force-directed layout
        G = nx.Graph()
        coords = adata.obsm['X_draw_graph_fa']
        
        # Add edges based on k-nearest neighbors in force-directed layout
        nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto').fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        
        for i in range(n_cells):
            for j in range(1, k+1):
                G.add_edge(i, indices[i, j])
        
        # Modify adjacency matrix based on principal graph
        for i in range(n_cells):
            for j in range(n_cells):
                if adj_matrix[i, j] == 1:
                    if not G.has_edge(i, j):
                        adj_matrix[i, j] = 0
    
    return adj_matrix

# R code:
# analyze_differential_expression <- function(cds) {
#   ...
# }

def analyze_differential_expression(adata):
    """
    Analyze differential expression using graph-based testing.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    AnnData
        Updated annotated data matrix with differential expression results
    """
    # Perform graph-based differential expression testing
    deg_results = graph_test(adata, neighbor_graph="knn", cores=1)
    
    # Filter for significant genes
    deg_genes = deg_results[deg_results['q_value'] < 0.05]
    
    # Store results in adata
    adata.uns['deg_results'] = deg_results
    adata.uns['deg_genes'] = deg_genes
    
    # Also perform standard differential expression analysis
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    return adata

# R code:
# load_and_preprocess_data <- function() {
#   ...
# }

def load_and_preprocess_data():
    """
    Load and preprocess the data.
    
    Returns:
    --------
    AnnData
        Preprocessed annotated data matrix
    """
    # Load sample data
    adata = sc.datasets.pbmc3k()
    
    # Replace infinite values with NaN and then with 0
    if scipy.sparse.issparse(adata.X):
        adata.X = csr_matrix(np.nan_to_num(adata.X.toarray(), nan=0.0, posinf=0.0, neginf=0.0))
    else:
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    return adata

# R code:
# preprocess_and_cluster <- function(adata) {
#   ...
# }

def preprocess_and_cluster(adata):
    """
    Preprocess and cluster the data.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    AnnData
        Processed and clustered annotated data matrix
    """
    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    
    # Scale the data
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    sc.tl.pca(adata, n_comps=50)
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Run UMAP
    sc.tl.umap(adata)
    
    # Run Leiden clustering
    sc.tl.leiden(adata, resolution=0.5)
    
    return adata

# R code:
# map_cell_types <- function(adata) {
#   ...
# }

def map_cell_types(adata):
    """
    Map cell types based on clustering results.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    AnnData
        Annotated data matrix with cell type annotations
    """
    # This is a simplified version - in a real analysis, you would use marker genes
    # to identify cell types based on the clusters
    
    # For demonstration, we'll just use the cluster labels as cell types
    adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
    
    return adata

# R code:
# find_gene_modules <- function(adata) {
#   ...
# }

def find_gene_modules(adata):
    """
    Find gene modules using correlation analysis.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    AnnData
        Annotated data matrix with gene module information
    """
    # This is a simplified version - in a real analysis, you would use more sophisticated
    # methods to identify gene modules
    
    # For demonstration, we'll just use correlation analysis
    sc.tl.score_genes(adata, gene_list=adata.var_names[:10], score_name='module_score')
    
    return adata

# R code:
# perform_trajectory_analysis <- function(adata) {
#   ...
# }

def perform_trajectory_analysis(adata):
    """
    Perform trajectory analysis using PAGA and diffusion maps.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    
    Returns:
    --------
    AnnData
        Annotated data matrix with trajectory information
    """
    # Run PAGA
    sc.tl.paga(adata)
    
    # Run diffusion map
    sc.tl.diffmap(adata)
    
    # Set root cell as the cell with minimum diffusion pseudotime component
    adata.uns['iroot'] = np.argmin(adata.obsm['X_diffmap'][:, 0])
    
    # Calculate pseudotime
    sc.tl.dpt(adata)
    
    return adata

# R code:
# plot_results <- function(adata) {
#   ...
# }

def plot_genes_in_pseudotime(adata, gene_list, color_by='embryo.time.bin', ncol=1, figsize=(8, 12)):
    """
    Plot gene expression along pseudotime.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    gene_list : list
        List of genes to plot
    color_by : str, optional (default: 'embryo.time.bin')
        Column name in adata.obs to color cells by
    ncol : int, optional (default: 1)
        Number of columns in the plot grid
    figsize : tuple, optional (default: (8, 12))
        Figure size
    """
    # Check if pseudotime exists
    if 'dpt_pseudotime' not in adata.obs:
        raise ValueError("No pseudotime found. Please run trajectory analysis first.")
    
    # Check if genes exist in adata
    missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
    if missing_genes:
        raise ValueError(f"Genes not found in adata: {missing_genes}")
    
    # Calculate number of rows needed
    nrow = (len(gene_list) + ncol - 1) // ncol
    
    # Create figure
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    if nrow == 1 and ncol == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot each gene
    for i, gene in enumerate(gene_list):
        ax = axes[i]
        
        # Get expression values
        if scipy.sparse.issparse(adata.X):
            expr = adata[:, gene].X.toarray().flatten()
        else:
            expr = adata[:, gene].X.flatten()
        
        # Get pseudotime values
        pseudotime = adata.obs['dpt_pseudotime'].values
        
        # Create scatter plot
        scatter = ax.scatter(pseudotime, expr + 1, c=adata.obs[color_by].cat.codes, 
                           cmap='tab20', alpha=0.7, s=10)
        
        # Add trend line using LOWESS
        smoothed = lowess(expr, pseudotime, frac=0.3)
        ax.plot(smoothed[:, 0], smoothed[:, 1] + 1, 'k-', linewidth=2)
        
        # Set labels and title
        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('Expression')
        ax.set_title(gene)
        
        # Set y-axis to log scale
        ax.set_yscale('log')
        ax.set_ylim(0.3, 100)  # Adjust limits as needed
        
        # Add colorbar if this is the first plot
        if i == 0:
            # Create custom legend
            if color_by in adata.obs:
                unique_categories = adata.obs[color_by].unique()
                handles = [plt.scatter([], [], c=plt.cm.tab20(j/len(unique_categories)), 
                                    label=cat, s=50) 
                          for j, cat in enumerate(sorted(unique_categories))]
                ax.legend(handles=handles, title=color_by, bbox_to_anchor=(1.05, 1), 
                         loc='upper left')
    
    # Hide empty subplots
    for i in range(len(gene_list), len(axes)):
        axes[i].set_visible(False)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.show()

def plot_results(adata):
    """
    Plot various results from the analysis.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot UMAP with cell types
    sc.pl.umap(adata, color='cell_type', title='Cell Types', show=False)
    plt.savefig(os.path.join(output_dir, 'umap_cell_types.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot differentially expressed genes
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    plt.savefig(os.path.join(output_dir, 'deg_genes.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot UMAP with module score
    sc.pl.umap(adata, color='module_score', title='Gene Module Score', show=False)
    plt.savefig(os.path.join(output_dir, 'umap_module_score.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot PAGA
    sc.pl.paga(adata, color='cell_type', show=False)
    plt.savefig(os.path.join(output_dir, 'paga.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot diffusion map
    sc.pl.diffmap(adata, color='cell_type', show=False)
    plt.savefig(os.path.join(output_dir, 'diffmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot genes in pseudotime if pseudotime exists
    if 'dpt_pseudotime' in adata.obs and 'deg_genes' in adata.uns:
        # Get top differentially expressed genes
        deg_genes = adata.uns['deg_genes']
        if len(deg_genes) > 0:
            top_genes = deg_genes.head(6)['gene_id'].tolist()
            plot_genes_in_pseudotime(adata, top_genes)
    
    print(f"Plots saved in: {output_dir}")

# R code:
# main <- function() {
#   ...
# }

def main():
    """
    Main function to run the entire analysis pipeline.
    """
    # Load and preprocess data
    adata = load_and_preprocess_data()
    
    # Process and cluster
    adata = preprocess_and_cluster(adata)
    
    # Map cell types
    adata = map_cell_types(adata)
    
    # Analyze differential expression
    adata = analyze_differential_expression(adata)
    
    # Find gene modules
    adata = find_gene_modules(adata)
    
    # Perform trajectory analysis
    adata = perform_trajectory_analysis(adata)
    
    # Plot results
    plot_results(adata)

if __name__ == "__main__":
    main()