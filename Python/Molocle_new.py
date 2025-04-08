#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, List, Dict, Tuple
import os
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
import umap
from sklearn.manifold import TSNE
import warnings
import leidenalg
import igraph
from scipy import stats
from scipy.sparse import csr_matrix
import networkx as nx
from sklearn.neighbors import NearestNeighbors
import random
import splines
from sklearn.preprocessing import LabelEncoder
import statsmodels.api as sm
from scipy.sparse.linalg import svds
from statsmodels.stats.multitest import multipletests




def new_cell_data_set(expression_data: Union[np.ndarray, scipy.sparse.spmatrix],
                     cell_metadata: Optional[pd.DataFrame] = None,
                     gene_metadata: Optional[pd.DataFrame] = None) -> sc.AnnData:
    """
    Creates a new cell data set object.
    
    Parameters:
    -----------
    expression_data : numpy.ndarray or scipy.sparse.spmatrix
        A matrix of expression values, where each row is a gene and each column is a cell
    cell_metadata : pandas.DataFrame, optional (default: None)
        A DataFrame containing cell metadata, where each row is a cell
    gene_metadata : pandas.DataFrame, optional (default: None)
        A DataFrame containing gene metadata, where each row is a gene
        
    Returns:
    --------
    scanpy.AnnData
        An AnnData object containing the expression data and metadata
    """
    # Create AnnData object
    adata = sc.AnnData(expression_data.T)  # Transpose to match AnnData format (cells x genes)
    
    # Add cell metadata if provided
    if cell_metadata is not None:
        # Ensure cell metadata has the correct index
        if not all(cell_metadata.index == adata.obs_names):
            # Create a mapping from cell_id to index
            cell_id_to_index = {f'cell_{i}': i for i in range(adata.n_obs)}
            
            # Create a new index for cell_metadata
            new_index = [cell_id_to_index.get(cell_id, i) for i, cell_id in enumerate(cell_metadata.index)]
            
            # Reindex cell_metadata
            cell_metadata = cell_metadata.reindex(new_index)
        
        # Add cell metadata to AnnData object
        for col in cell_metadata.columns:
            adata.obs[col] = cell_metadata[col].values
    
    # Add gene metadata if provided
    if gene_metadata is not None:
        # Ensure gene metadata has the correct index
        if not all(gene_metadata.index == adata.var_names):
            # Create a mapping from gene_id to index
            gene_id_to_index = {f'gene_{i}': i for i in range(adata.n_vars)}
            
            # Create a new index for gene_metadata
            new_index = [gene_id_to_index.get(gene_id, i) for i, gene_id in enumerate(gene_metadata.index)]
            
            # Reindex gene_metadata
            gene_metadata = gene_metadata.reindex(new_index)
        
        # Add gene metadata to AnnData object
        for col in gene_metadata.columns:
            adata.var[col] = gene_metadata[col].values
    
    return adata

def normalize_expr_data(adata: sc.AnnData,
                       norm_method: str = "log",
                       pseudo_count: Optional[float] = None) -> np.ndarray:
    """
    Helper function to normalize the expression data prior to dimensionality reduction.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    norm_method : str, optional (default: "log")
        Method to use for normalization. Options are "log", "size_only", or "none".
    pseudo_count : float, optional (default: None)
        The amount to add to expression values before normalization.
        If None, 1 is used for log normalization and 0 for size_only.
        
    Returns:
    --------
    numpy.ndarray
        Normalized expression matrix
    """
    if norm_method not in ["log", "size_only", "none"]:
        raise ValueError("norm_method must be one of 'log', 'size_only', or 'none'")
    
    # Get expression matrix
    if scipy.sparse.issparse(adata.X):
        FM = adata.X.toarray()
    else:
        FM = adata.X.copy()
    
    # Set default pseudo_count if not provided
    if pseudo_count is None:
        pseudo_count = 1 if norm_method == "log" else 0
    
    # Normalize by size factors if available
    if 'size_factors' in adata.obs:
        FM = FM / adata.obs['size_factors'].values[:, np.newaxis]
    
    # Apply normalization method
    if norm_method == "log":
        FM = np.log2(FM + pseudo_count)
    elif norm_method == "size_only":
        FM = FM + pseudo_count
    
    return FM

def tfidf(count_matrix: Union[np.ndarray, scipy.sparse.spmatrix],
          frequencies: bool = True,
          log_scale_tf: bool = True,
          scale_factor: float = 100000) -> Dict:
    """
    Calculate TF-IDF transformation of count matrix.
    
    Parameters:
    -----------
    count_matrix : numpy.ndarray or scipy.sparse.spmatrix
        Count matrix to transform
    frequencies : bool, optional (default: True)
        Whether to use term frequency method
    log_scale_tf : bool, optional (default: True)
        Whether to log scale the term frequencies
    scale_factor : float, optional (default: 100000)
        Scale factor for log scaling
        
    Returns:
    --------
    dict
        Dictionary containing transformed matrix and parameters
    """
    # Convert to sparse matrix if not already
    if not scipy.sparse.issparse(count_matrix):
        count_matrix = scipy.sparse.csr_matrix(count_matrix)
    
    # Calculate term frequencies
    if frequencies:
        col_sums = count_matrix.sum(axis=0).A1
        tf = count_matrix.multiply(1/col_sums)
    else:
        col_sums = None
        tf = count_matrix
    
    # Log scale if requested
    if log_scale_tf:
        if frequencies:
            tf.data = np.log1p(tf.data * scale_factor)
        else:
            tf.data = np.log1p(tf.data)
    
    # Calculate IDF
    num_cols = count_matrix.shape[1]
    row_sums = (count_matrix > 0).sum(axis=1).A1
    idf = np.log(1 + num_cols / row_sums)
    
    # Calculate TF-IDF
    tf_idf_counts = tf.multiply(idf[:, np.newaxis])
    
    return {
        'tf_idf_counts': tf_idf_counts,
        'frequencies': frequencies,
        'log_scale_tf': log_scale_tf,
        'scale_factor': scale_factor,
        'col_sums': col_sums,
        'row_sums': row_sums,
        'num_cols': num_cols
    }

def calculate_size_factors(cds: sc.AnnData) -> sc.AnnData:
    """
    Calculate size factors for normalization.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell dataset object
        
    Returns:
    --------
    scanpy.AnnData
        The cell dataset object with size factors added
    """
    # Calculate size factors using library size normalization
    sc.pp.normalize_total(cds, target_sum=1e4, key_added='size_factors')
    return cds

def preprocess_cds(cds: sc.AnnData,
                  method: str = "pca",
                  num_dim: int = 50,
                  norm_method: str = "log",
                  pseudo_count: Optional[float] = 1.0,
                  use_genes: Optional[List[str]] = None,
                  verbose: bool = False) -> sc.AnnData:
    """
    Preprocess the cell data set for downstream analysis.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set to preprocess
    method : str, optional (default: "pca")
        Method to use for preprocessing ("pca" or "lsi")
    num_dim : int, optional (default: 50)
        Number of dimensions to use
    norm_method : str, optional (default: "log")
        Method to use for normalization
    pseudo_count : float, optional (default: 1.0)
        Pseudo count to add before normalization
    use_genes : list of str, optional (default: None)
        List of genes to use for preprocessing
    verbose : bool, optional (default: False)
        Whether to print progress messages
        
    Returns:
    --------
    scanpy.AnnData
        Preprocessed cell data set
    """
    if method not in ["pca", "lsi"]:
        raise ValueError("method must be either 'pca' or 'lsi'")
    
    # Calculate size factors if not already present
    if 'size_factors' not in cds.obs:
        if verbose:
            print("Calculating size factors...")
        cds = calculate_size_factors(cds)
    
    # Store original data
    cds.layers['counts'] = cds.X.copy()
    
    # Normalize expression data
    if verbose:
        print("Normalizing expression data...")
    norm_expr = normalize_expr_data(cds, norm_method=norm_method, pseudo_count=pseudo_count)
    
    # Handle NaN values
    if verbose:
        print("Handling NaN values...")
    norm_expr = np.nan_to_num(norm_expr, nan=0.0)
    
    # Store normalized data
    cds.layers['normalized'] = norm_expr
    
    # Use specified genes if provided
    if use_genes is not None:
        if verbose:
            print(f"Using {len(use_genes)} specified genes...")
        gene_mask = cds.var_names.isin(use_genes)
        norm_expr = norm_expr[:, gene_mask]
    
    # Perform dimensionality reduction
    if verbose:
        print(f"Performing {method.upper()}...")
    
    if method == "pca":
        # Scale the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(norm_expr)
        
        # Handle NaN values in scaled data
        scaled_data = np.nan_to_num(scaled_data, nan=0.0)
        
        # Store scaled data
        cds.layers['scaled'] = scaled_data
        
        # Perform PCA
        pca = TruncatedSVD(n_components=min(num_dim, min(scaled_data.shape)-1))
        reduced = pca.fit_transform(scaled_data)
        
        # Store results
        cds.obsm['X_pca'] = reduced
        cds.uns['pca'] = {
            'variance_ratio': pca.explained_variance_ratio_,
            'variance': pca.explained_variance_,
            'components': pca.components_
        }
        
    else:  # LSI
        # Perform TF-IDF
        tfidf_result = tfidf(norm_expr)
        tf_idf_matrix = tfidf_result['tf_idf_counts']
        
        # Handle NaN values in TF-IDF matrix
        if scipy.sparse.issparse(tf_idf_matrix):
            tf_idf_matrix.data = np.nan_to_num(tf_idf_matrix.data, nan=0.0)
        else:
            tf_idf_matrix = np.nan_to_num(tf_idf_matrix, nan=0.0)
        
        # Store TF-IDF data
        cds.layers['tf_idf'] = tf_idf_matrix
        
        # Perform truncated SVD
        u, s, vt = svds(tf_idf_matrix, k=min(num_dim, min(tf_idf_matrix.shape)-1))
        reduced = u * s[:, np.newaxis]
        
        # Store results
        cds.obsm['X_lsi'] = reduced
        cds.uns['lsi'] = {
            'singular_values': s,
            'components': vt
        }
    
    # Store preprocessing parameters
    cds.uns['preprocess_params'] = {
        'method': method,
        'num_dim': num_dim,
        'norm_method': norm_method,
        'pseudo_count': pseudo_count
    }
    
    return cds

def reduce_dimension(adata: sc.AnnData,
                    max_components: int = 2,
                    reduction_method: str = "umap",
                    preprocess_method: Optional[str] = None,
                    umap_metric: str = "cosine",
                    umap_min_dist: float = 0.1,
                    umap_n_neighbors: int = 15,
                    umap_fast_sgd: bool = False,
                    umap_nn_method: str = "annoy",
                    verbose: bool = False,
                    cores: int = 1,
                    build_nn_index: bool = False,
                    **kwargs) -> sc.AnnData:
    """
    Compute a projection of a cell data set object into a lower dimensional space.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    max_components : int, optional (default: 2)
        The dimensionality of the reduced space
    reduction_method : str, optional (default: "umap")
        The algorithm to use for dimensionality reduction.
        Options are "umap", "tsne", "pca", "lsi", or "aligned"
    preprocess_method : str, optional (default: None)
        The preprocessing method used on the data.
        Options are "pca", "lsi", or "aligned"
    umap_metric : str, optional (default: "cosine")
        The distance metric to use for UMAP
    umap_min_dist : float, optional (default: 0.1)
        The minimum distance for UMAP
    umap_n_neighbors : int, optional (default: 15)
        The number of neighbors to use for UMAP
    umap_fast_sgd : bool, optional (default: False)
        Whether to use fast SGD for UMAP
    umap_nn_method : str, optional (default: "annoy")
        The nearest neighbor method to use for UMAP
    verbose : bool, optional (default: False)
        Whether to print progress messages
    cores : int, optional (default: 1)
        The number of cores to use
    build_nn_index : bool, optional (default: False)
        Whether to build a nearest neighbor index
    **kwargs : dict
        Additional arguments to pass to the dimensionality reduction function
        
    Returns:
    --------
    scanpy.AnnData
        The annotated data matrix with reduced dimensions
    """
    # Validate reduction method
    if reduction_method.lower() not in ["umap", "tsne", "pca", "lsi", "aligned"]:
        raise ValueError("reduction_method must be one of 'umap', 'tsne', 'pca', 'lsi', or 'aligned'")
    reduction_method = reduction_method.lower()
    
    # Determine preprocess method if not provided
    if preprocess_method is None:
        if "aligned" in adata.obsm:
            preprocess_method = "aligned"
            if verbose:
                print("No preprocess_method specified, using preprocess_method = 'aligned'")
        else:
            preprocess_method = "pca"
            if verbose:
                print("No preprocess_method specified, using preprocess_method = 'pca'")
    elif preprocess_method.lower() not in ["pca", "lsi", "aligned"]:
        raise ValueError("preprocess_method must be one of 'pca', 'lsi', or 'aligned'")
    preprocess_method = preprocess_method.lower()
    
    # Check if data has been preprocessed
    if f'X_{preprocess_method}' not in adata.obsm:
        # Try to find any preprocessed data
        available_methods = [key.split('_')[1] for key in adata.obsm.keys() if key.startswith('X_')]
        if available_methods:
            preprocess_method = available_methods[0]
            if verbose:
                print(f"Using available preprocessed data with method: {preprocess_method}")
        else:
            raise ValueError(f"No preprocessed data found. Please run preprocess_cds first.")
    
    # Get preprocessed matrix
    preprocess_mat = adata.obsm[f'X_{preprocess_method}']
    
    # Perform dimensionality reduction
    if reduction_method == "pca":
        if preprocess_method != "pca":
            raise ValueError("preprocess_method must be 'pca' when reduction_method = 'pca'")
        if verbose:
            print("Returning preprocessed PCA matrix")
        adata.obsm['X_pca_reduced'] = preprocess_mat[:, :max_components]
        
    elif reduction_method == "lsi":
        if preprocess_method != "lsi":
            raise ValueError("preprocess_method must be 'lsi' when reduction_method = 'lsi'")
        if verbose:
            print("Returning preprocessed LSI matrix")
        adata.obsm['X_lsi_reduced'] = preprocess_mat[:, :max_components]
        
    elif reduction_method == "aligned":
        if preprocess_method != "aligned":
            raise ValueError("preprocess_method must be 'aligned' when reduction_method = 'aligned'")
        if verbose:
            print("Returning preprocessed Aligned matrix")
        adata.obsm['X_aligned_reduced'] = preprocess_mat[:, :max_components]
        
    elif reduction_method == "tsne":
        if verbose:
            print("Reducing dimension by t-SNE...")
        
        # Set random seed for reproducibility
        np.random.seed(2016)
        
        # Perform t-SNE
        tsne = TSNE(n_components=max_components, 
                   metric='euclidean',
                   init='random',
                   random_state=2016,
                   n_jobs=cores,
                   **kwargs)
        tsne_data = tsne.fit_transform(preprocess_mat)
        
        # Store results
        adata.obsm['X_tsne'] = tsne_data
        adata.uns['tsne'] = {
            'params': {
                'n_components': max_components,
                'perplexity': tsne.perplexity,
                'early_exaggeration': tsne.early_exaggeration,
                'learning_rate': tsne.learning_rate,
                'n_iter': tsne.n_iter
            }
        }
        
    elif reduction_method == "umap":
        if verbose:
            print("Running Uniform Manifold Approximation and Projection")
        
        # Set random seed for reproducibility
        np.random.seed(2016)
        
        # Warn about non-deterministic results
        if umap_fast_sgd or cores > 1:
            warnings.warn("Note: reduce_dimension will produce slightly different output each time "
                         "you run it unless you set 'umap_fast_sgd = False' and 'cores = 1'")
        
        # Perform UMAP
        reducer = umap.UMAP(
            n_components=max_components,
            metric=umap_metric,
            min_dist=umap_min_dist,
            n_neighbors=umap_n_neighbors,
            low_memory=True if umap_fast_sgd else False,
            n_jobs=cores,
            verbose=verbose,
            **kwargs
        )
        
        # Fit and transform
        umap_data = reducer.fit_transform(preprocess_mat)
        
        # Store results
        adata.obsm['X_umap'] = umap_data
        adata.uns['umap'] = {
            'params': {
                'n_components': max_components,
                'metric': umap_metric,
                'min_dist': umap_min_dist,
                'n_neighbors': umap_n_neighbors,
                'fast_sgd': umap_fast_sgd,
                'nn_method': umap_nn_method
            }
        }
    
    # Clear out old graphs if they exist
    for key in ['principal_graph', 'principal_graph_aux', 'clusters']:
        if key in adata.uns and reduction_method in adata.uns[key]:
            del adata.uns[key][reduction_method]
    
    return adata

def cluster_cells_make_graph(data: np.ndarray,
                            weight: bool,
                            cell_names: List[str],
                            k: int = 20,
                            verbose: bool = False) -> Dict:
    """
    Create a k-nearest neighbor graph for cell clustering.
    
    Parameters:
    -----------
    data : numpy.ndarray
        The data matrix to use for graph creation
    weight : bool
        Whether to use weights in the graph
    cell_names : list of str
        List of cell names
    k : int, optional (default: 20)
        Number of nearest neighbors
    verbose : bool, optional (default: False)
        Whether to print progress messages
        
    Returns:
    --------
    dict
        Dictionary containing graph information
    """
    if verbose:
        print("Computing nearest neighbors...")
    
    # Create nearest neighbors object
    nn = NearestNeighbors(n_neighbors=k+1)  # +1 because a point is its own nearest neighbor
    nn.fit(data)
    
    # Find k-nearest neighbors
    distances, indices = nn.kneighbors()
    
    # Remove self-connections (first column)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    
    # Create adjacency matrix
    n_cells = len(cell_names)
    adjacency_matrix = np.zeros((n_cells, n_cells))
    
    # Fill adjacency matrix
    for i in range(n_cells):
        for j, dist in zip(indices[i], distances[i]):
            if weight:
                # Use Jaccard coefficient as weight
                intersection = len(set(indices[i]).intersection(indices[j]))
                union = len(set(indices[i]).union(indices[j]))
                weight_value = intersection / union if union > 0 else 0
            else:
                weight_value = 1
            
            adjacency_matrix[i, j] = weight_value
            adjacency_matrix[j, i] = weight_value
    
    return {
        'adjacency_matrix': adjacency_matrix,
        'distances': distances,
        'indices': indices,
        'cell_names': cell_names
    }

def louvain_clustering(data: np.ndarray,
                      pd: pd.DataFrame,
                      weight: bool = False,
                      k: int = 20,
                      louvain_iter: int = 1,
                      random_seed: int = 0,
                      verbose: bool = False) -> Dict:
    """
    Perform Louvain clustering.
    
    Parameters:
    -----------
    data : numpy.ndarray
        The data matrix
    pd : pandas.DataFrame
        Cell metadata
    weight : bool, optional (default: False)
        Whether to use weights in the graph
    k : int, optional (default: 20)
        Number of nearest neighbors
    louvain_iter : int, optional (default: 1)
        Number of iterations for Louvain clustering
    random_seed : int, optional (default: 0)
        Random seed
    verbose : bool, optional (default: False)
        Whether to print progress messages
        
    Returns:
    --------
    dict
        Dictionary containing clustering results
    """
    cell_names = pd.index.tolist()
    
    if not all(cell_names == pd.index.tolist()):
        raise ValueError("Phenotype and row name from the data doesn't match")
    
    # Create graph
    graph_result = cluster_cells_make_graph(
        data=data,
        weight=weight,
        cell_names=cell_names,
        k=k,
        verbose=verbose
    )
    
    if verbose:
        print("  Run louvain clustering...")
    
    # Convert NetworkX graph to igraph
    g_igraph = igraph.Graph()
    g_igraph.add_vertices(list(graph_result['g'].nodes()))
    g_igraph.add_edges([(u, v) for u, v in graph_result['g'].edges()])
    g_igraph.es['weight'] = [d['weight'] for u, v, d in graph_result['g'].edges(data=True)]
    
    # Set random seed
    if random_seed is not None:
        random.seed(random_seed)
    
    # Run Louvain clustering
    t_start = time.time()
    Qp = -1
    optim_res = None
    
    if louvain_iter < 1:
        warnings.warn("bad loop: louvain_iter is < 1")
    
    for iter in range(louvain_iter):
        if verbose:
            print(f"Running louvain iteration {iter+1}...")
        
        # Run Louvain algorithm
        Q = g_igraph.community_multilevel(weights='weight')
        
        if optim_res is None:
            Qp = Q.modularity
            optim_res = Q
        else:
            Qt = Q.modularity
            if Qt > Qp:
                optim_res = Q
                Qp = Qt
    
    t_end = time.time()
    
    if verbose:
        print(f'Maximal modularity is {Qp}')
        print(f"\nRun kNN based graph clustering DONE, totally takes {t_end - t_start} s.")
        print(f"  -Number of clusters: {len(set(optim_res.membership))}")
    
    # Set vertex names
    g_igraph.vs['name'] = [str(v) for v in g_igraph.vs]
    
    return {
        'g': g_igraph,
        'relations': graph_result['relations'],
        'distMatrix': graph_result['distMatrix'],
        'coord': None,
        'edge_links': None,
        'optim_res': optim_res
    }

def leiden_clustering(data: np.ndarray,
                     pd: pd.DataFrame,
                     weight: bool = False,
                     k: int = 20,
                     num_iter: int = 2,
                     resolution_parameter: Optional[float] = 0.0001,
                     random_seed: Optional[int] = None,
                     verbose: bool = False,
                     **kwargs) -> Dict:
    """
    Perform Leiden clustering on the data.
    
    Parameters:
    -----------
    data : numpy.ndarray
        The data matrix to cluster
    pd : pandas.DataFrame
        The pandas DataFrame containing cell metadata
    weight : bool, optional (default: False)
        Whether to use weights in the clustering
    k : int, optional (default: 20)
        Number of nearest neighbors
    num_iter : int, optional (default: 2)
        Number of iterations for the Leiden algorithm
    resolution_parameter : float, optional (default: 0.0001)
        Resolution parameter for the Leiden algorithm
    random_seed : int, optional (default: None)
        Random seed for reproducibility
    verbose : bool, optional (default: False)
        Whether to print progress messages
    **kwargs : dict
        Additional arguments to pass to the Leiden algorithm
        
    Returns:
    --------
    dict
        Dictionary containing clustering results
    """
    # Get cell names from pd DataFrame
    cell_names = pd.index.tolist()
    
    # Create k-nearest neighbor graph
    if verbose:
        print("Creating k-nearest neighbor graph...")
    
    knn_result = cluster_cells_make_graph(data, weight, cell_names, k, verbose)
    
    # Convert adjacency matrix to igraph
    if verbose:
        print("Converting adjacency matrix to igraph...")
    
    g = igraph.Graph.Weighted_Adjacency(knn_result['adjacency_matrix'].tolist(),
                                      mode=igraph.ADJ_UNDIRECTED,
                                      attr="weight")
    
    # Set random seed if provided
    if random_seed is not None:
        random.seed(random_seed)
    
    # Run Leiden algorithm
    if verbose:
        print("Running Leiden algorithm...")
    
    # Set default resolution parameter if not provided
    if resolution_parameter is None:
        resolution_parameter = 0.0001
    
    # Initialize partition dictionary
    partition = {}
    
    # Run Leiden clustering
    partition_result = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution_parameter,
        weights='weight',
        n_iterations=num_iter,
        seed=random_seed,
        **kwargs
    )
    
    # Convert membership to dictionary
    for i, cluster in enumerate(partition_result.membership):
        partition[cell_names[i]] = cluster
    
    # Return results
    return {
        'g': g,
        'optim_res': {
            'membership': partition,
            'modularity': partition_result.modularity,
            'resolution_parameter': resolution_parameter
        },
        'cell_membership': partition,
        'knn_result': knn_result
    }

def compute_partitions(g: igraph.Graph,
                      optim_res: Union[igraph.clustering.VertexClustering, Dict],
                      qval_thresh: float = 0.05,
                      verbose: bool = False) -> Dict:
    """
    Compute partitions from clustering results.
    
    Parameters:
    -----------
    g : igraph.Graph
        The graph object
    optim_res : igraph.clustering.VertexClustering or dict
        The clustering results
    qval_thresh : float, optional (default: 0.05)
        Q-value threshold for partitioning
    verbose : bool, optional (default: False)
        Whether to print progress messages
        
    Returns:
    --------
    dict
        Dictionary containing partition information
    """
    if verbose:
        print("Computing partitions...")
    
    # Get membership from optim_res
    if isinstance(optim_res, dict):
        membership = list(optim_res['membership'].values())
    else:
        membership = optim_res.membership
    
    # Get unique clusters
    unique_clusters = sorted(set(membership))
    num_clusters = len(unique_clusters)
    
    # Create cluster matrix
    cluster_mat = np.zeros((num_clusters, num_clusters))
    edges_per_module = np.zeros(num_clusters)
    
    # Count edges between and within modules
    for edge in g.es:
        source = membership[edge.source]
        target = membership[edge.target]
        weight = edge['weight'] if 'weight' in edge.attributes() else 1.0
        
        cluster_mat[source, target] += weight
        cluster_mat[target, source] += weight
        
        if source == target:
            edges_per_module[source] += weight
        else:
            edges_per_module[source] += weight / 2
            edges_per_module[target] += weight / 2
    
    # Convert to dense arrays
    edges_per_module = np.array(edges_per_module)
    total_edges = np.sum(edges_per_module)
    
    if total_edges == 0:
        if verbose:
            print("Warning: No edges found in the graph")
        return {
            'cluster_g': g,
            'cluster_mat': cluster_mat,
            'edges_per_module': edges_per_module,
            'total_edges': total_edges
        }
    
    # Compute theta matrix
    theta = np.outer(edges_per_module / total_edges, edges_per_module / total_edges)
    
    # Compute p-values
    pvals = np.zeros((num_clusters, num_clusters))
    for i in range(num_clusters):
        for j in range(i, num_clusters):
            if theta[i, j] > 0:
                # Use hypergeometric test
                k = cluster_mat[i, j]
                M = total_edges
                n = edges_per_module[i]
                N = edges_per_module[j]
                
                # Calculate p-value
                pval = 1 - stats.hypergeom.cdf(k-1, M, n, N)
                pvals[i, j] = pval
                pvals[j, i] = pval
    
    # Adjust p-values for multiple testing
    flat_pvals = pvals[np.triu_indices(num_clusters, k=1)]
    if len(flat_pvals) > 0:
        _, qvals, _, _ = multipletests(flat_pvals, method='fdr_bh')
        qval_mat = np.zeros_like(pvals)
        qval_mat[np.triu_indices(num_clusters, k=1)] = qvals
        qval_mat = qval_mat + qval_mat.T
    else:
        qval_mat = np.ones_like(pvals)
    
    # Create cluster graph
    cluster_g = igraph.Graph()
    cluster_g.add_vertices(num_clusters)
    
    # Add edges where q-value is below threshold
    edges = []
    for i in range(num_clusters):
        for j in range(i+1, num_clusters):
            if qval_mat[i, j] < qval_thresh:
                edges.append((i, j))
    
    cluster_g.add_edges(edges)
    
    return {
        'cluster_g': cluster_g,
        'cluster_mat': cluster_mat,
        'pvals': pvals,
        'qvals': qval_mat,
        'edges_per_module': edges_per_module,
        'total_edges': total_edges
    }

def cluster_cells(adata: sc.AnnData,
                 reduction_method: str = "umap",
                 k: int = 20,
                 cluster_method: str = "leiden",
                 num_iter: int = 2,
                 partition_qval: float = 0.05,
                 weight: bool = False,
                 resolution: Optional[float] = None,
                 random_seed: int = 42,
                 verbose: bool = False,
                 **kwargs) -> sc.AnnData:
    """
    Cluster cells using Louvain/Leiden community detection.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    reduction_method : str, optional (default: "umap")
        The dimensionality reduction method upon which to base clustering.
        Options are "umap", "tsne", "pca", "lsi", or "aligned"
    k : int, optional (default: 20)
        Number of nearest neighbors to use when creating the k-nearest neighbor graph
    cluster_method : str, optional (default: "leiden")
        The clustering method to use. Options are "louvain" or "leiden"
    num_iter : int, optional (default: 2)
        Number of iterations used for Louvain/Leiden clustering
    partition_qval : float, optional (default: 0.05)
        Q-value cutoff to determine when to partition
    weight : bool, optional (default: False)
        Whether to use Jaccard coefficients for two nearest neighbors as the weight
    resolution : float, optional (default: None)
        Parameter that controls the resolution of clustering
    random_seed : int, optional (default: 42)
        The seed used by the random number generator
    verbose : bool, optional (default: False)
        Whether to print progress messages
    **kwargs : dict
        Additional arguments passed to the clustering function
        
    Returns:
    --------
    scanpy.AnnData
        The annotated data matrix with cluster and partition information
    """
    # Validate reduction method
    if reduction_method.lower() not in ["umap", "tsne", "pca", "lsi", "aligned"]:
        raise ValueError("reduction_method must be one of 'umap', 'tsne', 'pca', 'lsi', or 'aligned'")
    reduction_method = reduction_method.lower()
    
    # Validate cluster method
    if cluster_method.lower() not in ["leiden", "louvain"]:
        raise ValueError("cluster_method must be one of 'leiden', 'louvain'")
    cluster_method = cluster_method.lower()
    
    # Check if data has been reduced
    reduced_key = f'X_{reduction_method}'
    if reduced_key not in adata.obsm:
        # Try to find any reduced dimension data
        available_methods = [key.split('_')[1] for key in adata.obsm.keys() if key.startswith('X_')]
        if available_methods:
            reduction_method = available_methods[0]
            reduced_key = f'X_{reduction_method}'
            if verbose:
                print(f"Using available reduced dimension data with method: {reduction_method}")
        else:
            raise ValueError(f"No reduced dimension data found. Please run reduce_dimension first.")
    
    # Get reduced dimension matrix
    reduced_dim_res = adata.obsm[reduced_key]
    
    # Set random seed
    if random_seed is not None:
        np.random.seed(random_seed)
        random.seed(random_seed)
    
    if verbose:
        print(f"Running {cluster_method} clustering algorithm...")
    
    # Perform clustering
    if cluster_method == 'louvain':
        cluster_result = louvain_clustering(
            data=reduced_dim_res,
            pd=adata.obs,
            weight=weight,
            k=k,
            louvain_iter=num_iter,
            random_seed=random_seed,
            verbose=verbose
        )
        
        # Get membership
        membership = cluster_result['optim_res'].membership
        
        # Compute partitions if more than one cluster
        if len(set(membership)) > 1:
            cluster_graph_res = compute_partitions(
                cluster_result['g'],
                cluster_result['optim_res'],
                partition_qval,
                verbose
            )
            
            # Get partition membership
            partition_membership = igraph.Graph.components(cluster_graph_res['cluster_g']).membership
            partitions = [partition_membership[m] for m in membership]
            partitions = pd.Series(partitions, index=adata.obs_names)
        else:
            partitions = pd.Series([1] * len(adata.obs_names), index=adata.obs_names)
        
        # Create clusters series
        clusters = pd.Series(membership, index=adata.obs_names)
        
        # Store results
        adata.uns['clusters'] = adata.uns.get('clusters', {})
        adata.uns['clusters'][reduction_method] = {
            'cluster_result': cluster_result,
            'partitions': partitions,
            'clusters': clusters
        }
        
    elif cluster_method == 'leiden':
        # Switch to leiden if resolution is provided and method is louvain
        if resolution is not None and cluster_method == 'louvain':
            print(f"Resolution can only be used when cluster_method is 'leiden'. "
                  f"Switching to leiden clustering.")
            cluster_method = 'leiden'
        
        cluster_result = leiden_clustering(
            data=reduced_dim_res,
            pd=adata.obs,
            weight=weight,
            k=k,
            num_iter=num_iter,
            resolution_parameter=resolution,
            random_seed=random_seed,
            verbose=verbose,
            **kwargs
        )
        
        # Get membership
        membership = list(cluster_result['optim_res']['membership'].values())
        
        # Compute partitions if more than one cluster
        if len(set(membership)) > 1:
            cluster_graph_res = compute_partitions(
                cluster_result['g'],
                cluster_result['optim_res'],
                partition_qval,
                verbose
            )
            
            # Get partition membership
            partition_membership = igraph.Graph.components(cluster_graph_res['cluster_g']).membership
            partitions = [partition_membership[m] for m in membership]
            partitions = pd.Series(partitions, index=adata.obs_names)
        else:
            partitions = pd.Series([1] * len(adata.obs_names), index=adata.obs_names)
        
        # Create clusters series
        clusters = pd.Series(membership, index=adata.obs_names)
        
        # Store results
        adata.uns['clusters'] = adata.uns.get('clusters', {})
        adata.uns['clusters'][reduction_method] = {
            'cluster_result': cluster_result,
            'partitions': partitions,
            'clusters': clusters
        }
    
    # Add cluster and partition information to obs
    adata.obs['cluster'] = adata.uns['clusters'][reduction_method]['clusters']
    adata.obs['partition'] = adata.uns['clusters'][reduction_method]['partitions']
    
    return adata

def plot_cells(adata: sc.AnnData,
               x: int = 1,
               y: int = 2,
               reduction_method: str = "umap",
               color_cells_by: str = "cluster",
               group_cells_by: str = "cluster",
               genes: Optional[Union[List[str], pd.DataFrame]] = None,
               show_trajectory_graph: bool = True,
               trajectory_graph_color: str = "grey28",
               trajectory_graph_segment_size: float = 0.75,
               norm_method: str = "log",
               label_cell_groups: bool = True,
               label_groups_by_cluster: bool = True,
               group_label_size: float = 2,
               labels_per_group: int = 1,
               label_branch_points: bool = True,
               label_roots: bool = True,
               label_leaves: bool = True,
               graph_label_size: float = 2,
               cell_size: float = 0.35,
               cell_stroke: float = None,
               alpha: float = 1,
               min_expr: float = 0.1,
               rasterize: bool = False,
               scale_to_range: bool = True,
               label_principal_points: bool = False) -> plt.Figure:
    """
    Plots the cells along with their trajectories.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    x : int, optional (default: 1)
        The column of reduced dimensions to plot on the horizontal axis
    y : int, optional (default: 2)
        The column of reduced dimensions to plot on the vertical axis
    reduction_method : str, optional (default: "umap")
        The lower dimensional space in which to plot cells.
        Must be one of "umap", "tsne", "pca", "lsi", or "aligned"
    color_cells_by : str, optional (default: "cluster")
        What to use for coloring the cells. Must be either the
        name of a column of adata.obs, or one of "clusters", "partitions", or
        "pseudotime".
    group_cells_by : str, optional (default: "cluster")
        How to group cells when labeling them. Must be either
        the name of a column of adata.obs, or one of "clusters" or "partitions".
        If a column in adata.obs, must be a categorical variable.
    genes : list of str or pandas.DataFrame, optional (default: None)
        Facet the plot, showing the expression of each gene in a facet
        panel. Must be either a list of gene ids (or short names), or a dataframe
        with two columns that groups the genes into modules that will be
        aggregated prior to plotting. If the latter, the first column must be gene
        ids, and the second must the group for each gene.
    show_trajectory_graph : bool, optional (default: True)
        Whether to render the principal graph for the trajectory.
        Requires that learn_graph() has been called on adata.
    trajectory_graph_color : str, optional (default: "grey28")
        The color to be used for plotting the trajectory graph.
    trajectory_graph_segment_size : float, optional (default: 0.75)
        The size of the line segments used for plotting the trajectory graph.
    norm_method : str, optional (default: "log")
        How to normalize gene expression scores prior to plotting them.
        Must be one of "log" or "size_only".
    label_cell_groups : bool, optional (default: True)
        Whether to label cells in each group (as specified
        by group_cells_by) according to the most frequently occurring label(s) (as
        specified by color_cells_by) in the group. If false, plot_cells() simply
        adds a traditional color legend.
    label_groups_by_cluster : bool, optional (default: True)
        Instead of labeling each cluster of cells,
        place each label once, at the centroid of all cells carrying that label.
    group_label_size : float, optional (default: 2)
        Font size to be used for cell group labels.
    labels_per_group : int, optional (default: 1)
        How many labels to plot for each group of cells.
        Defaults to 1, which plots only the most frequent label per group.
    label_branch_points : bool, optional (default: True)
        Whether to plot a label for each branch point in the principal graph.
    label_roots : bool, optional (default: True)
        Whether to plot a label for each root in the principal graph.
    label_leaves : bool, optional (default: True)
        Whether to plot a label for each leaf node in the principal graph.
    graph_label_size : float, optional (default: 2)
        How large to make the branch, root, and leaf labels.
    cell_size : float, optional (default: 0.35)
        The size of the point for each cell
    cell_stroke : float, optional (default: None)
        The stroke used for plotting each cell - default is 1/2 of the cell_size
    alpha : float, optional (default: 1)
        Alpha for the cells. Useful for reducing overplotting.
    min_expr : float, optional (default: 0.1)
        Minimum expression threshold for plotting genes
    rasterize : bool, optional (default: False)
        Whether to plot cells as a rastered bitmap. Requires the
        ggrastr package.
    scale_to_range : bool, optional (default: True)
        Logical indicating whether to scale expression to
        percent of maximum expression.
    label_principal_points : bool, optional (default: False)
        Logical indicating whether to label roots,
        leaves, and branch points with principal point names. This is useful for
        order_cells and choose_graph_segments in non-interactive mode.
        
    Returns:
    --------
    matplotlib.figure.Figure
        A matplotlib figure object
    """
    # Validate reduction method
    if reduction_method.lower() not in ["umap", "tsne", "pca", "lsi", "aligned"]:
        raise ValueError("reduction_method must be one of 'umap', 'tsne', 'pca', 'lsi', or 'aligned'")
    reduction_method = reduction_method.lower()
    
    # Check if data has been reduced
    if reduction_method not in adata.obsm:
        raise ValueError(f"No dimensionality reduction for {reduction_method} calculated. "
                        f"Please run reduce_dimension with reduction_method = {reduction_method} "
                        f"before running plot_cells")
    
    # Set cell stroke if not provided
    if cell_stroke is None:
        cell_stroke = cell_size / 2
    
    # Get reduced dimension matrix
    reduced_dim = adata.obsm[f'X_{reduction_method}']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot cells
    if genes is None:
        # Plot cells colored by the specified attribute
        if color_cells_by in ["clusters", "partitions"]:
            if color_cells_by not in adata.uns.get('clusters', {}).get(reduction_method, {}):
                raise ValueError(f"No {color_cells_by} found for {reduction_method}. "
                                f"Please run cluster_cells with reduction_method = {reduction_method} "
                                f"before running plot_cells")
            
            color_data = adata.uns['clusters'][reduction_method][color_cells_by]
            color_data = pd.Series(color_data, index=adata.obs_names)
        elif color_cells_by == "pseudotime":
            if "pseudotime" not in adata.obs:
                raise ValueError("No pseudotime found. "
                                f"Please run order_cells before running plot_cells")
            
            color_data = adata.obs["pseudotime"]
        else:
            if color_cells_by not in adata.obs:
                raise ValueError(f"Column {color_cells_by} not found in adata.obs")
            
            color_data = adata.obs[color_cells_by]
        
        # Create scatter plot
        scatter = ax.scatter(reduced_dim[:, x-1], reduced_dim[:, y-1], 
                            c=color_data, cmap='viridis', 
                            s=cell_size*100, alpha=alpha, 
                            edgecolors='black', linewidths=cell_stroke)
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax, label=color_cells_by)
        
        # Label cell groups if requested
        if label_cell_groups:
            if group_cells_by in ["clusters", "partitions"]:
                if group_cells_by not in adata.uns.get('clusters', {}).get(reduction_method, {}):
                    raise ValueError(f"No {group_cells_by} found for {reduction_method}. "
                                    f"Please run cluster_cells with reduction_method = {reduction_method} "
                                    f"before running plot_cells")
                
                group_data = adata.uns['clusters'][reduction_method][group_cells_by]
                group_data = pd.Series(group_data, index=adata.obs_names)
            else:
                if group_cells_by not in adata.obs:
                    raise ValueError(f"Column {group_cells_by} not found in adata.obs")
                
                group_data = adata.obs[group_cells_by]
            
            # Get unique groups
            unique_groups = group_data.unique()
            
            # For each group, find the centroid and add a label
            for group in unique_groups:
                group_indices = group_data == group
                group_coords = reduced_dim[group_indices]
                
                if len(group_coords) > 0:
                    centroid = np.mean(group_coords, axis=0)
                    
                    # Find the most common color_cells_by value in this group
                    if color_cells_by in ["clusters", "partitions"]:
                        color_values = color_data[group_indices]
                    elif color_cells_by == "pseudotime":
                        color_values = color_data[group_indices]
                    else:
                        color_values = color_data[group_indices]
                    
                    most_common = color_values.mode().iloc[0]
                    
                    # Add label
                    ax.text(centroid[x-1], centroid[y-1], str(most_common), 
                            fontsize=group_label_size*10, ha='center', va='center')
    else:
        # Plot gene expression
        if isinstance(genes, list):
            # Check if genes exist in adata
            for gene in genes:
                if gene not in adata.var_names:
                    raise ValueError(f"Gene {gene} not found in adata.var_names")
            
            # Get expression data
            if scipy.sparse.issparse(adata.X):
                expr_data = adata[:, genes].X.toarray()
            else:
                expr_data = adata[:, genes].X
            
            # Normalize expression data
            if norm_method == "log":
                expr_data = np.log1p(expr_data)
            elif norm_method == "size_only":
                pass
            else:
                raise ValueError("norm_method must be one of 'log' or 'size_only'")
            
            # Scale to range if requested
            if scale_to_range:
                expr_data = (expr_data - expr_data.min(axis=0)) / (expr_data.max(axis=0) - expr_data.min(axis=0))
            
            # Create subplots for each gene
            n_genes = len(genes)
            n_cols = min(3, n_genes)
            n_rows = (n_genes + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))
            axes = axes.flatten()
            
            for i, gene in enumerate(genes):
                ax = axes[i]
                
                # Create scatter plot
                scatter = ax.scatter(reduced_dim[:, x-1], reduced_dim[:, y-1], 
                                    c=expr_data[:, i], cmap='viridis', 
                                    s=cell_size*100, alpha=alpha, 
                                    edgecolors='black', linewidths=cell_stroke)
                
                # Add colorbar
                plt.colorbar(scatter, ax=ax, label=gene)
                
                # Set title
                ax.set_title(gene)
                
                # Remove unused subplots
                for j in range(i+1, len(axes)):
                    fig.delaxes(axes[j])
        else:
            # Group genes by module
            if not isinstance(genes, pd.DataFrame) or genes.shape[1] != 2:
                raise ValueError("genes must be either a list of gene ids or a dataframe with two columns")
            
            gene_ids = genes.iloc[:, 0].tolist()
            gene_modules = genes.iloc[:, 1].tolist()
            
            # Check if genes exist in adata
            for gene in gene_ids:
                if gene not in adata.var_names:
                    raise ValueError(f"Gene {gene} not found in adata.var_names")
            
            # Get unique modules
            unique_modules = np.unique(gene_modules)
            
            # Create subplots for each module
            n_modules = len(unique_modules)
            n_cols = min(3, n_modules)
            n_rows = (n_modules + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))
            axes = axes.flatten()
            
            for i, module in enumerate(unique_modules):
                ax = axes[i]
                
                # Get genes in this module
                module_genes = [gene for gene, m in zip(gene_ids, gene_modules) if m == module]
                
                # Get expression data
                if scipy.sparse.issparse(adata.X):
                    expr_data = adata[:, module_genes].X.toarray()
                else:
                    expr_data = adata[:, module_genes].X
                
                # Normalize expression data
                if norm_method == "log":
                    expr_data = np.log1p(expr_data)
                elif norm_method == "size_only":
                    pass
                else:
                    raise ValueError("norm_method must be one of 'log' or 'size_only'")
                
                # Scale to range if requested
                if scale_to_range:
                    expr_data = (expr_data - expr_data.min(axis=0)) / (expr_data.max(axis=0) - expr_data.min(axis=0))
                
                # Aggregate expression data
                agg_expr = np.mean(expr_data, axis=1)
                
                # Create scatter plot
                scatter = ax.scatter(reduced_dim[:, x-1], reduced_dim[:, y-1], 
                                    c=agg_expr, cmap='viridis', 
                                    s=cell_size*100, alpha=alpha, 
                                    edgecolors='black', linewidths=cell_stroke)
                
                # Add colorbar
                plt.colorbar(scatter, ax=ax, label=f"Module {module}")
                
                # Set title
                ax.set_title(f"Module {module}")
                
                # Remove unused subplots
                for j in range(i+1, len(axes)):
                    fig.delaxes(axes[j])
    
    # Plot trajectory graph if requested
    if show_trajectory_graph:
        if "principal_graph" not in adata.uns:
            warnings.warn("No principal graph found. Please run learn_graph before running plot_cells")
        else:
            # Get principal graph
            principal_graph = adata.uns["principal_graph"]
            
            # Plot edges
            for edge in principal_graph["edge_list"]:
                ax.plot([edge[0][x-1], edge[1][x-1]], 
                        [edge[0][y-1], edge[1][y-1]], 
                        color=trajectory_graph_color, 
                        linewidth=trajectory_graph_segment_size)
            
            # Plot nodes
            ax.scatter(principal_graph["node_positions"][:, x-1], 
                      principal_graph["node_positions"][:, y-1], 
                      color=trajectory_graph_color, 
                      s=cell_size*200)
            
            # Label branch points, roots, and leaves if requested
            if label_branch_points and "branch_points" in principal_graph:
                for point in principal_graph["branch_points"]:
                    ax.text(principal_graph["node_positions"][point, x-1], 
                           principal_graph["node_positions"][point, y-1], 
                           "B", fontsize=graph_label_size*10, 
                           ha='center', va='center')
            
            if label_roots and "roots" in principal_graph:
                for point in principal_graph["roots"]:
                    ax.text(principal_graph["node_positions"][point, x-1], 
                           principal_graph["node_positions"][point, y-1], 
                           "R", fontsize=graph_label_size*10, 
                           ha='center', va='center')
            
            if label_leaves and "leaves" in principal_graph:
                for point in principal_graph["leaves"]:
                    ax.text(principal_graph["node_positions"][point, x-1], 
                           principal_graph["node_positions"][point, y-1], 
                           "L", fontsize=graph_label_size*10, 
                           ha='center', va='center')
            
            # Label principal points if requested
            if label_principal_points and "principal_points" in principal_graph:
                for i, point in enumerate(principal_graph["principal_points"]):
                    ax.text(principal_graph["node_positions"][point, x-1], 
                           principal_graph["node_positions"][point, y-1], 
                           str(i), fontsize=graph_label_size*10, 
                           ha='center', va='center')
    
    # Set axis labels
    ax.set_xlabel(f"{reduction_method.upper()} {x}")
    ax.set_ylabel(f"{reduction_method.upper()} {y}")
    
    # Set title
    if genes is None:
        ax.set_title(f"Cells colored by {color_cells_by}")
    else:
        ax.set_title("Gene expression")
    
    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add grid
    ax.grid(False)
    
    # Tight layout
    plt.tight_layout()
    
    return fig

def graph_test(adata: sc.AnnData,
               neighbor_graph: str = "knn",
               reduction_method: str = "umap",
               k: int = 25,
               method: str = "moran_i",
               alternative: str = "greater",
               expression_family: str = "quasipoisson",
               cores: int = 1,
               verbose: bool = False,
               nn_control: Optional[Dict] = None) -> pd.DataFrame:
    """
    Test genes for differential expression based on the low dimensional
    embedding and the principal graph.
    
    We are often interested in finding genes that are differentially expressed 
    across a single-cell trajectory. This function introduces a new approach for 
    finding such genes that draws on a powerful technique in spatial correlation 
    analysis, the Moran's I test. Moran's I is a measure of multi-directional 
    and multi-dimensional spatial autocorrelation. The statistic tells you whether 
    cells at nearby positions on a trajectory will have similar (or dissimilar) 
    expression levels for the gene being tested. Although both Pearson correlation 
    and Moran's I ranges from -1 to 1, the interpretation of Moran's I is slightly 
    different: +1 means that nearby cells will have perfectly similar expression; 
    0 represents no correlation, and -1 means that neighboring cells will be 
    anti-correlated.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    neighbor_graph : str, optional (default: "knn")
        String indicating what neighbor graph to use.
        "principal_graph" and "knn" are supported. Default is "knn", but
        "principal_graph" is recommended for trajectory analysis.
    reduction_method : str, optional (default: "umap")
        The method used to reduce dimension.
        Currently only supported for "umap".
    k : int, optional (default: 25)
        Number of nearest neighbors used for building the kNN graph which
        is passed to knn2nb function during the Moran's I (Geary's C) test
        procedure.
    method : str, optional (default: "moran_i")
        A character string specifying the method (currently only
        'moran_i' is supported) for detecting significant genes showing
        correlation along the principal graph embedded in the low dimensional
        space.
    alternative : str, optional (default: "greater")
        A character string specifying the alternative hypothesis,
        must be one of greater (default), less or two.sided.
    expression_family : str, optional (default: "quasipoisson")
        A character string specifying the expression family
        function used for the test.
    cores : int, optional (default: 1)
        The number of cores to be used while testing each gene for
        differential expression.
    verbose : bool, optional (default: False)
        Whether to show spatial test (Moran's I) errors and warnings.
        Only valid for cores = 1.
    nn_control : dict, optional (default: None)
        An optional list of parameters used to make the nearest
        neighbor index. See the set_nn_control help for detailed information.
        
    Returns:
    --------
    pandas.DataFrame
        A data frame containing the p values and q-values from the Moran's I
        test on the parallel arrays of models.
    """
    # Validate parameters
    if neighbor_graph not in ["knn", "principal_graph"]:
        raise ValueError("neighbor_graph must be one of 'knn' or 'principal_graph'")
    
    if reduction_method.lower() != "umap":
        raise ValueError("reduction_method must be 'umap'")
    
    if method.lower() != "moran_i":
        raise ValueError("method must be 'moran_i'")
    
    if alternative not in ["greater", "less", "two.sided"]:
        raise ValueError("alternative must be one of 'greater', 'less', or 'two.sided'")
    
    if expression_family not in ["quasipoisson", "negbinomial", "poisson", "binomial"]:
        raise ValueError("expression_family must be one of 'quasipoisson', 'negbinomial', 'poisson', or 'binomial'")
    
    # Check if data has been reduced
    if reduction_method.lower() not in adata.obsm:
        raise ValueError(f"No dimensionality reduction for {reduction_method} calculated. "
                        f"Please run reduce_dimension with reduction_method = {reduction_method} "
                        f"before running graph_test")
    
    # Get reduced dimension matrix
    reduced_dim = adata.obsm[f'X_{reduction_method.lower()}']
    
    # Get expression data
    if scipy.sparse.issparse(adata.X):
        expr_data = adata.X.toarray()
    else:
        expr_data = adata.X.copy()
    
    # Create neighbor graph
    if neighbor_graph == "knn":
        # Create kNN graph
        nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto')
        nbrs.fit(reduced_dim)
        distances, indices = nbrs.kneighbors(reduced_dim)
        
        # Remove self from neighbors
        neighbor_matrix = indices[:, 1:]
        dist_matrix = distances[:, 1:]
        
        # Create spatial weights matrix
        spatial_weights = np.zeros((len(adata.obs_names), len(adata.obs_names)))
        for i in range(len(adata.obs_names)):
            for j, neighbor_idx in enumerate(neighbor_matrix[i]):
                spatial_weights[i, neighbor_idx] = 1.0 / (j + 1)  # Weight by rank
        
        # Normalize spatial weights
        row_sums = spatial_weights.sum(axis=1)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        spatial_weights = spatial_weights / row_sums[:, np.newaxis]
        
    else:  # principal_graph
        if "principal_graph" not in adata.uns:
            raise ValueError("No principal graph found. Please run learn_graph before running graph_test")
        
        # Get principal graph
        principal_graph = adata.uns["principal_graph"]
        
        # Create spatial weights matrix
        spatial_weights = np.zeros((len(adata.obs_names), len(adata.obs_names)))
        
        # For each cell, find its nearest node in the principal graph
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto')
        nbrs.fit(principal_graph["node_positions"])
        _, node_indices = nbrs.kneighbors(reduced_dim)
        node_indices = node_indices.flatten()
        
        # For each node, find cells that are closest to it
        for i in range(len(principal_graph["node_positions"])):
            cell_indices = np.where(node_indices == i)[0]
            
            # For each cell, connect it to other cells that are closest to the same node
            for cell_idx in cell_indices:
                for other_cell_idx in cell_indices:
                    if cell_idx != other_cell_idx:
                        # Calculate distance between cells
                        dist = np.linalg.norm(reduced_dim[cell_idx] - reduced_dim[other_cell_idx])
                        if dist > 0:
                            spatial_weights[cell_idx, other_cell_idx] = 1.0 / dist
        
        # Normalize spatial weights
        row_sums = spatial_weights.sum(axis=1)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        spatial_weights = spatial_weights / row_sums[:, np.newaxis]
    
    # Calculate Moran's I statistic for each gene
    n_genes = expr_data.shape[1]
    moran_i = np.zeros(n_genes)
    p_values = np.zeros(n_genes)
    z_scores = np.zeros(n_genes)
    
    # Function to calculate Moran's I for a single gene
    def calculate_moran_i(gene_expr):
        # Center the expression values
        centered_expr = gene_expr - np.mean(gene_expr)
        
        # Calculate Moran's I
        numerator = np.sum(spatial_weights * np.outer(centered_expr, centered_expr))
        denominator = np.sum(spatial_weights) * np.sum(centered_expr**2)
        
        if denominator == 0:
            return 0, 1, 0
        
        moran_i = (len(gene_expr) / (2 * np.sum(spatial_weights))) * (numerator / denominator)
        
        # Calculate expected value and variance under the null hypothesis
        n = len(gene_expr)
        s1 = np.sum(spatial_weights**2)
        s2 = np.sum((spatial_weights + spatial_weights.T)**2) / 2
        
        expected_i = -1 / (n - 1)
        
        # Calculate variance
        if n > 2:
            # Calculate moments of the expression data
            m2 = np.sum(centered_expr**2) / n
            m4 = np.sum(centered_expr**4) / n
            
            # Calculate variance
            var_i = ((n * ((n**2 - 3*n + 3)*s1 - n*s2 + 3*s1**2) - 
                     (k * ((n**2 - n)*s1 - 2*n*s2 + 3*s1**2)) / 4) / 
                    ((n-1)*(n-2)*(n-3)*s1**2) * (m4/m2**2) - 
                    ((n-1)**2 * expected_i**2) / ((n-1)*(n-2)*(n-3)))
        else:
            var_i = 0
        
        # Calculate z-score
        if var_i > 0:
            z_score = (moran_i - expected_i) / np.sqrt(var_i)
        else:
            z_score = 0
        
        # Calculate p-value
        if alternative == "greater":
            p_value = 1 - stats.norm.cdf(z_score)
        elif alternative == "less":
            p_value = stats.norm.cdf(z_score)
        else:  # two.sided
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        
        return moran_i, p_value, z_score
    
    # Calculate Moran's I for each gene
    if cores > 1:
        # Parallel processing
        from multiprocessing import Pool
        
        def process_gene(gene_idx):
            gene_expr = expr_data[:, gene_idx]
            return calculate_moran_i(gene_expr)
        
        with Pool(cores) as pool:
            results = pool.map(process_gene, range(n_genes))
        
        for i, (mi, pv, zs) in enumerate(results):
            moran_i[i] = mi
            p_values[i] = pv
            z_scores[i] = zs
    else:
        # Sequential processing
        for i in range(n_genes):
            gene_expr = expr_data[:, i]
            mi, pv, zs = calculate_moran_i(gene_expr)
            moran_i[i] = mi
            p_values[i] = pv
            z_scores[i] = zs
            
            if verbose and i % 100 == 0:
                print(f"Processed {i}/{n_genes} genes")
    
    # Calculate q-values (FDR-adjusted p-values)
    from statsmodels.stats.multitest import multipletests
    _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Create result dataframe
    result_df = pd.DataFrame({
        'gene_id': adata.var_names,
        'moran_i': moran_i,
        'p_value': p_values,
        'q_value': q_values,
        'z_score': z_scores
    })
    
    # Sort by p-value
    result_df = result_df.sort_values('p_value')
    
    return result_df

def find_gene_modules(adata: sc.AnnData,
                     reduction_method: str = "umap",
                     max_components: int = 2,
                     umap_metric: str = "cosine",
                     umap_min_dist: float = 0.1,
                     umap_n_neighbors: int = 15,
                     umap_fast_sgd: bool = False,
                     umap_nn_method: str = "annoy",
                     k: int = 20,
                     leiden_iter: int = 1,
                     partition_qval: float = 0.05,
                     weight: bool = False,
                     resolution: Optional[Union[float, List[float]]] = None,
                     random_seed: int = 0,
                     cores: int = 1,
                     verbose: bool = False,
                     preprocess_method: str = "pca",
                     nn_control: Optional[Dict] = None,
                     **kwargs) -> pd.DataFrame:
    """
    Cluster genes into modules that are co-expressed across cells.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    reduction_method : str, optional (default: "umap")
        The dimensionality reduction method used to generate the lower dimensional space 
        in which genes will be clustered. Currently only UMAP is supported.
    max_components : int, optional (default: 2)
        The number of dimensions in which to cluster genes into modules.
    umap_metric : str, optional (default: "cosine")
        Metric used by UMAP for measuring similarity between genes.
    umap_min_dist : float, optional (default: 0.1)
        Minimum distance parameter passed to UMAP.
    umap_n_neighbors : int, optional (default: 15)
        Number of nearest neighbors used by UMAP.
    umap_fast_sgd : bool, optional (default: False)
        Whether to allow UMAP to perform fast stochastic gradient descent.
        Setting False will result in slower, but deterministic behavior (if cores=1).
    umap_nn_method : str, optional (default: "annoy")
        The method used for nearest neighbor network construction during UMAP.
    k : int, optional (default: 20)
        Number of kNN used in creating the k nearest neighbor graph for Leiden clustering.
        The number of kNN is related to the resolution of the clustering result,
        bigger number of kNN gives low resolution and vice versa.
    leiden_iter : int, optional (default: 1)
        Integer number of iterations used for Leiden clustering.
        The clustering result with the largest modularity score is used as the final clustering result.
    partition_qval : float, optional (default: 0.05)
        Significance threshold used in Leiden community graph partitioning.
    weight : bool, optional (default: False)
        Whether to use Jaccard coefficient for two nearest neighbors as the weight
        used for Leiden clustering.
    resolution : float or list of float, optional (default: None)
        Resolution parameter passed to Leiden. Can be a list.
        If so, this method will evaluate modularity at each resolution and use the
        one with the highest value.
    random_seed : int, optional (default: 0)
        The seed used by the random number generator in Leiden.
    cores : int, optional (default: 1)
        Number of cores computer should use to execute function.
    verbose : bool, optional (default: False)
        Whether or not verbose output is printed.
    preprocess_method : str, optional (default: "pca")
        A string specifying the low-dimensional space to use for gene loadings,
        currently either "pca" or "lsi".
    nn_control : dict, optional (default: None)
        An optional list of parameters used to make the nearest neighbor index.
    **kwargs : dict
        Additional arguments passed to UMAP and Leiden analysis.
        
    Returns:
    --------
    pandas.DataFrame
        A dataframe with genes and the modules to which they are assigned.
    """
    # Validate parameters
    if reduction_method.lower() != "umap":
        raise ValueError("reduction_method must be 'umap'")
    
    if preprocess_method.lower() not in ["pca", "lsi"]:
        raise ValueError("preprocess_method must be one of 'pca' or 'lsi'")
    
    # Get expression data
    if scipy.sparse.issparse(adata.X):
        expr_data = adata.X.toarray()
    else:
        expr_data = adata.X.copy()
    
    # Transpose expression matrix to get genes as rows
    gene_expr = expr_data.T
    
    # Perform dimensionality reduction on genes
    if preprocess_method.lower() == "pca":
        # Perform PCA
        pca = TruncatedSVD(n_components=min(max_components, gene_expr.shape[1]-1))
        reduced_dim = pca.fit_transform(gene_expr)
    else:  # LSI
        # Calculate TF-IDF
        tfidf_result = tfidf(gene_expr)
        tf_idf_counts = tfidf_result['tf_idf_counts']
        
        # Perform SVD
        svd = TruncatedSVD(n_components=min(max_components, gene_expr.shape[1]-1))
        reduced_dim = svd.fit_transform(tf_idf_counts)
    
    # Perform UMAP
    reducer = umap.UMAP(
        n_components=max_components,
        metric=umap_metric,
        min_dist=umap_min_dist,
        n_neighbors=umap_n_neighbors,
        low_memory=True if umap_fast_sgd else False,
        n_jobs=cores,
        verbose=verbose,
        **kwargs
    )
    
    # Fit and transform
    gene_embeddings = reducer.fit_transform(reduced_dim)
    
    # Create k-nearest neighbors graph
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto')
    nbrs.fit(gene_embeddings)
    distances, indices = nbrs.kneighbors(gene_embeddings)
    
    # Remove self from neighbors
    neighbor_matrix = indices[:, 1:]
    dist_matrix = distances[:, 1:]
    
    # Create graph for Leiden clustering
    if weight:
        # Compute Jaccard coefficients
        links = jaccard_coeff(neighbor_matrix, weight=True)
    else:
        # Create unweighted links
        links = []
        for i in range(len(neighbor_matrix)):
            for j in neighbor_matrix[i]:
                links.append([i, j, 1.0])
        links = np.array(links)
    
    # Filter links
    links = links[links[:, 2] > 0]
    
    # Create relations dataframe
    relations = pd.DataFrame(links, columns=['from', 'to', 'weight'])
    
    # Create graph
    g = nx.Graph()
    for _, row in relations.iterrows():
        g.add_edge(int(row['from']), int(row['to']), weight=row['weight'])
    
    # Convert NetworkX graph to igraph
    g_igraph = igraph.Graph()
    g_igraph.add_vertices(list(g.nodes()))
    g_igraph.add_edges([(u, v) for u, v in g.edges()])
    g_igraph.es['weight'] = [d['weight'] for u, v, d in g.edges(data=True)]
    
    # Set random seed
    if random_seed is not None:
        random.seed(random_seed)
    
    # Run Leiden clustering
    if verbose:
        print("Running Leiden clustering...")
    
    # Handle resolution parameter
    if resolution is None:
        resolution = [0.0001]  # Default resolution
    elif not isinstance(resolution, (list, np.ndarray)):
        resolution = [resolution]
    
    best_modularity = -1
    best_result = None
    
    for cur_resolution in resolution:
        # Run Leiden algorithm
        partition = leidenalg.find_partition(
            g_igraph,
            leidenalg.CPMVertexPartition,
            resolution_parameter=cur_resolution,
            weights='weight',
            n_iterations=leiden_iter
        )
        
        if partition.modularity > best_modularity:
            best_result = partition
            best_modularity = partition.modularity
    
    # Create result dataframe
    result_df = pd.DataFrame({
        'gene_id': adata.var_names,
        'module': best_result.membership
    })
    
    # Sort by module
    result_df = result_df.sort_values('module')
    
    return result_df

def aggregate_gene_expression(adata: sc.AnnData,
                            gene_group_df: Optional[pd.DataFrame] = None,
                            cell_group_df: Optional[pd.DataFrame] = None,
                            norm_method: str = "log",
                            pseudocount: float = 1,
                            scale_agg_values: bool = True,
                            max_agg_value: float = 3,
                            min_agg_value: float = -3,
                            exclude_na: bool = True,
                            gene_agg_fun: str = "sum",
                            cell_agg_fun: str = "mean") -> np.ndarray:
    """
    Creates a matrix with aggregated expression values for arbitrary groups of genes.
    
    Parameters:
    -----------
    adata : scanpy.AnnData
        The annotated data matrix
    gene_group_df : pandas.DataFrame, optional (default: None)
        A dataframe in which the first column contains gene ids or short gene names 
        and the second contains groups. If None, genes are not grouped.
    cell_group_df : pandas.DataFrame, optional (default: None)
        A dataframe in which the first column contains cell ids and the second 
        contains groups. If None, cells are not grouped.
    norm_method : str, optional (default: "log")
        How to transform gene expression values before aggregating them.
        Options are "log", "binary", or "size_only".
    pseudocount : float, optional (default: 1)
        Value to add to expression prior to log transformation and aggregation.
    scale_agg_values : bool, optional (default: True)
        Whether to center and scale aggregated groups of genes.
    max_agg_value : float, optional (default: 3)
        If scale_agg_values is True, the maximum value the resulting Z scores can take.
    min_agg_value : float, optional (default: -3)
        If scale_agg_values is True, the minimum value the resulting Z scores can take.
    exclude_na : bool, optional (default: True)
        Whether to exclude NA values from the aggregated matrix.
    gene_agg_fun : str, optional (default: "sum")
        Function used for gene aggregation. Can be either "sum" or "mean".
    cell_agg_fun : str, optional (default: "mean")
        Function used for cell aggregation.
        
    Returns:
    --------
    numpy.ndarray
        A matrix of dimension NxM, where N is the number of gene groups and
        M is the number of cell groups.
    """
    # Get expression data
    if scipy.sparse.issparse(adata.X):
        expr_data = adata.X.toarray()
    else:
        expr_data = adata.X.copy()
    
    # Normalize expression data
    if norm_method == "log":
        expr_data = np.log2(expr_data + pseudocount)
    elif norm_method == "size_only":
        if 'size_factors' in adata.obs:
            expr_data = expr_data / adata.obs['size_factors'].values[:, np.newaxis]
        else:
            raise ValueError("Size factors not found in adata.obs")
    elif norm_method == "binary":
        expr_data = (expr_data > 0).astype(float)
    else:
        raise ValueError("norm_method must be one of 'log', 'binary', or 'size_only'")
    
    # Create gene groups if not provided
    if gene_group_df is None:
        gene_group_df = pd.DataFrame({
            'gene_id': adata.var_names,
            'group': range(len(adata.var_names))
        })
    
    # Create cell groups if not provided
    if cell_group_df is None:
        cell_group_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'group': range(len(adata.obs_names))
        })
    
    # Validate gene groups
    if not all(gene_group_df.iloc[:, 0].isin(adata.var_names)):
        raise ValueError("Some gene IDs in gene_group_df are not in adata.var_names")
    
    # Validate cell groups
    if not all(cell_group_df.iloc[:, 0].isin(adata.obs_names)):
        raise ValueError("Some cell IDs in cell_group_df are not in adata.obs_names")
    
    # Get unique groups
    gene_groups = gene_group_df.iloc[:, 1].unique()
    cell_groups = cell_group_df.iloc[:, 1].unique()
    
    # Initialize aggregated matrix
    agg_matrix = np.zeros((len(gene_groups), len(cell_groups)))
    
    # Aggregate expression values
    for i, gene_group in enumerate(gene_groups):
        # Get genes in this group
        group_genes = gene_group_df[gene_group_df.iloc[:, 1] == gene_group].iloc[:, 0]
        gene_indices = [list(adata.var_names).index(gene) for gene in group_genes]
        
        # Get expression values for these genes
        group_expr = expr_data[:, gene_indices]
        
        # Aggregate genes
        if gene_agg_fun == "sum":
            gene_agg = np.sum(group_expr, axis=1)
        elif gene_agg_fun == "mean":
            gene_agg = np.mean(group_expr, axis=1)
        else:
            raise ValueError("gene_agg_fun must be either 'sum' or 'mean'")
        
        for j, cell_group in enumerate(cell_groups):
            # Get cells in this group
            group_cells = cell_group_df[cell_group_df.iloc[:, 1] == cell_group].iloc[:, 0]
            cell_indices = [list(adata.obs_names).index(cell) for cell in group_cells]
            
            # Get aggregated values for these cells
            group_values = gene_agg[cell_indices]
            
            # Handle NA values
            if exclude_na:
                group_values = group_values[~np.isnan(group_values)]
            
            # Aggregate cells
            if cell_agg_fun == "mean":
                agg_matrix[i, j] = np.mean(group_values)
            elif cell_agg_fun == "sum":
                agg_matrix[i, j] = np.sum(group_values)
            else:
                raise ValueError("cell_agg_fun must be either 'mean' or 'sum'")
    
    # Scale aggregated values if requested
    if scale_agg_values:
        # Center and scale
        agg_matrix = (agg_matrix - np.mean(agg_matrix)) / np.std(agg_matrix)
        
        # Clip values
        agg_matrix = np.clip(agg_matrix, min_agg_value, max_agg_value)
    
    return agg_matrix

def learn_graph(cds: sc.AnnData,
                use_partition: bool = True,
                close_loop: bool = True,
                learn_graph_control: Optional[Dict] = None,
                euclidean_distance_ratio: float = 1,
                geodesic_distance_ratio: float = 1/2,
                minimal_branch_len: int = 10,
                orthogonal_proj_tip: bool = False,
                prune_graph: bool = True,
                scale: bool = False,
                ncenter: Optional[int] = None,
                num_paths: Optional[int] = None,
                tol: float = 1e-4,
                maxiter: int = 10,
                nn_method: str = "auto",
                k: int = 25,
                epsilon: float = 0,
                use_weights: bool = False,
                n_jobs: int = -1,
                verbose: bool = False) -> sc.AnnData:
    """
    Learn the principal graph from reduced dimension data.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set object
    use_partition : bool, optional (default: True)
        Whether to use partition information
    close_loop : bool, optional (default: True)
        Whether to close loops in the graph
    learn_graph_control : dict, optional (default: None)
        Control parameters for graph learning
    euclidean_distance_ratio : float, optional (default: 1)
        Ratio for Euclidean distance
    geodesic_distance_ratio : float, optional (default: 1/2)
        Ratio for geodesic distance
    minimal_branch_len : int, optional (default: 10)
        Minimum length for branches
    orthogonal_proj_tip : bool, optional (default: False)
        Whether to use orthogonal projection for tips
    prune_graph : bool, optional (default: True)
        Whether to prune the graph
    scale : bool, optional (default: False)
        Whether to scale the data
    ncenter : int, optional (default: None)
        Number of centers
    num_paths : int, optional (default: None)
        Number of paths
    tol : float, optional (default: 1e-4)
        Tolerance for convergence
    maxiter : int, optional (default: 10)
        Maximum number of iterations
    nn_method : str, optional (default: "auto")
        Method for nearest neighbor search
    k : int, optional (default: 25)
        Number of nearest neighbors
    epsilon : float, optional (default: 0)
        Epsilon parameter for graph learning
    use_weights : bool, optional (default: False)
        Whether to use weights
    n_jobs : int, optional (default: -1)
        Number of parallel jobs
    verbose : bool, optional (default: False)
        Whether to print progress messages
        
    Returns:
    --------
    scanpy.AnnData
        Updated cell data set with learned graph
    """
    # Check if reduced dimension data exists
    if "X_pca" not in cds.obsm and "X_umap" not in cds.obsm:
        raise ValueError("No reduced dimension data found. Please run reduce_dimension first.")
    
    # Get reduced dimension data
    if "X_umap" in cds.obsm:
        reduced_dim_data = cds.obsm["X_umap"]
    else:
        reduced_dim_data = cds.obsm["X_pca"]
    
    if verbose:
        print("Learning principal graph...")
    
    # Build k-nearest neighbors graph
    if verbose:
        print("Building k-nearest neighbors graph...")
    
    nbrs = NearestNeighbors(n_neighbors=k, algorithm=nn_method, n_jobs=n_jobs)
    nbrs.fit(reduced_dim_data)
    
    # Get distances and indices
    distances, indices = nbrs.kneighbors()
    
    # Create graph
    if verbose:
        print("Creating graph...")
    
    g = nx.Graph()
    
    # Add nodes
    for i in range(len(reduced_dim_data)):
        g.add_node(i, coord=reduced_dim_data[i])
    
    # Add edges
    for i in range(len(indices)):
        for j, dist in zip(indices[i][1:], distances[i][1:]):  # Skip self-connection
            if use_weights:
                g.add_edge(i, j, weight=np.exp(-dist * epsilon))
            else:
                g.add_edge(i, j)
    
    # Learn principal graph
    if verbose:
        print("Learning principal graph structure...")
    
    # Initialize centers if not provided
    if ncenter is None:
        ncenter = int(np.sqrt(len(reduced_dim_data)))
    
    # Initialize paths if not provided
    if num_paths is None:
        num_paths = ncenter
    
    # Convert to minimum spanning tree
    if prune_graph:
        if verbose:
            print("Converting to minimum spanning tree...")
        g = nx.minimum_spanning_tree(g)
    
    # Close loops if requested
    if close_loop:
        if verbose:
            print("Closing loops...")
        
        # Find nodes with degree 1 (tips)
        tips = [n for n in g.nodes() if g.degree(n) == 1]
        
        # Connect tips if they are close enough
        for i in range(len(tips)):
            for j in range(i+1, len(tips)):
                tip1, tip2 = tips[i], tips[j]
                dist = np.linalg.norm(reduced_dim_data[tip1] - reduced_dim_data[tip2])
                
                if dist < np.mean(distances) * euclidean_distance_ratio:
                    # Check if adding this edge creates a small cycle
                    path_length = nx.shortest_path_length(g, tip1, tip2)
                    if path_length > minimal_branch_len:
                        g.add_edge(tip1, tip2, weight=dist)
    
    # Store graph in AnnData object
    cds.uns['principal_graph'] = {
        'graph': g,
        'params': {
            'k': k,
            'epsilon': epsilon,
            'use_weights': use_weights,
            'ncenter': ncenter,
            'num_paths': num_paths
        }
    }
    
    return cds

def order_cells(cds: sc.AnnData,
                reduction_method: str = "UMAP",
                num_paths: Optional[int] = None,
                tol: float = 1e-4,
                maxiter: int = 10,
                fast_mode: bool = True,
                verbose: bool = False) -> sc.AnnData:
    """
    Order cells along the principal graph.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set object
    reduction_method : str, optional (default: "UMAP")
        The reduction method to use
    num_paths : int, optional (default: None)
        The number of paths to use
    tol : float, optional (default: 1e-4)
        The tolerance for convergence
    maxiter : int, optional (default: 10)
        The maximum number of iterations
    fast_mode : bool, optional (default: True)
        Whether to use a faster but less accurate ordering method
    verbose : bool, optional (default: False)
        Whether to print verbose output
        
    Returns:
    --------
    scanpy.AnnData
        Updated cell data set with cell ordering
    """
    if verbose:
        print("Ordering cells...")
    
    # Check if principal graph exists
    if 'principal_graph' not in cds.uns:
        raise ValueError("Principal graph not found. Please run learn_graph first.")
    
    # Get the graph from the stored dictionary
    g = cds.uns['principal_graph']['graph']
    
    # Use fast mode for small graphs or when fast_mode is True
    if fast_mode or len(g.nodes()) < 100:
        ordering = fast_cell_ordering(g, num_paths=num_paths)
    else:
        # Use approximate ordering for large graphs
        ordering = approximate_cell_ordering(g, num_paths=num_paths)
    
    # Store ordering in AnnData object
    cds.obs['cell_order'] = ordering
    
    if verbose:
        print("Cell ordering completed.")
    
    return cds

def fast_cell_ordering(g: nx.Graph, num_paths: Optional[int] = None) -> np.ndarray:
    """
    Fast cell ordering for small graphs.
    
    Parameters:
    -----------
    g : nx.Graph
        The principal graph
    num_paths : int, optional (default: None)
        The number of paths to use
        
    Returns:
    --------
    np.ndarray
        The cell ordering
    """
    # Initialize centers if not provided
    if num_paths is None:
        num_paths = min(3, int(np.sqrt(len(g.nodes()))))
    
    # Find terminal nodes (nodes with degree 1)
    terminal_nodes = [n for n in g.nodes() if g.degree(n) == 1]
    
    # If no terminal nodes, use all nodes
    if not terminal_nodes:
        terminal_nodes = list(g.nodes())
    
    # Select a subset of terminal nodes as centers
    if len(terminal_nodes) > num_paths:
        centers = np.random.choice(terminal_nodes, num_paths, replace=False)
    else:
        centers = terminal_nodes
    
    # Calculate shortest paths from centers to all nodes
    distances = np.zeros((len(g.nodes()), len(centers)))
    for i, center in enumerate(centers):
        paths = nx.single_source_shortest_path_length(g, center)
        for node, dist in paths.items():
            distances[node, i] = dist
    
    # Order cells based on their distances to centers
    ordering = np.argsort(np.min(distances, axis=1))
    
    return ordering

def root_nodes(cds: sc.AnnData, reduction_method: str = "UMAP") -> List[int]:
    """
    Find root nodes in the principal graph.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set object
    reduction_method : str, optional (default: "UMAP")
        The reduction method used
        
    Returns:
    --------
    list of int
        List of root node indices
    """
    # Check if principal graph exists
    if 'principal_graph' not in cds.uns:
        raise ValueError("Principal graph not found. Please run learn_graph first.")
    
    # Get the graph from the stored dictionary
    g = cds.uns['principal_graph']['graph']
    
    # Find nodes with degree 1 (potential roots)
    potential_roots = [n for n in g.nodes() if g.degree(n) == 1]
    
    # If no potential roots, use nodes with minimum degree
    if not potential_roots:
        min_degree = min(g.degree(n) for n in g.nodes())
        potential_roots = [n for n in g.nodes() if g.degree(n) == min_degree]
    
    # If multiple potential roots, select the one with the most descendants
    if len(potential_roots) > 1:
        # Count descendants for each potential root
        descendant_counts = {}
        for root in potential_roots:
            descendants = set()
            for node in g.nodes():
                if nx.has_path(g, root, node):
                    descendants.add(node)
            descendant_counts[root] = len(descendants)
        
        # Select the root with the most descendants
        root = max(descendant_counts, key=descendant_counts.get)
        return [root]
    
    return potential_roots

def leaf_nodes(cds: sc.AnnData, reduction_method: str = "UMAP") -> List[int]:
    """
    Find leaf nodes in the principal graph.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set object
    reduction_method : str, optional (default: "UMAP")
        The reduction method used
        
    Returns:
    --------
    list of int
        List of leaf node indices
    """
    # Check if principal graph exists
    if 'principal_graph' not in cds.uns:
        raise ValueError("Principal graph not found. Please run learn_graph first.")
    
    # Get the graph from the stored dictionary
    g = cds.uns['principal_graph']['graph']
    
    # Find nodes with degree 1 (leaves)
    leaves = [n for n in g.nodes() if g.degree(n) == 1]
    
    # Exclude root nodes
    roots = root_nodes(cds, reduction_method)
    leaves = [leaf for leaf in leaves if leaf not in roots]
    
    return leaves

def branch_nodes(cds: sc.AnnData, reduction_method: str = "UMAP") -> List[int]:
    """
    Find branch nodes in the principal graph.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set object
    reduction_method : str, optional (default: "UMAP")
        The reduction method used
        
    Returns:
    --------
    list of int
        List of branch node indices
    """
    # Check if principal graph exists
    if 'principal_graph' not in cds.uns:
        raise ValueError("Principal graph not found. Please run learn_graph first.")
    
    # Get the graph from the stored dictionary
    g = cds.uns['principal_graph']['graph']
    
    # Find nodes with degree > 2 (branches)
    branch_points = [node for node in g.nodes() if g.degree(node) > 2]
    
    # Exclude root nodes
    roots = root_nodes(cds, reduction_method)
    branch_points = [node for node in branch_points if node not in roots]
    
    return branch_points

def approximate_cell_ordering(g: nx.Graph, num_paths: Optional[int] = None) -> np.ndarray:
    """
    Approximate cell ordering for large graphs.
    
    Parameters:
    -----------
    g : nx.Graph
        The principal graph
    num_paths : int, optional (default: None)
        The number of paths to use
        
    Returns:
    --------
    np.ndarray
        The cell ordering
    """
    # Initialize centers if not provided
    if num_paths is None:
        num_paths = int(np.sqrt(len(g.nodes())))
    
    # Find terminal nodes (nodes with degree 1)
    terminal_nodes = [n for n in g.nodes() if g.degree(n) == 1]
    
    # If no terminal nodes, use all nodes
    if not terminal_nodes:
        terminal_nodes = list(g.nodes())
    
    # Select a subset of terminal nodes as centers
    if len(terminal_nodes) > num_paths:
        centers = np.random.choice(terminal_nodes, num_paths, replace=False)
    else:
        centers = terminal_nodes
    
    # Calculate shortest paths from centers to all nodes
    distances = np.zeros((len(g.nodes()), len(centers)))
    for i, center in enumerate(centers):
        paths = nx.single_source_shortest_path_length(g, center)
        for node, dist in paths.items():
            distances[node, i] = dist
    
    # Order cells based on their distances to centers
    ordering = np.argsort(np.min(distances, axis=1))
    
    return ordering

def exact_cell_ordering(g: nx.Graph,
                       num_paths: Optional[int] = None,
                       tol: float = 1e-4,
                       maxiter: int = 10) -> np.ndarray:
    """
    Exact cell ordering for small graphs.
    
    Parameters:
    -----------
    g : nx.Graph
        The principal graph
    num_paths : int, optional (default: None)
        The number of paths to use
    tol : float, optional (default: 1e-4)
        The tolerance for convergence
    maxiter : int, optional (default: 10)
        The maximum number of iterations
        
    Returns:
    --------
    np.ndarray
        The cell ordering
    """
    # Initialize centers if not provided
    if num_paths is None:
        num_paths = int(np.sqrt(len(g.nodes())))
    
    # Find terminal nodes (nodes with degree 1)
    terminal_nodes = [n for n in g.nodes() if g.degree(n) == 1]
    
    # If no terminal nodes, use all nodes
    if not terminal_nodes:
        terminal_nodes = list(g.nodes())
    
    # Select a subset of terminal nodes as centers
    if len(terminal_nodes) > num_paths:
        centers = np.random.choice(terminal_nodes, num_paths, replace=False)
    else:
        centers = terminal_nodes
    
    # Initialize cell ordering
    ordering = np.zeros(len(g.nodes()), dtype=int)
    ordered_nodes = set()
    
    # Order cells starting from each center
    for center in centers:
        if center in ordered_nodes:
            continue
        
        # Get nodes in the current path
        path_nodes = set()
        current = center
        while current is not None:
            path_nodes.add(current)
            # Find the next node in the path
            neighbors = list(g.neighbors(current))
            if not neighbors:
                break
            # Choose the neighbor that leads to the most unordered nodes
            best_neighbor = None
            best_score = -1
            for neighbor in neighbors:
                if neighbor in ordered_nodes:
                    continue
                # Count unordered nodes reachable from this neighbor
                score = len(set(nx.descendants(g, neighbor)) - ordered_nodes)
                if score > best_score:
                    best_score = score
                    best_neighbor = neighbor
            current = best_neighbor
        
        # Add nodes in this path to the ordering
        path_nodes = list(path_nodes)
        path_nodes.sort(key=lambda x: len(set(nx.descendants(g, x)) - ordered_nodes),
                       reverse=True)
        for node in path_nodes:
            if node not in ordered_nodes:
                ordering[len(ordered_nodes)] = node
                ordered_nodes.add(node)
    
    # Add any remaining unordered nodes
    remaining_nodes = set(g.nodes()) - ordered_nodes
    for node in remaining_nodes:
        ordering[len(ordered_nodes)] = node
        ordered_nodes.add(node)
    
    return ordering

def plot_genes_in_pseudotime(cds_subset: sc.AnnData,
                            min_expr: Optional[float] = None,
                            cell_size: float = 0.75,
                            nrow: Optional[int] = None,
                            ncol: int = 1,
                            panel_order: Optional[List[str]] = None,
                            color_cells_by: str = "pseudotime",
                            trend_formula: str = "~ splines::ns(pseudotime, df=3)",
                            label_by_short_name: bool = True,
                            vertical_jitter: Optional[float] = None,
                            horizontal_jitter: Optional[float] = None) -> plt.Figure:
    """
    Plots expression for one or more genes as a function of pseudotime.
    
    Parameters:
    -----------
    cds_subset : scanpy.AnnData
        Subset cell dataset including only the genes to be plotted
    min_expr : float, optional (default: None)
        The minimum (untransformed) expression level to plot
    cell_size : float, optional (default: 0.75)
        The size (in points) of each cell used in the plot
    nrow : int, optional (default: None)
        The number of rows used when laying out the panels for each gene's expression
    ncol : int, optional (default: 1)
        The number of columns used when laying out the panels for each gene's expression
    panel_order : list of str, optional (default: None)
        Vector of gene names indicating the order in which genes should be laid out
        (left-to-right, top-to-bottom). If label_by_short_name = True, use gene_short_name
        values, otherwise use feature IDs
    color_cells_by : str, optional (default: "pseudotime")
        The cell attribute to be used to color each cell
    trend_formula : str, optional (default: "~ splines::ns(pseudotime, df=3)")
        The model formula to be used for fitting the expression trend over pseudotime
    label_by_short_name : bool, optional (default: True)
        Label figure panels by gene_short_name (True) or feature ID (False)
    vertical_jitter : float, optional (default: None)
        A value passed to ggplot to jitter the points in the vertical dimension
    horizontal_jitter : float, optional (default: None)
        A value passed to ggplot to jitter the points in the horizontal dimension
        
    Returns:
    --------
    matplotlib.figure.Figure
        A matplotlib figure object containing the plot
    """
    # Check if pseudotime exists
    if 'pseudotime' not in cds_subset.obs:
        raise ValueError("No pseudotime information found. Please run order_cells first.")
    
    # Get gene names
    if label_by_short_name:
        gene_names = cds_subset.var['gene_short_name']
    else:
        gene_names = cds_subset.var_names
    
    # Filter genes based on minimum expression if specified
    if min_expr is not None:
        expr_matrix = cds_subset.X
        if scipy.sparse.issparse(expr_matrix):
            expr_matrix = expr_matrix.toarray()
        
        # Calculate mean expression for each gene
        mean_expr = np.mean(expr_matrix, axis=0)
        valid_genes = mean_expr >= min_expr
        gene_names = gene_names[valid_genes]
    
    # Set panel order if specified
    if panel_order is not None:
        if label_by_short_name:
            # Convert panel_order to gene indices
            gene_indices = [np.where(gene_names == name)[0][0] for name in panel_order]
        else:
            gene_indices = [np.where(cds_subset.var_names == name)[0][0] for name in panel_order]
    else:
        gene_indices = range(len(gene_names))
    
    # Calculate number of rows if not specified
    if nrow is None:
        nrow = int(np.ceil(len(gene_indices) / ncol))
    
    # Create figure and subplots
    fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow))
    if nrow == 1 and ncol == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot each gene
    for i, gene_idx in enumerate(gene_indices):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Get expression data for the gene
        if label_by_short_name:
            gene_name = gene_names[gene_idx]
            gene_idx = np.where(cds_subset.var['gene_short_name'] == gene_name)[0][0]
        else:
            gene_name = cds_subset.var_names[gene_idx]
        
        expr = cds_subset.X[:, gene_idx]
        if scipy.sparse.issparse(expr):
            expr = expr.toarray().flatten()
        
        # Get pseudotime values
        pseudotime = cds_subset.obs['pseudotime']
        
        # Create scatter plot
        if vertical_jitter is not None or horizontal_jitter is not None:
            x = pseudotime + np.random.normal(0, horizontal_jitter or 0, len(pseudotime))
            y = expr + np.random.normal(0, vertical_jitter or 0, len(expr))
        else:
            x = pseudotime
            y = expr
        
        scatter = ax.scatter(x, y, c=cds_subset.obs[color_cells_by],
                           s=cell_size, alpha=0.7)
        
        # Fit trend line
        if trend_formula is not None:
            # Create design matrix for spline
            if "splines::ns" in trend_formula:
                df = int(trend_formula.split("df=")[1].split(")")[0])
                spline_basis = splines.natural_spline(pseudotime, df=df)
                X = np.column_stack([np.ones_like(pseudotime), spline_basis])
            else:
                X = np.column_stack([np.ones_like(pseudotime), pseudotime])
            
            # Fit linear model
            coeffs = np.linalg.lstsq(X, expr, rcond=None)[0]
            
            # Generate trend line
            x_trend = np.linspace(min(pseudotime), max(pseudotime), 100)
            if "splines::ns" in trend_formula:
                spline_basis_trend = splines.natural_spline(x_trend, df=df)
                X_trend = np.column_stack([np.ones_like(x_trend), spline_basis_trend])
            else:
                X_trend = np.column_stack([np.ones_like(x_trend), x_trend])
            
            y_trend = X_trend @ coeffs
            
            # Plot trend line
            ax.plot(x_trend, y_trend, 'r-', linewidth=2)
        
        # Set labels and title
        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('Expression')
        ax.set_title(gene_name)
        
        # Add colorbar if color_cells_by is not pseudotime
        if color_cells_by != "pseudotime":
            plt.colorbar(scatter, ax=ax, label=color_cells_by)
    
    # Hide empty subplots
    for i in range(len(gene_indices), len(axes)):
        axes[i].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig

def align_cds(cds: sc.AnnData,
              preprocess_method: str = "PCA",
              alignment_group: Optional[str] = None,
              alignment_k: int = 20,
              residual_model_formula_str: Optional[str] = None,
              verbose: bool = False,
              build_nn_index: bool = False,
              nn_control: Optional[Dict] = None,
              **kwargs) -> sc.AnnData:
    """
    Align cells from different groups within a cell data set.
    
    Parameters:
    -----------
    cds : scanpy.AnnData
        The cell data set upon which to perform alignment
    preprocess_method : str, optional (default: "PCA")
        The low-dimensional space in which to perform alignment ("PCA" or "LSI")
    alignment_group : str, optional (default: None)
        Column name in cds.obs to use for aligning groups of cells
    alignment_k : int, optional (default: 20)
        The value of k used in mutual nearest neighbor alignment
    residual_model_formula_str : str, optional (default: None)
        Model formula for subtracting effects before dimensionality reduction
    verbose : bool, optional (default: False)
        Whether to emit verbose output during dimensionality reduction
    build_nn_index : bool, optional (default: False)
        Whether to build nearest neighbor index from aligned reduced matrix
    nn_control : dict, optional (default: None)
        Parameters for nearest neighbor index construction
    **kwargs : dict
        Additional arguments to pass to linear model fitting
        
    Returns:
    --------
    scanpy.AnnData
        Updated cell data set with aligned coordinates
    """
    # Validate input parameters
    if preprocess_method not in ["PCA", "LSI"]:
        raise ValueError("preprocess_method must be either 'PCA' or 'LSI'")
    
    if alignment_group is not None and alignment_group not in cds.obs.columns:
        raise ValueError(f"alignment_group '{alignment_group}' not found in cell metadata")
    
    # Get expression matrix
    if scipy.sparse.issparse(cds.X):
        expr_matrix = cds.X.toarray()
    else:
        expr_matrix = cds.X.copy()
    
    # Handle residual model if specified
    if residual_model_formula_str is not None:
        if verbose:
            print("Fitting and subtracting residual model...")
            
        # Create design matrix from formula
        model_matrix = pd.DataFrame(index=cds.obs.index)
        terms = residual_model_formula_str.replace("~", "").strip().split("+")
        for term in terms:
            term = term.strip()
            if term in cds.obs.columns:
                # Convert to numeric if possible
                try:
                    model_matrix[term] = pd.to_numeric(cds.obs[term])
                except:
                    model_matrix[term] = cds.obs[term]
        
        # Remove any columns with missing values
        model_matrix = model_matrix.dropna(axis=1)
        
        if model_matrix.empty:
            print("Warning: No valid terms found in residual model formula. Skipping residual model.")
        else:
            # Add intercept
            model_matrix['intercept'] = 1.0
            
            # Fit linear model and get residuals
            try:
                model = sm.OLS(expr_matrix, model_matrix)
                results = model.fit()
                expr_matrix = results.resid
            except Exception as e:
                print(f"Warning: Error fitting residual model: {str(e)}. Skipping residual model.")
    
    # Get number of components from preprocess_params if available
    n_components = min(50, expr_matrix.shape[1])  # Default value, but not more than features
    if 'preprocess_params' in cds.uns and 'num_dim' in cds.uns['preprocess_params']:
        n_components = min(cds.uns['preprocess_params']['num_dim'], expr_matrix.shape[1])
    
    # Perform dimensionality reduction
    if preprocess_method == "PCA":
        if verbose:
            print("Performing PCA...")
        pca = TruncatedSVD(n_components=n_components)
        reduced_matrix = pca.fit_transform(expr_matrix)
    else:  # LSI
        if verbose:
            print("Performing LSI...")
        # TF-IDF transformation
        tfidf_result = tfidf(expr_matrix)
        tf_idf_matrix = tfidf_result['tf_idf_counts']
        # SVD on TF-IDF matrix
        u, s, vt = svds(tf_idf_matrix, k=n_components)
        reduced_matrix = u * s[:, np.newaxis]
    
    # Perform alignment if alignment_group is specified
    if alignment_group is not None:
        if verbose:
            print(f"Aligning cells using {alignment_group}...")
            
        # Encode group labels
        le = LabelEncoder()
        group_labels = le.fit_transform(cds.obs[alignment_group])
        unique_groups = le.classes_
        
        # Initialize aligned coordinates matrix
        aligned_coords = np.zeros_like(reduced_matrix)
        
        # Perform mutual nearest neighbor alignment
        for i, group1 in enumerate(unique_groups):
            group1_mask = (group_labels == i)
            group1_coords = reduced_matrix[group1_mask]
            
            # Find mutual nearest neighbors with other groups
            for j, group2 in enumerate(unique_groups):
                if i == j:
                    continue
                    
                group2_mask = (group_labels == j)
                group2_coords = reduced_matrix[group2_mask]
                
                # Build k-nearest neighbors for both groups
                knn1 = NearestNeighbors(n_neighbors=min(alignment_k, len(group1_coords)))
                knn2 = NearestNeighbors(n_neighbors=min(alignment_k, len(group2_coords)))
                
                knn1.fit(group1_coords)
                knn2.fit(group2_coords)
                
                # Find mutual nearest neighbors
                dist1, idx1 = knn1.kneighbors(group2_coords)
                dist2, idx2 = knn2.kneighbors(group1_coords)
                
                # Compute correction vectors
                mutual_pairs = []
                for g2_idx in range(len(group2_coords)):
                    g1_neighbors = idx1[g2_idx]
                    for g1_idx in g1_neighbors:
                        if g2_idx in idx2[g1_idx]:
                            mutual_pairs.append((g1_idx, g2_idx))
                
                if mutual_pairs:
                    mutual_pairs = np.array(mutual_pairs)
                    correction_vectors = group2_coords[mutual_pairs[:, 1]] - group1_coords[mutual_pairs[:, 0]]
                    mean_correction = np.mean(correction_vectors, axis=0)
                    
                    # Apply correction to group1 coordinates
                    aligned_coords[group1_mask] = group1_coords + mean_correction
                else:
                    # If no mutual pairs found, keep original coordinates
                    aligned_coords[group1_mask] = group1_coords
        
        # Update reduced matrix with aligned coordinates
        reduced_matrix = aligned_coords
    
    # Store aligned coordinates in AnnData object
    cds.obsm[f"X_{preprocess_method.lower()}_aligned"] = reduced_matrix
    
    # Build nearest neighbor index if requested
    if build_nn_index:
        if verbose:
            print("Building nearest neighbor index...")
            
        nn_params = nn_control if nn_control is not None else {}
        nn = NearestNeighbors(**nn_params)
        nn.fit(reduced_matrix)
        cds.uns['neighbor_index'] = nn
    
    return cds

if __name__ == "__main__":
    # Example usage of the implemented functions
    
    # 1. Create a new cell dataset
    print("Creating a new cell dataset...")
    # In a real scenario, you would load your data here
    # For demonstration, we'll create a small synthetic dataset
    n_cells = 100
    n_genes = 50
    expression_matrix = np.random.poisson(lam=5, size=(n_genes, n_cells))
    
    # Create cell metadata
    cell_metadata = pd.DataFrame({
        'cell_id': [f'cell_{i}' for i in range(n_cells)],
        'batch': np.random.choice(['batch1', 'batch2'], size=n_cells),
        'cell_type': np.random.choice(['type1', 'type2', 'type3'], size=n_cells)
    })
    cell_metadata.set_index('cell_id', inplace=True)
    
    # Create gene metadata
    gene_metadata = pd.DataFrame({
        'gene_id': [f'gene_{i}' for i in range(n_genes)],
        'gene_short_name': [f'gene_{i}' for i in range(n_genes)]
    })
    gene_metadata.set_index('gene_id', inplace=True)
    
    # Create the cell dataset
    cds = new_cell_data_set(expression_matrix, cell_metadata, gene_metadata)
    print(f"Created cell dataset with {n_cells} cells and {n_genes} genes")
    
    # 2. Preprocess the cell dataset
    print("\nPreprocessing the cell dataset...")
    cds = preprocess_cds(cds, num_dim=10)
    print("Preprocessing complete")
    
    # 3. Align the cell dataset
    print("\nAligning the cell dataset...")
    cds = align_cds(cds, alignment_group="batch")
    print("Alignment complete")
    
    # 4. Reduce dimensionality
    print("\nReducing dimensionality...")
    cds = reduce_dimension(cds)
    print("Dimensionality reduction complete")
    
    # 5. Cluster cells
    print("\nClustering cells...")
    cds = cluster_cells(cds, resolution=0.1)
    print("Clustering complete")
    
    # 6. Learn graph
    print("\nLearning graph...")
    cds = learn_graph(cds)
    print("Graph learning complete")
    
    # 7. Order cells
    print("\nOrdering cells...")
    # For demonstration, we'll use the first cell as the root
    root_cells = [cds.obs_names[0]]
    cds = order_cells(cds, root_cells=root_cells)
    print("Cell ordering complete")
    
    # 8. Find gene modules
    print("\nFinding gene modules...")
    # For demonstration, we'll use all genes
    gene_module_df = find_gene_modules(cds, resolution=0.1)
    print(f"Found {len(gene_module_df['module'].unique())} gene modules")
    
    # 9. Aggregate gene expression
    print("\nAggregating gene expression...")
    # Create cell group dataframe
    cell_group_df = pd.DataFrame({
        'cell': cds.obs_names,
        'cell_group': cds.obs['cluster']
    })
    agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
    print(f"Aggregated expression matrix shape: {agg_mat.shape}")
    
    # 10. Plot genes in pseudotime
    print("\nPlotting genes in pseudotime...")
    # Select a subset of genes
    selected_genes = cds.var_names[:3]  # First 3 genes
    cds_subset = cds[:, selected_genes]
    
    # Create the plot
    fig = plot_genes_in_pseudotime(cds_subset, min_expr=0.1)
    plt.savefig('genes_in_pseudotime.png')
    print("Plot saved as 'genes_in_pseudotime.png'")
    
    print("\nDemo complete!")