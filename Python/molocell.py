#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from pathlib import Path
import os
import scipy.sparse
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import umap
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

def monocle_theme_opts():
    """Set the theme options for monocle plots"""
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'axes.edgecolor': 'black',
        'axes.linewidth': 0.25,
        'grid.color': 'none',
        'legend.frameon': False,
        'figure.figsize': [10, 10],
        'axes.grid': False,
        'axes.spines.top': False,
        'axes.spines.right': False
    })

class CellDataSet:
    def __init__(self, expression_matrix, cell_metadata, gene_metadata):
        # Ensure all cell indices match
        common_cells = expression_matrix.index.intersection(cell_metadata.index)
        
        self.expression_matrix = expression_matrix.loc[common_cells]
        self.cell_metadata = cell_metadata.loc[common_cells]
        self.gene_metadata = gene_metadata
        
        # Create colData structure
        self.colData = pd.DataFrame(index=common_cells)
        self.colData['cell_name'] = common_cells
        self.colData['assigned_cell_type'] = None  # Will be filled after clustering
        
        # Copy all columns from cell_metadata to colData
        for column in cell_metadata.columns:
            self.colData[column] = self.cell_metadata[column]
        
        # Create rowData structure
        self.rowData = pd.DataFrame(index=gene_metadata.index)
        self.rowData['gene_name'] = gene_metadata.index
        
        # Processing status
        self.preprocessed = False
        self.reduced = False
        self.clustered = False
        self.num_dim = None
        self.pca = None
        self.umap = None
        self.clusters = {}
        
        # Principal graph for trajectory
        self.principal_graph = {}
        self.principal_graph_aux = {}
        
        self.size_factors = None
        self.normalized_counts = None
        self.pseudotime = None
        self.reduced_dimensions = {}
        
    def get_colData(self):
        return self.colData
    
    def set_colData(self, col_name, values):
        self.colData[col_name] = values
        
    def filter_cells(self, cell_type_pattern, ignore_case=True):
        """
        Filter cells based on cell type pattern.
        
        Parameters:
        -----------
        cell_type_pattern : str
            Pattern to match in cell types
        ignore_case : bool
            Whether to ignore case when matching
            
        Returns:
        --------
        filtered_cds : CellDataSet
            Filtered CellDataSet object
        """
        cell_types = self.colData['assigned_cell_type']
        
        if ignore_case:
            mask = cell_types.str.contains(cell_type_pattern, case=False, na=False)
        else:
            mask = cell_types.str.contains(cell_type_pattern, na=False)
            
        # Convert mask to numpy array for indexing
        mask_array = mask.values
        
        # Print debug information
        print(f"\nDebug information:")
        print(f"Expression matrix shape: {self.expression_matrix.shape}")
        print(f"Number of cells in mask: {sum(mask_array)}")
        if hasattr(self, 'umap') and self.umap is not None:
            print(f"UMAP shape: {self.umap.shape}")
        
        # Filter expression matrix - transpose to match cell dimension
        filtered_matrix = self.expression_matrix[mask_array, :]
        
        # Filter metadata
        filtered_cell_metadata = self.cell_metadata[mask]
        
        # Create new CellDataSet with transposed matrix to maintain correct orientation
        filtered_cds = CellDataSet(filtered_matrix.T, filtered_cell_metadata, self.gene_metadata)
        
        # Copy processing status and results
        filtered_cds.preprocessed = self.preprocessed
        filtered_cds.reduced = self.reduced
        filtered_cds.clustered = self.clustered
        filtered_cds.num_dim = self.num_dim
        
        # Filter PCA, UMAP and clusters if they exist
        try:
            if hasattr(self, 'pca') and self.pca is not None:
                filtered_cds.pca = self.pca[mask_array]
            if hasattr(self, 'umap') and self.umap is not None:
                filtered_cds.umap = self.umap[mask_array]
            if hasattr(self, 'clusters') and self.clusters is not None:
                filtered_cds.clusters = self.clusters[mask_array]
        except Exception as e:
            print(f"\nError during dimensional data filtering: {str(e)}")
            # Continue without the dimensional data
            filtered_cds.pca = None
            filtered_cds.umap = None
            filtered_cds.clusters = None
            filtered_cds.reduced = False
            filtered_cds.clustered = False
            
        return filtered_cds
        
    def preprocess_cds(self, num_dim=50, method='PCA', norm_method='log', use_genes=None, 
                      pseudo_count=None, scaling=True, verbose=False):
        """
        Preprocess the CellDataSet object for dimensionality reduction.
        
        Parameters:
        -----------
        num_dim : int
            Number of dimensions to reduce to
        method : str
            Method for dimensionality reduction ('PCA' or 'LSI')
        norm_method : str
            Method for normalizing expression data ('log', 'size_only', or 'none')
        use_genes : list or None
            List of genes to use for dimensionality reduction
        pseudo_count : float or None
            Pseudo count to add before log transformation
        scaling : bool
            Whether to scale the data
        verbose : bool
            Whether to print verbose output
            
        Returns:
        --------
        self : CellDataSet
            The processed CellDataSet object
        """
        if method not in ['PCA', 'LSI']:
            raise ValueError("method must be one of 'PCA' or 'LSI'")
            
        if norm_method not in ['log', 'size_only', 'none']:
            raise ValueError("norm_method must be one of 'log', 'size_only', or 'none'")
            
        # Normalize expression data
        FM = self.normalize_expr_data(norm_method, pseudo_count)
        
        if use_genes is not None:
            if not all(gene in self.gene_metadata.index for gene in use_genes):
                raise ValueError("All genes in use_genes must be present in gene_metadata")
            FM = FM[use_genes, :]
            
        # Remove genes with zero variance
        fm_rowsums = np.sum(FM, axis=1)
        FM = FM[np.isfinite(fm_rowsums) & (fm_rowsums != 0), :]
        
        if method == 'PCA':
            # Initialize reduction metadata
            self.initialize_reduce_dim_metadata('PCA')
            self.initialize_reduce_dim_model_identity('PCA')
            
            if verbose:
                print("Remove noise by PCA ...")
                
            # Perform PCA
            from sklearn.decomposition import PCA
            pca = PCA(n_components=min(num_dim, min(FM.shape) - 1))
            preproc_res = pca.fit_transform(FM.T)
            
            # Store results
            self.pca = preproc_res
            self.num_dim = num_dim
            
            # Store PCA model information
            self.reduce_dim_aux = {
                'PCA': {
                    'model': {
                        'num_dim': num_dim,
                        'norm_method': norm_method,
                        'use_genes': use_genes,
                        'pseudo_count': pseudo_count,
                        'svd_v': pca.components_.T,
                        'svd_sdev': np.sqrt(pca.explained_variance_),
                        'svd_center': pca.mean_,
                        'svd_scale': None,
                        'prop_var_expl': pca.explained_variance_ratio_
                    }
                }
            }
            
        elif method == 'LSI':
            # Initialize reduction metadata
            self.initialize_reduce_dim_metadata('LSI')
            self.initialize_reduce_dim_model_identity('LSI')
            
            # Perform TF-IDF transformation
            tfidf_res = self.tfidf(FM)
            preproc_res = tfidf_res['tf_idf_counts']
            
            # Perform SVD
            from sklearn.decomposition import TruncatedSVD
            svd = TruncatedSVD(n_components=min(num_dim, min(FM.shape) - 1))
            preproc_res = svd.fit_transform(preproc_res.T)
            
            # Store results
            self.lsi = preproc_res
            self.num_dim = num_dim
            
            # Store LSI model information
            self.reduce_dim_aux = {
                'LSI': {
                    'model': {
                        'num_dim': num_dim,
                        'norm_method': norm_method,
                        'use_genes': use_genes,
                        'pseudo_count': pseudo_count,
                        'log_scale_tf': tfidf_res['log_scale_tf'],
                        'frequencies': tfidf_res['frequencies'],
                        'scale_factor': tfidf_res['scale_factor'],
                        'col_sums': tfidf_res['col_sums'],
                        'row_sums': tfidf_res['row_sums'],
                        'num_cols': tfidf_res['num_cols'],
                        'svd_v': svd.components_.T,
                        'svd_sdev': svd.singular_values_ / np.sqrt(max(1, FM.shape[1] - 1))
                    }
                }
            }
            
        self.preprocessed = True
        return self
        
    def normalize_expr_data(self, norm_method='log', pseudo_count=None):
        """
        Normalize expression data before dimensionality reduction.
        
        Parameters:
        -----------
        norm_method : str
            Method for normalizing expression data ('log', 'size_only', or 'none')
        pseudo_count : float or None
            Pseudo count to add before log transformation
            
        Returns:
        --------
        FM : ndarray
            Normalized expression matrix
        """
        FM = self.expression_matrix.values
        
        if pseudo_count is None:
            if norm_method == 'log':
                pseudo_count = 1
            else:
                pseudo_count = 0
                
        if norm_method == 'log':
            # Calculate size factors if not already present
            if self.size_factors is None:
                # Simple size factor calculation - sum of counts for each cell
                self.size_factors = np.sum(FM, axis=1)
                # Replace zero size factors with the mean of non-zero size factors
                non_zero_mean = np.mean(self.size_factors[self.size_factors > 0])
                self.size_factors[self.size_factors == 0] = non_zero_mean
                # Normalize size factors to have mean 1
                self.size_factors = self.size_factors / np.mean(self.size_factors)
                
            # Normalize by size factor
            FM = FM / self.size_factors.reshape(-1, 1)
            FM = FM + pseudo_count
            FM = np.log2(FM)
            
        elif norm_method == 'size_only':
            # Calculate size factors if not already present
            if self.size_factors is None:
                # Simple size factor calculation - sum of counts for each cell
                self.size_factors = np.sum(FM, axis=1)
                # Replace zero size factors with the mean of non-zero size factors
                non_zero_mean = np.mean(self.size_factors[self.size_factors > 0])
                self.size_factors[self.size_factors == 0] = non_zero_mean
                # Normalize size factors to have mean 1
                self.size_factors = self.size_factors / np.mean(self.size_factors)
                
            # Normalize by size factor
            FM = FM / self.size_factors.reshape(-1, 1)
            FM = FM + pseudo_count
            
        # Replace any remaining NaN or infinite values with 0
        FM = np.nan_to_num(FM, nan=0.0, posinf=0.0, neginf=0.0)
            
        # Store normalized counts
        self.normalized_counts = FM
            
        return FM
        
    def tfidf(self, count_matrix, frequencies=True, log_scale_tf=True, scale_factor=100000):
        """
        Perform TF-IDF transformation on count matrix.
        
        Parameters:
        -----------
        count_matrix : ndarray
            Count matrix to transform
        frequencies : bool
            Whether to use term frequency
        log_scale_tf : bool
            Whether to log scale term frequency
        scale_factor : float
            Scale factor for log scaling
            
        Returns:
        --------
        dict
            Dictionary containing transformed matrix and parameters
        """
        # Calculate term frequency
        if frequencies:
            col_sums = np.sum(count_matrix, axis=0)
            tf = count_matrix / col_sums[np.newaxis, :]
        else:
            col_sums = None
            tf = count_matrix
            
        # Log scale if requested
        if log_scale_tf:
            if frequencies:
                tf = np.log1p(tf * scale_factor)
            else:
                tf = np.log1p(tf)
                
        # Calculate inverse document frequency
        num_cols = count_matrix.shape[1]
        row_sums = np.sum(count_matrix > 0, axis=1)
        idf = np.log(1 + num_cols / row_sums)
        
        # Calculate TF-IDF
        tf_idf_counts = tf * idf[:, np.newaxis]
        
        return {
            'tf_idf_counts': tf_idf_counts,
            'frequencies': frequencies,
            'log_scale_tf': log_scale_tf,
            'scale_factor': scale_factor,
            'col_sums': col_sums,
            'row_sums': row_sums,
            'num_cols': num_cols
        }
        
    def initialize_reduce_dim_metadata(self, method):
        """
        Initialize reduction dimension metadata.
        
        Parameters:
        -----------
        method : str
            Method name
        """
        if not hasattr(self, 'reduce_dim_aux'):
            self.reduce_dim_aux = {}
        if method not in self.reduce_dim_aux:
            self.reduce_dim_aux[method] = {'model': {}}
            
    def initialize_reduce_dim_model_identity(self, method):
        """
        Initialize reduction dimension model identity.
        
        Parameters:
        -----------
        method : str
            Method name
        """
        if not hasattr(self, 'reduce_dim_model_identity'):
            self.reduce_dim_model_identity = {}
        if method not in self.reduce_dim_model_identity:
            self.reduce_dim_model_identity[method] = {}
        
    def reduce_dimension(self, method="UMAP", n_components=2, **kwargs):
        """
        Reduce dimensionality of the data.
        
        Parameters:
        -----------
        method : str
            Method to use for dimensionality reduction ('UMAP' or 'tSNE')
        n_components : int
            Number of dimensions to reduce to (should be 2 for UMAP1 and UMAP2)
        **kwargs : dict
            Additional arguments to pass to the reduction method
        """
        if self.normalized_counts is None:
            raise ValueError("Data must be normalized before dimension reduction")

        # Always use 2 components for UMAP
        if method == "UMAP":
            n_components = 2

        data = self.normalized_counts
        
        if method == "UMAP":
            # Calculate UMAP with specific parameters for better visualization
            reducer = umap.UMAP(
                n_components=2,  # Force 2D
                n_neighbors=30,
                min_dist=0.3,
                metric='euclidean',
                random_state=42
            )
            reduced_data = reducer.fit_transform(data)
            
            # Ensure the output is 2D
            if reduced_data.shape[1] != 2:
                raise ValueError("UMAP reduction failed to produce 2D output")
            
        elif method == "tSNE":
            reducer = TSNE(n_components=n_components, **kwargs)
            reduced_data = reducer.fit_transform(data)
        else:
            raise ValueError("Method must be either 'UMAP' or 'tSNE'")
        
        # Store the reduced dimensions
        self.reduced_dimensions[method] = reduced_data
        
        return self
        
    def align_cds(self, batch_key="patch_id", n_components=10, reduction_method="PCA", 
                 preprocess_method="PCA", preprocess_norm_method="log", 
                 preprocess_num_dim=50, preprocess_use_genes=None, 
                 preprocess_pseudo_count=None, preprocess_scaling=True, 
                 preprocess_verbose=False):
        """
        Align cells across batches using dimensionality reduction.
        
        Parameters:
        -----------
        batch_key : str
            Column name in cell_metadata containing batch information
        n_components : int
            Number of components to use for alignment
        reduction_method : str
            Method for dimensionality reduction ('PCA' or 'LSI')
        preprocess_method : str
            Method for preprocessing ('PCA' or 'LSI')
        preprocess_norm_method : str
            Method for normalizing expression data ('log', 'size_only', or 'none')
        preprocess_num_dim : int
            Number of dimensions to reduce to during preprocessing
        preprocess_use_genes : list or None
            List of genes to use for dimensionality reduction
        preprocess_pseudo_count : float or None
            Pseudo count to add before log transformation
        preprocess_scaling : bool
            Whether to scale the data
        preprocess_verbose : bool
            Whether to print verbose output
            
        Returns:
        --------
        self : CellDataSet
            The aligned CellDataSet object
        """
        # Check if batch_key exists in cell_metadata
        if batch_key not in self.cell_metadata.columns:
            raise ValueError(f"batch_key '{batch_key}' not found in cell_metadata")
            
        # Get unique batches
        batches = self.cell_metadata[batch_key].unique()
        if len(batches) < 2:
            print("Only one batch found, no alignment needed")
            return self
            
        # Preprocess data if not already done
        if not self.preprocessed:
            self.preprocess_cds(
                num_dim=preprocess_num_dim,
                method=preprocess_method,
                norm_method=preprocess_norm_method,
                use_genes=preprocess_use_genes,
                pseudo_count=preprocess_pseudo_count,
                scaling=preprocess_scaling,
                verbose=preprocess_verbose
            )
            
        # Get reduction coordinates
        if reduction_method == "PCA":
            if not hasattr(self, 'pca') or self.pca is None:
                raise ValueError("PCA reduction not found. Please run preprocess_cds with method='PCA' first.")
            coords = self.pca
        elif reduction_method == "LSI":
            if not hasattr(self, 'lsi') or self.lsi is None:
                raise ValueError("LSI reduction not found. Please run preprocess_cds with method='LSI' first.")
            coords = self.lsi
        else:
            raise ValueError("reduction_method must be one of 'PCA' or 'LSI'")
            
        # Initialize reduction metadata
        self.initialize_reduce_dim_metadata('Aligned')
        self.initialize_reduce_dim_model_identity('Aligned')
        
        # Perform alignment using sklearn's PCA
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(n_components, coords.shape[1]))
        aligned_coords = pca.fit_transform(coords)
        
        # Store aligned coordinates
        self.aligned = aligned_coords
        self.num_dim = n_components
        
        # Store alignment model information
        self.reduce_dim_aux['Aligned'] = {
            'model': {
                'num_dim': n_components,
                'reduction_method': reduction_method,
                'batch_key': batch_key,
                'svd_v': pca.components_.T,
                'svd_sdev': np.sqrt(pca.explained_variance_),
                'svd_center': pca.mean_,
                'svd_scale': None,
                'prop_var_expl': pca.explained_variance_ratio_
            }
        }
        
        # Set aligned flag
        self.aligned = True
        return self
        
    def cluster_cells(self, reduction_method="UMAP", k=20, cluster_method="leiden",
                     num_iter=2, resolution=None, random_seed=2016):
        """
        Cluster cells using the specified method.
        """
        # Initialize clusters dictionary if it doesn't exist
        if not hasattr(self, 'clusters'):
            self.clusters = {}
            
        # Get the reduced dimension data
        if reduction_method not in self.reduced_dimensions:
            raise ValueError(f"No {reduction_method} coordinates found. Please run reduce_dimension first.")
            
        data = self.reduced_dimensions[reduction_method]
            
        # Perform clustering
        if cluster_method == "leiden":
            # Build kNN graph
            from sklearn.neighbors import NearestNeighbors
            nbrs = NearestNeighbors(n_neighbors=k).fit(data)
            distances, indices = nbrs.kneighbors(data)
            
            # Convert to igraph
            import igraph as ig
            edges = []
            for i in range(len(indices)):
                for j in indices[i]:
                    if i != j:
                        edges.append((i, j))
                        
            g = ig.Graph(edges=edges, directed=False)
            g.simplify()  # Remove duplicate edges
            
            # Run Leiden clustering
            import leidenalg
            partition = leidenalg.find_partition(g, 
                                               leidenalg.RBConfigurationVertexPartition,
                                               resolution_parameter=resolution if resolution else 1.0,
                                               seed=random_seed)
                                               
            clusters = np.array(partition.membership)
            
        else:
            raise ValueError(f"Unsupported clustering method: {cluster_method}")
            
        # Store results
        self.clusters[reduction_method] = {
            'cluster': clusters,
            'method': cluster_method,
            'params': {
                'k': k,
                'num_iter': num_iter,
                'resolution': resolution,
                'random_seed': random_seed
            }
        }
        
        return self
        
    def learn_graph(self, close_loop=True, learn_graph_control=None, use_partition=True):
        """
        Learn the principal graph from the data.
        """
        # Initialize learn_graph_control if not provided
        if learn_graph_control is None:
            learn_graph_control = {
                'euclidean_distance_ratio': 1,
                'geodesic_distance_ratio': 1,
                'minimal_branch_len': 12,
                'orthogonal_proj_tip': True,
                'prune_graph': True,
                'scale': 0.9
            }
            
        # Get UMAP coordinates
        if "UMAP" not in self.reduced_dimensions:
            raise ValueError("No UMAP coordinates found. Please run reduce_dimension first.")
            
        data = self.reduced_dimensions["UMAP"]
        
        # Initialize principal graph dictionaries if they don't exist
        if not hasattr(self, 'principal_graph'):
            self.principal_graph = {}
            
        if not hasattr(self, 'principal_graph_aux'):
            self.principal_graph_aux = {}
            
        # Build kNN graph
        from sklearn.neighbors import NearestNeighbors
        k = 20  # You can make this a parameter if needed
        nbrs = NearestNeighbors(n_neighbors=k).fit(data)
        distances, indices = nbrs.kneighbors(data)
        
        # Convert to igraph
        import igraph as ig
        edges = []
        for i in range(len(indices)):
            for j in indices[i]:
                if i != j:
                    edges.append((i, j))
                    
        g = ig.Graph(edges=edges, directed=False)
        g.simplify()  # Remove duplicate edges
        
        # Find minimum spanning tree
        mst = g.spanning_tree(weights=None)
        
        # Store the graph
        self.principal_graph['UMAP'] = mst
        
        # Store auxiliary information
        self.principal_graph_aux['UMAP'] = {
            'control': learn_graph_control,
            'close_loop': close_loop,
            'use_partition': use_partition
        }
        
        return self
        
    def order_cells(self, reduction_method="UMAP", root_pr_nodes=None, root_cells=None, verbose=False):
        """
        Order cells along a trajectory.
        """
        from sklearn.metrics.pairwise import euclidean_distances
        
        # Validate reduction_method
        if reduction_method != "UMAP":
            raise ValueError("Currently only 'UMAP' is accepted as a reduction_method.")
            
        # Check if reduction exists
        if reduction_method not in self.reduced_dimensions:
            raise ValueError(f"No {reduction_method} coordinates found. Please run reduce_dimension first.")
                           
        # Check if clusters exist
        if not hasattr(self, 'clusters') or reduction_method not in self.clusters:
            raise ValueError(f"No cell clusters for {reduction_method} calculated. "
                           f"Please run cluster_cells with reduction_method = {reduction_method} "
                           f"and run learn_graph before running order_cells.")
                           
        # Check if principal graph exists
        if not hasattr(self, 'principal_graph') or reduction_method not in self.principal_graph:
            raise ValueError(f"No principal graph for {reduction_method} calculated. "
                           f"Please run learn_graph with reduction_method = {reduction_method} "
                           f"before running order_cells.")
                           
        # Check if principal graph is too large
        graph = self.principal_graph[reduction_method]
        if graph.vcount() >= 10000:  # Use vcount() instead of len()
            raise ValueError("Principal graph is too large. order_cells doesn't support "
                           "more than 10 thousand centroids.")
                           
        # Get the reduced dimension data
        data = self.reduced_dimensions[reduction_method]
                           
        # Validate root_pr_nodes if provided
        if root_pr_nodes is not None:
            # Check if nodes exist in the graph
            valid_nodes = [str(i) for i in range(graph.vcount())]
            if not all(node in valid_nodes for node in root_pr_nodes):
                raise ValueError("All provided root_pr_nodes must be present in the principal graph.")
                
        # Validate root_cells if provided
        if root_cells is not None:
            if not all(cell in self.colData.index for cell in root_cells):
                raise ValueError("All provided root_cells must be present in the cell data set.")
                
        # Check if either root_pr_nodes or root_cells is provided
        if root_cells is None and root_pr_nodes is None:
            # Use the first node as root
            root_pr_nodes = ["0"]
            
        elif root_cells is not None:
            # Find closest principal graph nodes for root cells
            if not hasattr(self, 'principal_graph_aux') or reduction_method not in self.principal_graph_aux:
                raise ValueError("Principal graph auxiliary data not found. "
                               "Please run learn_graph before running order_cells.")
                               
            # Calculate distances to all nodes
            cell_coords = data[root_cells]
            node_coords = np.array([data[i] for i in range(graph.vcount())])
            distances = euclidean_distances(cell_coords, node_coords)
            closest_nodes = np.argmin(distances, axis=1)
            root_pr_nodes = [str(node) for node in closest_nodes]
            root_pr_nodes = list(set(root_pr_nodes))  # Remove duplicates
            
        # Store root nodes
        if not hasattr(self, 'principal_graph_aux'):
            self.principal_graph_aux = {}
            
        if reduction_method not in self.principal_graph_aux:
            self.principal_graph_aux[reduction_method] = {}
            
        self.principal_graph_aux[reduction_method]['root_pr_nodes'] = root_pr_nodes
        
        # Calculate pseudotime using graph distances
        import numpy as np
        
        # Convert igraph to distance matrix
        distances = np.array(graph.distances())  # Use distances() instead of shortest_paths()
        
        # Get distances from root nodes
        root_indices = [int(node) for node in root_pr_nodes]
        root_distances = distances[root_indices]
        
        # For each cell, find its closest node in the graph
        cell_coords = data
        node_coords = np.array([data[i] for i in range(graph.vcount())])
        cell_to_node_distances = euclidean_distances(cell_coords, node_coords)
        closest_nodes = np.argmin(cell_to_node_distances, axis=1)
        
        # Assign pseudotime based on shortest path distance through the graph
        pseudotime = np.min(root_distances[:, closest_nodes], axis=0)
        
        # Normalize pseudotime to [0, 1]
        pseudotime = (pseudotime - np.min(pseudotime)) / (np.max(pseudotime) - np.min(pseudotime))
        
        # Store pseudotime
        self.principal_graph_aux[reduction_method]['pseudotime'] = pseudotime
        
        return self
        
    def plot_cells(self, reduction_method="UMAP", color_cells_by=None, group_cells_by=None,
                    show_trajectory_graph=False, show_backbone=False, cell_size=10, alpha=0.7):
        """
        Plot cells using dimensionality reduction coordinates.
        
        Parameters:
        -----------
        reduction_method : str
            The dimensionality reduction method to use for plotting
        color_cells_by : str
            Column name in colData to color cells by
        group_cells_by : str
            Column name in colData to group cells by
        show_trajectory_graph : bool
            Whether to show the trajectory graph
        show_backbone : bool
            Whether to show the backbone
        cell_size : int
            Size of cells in the plot
        alpha : float
            Transparency of cells
        """
        if reduction_method not in self.reduced_dimensions:
            raise ValueError(f"No {reduction_method} coordinates found. "
                            f"Please run reduce_dimension first.")
        
        coords = self.reduced_dimensions[reduction_method]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Color cells if specified
        if color_cells_by is not None:
            if color_cells_by not in self.colData.columns:
                raise ValueError(f"Column {color_cells_by} not found in colData")
            
            # Get colors for each cell
            colors = self.colData[color_cells_by].values
            unique_colors = np.unique(colors)
            color_map = plt.cm.tab20(np.linspace(0, 1, len(unique_colors)))
            color_dict = dict(zip(unique_colors, color_map))
            
            # Create color array with same length as coords
            cell_colors = np.array([color_dict[colors[i]] if i < len(colors) else color_dict[unique_colors[0]] 
                                  for i in range(len(coords))])
        else:
            cell_colors = 'blue'
        
        # Plot scatter
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=cell_colors, 
                            s=cell_size, alpha=alpha)
        
        # Show trajectory graph if requested
        if show_trajectory_graph and self.principal_graph is not None:
            graph = self.principal_graph[reduction_method]
            edges = graph.get_edgelist()
            for edge in edges:
                ax.plot([coords[edge[0], 0], coords[edge[1], 0]],
                       [coords[edge[0], 1], coords[edge[1], 1]],
                       'k-', alpha=0.5)
        
        # Show backbone if requested
        if show_backbone and self.principal_graph is not None:
            backbone = self.principal_graph['backbone']
            ax.plot(backbone[:, 0], backbone[:, 1], 'r-', linewidth=2)
        
        # Add legend
        if color_cells_by is not None:
            legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=color_dict[label],
                                        label=label, markersize=10)
                             for label in unique_colors]
            ax.legend(handles=legend_elements, title=color_cells_by)
        
        ax.set_title(f"{reduction_method} Plot")
        plt.show()
        return fig

def read_and_sample_data(file_path: str, sample_fraction: float = 0.1):
    """
    Read and sample data from an h5ad file.
    
    Parameters:
    -----------
    file_path : str
        Path to the h5ad file
    sample_fraction : float
        Fraction of cells to sample
        
    Returns:
    --------
    cds : CellDataSet
        Sampled CellDataSet object
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    adata = sc.read_h5ad(file_path)
    n_cells = int(adata.n_obs * sample_fraction)
    sc.pp.subsample(adata, n_obs=n_cells, random_state=42)
    
    cell_metadata = adata.obs.copy()
    gene_metadata = adata.var.copy()
    
    # Convert sparse matrix to dense if necessary
    if scipy.sparse.issparse(adata.X):
        expression_matrix = pd.DataFrame(adata.X.toarray(), 
                                       index=adata.obs_names,
                                       columns=adata.var_names)
    else:
        expression_matrix = pd.DataFrame(adata.X,
                                       index=adata.obs_names,
                                       columns=adata.var_names)
    
    # Create CellDataSet object
    cds = CellDataSet(expression_matrix, cell_metadata, gene_metadata)
    
    return cds

def load_data(file_path):
    """
    Load data from h5ad file and return AnnData object
    """
    try:
        # Read AnnData file
        adata = sc.read_h5ad(file_path)
        
        # Print basic information about the data
        print(f"Data shape: {adata.shape}")
        print(f"Number of cells: {adata.n_obs}")
        print(f"Number of genes: {adata.n_vars}")
        print("\nFirst 10 cell names:")
        print(adata.obs_names[:10].tolist())
        print("\nMetadata columns:")
        print(adata.obs.columns.tolist())
        print("\nFirst few rows of metadata:")
        print(adata.obs.head())
        
        return adata
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        return None

def reduce_dimension(adata, method='umap', n_components=2, perplexity=30):
    """
    Reduce data dimensions using UMAP or t-SNE
    """
    try:
        # Convert data to dense array if it's sparse
        if scipy.sparse.issparse(adata.X):
            data = adata.X.toarray()
        else:
            data = adata.X
            
        if method == 'umap':
            # Calculate UMAP
            reducer = umap.UMAP(n_components=n_components, random_state=42)
            embedding = reducer.fit_transform(data)
            
            # Store results in AnnData
            adata.obsm['X_umap'] = embedding
            
        elif method == 'tsne':
            # Calculate t-SNE
            tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
            embedding = tsne.fit_transform(data)
            
            # Store results in AnnData
            adata.obsm['X_tsne'] = embedding
            
        return adata
        
    except Exception as e:
        print(f"Error in dimension reduction: {str(e)}")
        return None

def cluster_cells(adata, method='kmeans', n_clusters=20):
    """
    Cluster cells using K-means or Louvain
    """
    try:
        # Convert data to dense array if it's sparse
        if scipy.sparse.issparse(adata.X):
            data = adata.X.toarray()
        else:
            data = adata.X
            
        if method == 'kmeans':
            # Use K-means
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            labels = kmeans.fit_predict(data)
            
            # Store labels in AnnData
            adata.obs['cluster'] = labels
            
        elif method == 'louvain':
            # Use Louvain algorithm
            sc.pp.neighbors(adata)
            sc.tl.louvain(adata)
            
        return adata
        
    except Exception as e:
        print(f"Error in clustering: {str(e)}")
        return None

def learn_graph(adata, n_neighbors=15, metric='euclidean'):
    """
    Learn graph from reduced data
    """
    try:
        # Convert data to dense array if it's sparse
        if scipy.sparse.issparse(adata.X):
            data = adata.X.toarray()
        else:
            data = adata.X
            
        # Use UMAP to build graph
        reducer = umap.UMAP(n_neighbors=n_neighbors, metric=metric)
        embedding = reducer.fit_transform(data)
        
        # Build kNN graph
        nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
        nbrs.fit(embedding)
        distances, indices = nbrs.kneighbors(embedding)
        
        # Build NetworkX graph
        G = nx.Graph()
        
        # Add edges
        for i in range(len(indices)):
            for j, dist in zip(indices[i], distances[i]):
                if i != j:
                    G.add_edge(i, j, weight=1/dist)
        
        # Store graph in adata
        if not hasattr(adata, 'uns'):
            adata.uns = {}
        adata.uns['graph'] = G
        
        return G
        
    except Exception as e:
        print(f"Error in graph learning: {str(e)}")
        return None

def plot_results(adata, method='umap', color_by='cluster', title=None, show_graph=True):
    """
    Plot dimension reduction and clustering results
    """
    try:
        plt.figure(figsize=(10, 8))
        
        if method == 'umap':
            coords = adata.obsm['X_umap']
        else:
            coords = adata.obsm['X_tsne']
            
        if color_by == 'cluster':
            # Get unique cluster labels
            unique_clusters = sorted(adata.obs['cluster'].unique())
            n_clusters = len(unique_clusters)
            
            # Create a color map for clusters
            colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
            
            # Create a scatter plot for each cluster
            for i, cluster in enumerate(unique_clusters):
                mask = adata.obs['cluster'] == cluster
                plt.scatter(coords[mask, 0], coords[mask, 1], 
                           label=f'Cluster {cluster}',
                           color=colors[i],
                           alpha=0.6)
            
            # Add legend with cluster names
            plt.legend(title='Clusters', loc='best')
            
        # Add graph if available and requested
        if show_graph and hasattr(adata, 'uns') and 'graph' in adata.uns:
            G = adata.uns['graph']
            # Get node positions from UMAP coordinates
            pos = {i: (coords[i, 0], coords[i, 1]) for i in range(len(coords))}
            # Draw edges
            for edge in G.edges():
                plt.plot([pos[edge[0]][0], pos[edge[1]][0]], 
                        [pos[edge[0]][1], pos[edge[1]][1]], 
                        'k-', alpha=0.2, linewidth=0.5)
            
        plt.title(title or f'Cell Clusters ({method.upper()})')
        plt.xlabel(f'{method.upper()} 1')
        plt.ylabel(f'{method.upper()} 2')
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error in plotting results: {str(e)}")

def main():
    """
    Main function following the Monocle3-like workflow
    """
    try:
        print("\n=== Step 1: Loading data ===")
        # Load all data first
        adata = load_data('Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5')
        if adata is None:
            raise ValueError("Failed to load data")
        print("✓ Data loaded successfully")
            
        print("\n=== Step 2: Creating new sample ===")
        sample_fraction = 0.1
        print(f"Sampling data with fraction: {sample_fraction}...")
        cds = read_and_sample_data('Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5', 
                                  sample_fraction=sample_fraction)
        print(f"✓ Created sample with {cds.expression_matrix.shape[0]} cells and {cds.expression_matrix.shape[1]} genes")
        
        print("\n=== Step 3: Preprocessing data ===")
        print("Starting preprocessing with PCA...")
        cds = cds.preprocess_cds(num_dim=50)
        if not cds.preprocessed:
            raise ValueError("Preprocessing failed")
        print("✓ Preprocessing completed")
        
        print("\n=== Step 4: Aligning data ===")
        if 'patch_id' in cds.cell_metadata.columns:
            print("Found patch_id in metadata, starting alignment...")
            unique_patches = cds.cell_metadata['patch_id'].unique()
            print(f"Found {len(unique_patches)} unique patches: {unique_patches}")
            cds = cds.align_cds(batch_key='patch_id')
            print("✓ Data alignment completed")
        else:
            print("⚠ No patch_id found in metadata, skipping alignment")
        
        print("\n=== Step 5: Reducing dimension ===")
        print("Starting UMAP dimension reduction...")
        cds = cds.reduce_dimension(method="UMAP")  # Will use 2D by default
        if "UMAP" not in cds.reduced_dimensions:
            raise ValueError("Dimension reduction failed")
        print("✓ Dimension reduction completed")
        print(f"UMAP coordinates shape: {cds.reduced_dimensions['UMAP'].shape}")
        
        print("\n=== Step 6: Clustering cells ===")
        print("Starting cell clustering...")
        cds = cds.cluster_cells(reduction_method="UMAP", k=20)
        if not hasattr(cds, 'clusters') or "UMAP" not in cds.clusters:
            raise ValueError("Clustering failed")
        n_clusters = len(np.unique(cds.clusters["UMAP"]['cluster']))
        print(f"✓ Clustering completed with {n_clusters} clusters")
        
        print("\n=== Step 7: Plotting results ===")
        print("Generating final plot...")
        cds.plot_cells(reduction_method="UMAP",
                      color_cells_by="cell_type",
                      show_trajectory_graph=False)  # Disable trajectory graph
        print("✓ Plot generated")
        
        print("\n=== Workflow completed successfully! ===")
        
        # Return the CellDataSet object for further analysis if needed
        return cds
        
    except Exception as e:
        print(f"\n❌ Error in workflow: {str(e)}")
        import traceback
        print("\nTraceback:")
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()