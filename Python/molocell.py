#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from pathlib import Path
import os
import scipy.sparse
import matplotlib.pyplot as plt
import seaborn as sns

def monocle_theme_opts():
    """Set the theme options for monocle plots"""
    plt.style.use('seaborn-white')
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.linewidth'] = 0.25
    plt.rcParams['grid.color'] = 'none'
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.key'] = 'white'

class CellDataSet:
    def __init__(self, expression_matrix, cell_metadata, gene_metadata):
        self.expression_matrix = expression_matrix
        self.cell_metadata = cell_metadata
        self.gene_metadata = gene_metadata
        
        # Create colData structure
        self.colData = pd.DataFrame(index=cell_metadata.index)
        self.colData['cell_name'] = cell_metadata.index
        self.colData['assigned_cell_type'] = None  # Will be filled after clustering
        
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
        self.clusters = None
        
    def get_colData(self):
        return self.colData
    
    def set_colData(self, col_name, values):
        self.colData[col_name] = values
        
    def filter_cells(self, cell_type_pattern, ignore_case=True):
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
        
    def plot_cells(self, x=1, y=2, reduction_method="UMAP", color_cells_by="cluster",
                  group_cells_by="cluster", genes=None, show_trajectory_graph=True,
                  trajectory_graph_color="grey28", trajectory_graph_segment_size=0.75,
                  norm_method="log", label_cell_groups=True, label_groups_by_cluster=True,
                  group_label_size=2, labels_per_group=1, label_branch_points=True,
                  label_roots=True, label_leaves=True, graph_label_size=2,
                  cell_size=0.35, cell_stroke=None, alpha=1, min_expr=0.1,
                  scale_to_range=True, label_principal_points=False):
        """
        Plot cells in a reduced dimension space.
        
        Parameters:
        -----------
        x : int
            The column of reduced dimension to plot on the horizontal axis
        y : int
            The column of reduced dimension to plot on the vertical axis
        reduction_method : str
            The reduction method to use (UMAP, tSNE, PCA, LSI, Aligned)
        color_cells_by : str
            What to use for coloring the cells (cluster, partition, pseudotime, or column name)
        group_cells_by : str
            How to group cells when labeling them
        genes : list or None
            List of genes to color cells by expression
        show_trajectory_graph : bool
            Whether to show the trajectory graph
        trajectory_graph_color : str
            Color of the trajectory graph
        trajectory_graph_segment_size : float
            Size of the trajectory graph segments
        norm_method : str
            How to normalize gene expression (log or size_only)
        label_cell_groups : bool
            Whether to label cell groups
        label_groups_by_cluster : bool
            Whether to label groups by cluster
        group_label_size : float
            Size of group labels
        labels_per_group : int
            Number of labels per group
        label_branch_points : bool
            Whether to label branch points
        label_roots : bool
            Whether to label roots
        label_leaves : bool
            Whether to label leaves
        graph_label_size : float
            Size of graph labels
        cell_size : float
            Size of cells in the plot
        cell_stroke : float or None
            Stroke width of cells
        alpha : float
            Transparency of cells
        min_expr : float
            Minimum expression threshold
        scale_to_range : bool
            Whether to scale expression to range
        label_principal_points : bool
            Whether to label principal points
            
        Returns:
        --------
        matplotlib.figure.Figure
            The plot figure
        """
        # Set cell stroke if not provided
        if cell_stroke is None:
            cell_stroke = cell_size / 2
            
        # Check if reduction exists
        if reduction_method == "UMAP" and self.umap is None:
            raise ValueError("UMAP reduction not found. Please run reduce_dimension first.")
        elif reduction_method == "PCA" and self.pca is None:
            raise ValueError("PCA reduction not found. Please run reduce_dimension first.")
            
        # Get reduction coordinates
        if reduction_method == "UMAP":
            coords = self.umap
        elif reduction_method == "PCA":
            coords = self.pca
            
        # Create plot
        plt.figure(figsize=(10, 10))
        
        # Plot cells
        if color_cells_by == "cluster" and self.clusters is not None:
            scatter = plt.scatter(coords[:, x-1], coords[:, y-1], 
                                c=self.clusters, cmap='viridis',
                                s=cell_size*100, alpha=alpha)
            plt.colorbar(scatter, label='Cluster')
        else:
            plt.scatter(coords[:, x-1], coords[:, y-1], 
                       s=cell_size*100, alpha=alpha)
            
        # Add labels
        if label_cell_groups:
            if group_cells_by == "cluster" and self.clusters is not None:
                for cluster in np.unique(self.clusters):
                    mask = self.clusters == cluster
                    center_x = np.mean(coords[mask, x-1])
                    center_y = np.mean(coords[mask, y-1])
                    plt.text(center_x, center_y, f'Cluster {cluster}',
                            fontsize=group_label_size*5)
                    
        plt.title(f'{reduction_method} Plot')
        plt.xlabel(f'{reduction_method} {x}')
        plt.ylabel(f'{reduction_method} {y}')
        
        # Apply monocle theme
        monocle_theme_opts()
        
        return plt.gcf()

def read_and_sample_data(file_path: str, sample_fraction: float = 0.01):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    adata = sc.read_h5ad(file_path)
    n_cells = int(adata.n_obs * sample_fraction)
    sc.pp.subsample(adata, n_obs=n_cells, random_state=42)
    
    cell_metadata = adata.obs.copy()
    gene_metadata = adata.var.copy()
    
    # Convert sparse matrix to dense if necessary
    if scipy.sparse.issparse(adata.X):
        expression_matrix = adata.X.toarray()
    else:
        expression_matrix = adata.X
    
    # Create CellDataSet object
    cds = CellDataSet(expression_matrix, cell_metadata, gene_metadata)
    
    return cds

if __name__ == "__main__":
    base_path = Path(__file__).parent.parent
    file_path = base_path / "Data" / "skin_TXK6Z4X_A1_processed" / "tmap" / "weighted_by_area_celltypist_cells_adata.h5"
    file_path = str(file_path)
    
    # Read and sample data
    cds = read_and_sample_data(file_path)
    
    # Print cell names
    print("First 10 cell names:")
    print(cds.colData['cell_name'].head(10).tolist())
    
    # Get a sample cell name
    sample_cell = cds.colData['cell_name'].iloc[0]
    print(f"\nUsing sample cell: {sample_cell}")
    
    # Create a test cell type assignment using the sample cell name
    cds.colData['assigned_cell_type'] = 'Other'
    cds.colData.loc[sample_cell, 'assigned_cell_type'] = 'TestCell'
    
    # Create dummy clusters and UMAP coordinates with correct dimensions
    n_cells = cds.expression_matrix.shape[0]  # Use first dimension for cells
    cds.clusters = np.random.randint(0, 5, size=n_cells)
    cds.umap = np.random.randn(n_cells, 2)
    cds.reduced = True
    
    print(f"\nCreated dummy data:")
    print(f"Number of cells: {n_cells}")
    print(f"UMAP shape: {cds.umap.shape}")
    print(f"Clusters shape: {cds.clusters.shape}")
    
    # Filter test cells
    test_cells = cds.filter_cells("TestCell")
    
    if test_cells.expression_matrix.shape[1] > 0:  # Check second dimension for cells after transpose
        # Plot test cells using the new plot_cells function
        fig = test_cells.plot_cells(color_cells_by="cluster")
        plt.savefig('test_cells_plot.png')
        plt.close()
        print(f"\nFound {test_cells.expression_matrix.shape[1]} test cells")
    else:
        print("No test cells found in the dataset") 