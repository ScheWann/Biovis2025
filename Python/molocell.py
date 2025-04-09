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
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

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
        
        # Principal graph for trajectory
        self.principal_graph = {}
        self.principal_graph_aux = {}
        
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
            raise ValueError(f"No dimensionality reduction for {reduction_method} calculated. "
                            f"Please run reduce_dimension with reduction_method = {reduction_method} "
                            f"before attempting to plot.")
        elif reduction_method == "PCA" and self.pca is None:
            raise ValueError(f"No dimensionality reduction for {reduction_method} calculated. "
                            f"Please run reduce_dimension with reduction_method = {reduction_method} "
                            f"before attempting to plot.")
            
        # Get reduction coordinates
        if reduction_method == "UMAP":
            coords = self.umap
        elif reduction_method == "PCA":
            coords = self.pca
            
        # Check if x and y are valid dimensions
        if coords.shape[1] < max(x, y):
            raise ValueError(f"x and/or y is too large. x and y must be dimensions in reduced dimension space.")
            
        # Check if color_cells_by is valid
        if color_cells_by is not None:
            valid_color_by = ["cluster", "partition", "pseudotime"]
            if color_cells_by not in valid_color_by and color_cells_by not in self.colData.columns:
                raise ValueError(f"color_cells_by must be one of {valid_color_by} or a column in the colData table.")
                
            # Check if pseudotime exists if color_cells_by is pseudotime
            if color_cells_by == "pseudotime":
                if not hasattr(self, 'pseudotime') or self.pseudotime is None:
                    raise ValueError(f"No pseudotime for {reduction_method} calculated. "
                                    f"Please run order_cells with reduction_method = {reduction_method} "
                                    f"before attempting to color by pseudotime.")
                    
        # Check if group_cells_by is valid
        if group_cells_by not in ["cluster", "partition"] and group_cells_by not in self.colData.columns:
            raise ValueError(f"group_cells_by must be one of ['cluster', 'partition'] or a column in the colData table.")
            
        # Check if trajectory graph exists
        if show_trajectory_graph and reduction_method not in self.principal_graph:
            print("No trajectory to plot. Has learn_graph() been called yet?")
            show_trajectory_graph = False
            
        # Check if principal points can be labeled
        if label_principal_points and reduction_method not in self.principal_graph:
            print("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
            label_principal_points = False
            
        # If label_principal_points is True, disable other labels
        if label_principal_points:
            label_branch_points = False
            label_leaves = False
            label_roots = False
            
        # Create data frame for plotting
        data_df = pd.DataFrame({
            'data_dim_1': coords[:, x-1],
            'data_dim_2': coords[:, y-1],
            'sample_name': self.colData.index
        })
        
        # Add colData to data_df
        for col in self.colData.columns:
            data_df[col] = self.colData[col].values
            
        # Add cell group based on group_cells_by
        if group_cells_by == "cluster" and self.clusters is not None:
            data_df['cell_group'] = self.clusters
        elif group_cells_by == "partition" and hasattr(self, 'partitions') and self.partitions is not None:
            data_df['cell_group'] = self.partitions
        elif group_cells_by in self.colData.columns:
            data_df['cell_group'] = self.colData[group_cells_by].values
            
        # Add cell color based on color_cells_by
        if color_cells_by == "cluster" and self.clusters is not None:
            data_df['cell_color'] = self.clusters
        elif color_cells_by == "partition" and hasattr(self, 'partitions') and self.partitions is not None:
            data_df['cell_color'] = self.partitions
        elif color_cells_by == "pseudotime" and hasattr(self, 'pseudotime') and self.pseudotime is not None:
            data_df['cell_color'] = self.pseudotime
        elif color_cells_by in self.colData.columns:
            data_df['cell_color'] = self.colData[color_cells_by].values
            
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Plot cells
        if genes is not None:
            # Handle gene expression coloring
            # This is a simplified version - in the full implementation, you would need to
            # calculate gene expression values and normalize them
            if norm_method == "size_only":
                # Simplified gene expression coloring
                scatter = ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                                    c=data_df['cell_color'], cmap='viridis',
                                    s=cell_size*100, alpha=alpha)
                plt.colorbar(scatter, label='Expression')
            else:
                # Log-transformed gene expression coloring
                scatter = ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                                    c=np.log10(data_df['cell_color'] + min_expr), cmap='viridis',
                                    s=cell_size*100, alpha=alpha)
                plt.colorbar(scatter, label='log10(Expression)')
        else:
            # Color by cell attribute
            if color_cells_by in ["cluster", "partition"]:
                if data_df['cell_color'] is None:
                    # Default coloring if no cluster/partition data
                    ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                              color='gray', s=cell_size*100, alpha=alpha)
                    print(f"{color_cells_by} not found in colData(cds), cells will not be colored")
                else:
                    # Color by cluster/partition
                    scatter = ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                                        c=data_df['cell_color'], cmap='viridis',
                                        s=cell_size*100, alpha=alpha)
                    plt.colorbar(scatter, label=color_cells_by)
            elif pd.api.types.is_numeric_dtype(data_df['cell_color']):
                # Color by numeric value
                scatter = ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                                    c=data_df['cell_color'], cmap='viridis',
                                    s=cell_size*100, alpha=alpha)
                plt.colorbar(scatter, label=color_cells_by)
            else:
                # Color by categorical value
                scatter = ax.scatter(data_df['data_dim_1'], data_df['data_dim_2'], 
                                    c=pd.Categorical(data_df['cell_color']).codes, cmap='viridis',
                                    s=cell_size*100, alpha=alpha)
                plt.colorbar(scatter, label=color_cells_by)
                
        # Add trajectory graph if available
        if show_trajectory_graph and reduction_method in self.principal_graph:
            # This is a simplified version - in the full implementation, you would need to
            # extract the graph edges and nodes from the principal graph
            pass
            
        # Add labels for cell groups
        if label_cell_groups and data_df['cell_color'] is not None:
            if group_cells_by == "cluster" and data_df['cell_group'] is not None:
                # Group by cluster and color
                for cluster in data_df['cell_group'].unique():
                    cluster_data = data_df[data_df['cell_group'] == cluster]
                    for color in cluster_data['cell_color'].unique():
                        color_data = cluster_data[cluster_data['cell_color'] == color]
                        if len(color_data) > 0:
                            center_x = color_data['data_dim_1'].median()
                            center_y = color_data['data_dim_2'].median()
                            ax.text(center_x, center_y, f'{color}', 
                                    fontsize=group_label_size*5)
            else:
                # Group by color only
                for color in data_df['cell_color'].unique():
                    color_data = data_df[data_df['cell_color'] == color]
                    if len(color_data) > 0:
                        center_x = color_data['data_dim_1'].median()
                        center_y = color_data['data_dim_2'].median()
                        ax.text(center_x, center_y, f'{color}', 
                                fontsize=group_label_size*5)
                
        # Set title and labels
        ax.set_title(f'{reduction_method} Plot')
        ax.set_xlabel(f'{reduction_method} {x}')
        ax.set_ylabel(f'{reduction_method} {y}')
        
        # Apply monocle theme
        monocle_theme_opts()
        
        return fig

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
    
    # Create cell type assignments for all cells
    # Assign random cell types to all cells for demonstration
    cell_types = ['Type_A', 'Type_B', 'Type_C', 'Type_D', 'Type_E']
    cds.colData['assigned_cell_type'] = np.random.choice(cell_types, size=len(cds.colData))
    
    # Create dummy clusters and UMAP coordinates with correct dimensions
    n_cells = cds.expression_matrix.shape[0]  # Use first dimension for cells
    cds.clusters = np.random.randint(0, 5, size=n_cells)
    cds.umap = np.random.randn(n_cells, 2)
    cds.reduced = True
    
    print(f"\nCreated dummy data:")
    print(f"Number of cells: {n_cells}")
    print(f"UMAP shape: {cds.umap.shape}")
    print(f"Clusters shape: {cds.clusters.shape}")
    
    # Plot all cells
    print("\nPlotting all cells...")
    fig = cds.plot_cells(color_cells_by="assigned_cell_type")
    plt.savefig('all_cells_plot.png')
    plt.close()
    
    # Also plot cells colored by cluster
    print("Plotting cells colored by cluster...")
    fig = cds.plot_cells(color_cells_by="cluster")
    plt.savefig('all_cells_cluster_plot.png')
    plt.close()
    
    print(f"\nPlots saved as 'all_cells_plot.png' and 'all_cells_cluster_plot.png'") 