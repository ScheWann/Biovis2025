import scanpy as sc
import anndata
import os
from scipy import sparse
import numpy as np
import pandas as pd
import json
import sys
from typing import Dict, Any
sc.settings.verbosity = 3
sc.logging.print_version_and_date()
from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import umap
import scanpy as sc
import os
import sys
sys.executable
import sympy
import pandas as pd
import scipy.stats as stats
from sklearn.decomposition import PCA

# Set the data path directly
data_path = "D:/Research/vis2025/siyuan/Biovis2025/Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5"

try:
    print(f"Reading data from: {data_path}")
    # Read the AnnData file
    adata = sc.read_h5ad(data_path)
    print(f"Successfully loaded data with shape: {adata.shape}")
    
    # Sample 1% of cells
    n_cells = int(len(adata) * 0.01)
    print(f"Sampling {n_cells} cells from {len(adata)} total cells...")
    
    adata_sampled = adata[np.random.choice(adata.shape[0], n_cells, replace=False), :]
    
    # Perform PCA with fewer components
    print("Performing PCA...")
    sc.pp.pca(adata_sampled, n_comps=20, svd_solver='randomized')
    
    # Run UMAP with optimized parameters
    print("Running UMAP...")
    sc.pp.neighbors(adata_sampled, n_neighbors=5, use_rep='X_pca')
    sc.tl.umap(adata_sampled)
    
    # Get UMAP coordinates
    umap_coords = adata_sampled.obsm['X_umap']
    
    # Get cell types
    cell_types = adata_sampled.obs['cell_type'].tolist()
    
    # Get gene expression data for top marker genes
    genes_of_interest = [
        'CXCL8', 'S100A9', 'SAT1', 'FTH1', 'S100A8', 'G0S2', 
        'TMSB4X', 'LITAF', 'COL1A1', 'MT-CYB'
    ]
    gene_expression = {}
    print("Processing gene expression data...")
    for gene in genes_of_interest:
        if gene in adata_sampled.var_names:
            gene_expression[gene] = adata_sampled[:, gene].X.toarray().flatten().tolist()
    
    # Prepare data for D3.js
    d3_data = {
        'umap_coordinates': {
            'x': umap_coords[:, 0].tolist(),
            'y': umap_coords[:, 1].tolist()
        },
        'cell_types': cell_types,
        'gene_expression': gene_expression,
        'metadata': {
            'total_cells': len(adata_sampled),
            'unique_cell_types': list(set(cell_types))
        }
    }
    
    print("Processing completed successfully!")
    print(json.dumps(d3_data))
    
except FileNotFoundError:
    print(f"Error: Data file not found: {data_path}", file=sys.stderr)
    error_data = {
        'error': f"Data file not found: {data_path}",
        'error_type': 'FileNotFoundError',
        'data_path': data_path
    }
    print(json.dumps(error_data))
    sys.exit(1)
except Exception as e:
    print(f"Error occurred: {str(e)}", file=sys.stderr)
    error_data = {
        'error': str(e),
        'error_type': type(e).__name__,
        'data_path': data_path
    }
    print(json.dumps(error_data))
    sys.exit(1)
