# Install required packages
# !pip install deaplog scanpy numpy pandas matplotlib seaborn

# Import necessary libraries
import deaplog
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white')

# Load sample data
adata = sc.datasets.pbmc3k()

# Basic preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.scale(adata)

# Run PCA
sc.tl.pca(adata)

# Compute neighborhood graph
sc.neighbors(adata)

# Run clustering
sc.tl.leiden(adata)

# Run UMAP
sc.tl.umap(adata)

# Create rdata (normalized counts)
rdata = adata.copy()
sc.pp.normalize_total(rdata)
sc.pp.log1p(rdata)
sc.pp.scale(rdata)

# Find genes differentially expressed in single cell types
markers_unique = deaplog.get_DEG_uniq(
    rdata, 
    adata,
    group_key='leiden',
    power=11,
    ratio=0.2,
    p_threshold=0.01,
    q_threshold=0.05
)

# Display results for unique markers
print("Number of unique marker genes found:", len(markers_unique))
print("\nFirst few marker genes:")
print(markers_unique.head())

# Find genes differentially expressed in multiple cell types
markers_multi = deaplog.get_DEG_multi(
    rdata, 
    adata,
    group_key='leiden',
    power=11,
    ratio=0.2,
    p_threshold=0.01,
    q_threshold=0.05
)

# Display results for multi-markers
print("\nNumber of multi-marker genes found:", len(markers_multi))
print("\nFirst few multi-marker genes:")
print(markers_multi.head())

# Calculate gene locations and pseudotime
genes_location = deaplog.get_genes_location_pseudotime(
    rdata, 
    adata,
    group_key='leiden',
    power=11,
    gene_matrix=markers_unique,
    obsm='X_umap'
)

# Display results for gene locations
print("\nGene locations calculated successfully")
print("\nFirst few gene locations:")
print(genes_location.head())

# Visualization
# Plot UMAP with cell clusters
sc.pl.umap(adata, color='leiden', title='Cell Clusters', show=True)

# Plot expression of top marker genes
top_genes = markers_unique.index[:5]
sc.pl.umap(adata, color=top_genes, ncols=2, title='Top Marker Genes', show=True) 