#!/path/to/deaplog_env/bin/python

import scanpy as sc
import numpy as np
import pandas as pd
from deaplog.deaplog import get_DEG_uniq, get_DEG_multi
import matplotlib.pyplot as plt

def load_and_preprocess_data(subsample_ratio=0.001):
    """Load and preprocess the data"""
    print("Loading and preprocessing data...")
    
    # Load data
    data_path = "Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5"
    print(f"Loading data from: {data_path}")
    adata = sc.read_h5ad(data_path)
    umi_counts = adata.to_df()
    
    # Apply cell type to adata
    adata.obs["cell_type"] = adata.obs["cell_type"]
    
    # Subsample if needed
    if subsample_ratio < 1:
        n_cells = int(adata.n_obs * subsample_ratio)
        print(f"Subsampling to {n_cells} cells...")
        sc.pp.subsample(adata, n_obs=n_cells, random_state=42)
    
    print("Preprocessing data...")
    # Basic preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, n_comps=50)
    
    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    
    # Try to run clustering
    try:
        print("Running Leiden clustering...")
        sc.tl.leiden(adata, resolution=0.8)
    except BaseException as e:
        print(f"Warning: Leiden clustering failed: {str(e)}")
        print("Trying Louvain clustering instead...")
        try:
            sc.tl.louvain(adata, resolution=0.8)
            adata.obs['leiden'] = adata.obs['louvain']
        except Exception as e:
            print(f"Warning: Louvain clustering also failed: {str(e)}")
            print("Using random cluster labels for testing...")
            adata.obs['leiden'] = pd.Categorical(
                np.random.randint(0, 5, size=adata.n_obs).astype(str)
            )
    
    # Compute UMAP
    sc.tl.umap(adata)
    
    return adata

def run_differential_expression_analysis(adata):
    """Run differential expression analysis"""
    print("\nRunning differential expression analysis...")
    
    # Print expression data statistics
    print("\nExpression data statistics:")
    if isinstance(adata.X, np.ndarray):
        data_array = adata.X
    else:
        data_array = adata.X.toarray()
    print(f"Mean expression: {np.mean(data_array)}")
    print(f"Max expression: {np.max(data_array)}")
    print(f"Number of non-zero values: {np.count_nonzero(data_array)}")
    
    # Find unique marker genes
    print("\nFinding unique marker genes...")
    markers_unique = get_DEG_uniq(adata, adata, group_key='leiden', power=5, ratio=0.05, p_threshold=0.1, q_threshold=0.2)
    
    # Find shared marker genes
    print("\nFinding shared marker genes...")
    markers_multi = get_DEG_multi(adata, adata, group_key='leiden', power=5, ratio=0.05, p_threshold=0.1, q_threshold=0.2)
    
    # Display results
    print("\nAnalysis Results:")
    print(f"Found {len(markers_unique)} unique marker genes")
    print(f"Found {len(markers_multi)} shared marker genes")
    
    if not markers_unique.empty:
        print("\nTop 5 unique marker genes by score:")
        print(markers_unique.nlargest(5, 'score')[['gene_name', 'cell_type', 'score']])
    
    return markers_unique, markers_multi

def visualize_results(adata, markers_unique, markers_shared):
    """Visualize the results"""
    print("\nVisualizing results...")
    
    # Plot UMAP with cell clusters
    print("\nPlotting UMAP with cell clusters...")
    sc.settings.set_figure_params(dpi=100, frameon=False)
    
    # Ensure leiden is categorical
    adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'])
    
    # Plot clusters
    plt.figure(figsize=(8, 6))
    sc.pl.umap(adata, color='leiden', title='Cell Clusters', show=False)
    plt.savefig('Python/figures/umap_clusters.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    # Plot UMAP based on cell types
    print("\nPlotting UMAP based on cell types...")
    plt.figure(figsize=(8, 6))
    sc.pl.umap(adata, color='cell_type', title='Cell Types', show=False)
    plt.savefig('Python/figures/umap_cell_types.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    # Plot top marker genes
    if not markers_unique.empty:
        print("\nPlotting top marker genes...")
        top_genes = markers_unique.nlargest(5, 'score')['gene_name'].tolist()
        print(f"\nTop marker genes to visualize: {', '.join(top_genes)}")
        
        # Plot marker genes
        plt.figure(figsize=(15, 10))
        sc.pl.umap(adata, color=top_genes, ncols=3, show=False)
        plt.savefig('Python/figures/umap_marker_genes.png', bbox_inches='tight', dpi=300)
        plt.close()
    else:
        print("\nNo marker genes found to visualize.")

def main():
    """Main function to run the analysis"""
    # Load and preprocess data
    adata = load_and_preprocess_data(subsample_ratio=0.001)
    
    # Run differential expression analysis
    markers_unique, markers_multi = run_differential_expression_analysis(adata)
    
    # Visualize results
    visualize_results(adata, markers_unique, markers_multi)
    
    print("\nAnalysis completed!")

if __name__ == "__main__":
    main()