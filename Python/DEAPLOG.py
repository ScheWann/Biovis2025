import scanpy as sc
import anndata
import os
from scipy import sparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import umap
import scipy.stats as stats
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit

def sample_data(adata, sample_ratio=0.01, random_state=42):
    """
    Sample a subset of cells from the dataset
    
    Parameters:
    -----------
    adata : AnnData
        Input dataset
    sample_ratio : float
        Ratio of cells to sample (default: 0.01 = 1%)
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    AnnData
        Sampled dataset
    """
    print(f"Sampling {sample_ratio*100}% of cells...")
    n_cells = int(adata.n_obs * sample_ratio)
    np.random.seed(random_state)
    sampled_indices = np.random.choice(adata.n_obs, n_cells, replace=False)
    return adata[sampled_indices]

def preprocess_data(adata, target_sum=1e4):
    """
    Preprocess the data with normalization and scaling
    """
    print("Preprocessing data...")
    # Basic preprocessing
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)
    
    # Compute PCA and UMAP
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    
    return adata

def power_law(x, a, b):
    """Power law function for curve fitting"""
    return a * np.power(x, b)

def get_highly_cells_for_each_gene(rdata_df, gene, power=11):
    """Get cells with high expression for a given gene"""
    mean_expr = rdata_df[gene].mean()
    std_expr = rdata_df[gene].std()
    threshold = mean_expr + (power * std_expr)
    highly_expressed_cells = rdata_df[rdata_df[gene] > threshold].index.tolist()
    return highly_expressed_cells

def Fisher_test_for_each_gene(rdata_df, cell_sets, num_allCells, gene, power=8):
    """Perform Fisher's exact test for a given gene"""
    gene_data = rdata_df[gene]
    
    gene_ratio = {}
    gene_pv = {}
    gene_qv = {}
    gene_score = {}
    gene_means = {}
    
    for ct, cells in cell_sets.items():
        high_expression_cells = gene_data[cells] > 0
        
        num_high_expr = sum(high_expression_cells)
        num_low_expr = len(high_expression_cells) - num_high_expr
        num_cells_in_type = len(cells)
        num_high_expr_type = sum(high_expression_cells)
        num_low_expr_type = num_cells_in_type - num_high_expr_type
        
        contingency_table = [
            [num_high_expr_type, num_low_expr_type],
            [num_high_expr, num_low_expr]
        ]
        
        if any(val < 0 for row in contingency_table for val in row):
            continue
        
        _, p_value = stats.fisher_exact(contingency_table)
        
        gene_ratio[ct] = num_high_expr / len(cells)
        gene_pv[ct] = p_value
        gene_qv[ct] = p_value
        gene_score[ct] = p_value
        gene_means[ct] = gene_data[cells].mean()
    
    return gene_ratio, gene_pv, gene_qv, gene_score, gene_means

def analyze_gene_expression(adata, genes_to_analyze):
    """Analyze expression patterns of specific genes"""
    print("\nAnalyzing gene expression patterns...")
    rdata_df = adata.to_df()
    
    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    # Keep track of genes with enough expression
    genes_with_expression = []
    
    for gene in genes_to_analyze:
        print(f"\nAnalyzing gene: {gene}")
        g_highly_cells = get_highly_cells_for_each_gene(rdata_df=rdata_df, gene=gene, power=11)
        
        # Extract gene expression data and filter out values <= 0
        gene_list = rdata_df[gene]
        gene_filter = gene_list.loc[gene_list > 0]
        gene_sort_df_filter = gene_filter.sort_values(ascending=False)
        
        if len(gene_sort_df_filter) > 1:  # Only fit if we have more than one point
            genes_with_expression.append(gene)
            # Prepare the data for curve fitting
            x_data = np.arange(1, len(gene_sort_df_filter) + 1)
            y_data = gene_sort_df_filter.values
            
            try:
                # Initial parameter guess
                p0 = [np.max(y_data), -1.0]
                
                # Fit the curve to the data using the power law
                popt, _ = curve_fit(power_law, x_data, y_data, p0=p0, maxfev=5000)
                
                # Plot results
                plt.figure(figsize=(10, 6))
                plt.plot(x_data, y_data, label='Data', color='blue')
                plt.plot(x_data, power_law(x_data, *popt), 
                        label=f'Fitted curve: a={popt[0]:.2f}, b={popt[1]:.2f}', color='red')
                plt.xlabel('Index of ranked cells')
                plt.ylabel(f'Expression level of {gene}')
                if len(g_highly_cells) > 0:
                    plt.vlines(x=len(g_highly_cells), ymin=0, ymax=gene_sort_df_filter.max() - 0.5,
                             colors='black', linestyles='dashed')
                    plt.text(x=len(g_highly_cells), y=gene_sort_df_filter.max() - 0.5,
                            s=f'x={len(g_highly_cells)}')
                plt.grid(True)
                plt.legend()
                plt.savefig(f'figures/gene_expression_{gene}.png', bbox_inches='tight', dpi=300)
                plt.close()
            except RuntimeError as e:
                print(f"Could not fit curve for gene {gene}: {str(e)}")
        else:
            print(f"Not enough non-zero expression values for gene {gene}")
        
        # Add highly expressed cells information to adata
        adata.obs[f'{gene}HighlyCells'] = 'False'
        adata.obs.loc[g_highly_cells, f'{gene}HighlyCells'] = 'True'
        adata.obs[f'{gene}HighlyCells'] = pd.Categorical(
            adata.obs[f'{gene}HighlyCells'],
            categories=['True', 'False'],
            ordered=True
        )
        adata.uns[f'{gene}HighlyCells_colors'] = ['#ff3366', '#cccccc']
    
    return genes_with_expression

def run_differential_expression_analysis(adata):
    """Run differential expression analysis"""
    print("\nRunning differential expression analysis...")
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=10)
    highly_variable_genes = adata.var_names[adata.var['highly_variable']].tolist()
    print("\nHighly variable genes:", highly_variable_genes)
    
    # Prepare cell type data
    num_allCells = len(adata.obs_names)
    print('Structuring the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs['cell_type_id'])
    cell_type_index = pd.Categorical(adata.obs['cell_type_id']).categories
    cell_sets = {ct: list(adata_cell_type_df.loc[adata_cell_type_df['cell_type_id'] == ct, :].index) 
                 for ct in cell_type_index}
    
    # Initialize results DataFrame
    all_gene_Fisher_r = pd.DataFrame(index=cell_type_index)
    
    # Analyze each gene
    rdata_df = adata.to_df()
    for gene in highly_variable_genes:
        gene_ratio, gene_pv, gene_qv, gene_score, gene_means = Fisher_test_for_each_gene(
            rdata_df=rdata_df,
            cell_sets=cell_sets,
            num_allCells=num_allCells,
            gene=gene,
            power=8
        )
        
        # Store results
        Fisher_r = pd.DataFrame({
            f'{gene}_ratio': gene_ratio,
            f'{gene}_pValue': gene_pv,
            f'{gene}_qValue': gene_qv,
            f'{gene}_score': gene_score
        })
        
        all_gene_Fisher_r = all_gene_Fisher_r.join(Fisher_r)
    
    # Save results
    os.makedirs('./write', exist_ok=True)
    all_gene_Fisher_r.to_csv('./write/adata_Genes_HighlyCells_Fisher_r.csv')
    
    return all_gene_Fisher_r, highly_variable_genes

def visualize_results(adata, genes_to_analyze):
    """Visualize analysis results"""
    print("\nVisualizing results...")
    
    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    # Set figure parameters
    sc.settings.set_figure_params(dpi=100, frameon=False)
    
    # Plot UMAP with cell types
    plt.figure(figsize=(8, 6))
    sc.pl.umap(adata, color='cell_type_id', title='Cell Types', show=False)
    plt.savefig('figures/umap_cell_types.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    if genes_to_analyze:
        # Plot gene expression on UMAP
        plt.figure(figsize=(15, 10))
        sc.pl.umap(adata, color=genes_to_analyze[:5], ncols=3, show=False)
        plt.savefig('figures/umap_gene_expression.png', bbox_inches='tight', dpi=300)
        plt.close()
        
        # Create clustermap
        adata_markers_df = adata[:, genes_to_analyze].to_df()
        adata_markers_df['cell_type_id'] = adata.obs['cell_type_id']
        adata_markers_df = adata_markers_df.sort_values(by='cell_type_id')
        
        # Prepare colors for cell types
        if 'cell_type_id_colors' not in adata.uns:
            louvain_tag = adata.obs['cell_type_id'].cat.categories
            louvain_colors = sc.pl.palettes.vega_10
            if len(louvain_colors) < len(louvain_tag):
                louvain_colors = louvain_colors * (len(louvain_tag) // len(louvain_colors)) + \
                                louvain_colors[:len(louvain_tag) % len(louvain_colors)]
            adata.uns['cell_type_id_colors'] = louvain_colors[:len(louvain_tag)]
        
        cluster_colors = pd.Series(
            adata.uns['cell_type_id_colors'],
            index=adata.obs['cell_type_id'].cat.categories
        ).loc[adata_markers_df['cell_type_id']]
        
        # Create and save clustermap
        sns.clustermap(
            adata_markers_df.drop('cell_type_id', axis=1).transpose(),
            metric="correlation",
            cmap='viridis',
            row_cluster=False,
            col_cluster=False,
            robust=True,
            xticklabels=False,
            figsize=(8, 8),
            yticklabels=True,
            z_score=0,
            col_colors=[cluster_colors]
        ).savefig('figures/gene_expression_heatmap.png', bbox_inches='tight', dpi=300)
    else:
        print("No genes with sufficient expression to visualize")

def main():
    """Main function to run the analysis"""
    # Load data
    data_path = "Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5"
    print(f"Loading data from: {data_path}")
    adata = sc.read_h5ad(data_path)
    
    # Filter unwanted cell types
    adata = adata[~adata.obs['cell_type'].isin(['Endothelial cell', 'Eosinophil/Basophil/Mast'])]
    
    # Create cell type mapping
    category_mapping = {
        'cd19+cd20+ b': 'B', 'cd4+ t cells': 'cycl-B', 'cd8+ t cells': 'EarlyLymp',
        'cdc': 'EarlyMyel', 'cms1': 'Endo', 'cms2': 'E/B/M', 'cms3': 'Ery',
        'enteric glial cells': 'cycl-Ery', 'gamma delta t cells': 'HSC/MPP',
        'iga+ plasma': 'Late-Ery', 'igg+ plasma': 'MPP', 'intermediate': 'Lymp-MPP',
        'lymphatic ecs': 'Myel-MPP', 'mast cells': 'Mega', 'mature enterocytes type 1': 'ProMega',
        'mature enterocytes type 2': 'Mono', 'myofibroblasts': 'neutro-Myel', 'nk cells': 'cycl-Myel',
        'pericytes': 'mono-Myel', 'pro-inflammatory': 'PreDC', 'proliferating': 'cycl-PreDC',
        'regulatory t cells': 'B', 'smooth muscle cells': 'B', 'spp1+': 'B', 'stalk-like ecs': 'B',
        'stem-like/ta': 'B', 'stromal 1': 'B', 'stromal 2': 'B', 'stromal 3': 'B',
        't follicular helper cells': 'B', 't helper 17 cells': 'B', 'tip-like ecs': 'B', 'unknown': 'B'
    }
    
    # Apply cell type mapping
    adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'])
    adata.obs['cell_type_id'] = pd.Categorical(adata.obs['cell_type'].map(category_mapping))
    
    # Sample data
    adata = sample_data(adata, sample_ratio=0.01)
    
    # Preprocess data
    adata = preprocess_data(adata)
    
    # Run differential expression analysis
    fisher_results, genes_to_analyze = run_differential_expression_analysis(adata)
    
    # Analyze gene expression patterns
    analyze_gene_expression(adata, genes_to_analyze)
    
    # Visualize results
    visualize_results(adata, genes_to_analyze)
    
    print("\nAnalysis completed!")

if __name__ == "__main__":
    main()
