"""
Monocle3 Pseudotime Analysis in Python
=====================================

This module provides pseudotime analysis using the R package 'monocle3' 
through rpy2 interface. It analyzes cell trajectories and identifies
significantly changing genes along pseudotime paths.

Requirements:
- rpy2
- scanpy
- pandas
- numpy
- R with monocle3 package installed

Author: Generated for Biovis2025 project
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import List, Dict, Any, Optional, Tuple
import warnings
import logging
from scipy.sparse import issparse

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    from rpy2 import rinterface_lib
    rpy2_available = True
except ImportError:
    logger.warning("rpy2 not available. Please install rpy2 to use monocle3 functionality.")
    rpy2_available = False


class Monocle3PseudotimeAnalyzer:
    """
    A class for performing pseudotime analysis using monocle3 in Python.
    """
    
    def __init__(self):
        """Initialize the Monocle3PseudotimeAnalyzer."""
        if not rpy2_available:
            raise ImportError("rpy2 is required but not available. Please install rpy2.")
        
        # Set up converters (don't activate globally to avoid conflicts)
        self.pandas_converter = pandas2ri.converter
        self.numpy_converter = numpy2ri.converter
        
        # Import R packages
        try:
            self.monocle3 = importr('monocle3')
            self.base = importr('base')
            self.utils = importr('utils')
            self.stats = importr('stats')
            self.singlecellexperiment = importr('SingleCellExperiment')
            logger.info("Successfully imported R packages")
        except Exception as e:
            logger.error(f"Failed to import R packages: {e}")
            raise ImportError("Required R packages not available. Please install monocle3 in R.")
    
    def convert_adata_to_cds(self, adata: ad.AnnData, cell_ids: Optional[List[str]] = None) -> Any:
        """
        Convert AnnData object to Monocle3 Cell Data Set (CDS).
        
        Parameters:
        -----------
        adata : ad.AnnData
            Annotated data matrix with cells as observations and genes as variables
        cell_ids : List[str], optional
            List of cell IDs to subset the analysis. If None, use all cells.
            
        Returns:
        --------
        cds : R object
            Monocle3 Cell Data Set object
        """
        try:
            # Subset data if cell_ids provided
            if cell_ids is not None:
                # Convert cell_ids to indices
                cell_indices = []
                for cell_id in cell_ids:
                    if cell_id in adata.obs.index:
                        cell_indices.append(adata.obs.index.get_loc(cell_id))
                    else:
                        logger.warning(f"Cell ID {cell_id} not found in adata")
                
                if not cell_indices:
                    raise ValueError("No valid cell IDs found in the data")
                
                # Subset adata
                adata_subset = adata[cell_indices, :].copy()
            else:
                adata_subset = adata.copy()
            
            logger.info(f"Creating CDS with {adata_subset.n_obs} cells and {adata_subset.n_vars} genes")
            
            # Prepare expression matrix (genes x cells for monocle3)
            if issparse(adata_subset.X):
                expr_matrix = adata_subset.X.T.toarray()
            else:
                expr_matrix = adata_subset.X.T
            
            # Convert to R matrix
            with localconverter(ro.default_converter + self.numpy_converter):
                r_expr_matrix = ro.conversion.py2rpy(expr_matrix)
            
            # Prepare gene metadata
            gene_metadata = pd.DataFrame({
                'gene_short_name': adata_subset.var.index.tolist(),
                'gene_id': adata_subset.var.index.tolist()
            })
            gene_metadata.index = adata_subset.var.index
            
            # Prepare cell metadata  
            cell_metadata = adata_subset.obs.copy()
            cell_metadata['cell_id'] = adata_subset.obs.index.tolist()
            
            # Convert metadata to R data frames
            with localconverter(ro.default_converter + self.pandas_converter):
                r_gene_metadata = ro.conversion.py2rpy(gene_metadata)
                r_cell_metadata = ro.conversion.py2rpy(cell_metadata)
            
            # Set row and column names
            ro.r.assign('expr_matrix', r_expr_matrix)
            ro.r.assign('gene_metadata', r_gene_metadata)
            ro.r.assign('cell_metadata', r_cell_metadata)
            
            # Create CDS object in R
            ro.r('''
                rownames(expr_matrix) <- gene_metadata$gene_id
                colnames(expr_matrix) <- cell_metadata$cell_id
                rownames(gene_metadata) <- gene_metadata$gene_id
                rownames(cell_metadata) <- cell_metadata$cell_id
                
                cds <- new_cell_data_set(expr_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_metadata)
            ''')
            
            cds = ro.r('cds')
            logger.info("Successfully created Monocle3 CDS object")
            return cds
            
        except Exception as e:
            logger.error(f"Error converting AnnData to CDS: {e}")
            raise
    
    def preprocess_cds(self, cds: Any, num_dim: int = 50) -> Any:
        """
        Preprocess the CDS object for trajectory analysis.
        
        Parameters:
        -----------
        cds : R object
            Monocle3 Cell Data Set object
        num_dim : int
            Number of dimensions for PCA
            
        Returns:
        --------
        cds : R object
            Preprocessed CDS object
        """
        try:
            logger.info("Preprocessing CDS...")
            ro.r.assign('cds', cds)
            ro.r.assign('num_dim', num_dim)
            
            ro.r('''
                # Preprocess the data
                cds <- preprocess_cds(cds, num_dim = num_dim)
                
                # Reduce dimensions
                cds <- reduce_dimension(cds, reduction_method = "UMAP")
                
                # Cluster cells
                cds <- cluster_cells(cds)
            ''')
            
            preprocessed_cds = ro.r('cds')
            logger.info("CDS preprocessing completed")
            return preprocessed_cds
            
        except Exception as e:
            logger.error(f"Error preprocessing CDS: {e}")
            raise
    
    def learn_trajectory(self, cds: Any, use_partition: bool = True) -> Any:
        """
        Learn the trajectory structure using monocle3.
        
        Parameters:
        -----------
        cds : R object
            Preprocessed CDS object
        use_partition : bool
            Whether to use cell partitions for trajectory learning
            
        Returns:
        --------
        cds : R object
            CDS object with learned trajectory
        """
        try:
            logger.info("Learning trajectory...")
            ro.r.assign('cds', cds)
            ro.r.assign('use_partition', use_partition)
            
            ro.r('''
                # Learn the trajectory graph
                cds <- learn_graph(cds, use_partition = use_partition)
                
                # Order cells in pseudotime (using automatic root selection)
                cds <- order_cells(cds)
            ''')
            
            trajectory_cds = ro.r('cds')
            logger.info("Trajectory learning completed")
            return trajectory_cds
            
        except Exception as e:
            logger.error(f"Error learning trajectory: {e}")
            raise
    
    def extract_pseudotime_data(self, cds: Any) -> Dict[str, Any]:
        """
        Extract pseudotime and trajectory information from CDS.
        
        Parameters:
        -----------
        cds : R object
            CDS object with learned trajectory and pseudotime
            
        Returns:
        --------
        Dict containing:
            - cell_pseudotime: Dict mapping cell IDs to pseudotime values
            - partitions: Dict mapping cell IDs to partition assignments
            - clusters: Dict mapping cell IDs to cluster assignments
        """
        try:
            logger.info("Extracting pseudotime data...")
            ro.r.assign('cds', cds)
            
            # Extract pseudotime, partitions, and clusters
            ro.r('''
                pseudotime_data <- pseudotime(cds)
                partition_data <- partitions(cds)
                cluster_data <- clusters(cds)
                cell_ids <- colnames(cds)
            ''')
            
            # Convert R data to Python
            with localconverter(ro.default_converter + self.numpy_converter):
                pseudotime_values = ro.conversion.rpy2py(ro.r('pseudotime_data'))
                partition_values = ro.conversion.rpy2py(ro.r('partition_data'))
                cluster_values = ro.conversion.rpy2py(ro.r('cluster_data'))
                cell_ids = ro.conversion.rpy2py(ro.r('cell_ids'))
            
            # Create dictionaries
            cell_pseudotime = dict(zip(cell_ids, pseudotime_values))
            partitions = dict(zip(cell_ids, partition_values))
            clusters = dict(zip(cell_ids, cluster_values))
            
            logger.info(f"Extracted pseudotime data for {len(cell_ids)} cells")
            
            return {
                'cell_pseudotime': cell_pseudotime,
                'partitions': partitions,
                'clusters': clusters
            }
            
        except Exception as e:
            logger.error(f"Error extracting pseudotime data: {e}")
            raise
    
    def find_trajectory_paths(self, pseudotime_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Identify distinct trajectory paths based on partitions.
        
        Parameters:
        -----------
        pseudotime_data : Dict
            Dictionary containing pseudotime, partition, and cluster data
            
        Returns:
        --------
        List of dictionaries, each representing a trajectory path with:
            - path_id: Unique identifier for the path
            - cells: List of cell IDs in this path
            - pseudotime_range: Tuple of (min_time, max_time)
            - n_cells: Number of cells in the path
        """
        try:
            partitions = pseudotime_data['partitions']
            cell_pseudotime = pseudotime_data['cell_pseudotime']
            
            # Group cells by partition
            partition_groups = {}
            for cell_id, partition in partitions.items():
                if partition not in partition_groups:
                    partition_groups[partition] = []
                partition_groups[partition].append(cell_id)
            
            # Create path information
            paths = []
            for partition_id, cells in partition_groups.items():
                # Filter out cells with infinite pseudotime
                valid_cells = [cell for cell in cells 
                             if not np.isinf(cell_pseudotime.get(cell, np.inf))]
                
                if len(valid_cells) > 0:
                    pseudotimes = [cell_pseudotime[cell] for cell in valid_cells]
                    path_info = {
                        'path_id': f"path_{partition_id}",
                        'cells': valid_cells,
                        'pseudotime_range': (min(pseudotimes), max(pseudotimes)),
                        'n_cells': len(valid_cells),
                        'partition_id': partition_id
                    }
                    paths.append(path_info)
            
            logger.info(f"Identified {len(paths)} trajectory paths")
            return paths
            
        except Exception as e:
            logger.error(f"Error finding trajectory paths: {e}")
            raise
    
    def find_differential_genes(self, cds: Any, q_value_threshold: float = 0.05) -> pd.DataFrame:
        """
        Find genes that change significantly along pseudotime.
        
        Parameters:
        -----------
        cds : R object
            CDS object with trajectory and pseudotime
        q_value_threshold : float
            Q-value threshold for significance
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with differential gene test results
        """
        try:
            logger.info("Finding differential genes...")
            ro.r.assign('cds', cds)
            ro.r.assign('q_threshold', q_value_threshold)
            
            ro.r('''
                # Test for genes that change along pseudotime
                diff_test_res <- graph_test(cds, neighbor_graph="principal_graph")
                
                # Filter significant genes
                sig_genes <- subset(diff_test_res, q_value < q_threshold)
                sig_genes <- sig_genes[order(sig_genes$q_value),]
            ''')
            
            # Convert results to pandas DataFrame
            with localconverter(ro.default_converter + self.pandas_converter):
                diff_genes_df = ro.conversion.rpy2py(ro.r('sig_genes'))
            
            logger.info(f"Found {len(diff_genes_df)} significantly changing genes")
            return diff_genes_df
            
        except Exception as e:
            logger.error(f"Error finding differential genes: {e}")
            raise
    
    def analyze_gene_expression_along_paths(self, 
                                          adata: ad.AnnData,
                                          paths: List[Dict[str, Any]],
                                          pseudotime_data: Dict[str, Any],
                                          diff_genes_df: pd.DataFrame,
                                          n_top_genes: int = 50,
                                          n_bins: int = 10) -> Dict[str, Any]:
        """
        Analyze gene expression changes along each trajectory path.
        
        Parameters:
        -----------
        adata : ad.AnnData
            Original annotated data
        paths : List[Dict]
            List of trajectory paths
        pseudotime_data : Dict
            Pseudotime data for cells
        diff_genes_df : pd.DataFrame
            Significantly changing genes
        n_top_genes : int
            Number of top genes to analyze per path
        n_bins : int
            Number of pseudotime bins for expression analysis
            
        Returns:
        --------
        Dict containing path-wise gene expression analysis
        """
        try:
            results = {}
            cell_pseudotime = pseudotime_data['cell_pseudotime']
            
            # Get top significant genes
            top_genes = diff_genes_df.head(n_top_genes)['gene_short_name'].tolist()
            
            for path in paths:
                path_id = path['path_id']
                path_cells = path['cells']
                
                logger.info(f"Analyzing {path_id} with {len(path_cells)} cells")
                
                # Get pseudotime values for cells in this path
                path_pseudotimes = np.array([cell_pseudotime[cell] for cell in path_cells])
                
                # Create pseudotime bins
                pseudotime_min, pseudotime_max = path_pseudotimes.min(), path_pseudotimes.max()
                bin_edges = np.linspace(pseudotime_min, pseudotime_max, n_bins + 1)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                
                # Analyze gene expression in each bin
                path_gene_data = []
                
                for gene in top_genes:
                    if gene not in adata.var.index:
                        continue
                    
                    gene_idx = adata.var.index.get_loc(gene)
                    
                    # Get cell indices in adata
                    cell_indices = [adata.obs.index.get_loc(cell) for cell in path_cells 
                                  if cell in adata.obs.index]
                    
                    if len(cell_indices) == 0:
                        continue
                    
                    # Get gene expression for path cells
                    if issparse(adata.X):
                        gene_expr = adata.X[cell_indices, gene_idx].toarray().flatten()
                    else:
                        gene_expr = adata.X[cell_indices, gene_idx]
                    
                    # Calculate expression in each pseudotime bin
                    bin_expressions = []
                    bin_cell_counts = []
                    
                    for i in range(n_bins):
                        mask = (path_pseudotimes >= bin_edges[i]) & (path_pseudotimes < bin_edges[i + 1])
                        if i == n_bins - 1:  # Include the last point in the final bin
                            mask = (path_pseudotimes >= bin_edges[i]) & (path_pseudotimes <= bin_edges[i + 1])
                        
                        if mask.sum() > 0:
                            mean_expr = float(gene_expr[mask].mean())
                            cell_count = int(mask.sum())
                        else:
                            mean_expr = 0.0
                            cell_count = 0
                        
                        bin_expressions.append(mean_expr)
                        bin_cell_counts.append(cell_count)
                    
                    # Get gene statistics from differential analysis
                    gene_stats = diff_genes_df[diff_genes_df['gene_short_name'] == gene]
                    if len(gene_stats) > 0:
                        morans_i = float(gene_stats.iloc[0]['morans_I']) if 'morans_I' in gene_stats.columns else None
                        q_value = float(gene_stats.iloc[0]['q_value']) if 'q_value' in gene_stats.columns else None
                    else:
                        morans_i = None
                        q_value = None
                    
                    path_gene_data.append({
                        'gene': gene,
                        'morans_I': morans_i,
                        'q_value': q_value,
                        'pseudotime_bins': [float(x) for x in bin_centers],
                        'expression_values': bin_expressions,
                        'bin_cell_counts': bin_cell_counts,
                        'total_expression_change': float(max(bin_expressions) - min(bin_expressions))
                    })
                
                results[path_id] = {
                    'path_info': path,
                    'n_cells': len(path_cells),
                    'pseudotime_range': (float(pseudotime_min), float(pseudotime_max)),
                    'genes': path_gene_data
                }
            
            logger.info(f"Completed analysis for {len(results)} paths")
            return results
            
        except Exception as e:
            logger.error(f"Error analyzing gene expression along paths: {e}")
            raise
    
    def run_full_analysis(self, 
                         adata: ad.AnnData,
                         cell_ids: Optional[List[str]] = None,
                         num_dim: int = 50,
                         q_value_threshold: float = 0.05,
                         n_top_genes: int = 50,
                         n_bins: int = 10) -> Dict[str, Any]:
        """
        Run complete monocle3 pseudotime analysis pipeline.
        
        Parameters:
        -----------
        adata : ad.AnnData
            Annotated data matrix
        cell_ids : List[str], optional
            Specific cell IDs to analyze
        num_dim : int
            Number of PCA dimensions
        q_value_threshold : float
            Significance threshold for differential genes
        n_top_genes : int
            Number of top genes to analyze per path
        n_bins : int
            Number of pseudotime bins
            
        Returns:
        --------
        Dict containing complete analysis results:
            - paths: Trajectory path information
            - pseudotime_data: Cell pseudotime assignments
            - differential_genes: Significantly changing genes
            - path_analysis: Gene expression analysis for each path
        """
        try:
            logger.info("Starting full monocle3 pseudotime analysis...")
            
            # Step 1: Convert to CDS
            cds = self.convert_adata_to_cds(adata, cell_ids)
            
            # Step 2: Preprocess
            cds = self.preprocess_cds(cds, num_dim)
            
            # Step 3: Learn trajectory
            cds = self.learn_trajectory(cds)
            
            # Step 4: Extract pseudotime data
            pseudotime_data = self.extract_pseudotime_data(cds)
            
            # Step 5: Find trajectory paths
            paths = self.find_trajectory_paths(pseudotime_data)
            
            # Step 6: Find differential genes
            diff_genes_df = self.find_differential_genes(cds, q_value_threshold)
            
            # Step 7: Analyze gene expression along paths
            path_analysis = self.analyze_gene_expression_along_paths(
                adata, paths, pseudotime_data, diff_genes_df, n_top_genes, n_bins
            )
            
            logger.info("Full analysis completed successfully!")
            
            return {
                'paths': paths,
                'pseudotime_data': pseudotime_data,
                'differential_genes': diff_genes_df,
                'path_analysis': path_analysis,
                'analysis_parameters': {
                    'num_dim': num_dim,
                    'q_value_threshold': q_value_threshold,
                    'n_top_genes': n_top_genes,
                    'n_bins': n_bins,
                    'n_input_cells': len(cell_ids) if cell_ids else adata.n_obs
                }
            }
            
        except Exception as e:
            logger.error(f"Error in full analysis: {e}")
            raise


def run_monocle3_analysis(adata: ad.AnnData, 
                         cell_ids: Optional[List[str]] = None,
                         **kwargs) -> Dict[str, Any]:
    """
    Convenience function to run monocle3 pseudotime analysis.
    
    Parameters:
    -----------
    adata : ad.AnnData
        Annotated data matrix with cells as observations and genes as variables
    cell_ids : List[str], optional
        List of specific cell IDs to analyze. If None, analyze all cells.
    **kwargs : dict
        Additional parameters for the analysis
        
    Returns:
    --------
    Dict containing complete analysis results
    """
    analyzer = Monocle3PseudotimeAnalyzer()
    return analyzer.run_full_analysis(adata, cell_ids, **kwargs)


# Example usage
if __name__ == "__main__":
    # Example of how to use the analyzer
    print("Monocle3 Pseudotime Analyzer")
    print("============================")
    print("This module provides pseudotime analysis using R's monocle3 package.")
    print("\nExample usage:")
    print("""
    import scanpy as sc
    from monocle3_pseudotime import run_monocle3_analysis
    
    # Load your data
    adata = sc.read_h5ad('your_data.h5ad')
    
    # Define cell IDs to analyze (optional)
    cell_ids = ['cell_1', 'cell_2', 'cell_3', ...]
    
    # Run analysis
    results = run_monocle3_analysis(adata, cell_ids)
    
    # Access results
    paths = results['paths']
    pseudotime_data = results['pseudotime_data']
    differential_genes = results['differential_genes']
    path_analysis = results['path_analysis']
    """)