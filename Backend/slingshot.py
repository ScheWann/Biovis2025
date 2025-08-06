"""
Slingshot Trajectory and Gene Analysis Module

This module contains functions extracted from the slingshot_real.ipynb notebook
for trajectory inference using Slingshot and gene expression analysis along trajectories.
Visualization functions have been removed for backend compatibility.
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import warnings
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import tempfile
import os
from scipy import stats

warnings.filterwarnings("ignore")


def run_slingshot_via_rpy2_improved(
    adata,
    cluster_key="leiden",
    embedding_key="X_umap",
    start_cluster=None,
    end_clusters=None,
):
    """
    Run Slingshot trajectory inference using rpy2 with improved error handling.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix
    cluster_key : str
        The key in adata.obs for cluster labels
    embedding_key : str
        The key in adata.obsm for reduced dimension coordinates (e.g., UMAP)
    start_cluster : str, optional
        Starting cluster for trajectory inference
    end_clusters : list, optional
        End clusters for trajectory inference

    Returns:
    --------
    AnnData
        Updated adata with pseudotime and weight information
    """
    try:
        # Ensure conversion context is properly set
        with localconverter(ro.default_converter + pandas2ri.converter):
            # Import R packages
            base = importr("base")
            utils = importr("utils")

            # Create temporary directory for data exchange
            with tempfile.TemporaryDirectory() as temp_dir:
                print("Creating temporary files...")

                # Export data to CSV files
                if sparse.issparse(adata.X):
                    expr_df = pd.DataFrame(adata.X.toarray())
                else:
                    expr_df = pd.DataFrame(adata.X)

                umap_df = pd.DataFrame(
                    adata.obsm[embedding_key], columns=["UMAP1", "UMAP2"]
                )
                clusters_df = pd.DataFrame({"clusters": adata.obs[cluster_key].astype(str)})

                expr_file = os.path.join(temp_dir, "expr.csv")
                umap_file = os.path.join(temp_dir, "umap.csv")
                clusters_file = os.path.join(temp_dir, "clusters.csv")

                expr_df.to_csv(expr_file, index=False)
                umap_df.to_csv(umap_file, index=False)
                clusters_df.to_csv(clusters_file, index=False)

                print("Reading data and running analysis in R...")

                # Build R command
                r_cmd = f"""
                # Import and load packages
                if (!requireNamespace("BiocManager", quietly = TRUE))
                    install.packages("BiocManager")
                
                required_packages <- c("slingshot", "SingleCellExperiment")
                for (pkg in required_packages) {{
                    if (!requireNamespace(pkg, quietly = TRUE)) {{
                        BiocManager::install(pkg)
                    }}
                }}
                
                library(slingshot)
                library(SingleCellExperiment)
                
                # Load data
                expr_matrix <- as.matrix(read.csv("{expr_file}"))
                umap_coords <- as.matrix(read.csv("{umap_file}"))
                clusters <- read.csv("{clusters_file}")$clusters
                
                # Transpose gene expression matrix (gene x cell)
                expr_matrix <- t(expr_matrix)

                # Create SingleCellExperiment object
                sce <- SingleCellExperiment(
                    assays = list(counts = expr_matrix)
                )
                
                # Add UMAP and cluster info
                reducedDims(sce) <- list(UMAP = umap_coords)
                colData(sce)$clusters <- clusters
                
                # Run Slingshot
                """

                # Add start and end cluster parameters if provided
                if start_cluster is not None:
                    r_cmd += f'start_clus <- "{start_cluster}"\n'
                    r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP", start.clus = start_clus)\n'
                else:
                    r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP")\n'

                r_cmd += f"""

                # Get results
                pseudotimes <- slingPseudotime(sce)
                weights <- slingCurveWeights(sce)

                write.csv(pseudotimes, "{temp_dir}/pseudotimes.csv")
                write.csv(weights, "{temp_dir}/weights.csv")
                
                # Return number of trajectories
                n_lineages <- ncol(pseudotimes)
                cat("Found", n_lineages, "trajectories\\n")
                """

                # Execute R command
                ro.r(r_cmd)

                print("Reading results...")

                pseudotimes_file = os.path.join(temp_dir, "pseudotimes.csv")
                weights_file = os.path.join(temp_dir, "weights.csv")

                if os.path.exists(pseudotimes_file):
                    pseudotimes_df = pd.read_csv(pseudotimes_file, index_col=0)
                    weights_df = pd.read_csv(weights_file, index_col=0)

                    # Add results to adata with embedding_key as unique identifier
                    for i, col in enumerate(pseudotimes_df.columns):
                        adata.obs[f"slingshot_pseudotime_{embedding_key}_{i+1}"] = pseudotimes_df.iloc[
                            :, i
                        ].values

                    for i, col in enumerate(weights_df.columns):
                        adata.obs[f"slingshot_weight_{embedding_key}_{i+1}"] = weights_df.iloc[:, i].values

                    print(
                        f"Slingshot analysis completed! Found {len(pseudotimes_df.columns)} trajectories"
                    )
                    return adata
                else:
                    print("Could not find result files")
                    return None

    except ImportError:
        print("Error: Please install rpy2 package")
        print("Run: pip install rpy2")
        print("And install Slingshot in R: BiocManager::install('slingshot')")
        return None
    except Exception as e:
        print(f"Error running Slingshot: {e}")
        return None


def analyze_trajectory_cluster_transitions(
    adata, cluster_key="leiden_subset", embedding_key="X_umap_subset"
):
    """
    Analyze cluster transition patterns in Slingshot trajectories.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix with Slingshot results
    cluster_key : str
        The key in adata.obs for cluster labels
    embedding_key : str
        The key in adata.obsm for reduced dimension coordinates

    Returns:
    --------
    dict
        Dictionary containing trajectory analysis results
    """
    # Get pseudotime columns with embedding key
    pseudotime_cols = [
        col for col in adata.obs.columns if col.startswith(f"slingshot_pseudotime_{embedding_key}_")
    ]

    if not pseudotime_cols:
        print("No Slingshot pseudotime results found")
        return {}

    print(
        f"Analyzing cluster transition patterns for {len(pseudotime_cols)} trajectories"
    )
    print("=" * 60)

    # Analyze cluster transitions for each trajectory
    trajectory_analysis = {}

    for i, pt_col in enumerate(pseudotime_cols):
        lineage_num = i + 1
        print(f"Trajectory {lineage_num} ({pt_col}):")
        print("-" * 40)

        # Get cells on this trajectory
        valid_mask = ~np.isnan(adata.obs[pt_col])
        if valid_mask.sum() == 0:
            print("No valid cells")
            continue

        # Extract data for this trajectory
        trajectory_data = adata.obs[valid_mask].copy()
        trajectory_data = trajectory_data.sort_values(pt_col)

        # Cluster distribution
        cluster_counts = trajectory_data[cluster_key].value_counts().sort_index()
        print(f"  Involved clusters: {', '.join(cluster_counts.index.astype(str))}")
        print(f"  Cell count distribution: {dict(cluster_counts)}")

        # Calculate average pseudotime for each cluster
        cluster_pseudotime = trajectory_data.groupby(cluster_key)[pt_col].agg(
            ["mean", "std", "min", "max", "count"]
        )

        # Exclude clusters with NaN average pseudotime
        valid_clusters_mask = ~np.isnan(cluster_pseudotime["mean"])
        cluster_pseudotime_filtered = cluster_pseudotime[valid_clusters_mask]

        if len(cluster_pseudotime_filtered) == 0:
            print("  Warning: No valid cluster data found for this trajectory")
            continue

        # Sort by average pseudotime
        cluster_pseudotime_filtered = cluster_pseudotime_filtered.sort_values("mean")

        print(f"  Valid clusters sorted by pseudotime:")
        for cluster_id, row in cluster_pseudotime_filtered.iterrows():
            std_str = f"{row['std']:.2f}" if not np.isnan(row["std"]) else "nan"
            print(
                f"Cluster {cluster_id}: Average pseudotime {row['mean']:.2f} ± {std_str} "
                f"(Range: {row['min']:.2f}-{row['max']:.2f}, Cell count: {row['count']})"
            )

        # Check for excluded clusters
        excluded_clusters = cluster_pseudotime[~valid_clusters_mask]
        if len(excluded_clusters) > 0:
            print(
                f"Excluded clusters (NaN average pseudotime): {', '.join([str(c) for c in excluded_clusters.index])}"
            )

        # Infer transition order (only include valid clusters)
        ordered_clusters = cluster_pseudotime_filtered.index.tolist()
        if len(ordered_clusters) > 1:
            transitions = " → ".join([str(c) for c in ordered_clusters])
            print(f"Inferred trajectory: {transitions}")
        elif len(ordered_clusters) == 1:
            print(f"Single cluster trajectory: {ordered_clusters[0]}")

        # Save analysis results
        trajectory_analysis[f"lineage_{lineage_num}"] = {
            "clusters_involved": ordered_clusters,
            "transition_path": ordered_clusters,
            "cluster_stats": cluster_pseudotime_filtered,
            "excluded_clusters": (
                excluded_clusters.index.tolist() if len(excluded_clusters) > 0 else []
            ),
            "total_cells": valid_mask.sum(),
            "valid_clusters_count": len(ordered_clusters),
        }

    return trajectory_analysis


def analyze_trajectory_relationships(trajectory_analysis):
    """
    Analyze relationships between trajectories to find subsets and branching points.

    Parameters:
    -----------
    trajectory_analysis : dict
        Dictionary containing trajectory analysis results

    Returns:
    --------
    list
        List of relationships between trajectories
    """
    print("=" * 60)

    relationships = []
    all_lineages = list(trajectory_analysis.keys())

    # Check each pair of trajectories for subset or branching relationships
    for i, lineage1 in enumerate(all_lineages):
        for j, lineage2 in enumerate(all_lineages):
            if i >= j:
                continue

            path1 = trajectory_analysis[lineage1]["clusters_involved"]
            path2 = trajectory_analysis[lineage2]["clusters_involved"]

            # Check for subset relationship
            if is_subpath(path1, path2):
                relationships.append(
                    {
                        "type": "subset",
                        "shorter": lineage1,
                        "longer": lineage2,
                        "shorter_path": path1,
                        "longer_path": path2,
                        "divergence_point": len(path1),
                        "extension": path2[len(path1) :],
                    }
                )
            elif is_subpath(path2, path1):
                relationships.append(
                    {
                        "type": "subset",
                        "shorter": lineage2,
                        "longer": lineage1,
                        "shorter_path": path2,
                        "longer_path": path1,
                        "divergence_point": len(path2),
                        "extension": path1[len(path2) :],
                    }
                )
            else:
                # Check for common prefix
                common_prefix = find_common_prefix(path1, path2)
                if len(common_prefix) > 1:  # At least 2 common steps
                    relationships.append(
                        {
                            "type": "branching",
                            "lineage1": lineage1,
                            "lineage2": lineage2,
                            "path1": path1,
                            "path2": path2,
                            "common_prefix": common_prefix,
                            "branch1": path1[len(common_prefix) :],
                            "branch2": path2[len(common_prefix) :],
                            "divergence_point": len(common_prefix),
                        }
                    )

    # Show results
    if not relationships:
        print("No relationships found between trajectories.")
        return relationships

    subset_relations = [r for r in relationships if r["type"] == "subset"]
    branching_relations = [r for r in relationships if r["type"] == "branching"]

    # Show subset relationships
    if subset_relations:
        print("\n📦 Found subset relationships:")
        for rel in subset_relations:
            shorter_num = rel["shorter"].split("_")[1]
            longer_num = rel["longer"].split("_")[1]
            print(f"Trajectory {shorter_num} ⊆ Trajectory {longer_num}")
            print(f"Shorter path: {' → '.join(map(str, rel['shorter_path']))}")
            print(f"Longer path: {' → '.join(map(str, rel['longer_path']))}")
            print(
                f"Divergence point: Step {rel['divergence_point']} (Cluster {rel['shorter_path'][-1]})"
            )
            print(f"Extension: {' → '.join(map(str, rel['extension']))}")

            # Analyze cell counts
            shorter_cells = trajectory_analysis[rel["shorter"]]["total_cells"]
            longer_cells = trajectory_analysis[rel["longer"]]["total_cells"]
            print(
                f"    Cell counts: Shorter ({shorter_cells}) vs Longer ({longer_cells})"
            )

    # Show branching relationships
    if branching_relations:
        print("Found branching relationships:")
        for rel in branching_relations:
            lineage1_num = rel["lineage1"].split("_")[1]
            lineage2_num = rel["lineage2"].split("_")[1]
            print(f"Trajectory {lineage1_num} ↔ Trajectory {lineage2_num}")
            print(f"Common prefix: {' → '.join(map(str, rel['common_prefix']))}")
            print(f"Branch 1: {' → '.join(map(str, rel['branch1']))}")
            print(f"Branch 2: {' → '.join(map(str, rel['branch2']))}")
            print(
                f"Divergence point: Step {rel['divergence_point']} (Cluster {rel['common_prefix'][-1]})"
            )

    return relationships


def is_subpath(shorter_path, longer_path):
    """Check if shorter_path is a prefix of longer_path."""
    if len(shorter_path) > len(longer_path):
        return False
    return shorter_path == longer_path[: len(shorter_path)]


def find_common_prefix(path1, path2):
    """Find the common prefix between two paths."""
    common = []
    for i in range(min(len(path1), len(path2))):
        if path1[i] == path2[i]:
            common.append(path1[i])
        else:
            break
    return common


def suggest_trajectory_merging(trajectory_analysis, relationships):
    """
    Based on subset relationships, suggest trajectory merging.

    Parameters:
    -----------
    trajectory_analysis : dict
        Dictionary containing trajectory analysis results
    relationships : list
        List of relationships between trajectories
    """
    print("=" * 50)

    subset_relations = [r for r in relationships if r["type"] == "subset"]

    if not subset_relations:
        print("No subset relationships found, no merging suggestions.")
        return

    for rel in subset_relations:
        shorter_num = rel["shorter"].split("_")[1]
        longer_num = rel["longer"].split("_")[1]

        print(f"Suggest merging Trajectory {shorter_num} and Trajectory {longer_num}:")
        print(
            f"Reason: Trajectory {shorter_num} is a complete subset of Trajectory {longer_num}"
        )

        # Analyze biological significance based on cell counts
        shorter_cells = trajectory_analysis[rel["shorter"]]["total_cells"]
        longer_cells = trajectory_analysis[rel["longer"]]["total_cells"]

        if shorter_cells > longer_cells * 0.8:
            print(f"Biological explanation: Possible early differentiation stalling")
            print(
                f"- Most cells remain at cluster {rel['shorter_path'][-1]} (differentiation point)"
            )
            print(
                f"- Some cells continue to differentiate into {' → '.join(map(str, rel['extension']))}"
            )
        else:
            print(f"Biological explanation: Possible branching differentiation")
            print(
                f"- Main differentiation path: {' → '.join(map(str, rel['longer_path']))}"
            )
            print(f"- Some cells terminate early at cluster {rel['shorter_path'][-1]}")

        print(f"Unified path after merging: {' → '.join(map(str, rel['longer_path']))}")
        print(
            f"Key differentiation point: Cluster {rel['shorter_path'][-1]} → Cluster {rel['extension'][0] if rel['extension'] else 'N/A'}"
        )


def merge_subset_trajectories(
    adata, trajectory_analysis, relationships, merge_strategy="keep_longer"
):
    """
    Merge trajectories that have subset relationships.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix
    trajectory_analysis : dict
        Dictionary containing trajectory analysis results
    relationships : list
        List of relationships between trajectories
    merge_strategy : str
        Strategy for merging ('keep_longer', 'keep_shorter', 'combine')

    Returns:
    --------
    dict
        Updated trajectory analysis after merging
    """
    print("Executing trajectory merging based on subset relationships")
    print("=" * 50)

    subset_relations = [r for r in relationships if r["type"] == "subset"]

    if not subset_relations:
        print("No subset relationships found, no merging needed.")
        return trajectory_analysis

    merged_analysis = trajectory_analysis.copy()

    for rel in subset_relations:
        shorter_key = rel["shorter"]
        longer_key = rel["longer"]
        shorter_num = shorter_key.split("_")[1]
        longer_num = longer_key.split("_")[1]

        print(f"Merged trajectory {shorter_num} into {longer_num}")

        if merge_strategy == "keep_longer":
            # Keep the longer trajectory and remove the shorter one
            print(
                f"Strategy: Keep complete trajectory {longer_num}, remove subset trajectory {shorter_num}"
            )
            del merged_analysis[shorter_key]

            # Update the longer trajectory with merger info
            merged_analysis[longer_key]["merger_info"] = {
                "merged_from": shorter_key,
                "divergence_point": rel["divergence_point"],
                "divergence_cluster": rel["shorter_path"][-1],
                "extension_path": rel["extension"],
                "cells_before_divergence": trajectory_analysis[shorter_key][
                    "total_cells"
                ],
                "cells_after_divergence": trajectory_analysis[longer_key]["total_cells"]
                - trajectory_analysis[shorter_key]["total_cells"],
            }

        elif merge_strategy == "keep_shorter":
            # Keep the shorter trajectory as the main path and mark the longer trajectory as an extension
            print(
                f"Strategy: Keep main trajectory {shorter_num}, mark trajectory {longer_num} as extension"
            )
            merged_analysis[shorter_key]["extension_info"] = {
                "extended_in": longer_key,
                "extension_path": rel["extension"],
                "total_extended_cells": trajectory_analysis[longer_key]["total_cells"],
            }
            del merged_analysis[longer_key]

        elif merge_strategy == "combine":
            # Combine the paths into a new trajectory
            combined_key = f"combined_{shorter_num}_{longer_num}"
            merged_analysis[combined_key] = {
                "clusters_involved": rel["longer_path"],
                "transition_path": rel["longer_path"],
                "main_branch": rel["shorter_path"],
                "extension_branch": rel["extension"],
                "divergence_point": rel["divergence_point"],
                "divergence_cluster": rel["shorter_path"][-1],
                "main_branch_cells": trajectory_analysis[shorter_key]["total_cells"],
                "extension_cells": trajectory_analysis[longer_key]["total_cells"]
                - trajectory_analysis[shorter_key]["total_cells"],
                "total_cells": trajectory_analysis[longer_key]["total_cells"],
                "original_trajectories": [shorter_key, longer_key],
            }
            # Delete the original trajectories
            del merged_analysis[shorter_key]
            del merged_analysis[longer_key]

    return merged_analysis


def analyze_gene_expression_along_trajectories(
    adata, gene_names, trajectory_analysis=None, use_merged=True, embedding_key="X_umap"
):
    """
    Analyze the expression of specified genes along trajectories based on pseudotime.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix with gene expression and pseudotime data
    gene_names : str or list
        Gene name(s) to analyze
    trajectory_analysis : dict, optional
        Dictionary containing trajectory analysis results
    use_merged : bool
        Whether to use merged trajectory analysis if available
    embedding_key : str
        The embedding key used for the trajectory analysis

    Returns:
    --------
    dict
        Dictionary containing gene expression analysis results for each trajectory
    """
    # Input gene_names can be a single gene or a list of genes
    if isinstance(gene_names, str):
        gene_names = [gene_names]

    # Use provided trajectory analysis or look for existing ones
    if trajectory_analysis is None:
        print("No trajectory analysis provided. Please run trajectory analysis first.")
        return {}

    print(f"Gene Analysis: {', '.join(gene_names)}")
    print(f"Number of Trajectories: {len(trajectory_analysis)}")
    print("=" * 60)

    # Check for available genes
    available_genes = []
    missing_genes = []

    for gene in gene_names:
        if gene in adata.var_names:
            available_genes.append(gene)
        else:
            missing_genes.append(gene)

    if missing_genes:
        print(f"Genes not found: {', '.join(missing_genes)}")

        # Try to find similar genes in highly variable genes
        for missing_gene in missing_genes:
            # Search in highly variable genes
            if hasattr(adata.var, "highly_variable"):
                hvg_genes = adata.var_names[adata.var.highly_variable]
                matches = [g for g in hvg_genes if missing_gene.lower() in g.lower()]
                if matches:
                    print(f"'{missing_gene}' possible matches: {matches[:5]}")

    if not available_genes:
        print("No available genes for analysis")
        return {}

    print(f"Available genes for analysis: {', '.join(available_genes)}")

    # Get all pseudotime columns with embedding key
    pseudotime_cols = [
        col for col in adata.obs.columns if col.startswith(f"slingshot_pseudotime_{embedding_key}_")
    ]

    # Create analysis for each gene
    gene_results = {}

    for gene in available_genes:
        print(f"Analyzing: {gene}")
        print("-" * 50)

        # Get gene expression data
        gene_idx = adata.var_names.get_loc(gene)
        if hasattr(adata.X, "toarray"):
            gene_expression = adata.X[:, gene_idx].toarray().flatten()
        else:
            gene_expression = adata.X[:, gene_idx]

        # Analyze each trajectory
        trajectory_data = {}

        for traj_key, traj_info in trajectory_analysis.items():
            if "clusters_involved" not in traj_info:
                continue

            traj_num = traj_key.split("_")[-1] if "_" in traj_key else traj_key

            # Find the corresponding pseudotime column with embedding key
            pt_col = None
            for col in pseudotime_cols:
                if col.endswith(f"_{traj_num}"):
                    pt_col = col
                    break

            if pt_col is None:
                print(f"No pseudotime data found for trajectory {traj_num}")
                continue

            # Get valid cells for this trajectory
            valid_mask = ~np.isnan(adata.obs[pt_col])
            if valid_mask.sum() == 0:
                continue

            # Extract pseudotime and gene expression for valid cells
            pseudotime = adata.obs[pt_col][valid_mask].values
            expression = gene_expression[valid_mask]

            # Calculate Spearman correlation
            correlation, p_value = stats.spearmanr(pseudotime, expression)

            trajectory_data[traj_key] = {
                "pseudotime": pseudotime,
                "expression": expression,
                "correlation": correlation,
                "p_value": p_value,
                "n_cells": len(pseudotime),
                "traj_name": f"Trajectory{traj_num}",
            }

            print(
                f"{trajectory_data[traj_key]['traj_name']}: "
                f"Correlation={correlation:.3f}, p={p_value:.3e}, Cell Count={len(pseudotime)}"
            )

        gene_results[gene] = trajectory_data

    return gene_results


def show_highly_variable_genes(adata, n_genes=20):
    """
    Show available highly variable genes in the dataset.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix
    n_genes : int
        Number of genes to display
    """
    print(f"Available genes: {n_genes}")
    print("=" * 50)

    if hasattr(adata.var, "highly_variable"):
        hvg_genes = adata.var_names[adata.var.highly_variable]
        print(f"Total {len(hvg_genes)} highly variable genes found.")
        print(f"Top {n_genes} genes:")
        for i, gene in enumerate(hvg_genes[:n_genes]):
            print(f"  {i+1:2d}. {gene}")

        if len(hvg_genes) > n_genes:
            print(f"... and {len(hvg_genes) - n_genes} more genes")
    else:
        print("No highly variable genes found in the dataset.")
        print(f"Total gene count: {adata.n_vars}")
        print(f"Top {n_genes} genes:")
        for i, gene in enumerate(adata.var_names[:n_genes]):
            print(f"  {i+1:2d}. {gene}")
