"""
Slingshot Trajectory and Gene Analysis Module
"""

import pandas as pd
import numpy as np
import scanpy as sc
import subprocess
from scipy import sparse
import warnings
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import tempfile
import os
from scipy import stats
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")



def check_r_availability():
    """
    Check if R and required packages are available for Slingshot analysis.
    """
    availability = {
        "rpy2": False,
        "r_installed": False,
        "slingshot_package": False,
        "singlecellexperiment_package": False,
        "errors": []
    }
    
    # Check RPy2
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        availability["rpy2"] = True
    except ImportError:
        availability["errors"].append("RPy2 not installed. Install with: pip install rpy2")
        return availability

    with localconverter(ro.default_converter + pandas2ri.converter):
        # Check R installation
        try:
            base = importr("base")
            availability["r_installed"] = True
        except Exception as e:
            availability["errors"].append(f"R not available: {e}")
            return availability
        
        # Check BiocManager
        try:
            biocmanager = importr("BiocManager")
            availability["biocmanager"] = True
        except Exception as e:
            availability["errors"].append(f"BiocManager not available: {e}")
        
        # Check Slingshot package
        try:
            slingshot = importr("slingshot")
            availability["slingshot_package"] = True
        except Exception as e:
            availability["errors"].append(f"Slingshot package not available: {e}")
        
        # Check SingleCellExperiment package
        try:
            sce = importr("SingleCellExperiment")
            availability["singlecellexperiment_package"] = True
        except Exception as e:
            availability["errors"].append(f"SingleCellExperiment package not available: {e}")
    
    return availability


def run_slingshot(
    adata,
    cluster_key="leiden",
    embedding_key="X_umap",
    start_cluster=None,
    end_clusters=None,
):
    """
    Run Slingshot analysis with improved error handling and fallback mechanisms.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with Slingshot results
    cluster_key : str
        The key in adata.obs for cluster labels
    embedding_key : str
        The key in adata.obsm for reduced dimension coordinates
    start_cluster : str
        The starting cluster for trajectory analysis
    end_clusters : list of str
        The end clusters for trajectory analysis

    Returns:
    -------
    AnnData : The annotated data matrix with Slingshot results
    """
    # Check R availability first
    availability = check_r_availability()

    if not availability["rpy2"] or not availability["r_installed"]:
        print("R/RPy2 not available, trying subprocess fallback...")
        return run_slingshot_by_subprocess(adata, cluster_key, embedding_key, start_cluster, end_clusters)
    
    # Try RPy2 first, then fallback to subprocess if needed
    result = run_slingshot_by_rpy2(adata, cluster_key, embedding_key, start_cluster, end_clusters)
    if result is not None:
        return result
    
    print("RPy2 method failed, trying subprocess fallback...")
    return run_slingshot_by_subprocess(adata, cluster_key, embedding_key, start_cluster, end_clusters)


def run_slingshot_by_rpy2(
    adata,
    cluster_key="leiden",
    embedding_key="X_umap",
    start_cluster=None,
    end_clusters=None,
):
    """
    Run Slingshot using RPy2 with proper context management.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with Slingshot results
    cluster_key : str
        The key in adata.obs for cluster labels
    embedding_key : str
        The key in adata.obsm for reduced dimension coordinates
    start_cluster : str
        The starting cluster for trajectory analysis
    end_clusters : list of str
        The end clusters for trajectory analysis

    Returns:
    -------
    AnnData : The annotated data matrix with Slingshot results
    """
    try:
        with localconverter(ro.default_converter + pandas2ri.converter):
            # Import R packages with error handling
            try:
                base = importr("base")
                utils = importr("utils")
                print("Successfully imported R packages.")
            except Exception as import_error:
                print(f"Error importing R packages: {import_error}")
                return None

            # Creating temporary directory
            with tempfile.TemporaryDirectory() as temp_dir:
                print("Creating temporary files...")

                # Export to CSV file with proper error handling
                try:
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
                except Exception as file_error:
                    print(f"Error creating temporary files: {file_error}")
                    return None

                print("Reading data and running analysis in R...")

                # Build R command with better error handling
                r_cmd = f"""
                # Error handling wrapper
                tryCatch({{
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
                    expr_matrix <- as.matrix(read.csv(\"{expr_file}\"))
                    umap_coords <- as.matrix(read.csv(\"{umap_file}\"))
                    clusters <- read.csv(\"{clusters_file}\")$clusters
                    
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

                # Start and End cluster parameters (Optional)
                if start_cluster is not None:
                    r_cmd += f'start_clus <- \"{start_cluster}\"\n'
                    r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP", start.clus = start_clus)\n'
                else:
                    r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP")\n'

                r_cmd += f"""
                    # Get result
                    pseudotimes <- slingPseudotime(sce)
                    lineages <- slingLineages(sce)

                    write.csv(pseudotimes, "{temp_dir}/pseudotimes.csv")

                    # lineage list
                    lineages_df <- data.frame(
                        Lineage = rep(names(lineages), lengths(lineages)),
                        Cluster = unlist(lineages)
                    )
                    write.csv(lineages_df, "{temp_dir}/lineages.csv", row.names=FALSE)
                    write.csv(lineages_df, "./backup/lineages.csv", row.names=FALSE)
                    
                    # Return number of trajectories
                    n_lineages <- ncol(pseudotimes)
                    cat("Found", n_lineages, "trajectories\\n")
                    
                }}, error = function(e) {{
                    cat("R Error:", conditionMessage(e), "\\n")
                }})
                """

                # Execute R with proper context management
                try:
                    ro.r(r_cmd)
                except Exception as r_exec_error:
                    print(f"Error executing R command: {r_exec_error}")
                    return None

                print("Reading results...")

                pseudotimes_file = os.path.join(temp_dir, "pseudotimes.csv")
                lineages_file = os.path.join(temp_dir, "lineages.csv")

                if os.path.exists(pseudotimes_file):
                    try:
                        pseudotimes_df = pd.read_csv(pseudotimes_file, index_col=0)
                        lineages_df = pd.read_csv(lineages_file, index_col=0)

                        # Check if results are empty (indicating R error)
                        if pseudotimes_df.empty or lineages_df.empty:
                            print("R analysis failed - empty results returned")
                            return None

                        # Add to adata
                        for i, col in enumerate(pseudotimes_df.columns):
                            adata.obs[f"slingshot_pseudotime_{embedding_key}_{col}"] = pseudotimes_df.iloc[:, i].values

                        # Store lineage information as metadata instead of obs
                        lineage_paths = lineages_df.groupby("Lineage")["Cluster"].apply(list).to_dict()
                        adata.obs[f"slingshot_lineages_{embedding_key}"] = lineage_paths
                        
                        print(
                            f"Slingshot analysis completed! Found {len(pseudotimes_df.columns)} trajectories"
                        )
                        return adata
                    except Exception as read_error:
                        print(f"Error reading result files: {read_error}")
                        return None
                else:
                    print("Could not find result files")
                    return None

    except ImportError:
        print("Error: Please install rpy2 package")
        print("Run: pip install rpy2")
        print("And install Slingshot in R: BiocManager::install('slingshot')")
        return None
    except Exception as e:
        print(f"Error running Slingshot with RPy2: {e}")
        return None


def run_slingshot_by_subprocess(
    adata,
    cluster_key="leiden",
    embedding_key="X_umap",
    start_cluster=None,
    end_clusters=None,
):
    """
    Run Slingshot using subprocess to call R directly.
    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix with Slingshot results
    cluster_key : str
        The key in adata.obs for cluster labels
    embedding_key : str
        The key in adata.obsm for reduced dimension coordinates
    start_cluster : str
        The starting cluster for trajectory analysis
    end_clusters : list of str
        The end clusters for trajectory analysis
        
    Returns:
    --------
    AnnData : The annotated data matrix with Slingshot results
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            print("Creating temporary files for subprocess R call...")

            # Export to CSV file
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

            # Create R script
            r_script = f"""
                # Install and load required packages
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
            """

            if start_cluster is not None:
                r_script += f'''
                    start_clus <- "{start_cluster}"
                    sce <- slingshot(sce, clusterLabels = "clusters", omega = TRUE, reducedDim = "UMAP", start.clus = start_clus)
                '''
            else:
                r_script += '''
                    sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP")
                '''

            r_script += f'''
                    # Get results
                    pseudotimes <- slingPseudotime(sce)
                    lineages <- slingLineages(sce)

                    # Save results
                    write.csv(pseudotimes, "{temp_dir}/pseudotimes.csv")

                    # lineage list
                    lineages_df <- data.frame(
                        Lineage = rep(names(lineages), lengths(lineages)),
                        Cluster = unlist(lineages)
                    )
                    write.csv(lineages_df, "{temp_dir}/lineages.csv", row.names=FALSE)

                    # Print number of trajectories
                    n_lineages <- ncol(pseudotimes)
                    cat("Found", n_lineages, "trajectories\\n")
                '''

            # Write R script to file
            r_script_file = os.path.join(temp_dir, "slingshot_script.R")
            with open(r_script_file, 'w') as f:
                f.write(r_script)

            print("Executing R script via subprocess...")
            
            # Run R script
            try:
                result = subprocess.run(
                    ['Rscript', r_script_file],
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minute timeout
                )
                
                if result.returncode != 0:
                    print(f"R script failed with return code {result.returncode}")
                    print(f"STDOUT: {result.stdout}")
                    print(f"STDERR: {result.stderr}")
                    return None
                    
                print(f"R script output: {result.stdout}")
                
            except subprocess.TimeoutExpired:
                print("R script timed out after 5 minutes")
                return None
            except FileNotFoundError:
                print("Rscript not found. Please ensure R is installed and Rscript is in PATH")
                return None

            # Read results
            pseudotimes_file = os.path.join(temp_dir, "pseudotimes.csv")
            lineages_file = os.path.join(temp_dir, "lineages.csv")

            if os.path.exists(pseudotimes_file):
                try:
                    pseudotimes_df = pd.read_csv(pseudotimes_file, index_col=0)
                    lineages_df = pd.read_csv(lineages_file, index_col=0)

                    # Check if results are empty
                    if pseudotimes_df.empty or lineages_df.empty:
                        print("R analysis failed - empty results returned")
                        return None

                    # Add to adata
                    for i, col in enumerate(pseudotimes_df.columns):
                        adata.obs[f"slingshot_pseudotime_{embedding_key}_{col}"] = pseudotimes_df.iloc[:, i].values

                    lineage_paths = lineages_df.groupby("Lineage")["Cluster"].apply(list).to_dict()
                    adata.uns[f"slingshot_lineages_{embedding_key}"] = lineage_paths

                    print(
                        f"Slingshot analysis completed via subprocess! Found {len(pseudotimes_df.columns)} trajectories"
                    )
                    return adata
                except Exception as read_error:
                    print(f"Error reading result files: {read_error}")
                    return None
            else:
                print("Could not find result files")
                return None

    except Exception as e:
        print(f"Error running Slingshot with subprocess: {e}")
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

        if merge_strategy == 'keep_longer':
            if shorter_key in merged_analysis:
                del merged_analysis[shorter_key]

            if longer_key in merged_analysis:
                merged_analysis[longer_key]['merger_info'] = {
                    'merged_from': shorter_key,
                    'divergence_point': rel['divergence_point'],
                    'divergence_cluster': rel['shorter_path'][-1],
                    'extension_path': rel['extension'],
                    'cells_before_divergence': trajectory_analysis[shorter_key]['total_cells'],
                    'cells_after_divergence': trajectory_analysis[longer_key]['total_cells'] - trajectory_analysis[shorter_key]['total_cells']
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


def direct_slingshot_analysis(
    adata,
    start_cluster,
    cluster_key="leiden",
    embedding_key="X_umap",
    end_clusters=None,
    **kwargs
):
    """
    Direct Slingshot Analysis with a specified start cluster.

    Parameters:
    -----------
    adata : AnnData
        Single-cell data
    start_cluster : str or int
        The cluster ID to use as the starting point for trajectory analysis
    cluster_key : str
        Clustering information key in adata.obs
    embedding_key : str
        Embedding coordinates key in adata.obsm
    end_clusters : list, optional
        List of cluster IDs to use as end points. If None, Slingshot will
        automatically determine end points.
    **kwargs : dict
        Additional arguments passed to run_slingshot_via_rpy2_improved
        
    Returns:
    --------
    dict : containing analysis results and updated adata
    """
    
    print("=====Direct Slingshot Analysis=====")
    print("=" * 50)
    print(f"Start cluster: {start_cluster}")
    print(f"Cluster key: {cluster_key}")
    print(f"Embedding key: {embedding_key}")

    if end_clusters:
        print(f"End clusters: {end_clusters}")
    else:
        print("End clusters: Auto-determined by Slingshot")
    print("=" * 50)
    
    # Validate start_cluster exists in the data
    unique_clusters = adata.obs[cluster_key].unique()
    unique_clusters = [str(c) for c in unique_clusters if pd.notna(c)]
    
    if str(start_cluster) not in unique_clusters:
        print(f"Error: Start cluster '{start_cluster}' not found in data.")
        print(f"Available clusters: {unique_clusters}")
        return None
    
    # Check cluster size
    cluster_mask = adata.obs[cluster_key] == start_cluster
    cluster_size = cluster_mask.sum()
    print(f"Start cluster size: {cluster_size} cells")
    
    if cluster_size < 5:
        print(f"Warning: Start cluster has only {cluster_size} cells, which may be too small for reliable trajectory analysis.")
    
    # Validate embedding exists
    if embedding_key not in adata.obsm:
        print(f"Error: Embedding '{embedding_key}' not found in adata.obsm")
        print(f"Available embeddings: {list(adata.obsm.keys())}")
        return None
    
    # Run Slingshot analysis directly
    print("Running Slingshot analysis...")
    try:
        result_adata = run_slingshot(
            adata.copy(),
            cluster_key=cluster_key,
            embedding_key=embedding_key,
            start_cluster=str(start_cluster),
            end_clusters=end_clusters,
            **kwargs
        )
        
        if result_adata is None:
            print("Slingshot analysis failed.")
            return None
        
        # Get pseudotime columns
        pt_cols = [col for col in result_adata.obs.columns if col.startswith(f"slingshot_pseudotime_{embedding_key}")]
        # weight_cols = [col for col in result_adata.obs.columns if col.startswith("slingshot_weight")]
        lineages_cols = [col for col in result_adata.obs.columns if col.startswith(f"slingshot_lineages_{embedding_key}")]
        
        print(f"Analysis completed successfully!")
        print(f"Found {len(pt_cols)} trajectories")
        
        # Basic trajectory information
        trajectory_info = {}
        for i, pt_col in enumerate(pt_cols):
            traj_name = f"Trajectory_{i+1}"
            valid_mask = ~np.isnan(result_adata.obs[pt_col])
            valid_cells = valid_mask.sum()
            
            trajectory_info[traj_name] = {
                "pseudotime_column": pt_col,
                # "weight_column": weight_cols[i] if i < len(weight_cols) else None,
                "valid_cells": valid_cells,
                "total_cells": len(result_adata),
                "coverage": valid_cells / len(result_adata)
            }
            
            print(f"  {traj_name}: {valid_cells} cells ({trajectory_info[traj_name]['coverage']:.1%} coverage)")
        
        # Analyze cluster transitions using lineages data
        print("\nAnalyzing cluster transitions...")
        cluster_transitions = {}
        
        # Get lineages data from adata.uns
        lineages_key = f"slingshot_lineages_{embedding_key}"
        if lineages_key in result_adata.obs:
            lineage_paths = result_adata.obs[lineages_key]
            print(f"Found lineages data: {lineage_paths}")
            
            # Match trajectories with lineage paths
            for i, (traj_name, traj_info) in enumerate(trajectory_info.items()):
                pt_col = traj_info["pseudotime_column"]
                valid_mask = ~np.isnan(result_adata.obs[pt_col])
                
                if valid_mask.sum() == 0:
                    continue
                
                # Find corresponding lineage path
                lineage_key = f"Lineage{i+1}"  # Assuming lineages are numbered starting from 1
                if lineage_key in lineage_paths:
                    cluster_path = lineage_paths[lineage_key]
                    
                    # Calculate cluster statistics for validation
                    traj_data = result_adata.obs[valid_mask].copy()
                    traj_data = traj_data.sort_values(pt_col)
                    
                    cluster_stats = traj_data.groupby(cluster_key)[pt_col].agg([
                        "mean", "std", "min", "max", "count"
                    ]).sort_values("mean")
                    
                    cluster_transitions[traj_name] = {
                        "ordered_clusters": cluster_path,
                        "cluster_statistics": cluster_stats.to_dict() if not cluster_stats.empty else {},
                        "transition_path": " → ".join([str(c) for c in cluster_path])
                    }
                    
                    print(f"  {traj_name}: {cluster_transitions[traj_name]['transition_path']}")
                else:
                    print(f"  Warning: No lineage path found for {traj_name}")
        else:
            print(f"  Warning: No lineages data found in adata.uns['{lineages_key}']")
            # Fallback to original method if lineages data is not available
            for traj_name, traj_info in trajectory_info.items():
                pt_col = traj_info["pseudotime_column"]
                valid_mask = ~np.isnan(result_adata.obs[pt_col])
                
                if valid_mask.sum() == 0:
                    continue
                
                # Get trajectory data
                traj_data = result_adata.obs[valid_mask].copy()
                traj_data = traj_data.sort_values(pt_col)
                
                # Calculate cluster statistics along trajectory
                cluster_stats = traj_data.groupby(cluster_key)[pt_col].agg([
                    "mean", "std", "min", "max", "count"
                ]).sort_values("mean")
                
                # Filter out clusters with NaN mean pseudotime
                valid_clusters = cluster_stats[~np.isnan(cluster_stats["mean"])]
                
                if len(valid_clusters) > 0:
                    cluster_transitions[traj_name] = {
                        "ordered_clusters": valid_clusters.index.tolist(),
                        "cluster_statistics": valid_clusters.to_dict(),
                        "transition_path": " → ".join([str(c) for c in valid_clusters.index])
                    }
                    
                    print(f"  {traj_name}: {cluster_transitions[traj_name]['transition_path']}")
        
        return {
            "adata": result_adata,
            "start_cluster": start_cluster,
            "trajectory_info": trajectory_info,
            "cluster_transitions": cluster_transitions,
            "analysis_type": "direct",
            "parameters": {
                "cluster_key": cluster_key,
                "embedding_key": embedding_key,
                "end_clusters": end_clusters,
                **kwargs
            }
        }
        
    except Exception as e:
        print(f"Error during Slingshot analysis: {e}")
        return None