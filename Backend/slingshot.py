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
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

SKIN_GENE_SIGNATURES = {
    # Epidermis-related genes
    "Epidermis": [
        "KRT5",
        "KRT14",
        "KRT15",
        "KRT19",  # Basal layer keratin
        "TP63",
        "TP73",  # Epidermal stem cell transcription factor
        "ITGA6",
        "ITGB4",  # Epidermal stem cell marker
        "COL17A1",
        "LAMA5",  # Basement membrane component
    ],
    # Dermis-related genes
    "Dermis": [
        "COL1A1",
        "COL1A2",
        "COL3A1",  # Collagen
        "VIM",
        "ACTA2",  # Fibroblast marker
        "PDGFRA",
        "PDGFRB",  # Fibroblast receptor
        "ELN",
        "FBN1",  # Elastic fiber
    ],
    # Hair follicle-related genes
    "Hair_Follicle": [
        "LGR5",
        "LGALS7",  # Hair follicle stem cells
        "SOX9",
        "NFATC1",  # Hair follicle development
        "BMP2",
        "BMP4",  # Hair follicle morphogenesis
        "WNT3",
        "WNT10A",  # Wnt signaling pathway
    ],
    # Sebaceous gland-related genes
    "Sebaceous": [
        "BLIMP1",
        "CEBPA",  # Sebaceous gland differentiation
        "PPARG",
        "SREBF1",  # Lipid metabolism
        "FASN",
        "SCD",  # Fatty acid synthesis
    ],
    # Immune/inflammation-related genes
    "Immune": [
        "CD68",
        "CD163",  # Macrophages
        "CD3E",
        "CD8A",  # T cells
        "IL1B",
        "TNF",  # Inflammatory factors
        "PTPRC",  # Common leukocyte antigen
    ],
    # Vascular-related genes
    "Vascular": [
        "PECAM1",
        "CDH5",  # Endothelial cells
        "ACTA2",
        "TAGLN",  # Smooth muscle cells
        "VEGFA",
        "ANGPT1",  # Angiogenesis
    ],
}


def intelligent_slingshot_analysis(
    adata,
    gene_signatures=None,
    cluster_key="leiden",
    embedding_key="X_umap",
    expected_trajectories=None,
    min_cluster_size=10,
    plot_analysis=True,
    auto_infer_trajectories=True,
):
    """
    基于基因表达的智能Slingshot轨迹分析

    Parameters:
    -----------
    adata : AnnData
        single-cell data
    gene_signatures : dict
        gene signature dictionary, if None, use the default SKIN_GENE_SIGNATURES
    cluster_key : str
        clustering information key
    embedding_key : str
        embedding coordinates key
    expected_trajectories : list, optional
        expected trajectory directions, e.g. [('Hair_Follicle', 'Epidermis'), ('Epidermis', 'Dermis')]
        If None and auto_infer_trajectories=True, will infer automatically
    min_cluster_size : int
        minimum cluster size threshold
    plot_analysis : bool
        whether to plot analysis results
    auto_infer_trajectories : bool
        whether to automatically infer trajectory directions (based on gene expression)

    Returns:
    --------
    dict : containing analysis results and adjusted adata
    """

    if gene_signatures is None:
        gene_signatures = SKIN_GENE_SIGNATURES

    print("=====Starting Intelligent Slingshot Analysis=====")
    print("=" * 80)

    # 1. Analyzing gene expression patterns
    print("1. Analyzing gene expression patterns...")
    gene_analysis = analyze_cluster_gene_signatures(
        adata, gene_signatures, cluster_key, min_cluster_size
    )

    if not gene_analysis["valid_clusters"]:
        print("No enough valid clusters found.")
        return None

    # 2. Determining the best trajectory starting points
    print("2. Determining the best trajectory starting points...")

    # If auto-inference is needed and no expected trajectories are provided
    if auto_infer_trajectories and expected_trajectories is None:
        print(
            " Automatically inferring trajectory directions based on gene expression..."
        )
        inferred_trajectories = auto_infer_developmental_trajectories(
            adata, gene_analysis, cluster_key, embedding_key
        )
        expected_trajectories = inferred_trajectories["trajectories"]

        print(f"There are {len(expected_trajectories)} potential trajectories:")
        for i, (start, end) in enumerate(expected_trajectories):
            print(f"Trajectory {i+1}: {start} → {end}")

    trajectory_plan = determine_trajectory_starting_points(
        adata, gene_analysis, expected_trajectories, cluster_key
    )

    # 3. Running optimized Slingshot analysis
    print("3. Running optimized Slingshot analysis")
    optimized_results = {}

    for i, start_info in enumerate(trajectory_plan["starting_points"]):
        start_cluster = start_info["cluster_id"]
        rationale = start_info["rationale"]

        print(f"Trajectory {i+1}: Starting Cluster {start_cluster}")
        print(f"Rationale: {rationale}")

        # Run Slingshot
        adata_result = run_slingshot_via_rpy2_improved(
            adata.copy(),
            cluster_key=cluster_key,
            embedding_key=embedding_key,
            start_cluster=str(start_cluster),
        )

        if adata_result is not None:
            optimized_results[f"start_cluster_{start_cluster}"] = adata_result

    # If no optimized results found, run default analysis (No specified starting points)
    if not optimized_results:
        print("Running default Slingshot analysis...(No specified starting points)")
        default_result = run_slingshot_via_rpy2_improved(
            adata.copy(), cluster_key=cluster_key, embedding_key=embedding_key
        )
        if default_result is not None:
            optimized_results["default"] = default_result

    # 4. Evaluating and selecting the best trajectory result
    print("4. Evaluating and selecting the best trajectory result")
    best_result = select_best_trajectory_result(
        optimized_results, gene_signatures, cluster_key
    )

    if best_result is None:
        print("No valid trajectory result found.")
        return None

    # 5. Validating and visualizing
    print("5. Validating trajectory biological significance")
    validation_results = validate_trajectory_biology(
        best_result["adata"], gene_signatures, cluster_key
    )

    # Plotting the comprehensive analysis results
    if plot_analysis:
        plot_intelligent_analysis_summary(
            best_result["adata"],
            gene_analysis,
            validation_results,
            cluster_key,
            embedding_key,
        )

    return {
        "final_adata": best_result["adata"],
        "gene_analysis": gene_analysis,
        "trajectory_plan": trajectory_plan,
        "all_results": optimized_results,
        "best_result_info": best_result,
        "validation": validation_results,
    }


def analyze_cluster_gene_signatures(
    adata, gene_signatures, cluster_key, min_cluster_size
):
    """Analyzing cluster gene signatures"""

    print("Calculating cluster gene signature scores...")

    cluster_gene_scores = {}
    cluster_info = {}

    for cluster_id in adata.obs[cluster_key].unique():
        if pd.isna(cluster_id):
            continue

        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_size = cluster_mask.sum()

        if cluster_size < min_cluster_size:
            print(f"Skipping cluster {cluster_id} due to insufficient cell count.")
            continue

        cluster_cells = adata[cluster_mask]
        cluster_scores = {}

        # Calculating gene signature scores
        for cell_type, genes in gene_signatures.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if len(available_genes) == 0:
                cluster_scores[cell_type] = 0
                continue

            if hasattr(cluster_cells.X, "toarray"):
                gene_expr = cluster_cells[:, available_genes].X.toarray()
            else:
                gene_expr = cluster_cells[:, available_genes].X

            # Using mean expression as signature score
            signature_score = np.mean(gene_expr)
            cluster_scores[cell_type] = signature_score

        cluster_gene_scores[cluster_id] = cluster_scores
        cluster_info[cluster_id] = {
            "size": cluster_size,
            "dominant_signature": max(cluster_scores, key=cluster_scores.get),
            "max_score": max(cluster_scores.values()),
            "scores": cluster_scores,
        }

        # Output cluster features
        dominant = cluster_info[cluster_id]["dominant_signature"]
        max_score = cluster_info[cluster_id]["max_score"]
        print(
            f"Cluster {cluster_id:2}: {cluster_size:3} cells | Dominant Feature: {dominant:12} (Score: {max_score:.3f})"
        )

    return {
        "cluster_scores": cluster_gene_scores,
        "cluster_info": cluster_info,
        "valid_clusters": list(cluster_info.keys()),
    }


def auto_infer_developmental_trajectories(
    adata, gene_analysis, cluster_key, embedding_key
):
    """
    Automatically inferring developmental trajectories based on gene expression and spatial proximity
    """

    cluster_info = gene_analysis["cluster_info"]
    cluster_scores = gene_analysis["cluster_scores"]

    print("Analyzing cluster cell type primitiveness and differentiation potential...")

    # 1. Calculating cell type primitiveness scores
    primitiveness_scores = calculate_cluster_primitiveness(cluster_info)

    # 2. Calculating transcriptional similarity between clusters
    transcriptional_distances = calculate_cluster_transcriptional_distances(
        adata, cluster_key, list(cluster_info.keys())
    )

    # 3. Calculating spatial adjacency (based on UMAP)
    spatial_adjacency = calculate_cluster_spatial_adjacency(
        adata, cluster_key, embedding_key, list(cluster_info.keys())
    )

    # 4. Inferring potential trajectories
    potential_trajectories = infer_trajectories_from_analysis(
        cluster_info, primitiveness_scores, transcriptional_distances, spatial_adjacency
    )

    return {
        "trajectories": potential_trajectories,
        "primitiveness_scores": primitiveness_scores,
        "transcriptional_distances": transcriptional_distances,
        "spatial_adjacency": spatial_adjacency,
        "analysis_details": {
            "method": "gene_expression_based",
            "n_trajectories": len(potential_trajectories),
        },
    }


def calculate_cluster_primitiveness(cluster_info):
    """Calculating the primitiveness scores for each cluster"""

    # Defining the developmental hierarchy of cell types
    developmental_hierarchy = {
        "Hair_Follicle": 1.0,  # Hair follicle stem cells - most primitive
        "Sebaceous": 0.8,  # Sebaceous gland progenitors
        "Epidermis": 0.6,  # Epidermal cells
        "Dermis": 0.4,  # Dermal fibroblasts
        "Vascular": 0.2,  # Vascular cells
        "Immune": 0.1,  # Immune cells - most differentiated
    }

    primitiveness_scores = {}

    for cluster_id, info in cluster_info.items():
        dominant_type = info["dominant_signature"]
        base_score = developmental_hierarchy.get(dominant_type, 0.5)

        signature_strength = info["max_score"]

        # Original Score: Primitiveness Score = Base Developmental Level * Gene Signature Strength
        primitiveness_score = base_score * min(signature_strength, 1.0)

        primitiveness_scores[cluster_id] = {
            "score": primitiveness_score,
            "cell_type": dominant_type,
            "signature_strength": signature_strength,
        }

        print(
            f"Cluster {cluster_id}: {dominant_type:12} | Primitiveness: {primitiveness_score:.3f}"
        )

    return primitiveness_scores


def calculate_cluster_transcriptional_distances(adata, cluster_key, valid_clusters):
    """
    Calculating the transcriptional distances between clusters
    """

    print("   📊 Calculating transcriptional similarity between clusters...")

    # Calculating the average gene expression for each cluster
    cluster_profiles = {}

    for cluster_id in valid_clusters:
        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_cells = adata[cluster_mask]

        if hasattr(cluster_cells.X, "toarray"):
            cluster_expr = np.mean(cluster_cells.X.toarray(), axis=0)
        else:
            cluster_expr = np.mean(cluster_cells.X, axis=0)

        cluster_profiles[cluster_id] = cluster_expr

    # Calculating the distance matrix
    cluster_ids = list(cluster_profiles.keys())
    n_clusters = len(cluster_ids)
    distance_matrix = np.zeros((n_clusters, n_clusters))

    for i, cluster1 in enumerate(cluster_ids):
        for j, cluster2 in enumerate(cluster_ids):
            if i != j:
                # Using Pearson correlation coefficient distance
                correlation = np.corrcoef(
                    cluster_profiles[cluster1], cluster_profiles[cluster2]
                )[0, 1]
                distance = 1 - abs(correlation)  # Distance = 1 - |Correlation|
                distance_matrix[i, j] = distance

    # Converting to dictionary format
    distance_dict = {}
    for i, cluster1 in enumerate(cluster_ids):
        distance_dict[cluster1] = {}
        for j, cluster2 in enumerate(cluster_ids):
            distance_dict[cluster1][cluster2] = distance_matrix[i, j]

    return distance_dict


def calculate_cluster_spatial_adjacency(
    adata, cluster_key, embedding_key, valid_clusters
):
    """Calculating the spatial adjacency of clusters"""

    print("   🗺️  Calculating spatial adjacency...")
    if embedding_key not in adata.obsm:
        print(
            "      Warning: Embedding coordinates not found, skipping spatial analysis"
        )
        return {}

    # Calculating the center point for each cluster
    cluster_centers = {}
    embedding = adata.obsm[embedding_key]

    for cluster_id in valid_clusters:
        cluster_mask = adata.obs[cluster_key] == cluster_id
        cluster_coords = embedding[cluster_mask]
        center = np.mean(cluster_coords, axis=0)
        cluster_centers[cluster_id] = center

    # Calculate the distance between cluster centers
    adjacency_dict = {}
    cluster_ids = list(cluster_centers.keys())

    for cluster1 in cluster_ids:
        adjacency_dict[cluster1] = {}
        center1 = cluster_centers[cluster1]

        for cluster2 in cluster_ids:
            if cluster1 != cluster2:
                center2 = cluster_centers[cluster2]
                distance = np.linalg.norm(center1 - center2)
                # Convert distance to adjacency (smaller distance means higher adjacency)
                adjacency = 1.0 / (1.0 + distance)
                adjacency_dict[cluster1][cluster2] = adjacency
            else:
                adjacency_dict[cluster1][cluster2] = 1.0

    return adjacency_dict


def infer_trajectories_from_analysis(
    cluster_info, primitiveness_scores, transcriptional_distances, spatial_adjacency
):
    """
    Based on multi-omics analysis to infer developmental trajectories
    """

    print("Integrating analysis results to infer trajectories...")

    trajectories = []
    cluster_ids = list(cluster_info.keys())

    # Sort clusters by primitiveness
    sorted_clusters = sorted(
        cluster_ids, key=lambda x: primitiveness_scores[x]["score"], reverse=True
    )

    # Strategy 1: Start from the most primitive clusters and look for potential differentiation paths
    for i, start_cluster in enumerate(
        sorted_clusters[:3]
    ):  # Consider the top 3 most primitive clusters
        start_type = cluster_info[start_cluster]["dominant_signature"]
        start_primitiveness = primitiveness_scores[start_cluster]["score"]

        # Find potential end clusters
        potential_ends = []

        for end_cluster in cluster_ids:
            if end_cluster == start_cluster:
                continue

            end_type = cluster_info[end_cluster]["dominant_signature"]
            end_primitiveness = primitiveness_scores[end_cluster]["score"]

            # Condition 1: The end cluster should be more differentiated than the start cluster (lower primitiveness)
            if end_primitiveness >= start_primitiveness:
                continue

            # Condition 2: Check biological plausibility
            if not is_biologically_plausible_transition(start_type, end_type):
                continue

            # Condition 3: Calculate overall connection strength
            transcriptional_similarity = 1 - transcriptional_distances.get(
                start_cluster, {}
            ).get(end_cluster, 1.0)
            spatial_proximity = spatial_adjacency.get(start_cluster, {}).get(
                end_cluster, 0.0
            )
            primitiveness_gradient = start_primitiveness - end_primitiveness

            connection_score = (
                transcriptional_similarity * 0.4
                + spatial_proximity * 0.3
                + primitiveness_gradient * 0.3
            )

            potential_ends.append(
                {"cluster": end_cluster, "type": end_type, "score": connection_score}
            )

        # Select the best end cluster
        if potential_ends:
            potential_ends.sort(key=lambda x: x["score"], reverse=True)
            best_end = potential_ends[0]

            if best_end["score"] > 0.2:
                trajectories.append((start_type, best_end["type"]))
                print(
                    f"      Found trajectory: {start_type} → {best_end['type']} (Score: {best_end['score']:.3f})"
                )

    # Remove duplicates and limit the number
    trajectories = list(set(trajectories))[:3]  # At most 3 trajectories

    if not trajectories:
        print("      ⚠️  No clear trajectories found, using default strategy")
        # Backup strategy: Select the most primitive and most differentiated cell types
        if len(sorted_clusters) >= 2:
            start_type = cluster_info[sorted_clusters[0]]["dominant_signature"]
            end_type = cluster_info[sorted_clusters[-1]]["dominant_signature"]
            trajectories = [(start_type, end_type)]

    return trajectories


def is_biologically_plausible_transition(start_type, end_type):
    """Check if the transition between two cell types is biologically plausible"""

    # Define known biological transition relationships
    plausible_transitions = {
        "Hair_Follicle": ["Epidermis", "Sebaceous"],
        "Sebaceous": ["Epidermis"],
        "Epidermis": ["Dermis"],
        "Dermis": ["Vascular"],
        "Vascular": [],
        "Immune": [],  # Immune cells are usually not precursors to other cell types
    }

    # Check if the transition is in the known list
    allowed_targets = plausible_transitions.get(start_type, [])
    return (
        end_type in allowed_targets or len(allowed_targets) == 0
    )  # If not defined, allow transition


def determine_trajectory_starting_points(
    adata, gene_analysis, expected_trajectories, cluster_key
):
    """Determine trajectory starting points"""

    cluster_info = gene_analysis["cluster_info"]

    # Define the primitiveness hierarchy of cell types (lower values are more primitive)
    primitiveness_hierarchy = {
        "Hair_Follicle": 1,  # Hair follicle stem cells are the most primitive
        "Sebaceous": 2,  # Sebaceous stem cells
        "Epidermis": 3,  # Epidermal cells
        "Dermis": 4,  # Dermal fibroblasts
        "Vascular": 5,  # Vascular cells
        "Immune": 6,  # Immune cells (usually not developmental starting points)
    }

    starting_points = []

    # Strategy 1: If there are expected trajectories, use the expected starting points
    if expected_trajectories:
        print("Determine starting points based on expected trajectories:")
        for start_type, end_type in expected_trajectories:
            # Find the cluster that best matches the starting type
            best_cluster = None
            best_score = -1

            for cluster_id, info in cluster_info.items():
                if info["dominant_signature"] == start_type:
                    score = info["scores"][start_type]
                    if score > best_score:
                        best_score = score
                        best_cluster = cluster_id

            if best_cluster is not None:
                starting_points.append(
                    {
                        "cluster_id": best_cluster,
                        "rationale": f"Expected starting point for {start_type}→{end_type}",
                        "confidence": "high",
                    }
                )
                print(f"{start_type}→{end_type}: Choosing cluster {best_cluster}")

    # Strategy 2: Automatically select starting points based on primitiveness
    if not starting_points:
        print("Automatically selecting starting points based on cell primitiveness:")

        # Sort clusters by primitiveness
        cluster_primitiveness = []
        for cluster_id, info in cluster_info.items():
            dominant_type = info["dominant_signature"]
            primitiveness_score = primitiveness_hierarchy.get(dominant_type, 10)

            cluster_primitiveness.append(
                {
                    "cluster_id": cluster_id,
                    "type": dominant_type,
                    "primitiveness": primitiveness_score,
                    "expression_score": info["max_score"],
                    "size": info["size"],
                }
            )

        # Sort by: primitiveness > expression score > cluster size
        cluster_primitiveness.sort(
            key=lambda x: (x["primitiveness"], -x["expression_score"], -x["size"])
        )

        # Select the most primitive 1-2 clusters as starting points
        for i, cluster_data in enumerate(cluster_primitiveness[:2]):
            starting_points.append(
                {
                    "cluster_id": cluster_data["cluster_id"],
                    "rationale": f"The most primitive cell type ({cluster_data['type']}, primitiveness: {cluster_data['primitiveness']})",
                    "confidence": "medium" if i == 0 else "low",
                }
            )
            print(f"Choosing {cluster_data['cluster_id']}: {cluster_data['type']}")

    # Strategy 3: If all else fails, select the largest cluster
    if not starting_points:
        print("Backup strategy: Select the largest cluster")
        largest_cluster = max(
            cluster_info.keys(), key=lambda x: cluster_info[x]["size"]
        )
        starting_points.append(
            {
                "cluster_id": largest_cluster,
                "rationale": f"max size {cluster_info[largest_cluster]['size']}",
                "confidence": "low",
            }
        )
        print(f"Choose Cluster {largest_cluster}")

    return {
        "starting_points": starting_points,
        "primitiveness_hierarchy": primitiveness_hierarchy,
        "strategy_used": "expected" if expected_trajectories else "automatic",
    }


def select_best_trajectory_result(results_dict, gene_signatures, cluster_key):
    """Choose the best trajectory from multiple results"""

    if not results_dict:
        return None

    print("Evaluating the quality of each trajectory result:")

    result_scores = {}

    for result_name, adata_result in results_dict.items():
        # Get pseudotime columns
        pt_cols = [
            col
            for col in adata_result.obs.columns
            if col.startswith("slingshot_pseudotime")
        ]

        if not pt_cols:
            continue

        total_score = 0
        trajectory_count = len(pt_cols)

        for pt_col in pt_cols:
            # Calculate the quality score for this trajectory
            valid_mask = ~np.isnan(adata_result.obs[pt_col])
            if valid_mask.sum() < 10:
                continue

            # Calculate gene expression correlation score
            correlation_score = calculate_trajectory_gene_correlation_score(
                adata_result, pt_col, gene_signatures
            )

            # Calculate cell count-based score
            coverage_score = valid_mask.sum() / len(adata_result)

            # Combine scores
            trajectory_score = correlation_score * 0.7 + coverage_score * 0.3
            total_score += trajectory_score

        avg_score = total_score / trajectory_count if trajectory_count > 0 else 0
        result_scores[result_name] = {
            "score": avg_score,
            "trajectory_count": trajectory_count,
            "adata": adata_result,
        }

        print(f"{result_name}: Score {avg_score:.3f} ({trajectory_count} trajectories)")

    # Select the best result
    if result_scores:
        best_name = max(result_scores, key=lambda x: result_scores[x]["score"])
        best_result = result_scores[best_name]
        print(f"Best result: {best_name} (Score: {best_result['score']:.3f})")
        return {
            "name": best_name,
            "adata": best_result["adata"],
            "score": best_result["score"],
        }

    return None


def calculate_trajectory_gene_correlation_score(adata, pt_col, gene_signatures):
    """Calculate the correlation score between trajectory and gene expression"""

    valid_mask = ~np.isnan(adata.obs[pt_col])
    if valid_mask.sum() < 10:
        return 0

    valid_adata = adata[valid_mask]
    pseudotime = valid_adata.obs[pt_col].values

    correlation_scores = []

    for cell_type, genes in gene_signatures.items():
        available_genes = [g for g in genes if g in adata.var_names]
        if len(available_genes) == 0:
            continue

        # Calculate gene signature score
        if hasattr(valid_adata.X, "toarray"):
            gene_expr = valid_adata[:, available_genes].X.toarray()
        else:
            gene_expr = valid_adata[:, available_genes].X

        signature_score = np.mean(gene_expr, axis=1)

        # Calculate correlation
        correlation, p_value = spearmanr(pseudotime, signature_score)

        if p_value < 0.05:
            correlation_scores.append(abs(correlation))

    return np.mean(correlation_scores) if correlation_scores else 0


def validate_trajectory_biology(adata, gene_signatures, cluster_key):
    """Verify the biological relevance of trajectories"""

    pt_cols = [
        col for col in adata.obs.columns if col.startswith("slingshot_pseudotime")
    ]

    validation_results = {}

    for pt_col in pt_cols:
        traj_name = pt_col.replace("slingshot_pseudotime_", "Trajectory_")

        valid_mask = ~np.isnan(adata.obs[pt_col])
        if valid_mask.sum() < 10:
            continue

        # Analysis the trajectory direction
        correlations = {}
        for cell_type, genes in gene_signatures.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if len(available_genes) == 0:
                continue

            valid_adata = adata[valid_mask]
            pseudotime = valid_adata.obs[pt_col].values

            if hasattr(valid_adata.X, "toarray"):
                gene_expr = valid_adata[:, available_genes].X.toarray()
            else:
                gene_expr = valid_adata[:, available_genes].X

            signature_score = np.mean(gene_expr, axis=1)
            correlation, p_value = spearmanr(pseudotime, signature_score)

            correlations[cell_type] = {
                "correlation": correlation,
                "p_value": p_value,
                "significant": p_value < 0.05,
            }

        # Infer trajectory direction
        start_features = [
            ct
            for ct, data in correlations.items()
            if data["correlation"] < -0.3 and data["significant"]
        ]
        end_features = [
            ct
            for ct, data in correlations.items()
            if data["correlation"] > 0.3 and data["significant"]
        ]

        validation_results[traj_name] = {
            "correlations": correlations,
            "inferred_start": start_features,
            "inferred_end": end_features,
            "direction": (
                f"{start_features} → {end_features}"
                if start_features and end_features
                else "N/A"
            ),
            "valid_cells": valid_mask.sum(),
        }

        print(
            f"{traj_name}: {validation_results[traj_name]['direction']} ({valid_mask.sum()} cells)"
        )

    return validation_results


# Improved Slingshot Function
def run_slingshot_via_rpy2_improved(
    adata,
    cluster_key="leiden",
    embedding_key="X_umap",
    start_cluster=None,
    end_clusters=None,
):
    try:
        from rpy2 import robjects as ro
        from rpy2.robjects.packages import importr
        import tempfile
        from scipy import sparse

        # Import R
        base = importr("base")
        utils = importr("utils")

        # Creating temporary
        with tempfile.TemporaryDirectory() as temp_dir:
            print("Creating temporary files...")

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

            # Start and End cluster parameters (Optional)
            if start_cluster is not None:
                r_cmd += f'start_clus <- "{start_cluster}"\n'
                r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP", start.clus = start_clus)\n'
            else:
                r_cmd += 'sce <- slingshot(sce, clusterLabels = "clusters", reducedDim = "UMAP")\n'

            r_cmd += f"""

            # Get result
            pseudotimes <- slingPseudotime(sce)
            weights <- slingCurveWeights(sce)

            write.csv(pseudotimes, "{temp_dir}/pseudotimes.csv")
            write.csv(weights, "{temp_dir}/weights.csv")
            
            # Return number of trajectories
            n_lineages <- ncol(pseudotimes)
            cat("Found", n_lineages, "trajectories\\n")
            """

            # Execute R
            ro.r(r_cmd)

            print("Reading results...")

            pseudotimes_file = os.path.join(temp_dir, "pseudotimes.csv")
            weights_file = os.path.join(temp_dir, "weights.csv")

            if os.path.exists(pseudotimes_file):
                pseudotimes_df = pd.read_csv(pseudotimes_file, index_col=0)
                weights_df = pd.read_csv(weights_file, index_col=0)

                # Add to adata
                for i, col in enumerate(pseudotimes_df.columns):
                    adata.obs[f"slingshot_pseudotime_{i+1}"] = pseudotimes_df.iloc[
                        :, i
                    ].values

                for i, col in enumerate(weights_df.columns):
                    adata.obs[f"slingshot_weight_{i+1}"] = weights_df.iloc[:, i].values

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


def smart_slingshot_skin_analysis(adata, expected_trajectories=None, auto_infer=True, **kwargs):
    # If automatic inference is enabled and no expected trajectories are provided, set to None for automatic analysis
    if auto_infer and expected_trajectories is None:
        print("Automatic trajectory inference enabled.")
        expected_trajectories = None
    elif not auto_infer and expected_trajectories is None:
        expected_trajectories = [
            ("Hair_Follicle", "Epidermis"),
            ("Epidermis", "Dermis"),
            ("Sebaceous", "Epidermis"),
        ]
        print("Default trajectories set.")

    return intelligent_slingshot_analysis(
        adata,
        gene_signatures=SKIN_GENE_SIGNATURES,
        expected_trajectories=expected_trajectories,
        auto_infer_trajectories=auto_infer,
        **kwargs,
    )


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
