import numpy as np
import pandas as pd
import json
import io
from collections import defaultdict
from scipy.optimize import fsolve
import scanpy as sc
from PIL import Image
import gseapy as gp
from scipy.sparse import issparse
import networkx as nx

# Disable the PIL image limit entirely
Image.MAX_IMAGE_PIXELS = None

# Global AnnData cache to avoid repeated loading
ADATA_CACHE = {}

JSON_PATH = "./samples_list.json"
"""
    Load sample list from a JSON file.
"""
with open(JSON_PATH, "r") as f:
    SAMPLES = json.load(f)


def load_adata_to_cache(sample_ids):
    """
    Load AnnData objects for the given sample IDs into the global cache.
    This should be called once when samples are confirmed.
    """
    global ADATA_CACHE
    
    for sample_id in sample_ids:
        if sample_id in SAMPLES and sample_id not in ADATA_CACHE:
            sample_info = SAMPLES[sample_id]
            if "adata_path" in sample_info:
                try:
                    adata = sc.read_h5ad(sample_info["adata_path"])
                    ADATA_CACHE[sample_id] = adata
                    print(f"Loaded AnnData for sample {sample_id}")
                except Exception as e:
                    print(f"Error loading AnnData for sample {sample_id}: {e}")
                    ADATA_CACHE[sample_id] = None
            elif "gene_expression_path" in sample_info:
                # Handle legacy format with gene_expression_path
                print(f"Sample {sample_id} uses legacy format (gene_expression_path), skipping cache loading")
                ADATA_CACHE[sample_id] = None
            else:
                print(f"No AnnData path found for sample {sample_id}")
                ADATA_CACHE[sample_id] = None


def get_cached_adata(sample_id):
    """
    Get AnnData object from cache for the given sample ID.
    """
    if sample_id not in ADATA_CACHE:
        raise ValueError(f"AnnData for sample {sample_id} not found in cache. Call load_adata_to_cache() first.")
    
    if ADATA_CACHE[sample_id] is None:
        raise ValueError(f"AnnData for sample {sample_id} failed to load.")
    
    return ADATA_CACHE[sample_id]


def clear_adata_cache():
    """
    Clear the global AnnData cache to free memory.
    """
    global ADATA_CACHE
    ADATA_CACHE.clear()
    print("AnnData cache cleared")


def get_samples_option():
    """
    Return a list of tissue samples for selector, group by example data and upload data.
    """
    groups = defaultdict(list)

    for sample in SAMPLES.values():
        cell_scale_group = sample.get("group", "default")
        groups[cell_scale_group].append(
            {"value": sample["id"], "label": sample["name"]}
        )

    return [
        {"label": group, "title": group, "options": options}
        for group, options in groups.items()
    ]


def get_hires_image_size(sample_ids):
    """
    Get the size of high-resolution images for the given sample IDs.
    """
    tissue_image_size = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            try:
                image = Image.open(SAMPLES[sample_id]["image_tif_path"])
                tissue_image_size[sample_id] = image.size
                image.close()  # Close the image to free memory
            except Exception as e:
                print(f"Error loading image for sample {sample_id}: {e}")
                tissue_image_size[sample_id] = (0, 0)

    return tissue_image_size


def get_coordinates(sample_ids):
    """
    Get cell coordinates for the given sample IDs.
    """
    cell_coordinate_result = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            sample_info = SAMPLES[sample_id]
            
            # Check if we have adata_path (new format)
            if "adata_path" in sample_info:
                try:
                    adata = get_cached_adata(sample_id)
                    scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_0.5_mpp_150_buffer_scalef']
                    # Use array_row and array_col from obs metadata
                    if 'array_row' in adata.obs and 'array_col' in adata.obs:
                        coords_df = pd.DataFrame({
                            'cell_x': adata.obsm['spatial_cropped_150_buffer'][:, 0] * scalef,
                            'cell_y': adata.obsm['spatial_cropped_150_buffer'][:, 1] * scalef
                        }, index=adata.obs_names)
                    else:
                        print(f"Warning: No spatial coordinates found for sample {sample_id}")
                        coords_df = pd.DataFrame(columns=['cell_x', 'cell_y'])
                    
                    # Add ID column and reset index
                    coords_df = coords_df.reset_index().rename(columns={'index': 'id'})
                    cell_coordinate_result[sample_id] = coords_df.to_dict(orient="records")
                    
                except Exception as e:
                    print(f"Error loading adata coordinates for sample {sample_id}: {e}")
                    cell_coordinate_result[sample_id] = []
            else:
                print(f"Warning: No coordinate data available for sample {sample_id}")
                cell_coordinate_result[sample_id] = []

    return cell_coordinate_result


def get_gene_list(sample_ids):
    """
    Get the list of genes for the given sample IDs.
    """
    sample_gene_dict = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            sample_info = SAMPLES[sample_id]
            
            if "adata_path" in sample_info:
                try:
                    adata = get_cached_adata(sample_id)
                    gene_list = adata.var_names.tolist()
                    sample_gene_dict[sample_id] = gene_list
                except Exception as e:
                    print(f"Error loading adata gene list for sample {sample_id}: {e}")
                    sample_gene_dict[sample_id] = []
            else:
                print(f"Warning: No gene list data available for sample {sample_id}")
                sample_gene_dict[sample_id] = []

    return sample_gene_dict


# get kosara data
def get_kosara_data(sample_ids, gene_list, cell_list):
    radius = 5
    d = np.sqrt(2) * radius

    special_r, special_d = np.sqrt(2) * radius, np.sqrt(2) * radius
    special_angle1 = np.degrees(
        np.arccos((special_r**2 + d**2 - radius**2) / (2 * special_r * d))
    )
    special_angle2 = np.degrees(
        np.arccos((radius**2 + d**2 - special_r**2) / (2 * radius * d))
    )
    special_result = (
        special_angle1 * special_r**2
        + special_angle2 * radius**2
        - d * special_r * np.sin(special_angle1)
    )
    special_value = special_result / (np.pi * radius**2)

    def equation(r, d, a, radius_val, special_value):
        angle1 = np.arccos((r**2 + d**2 - radius_val**2) / (2 * r * d))
        angle2 = np.arccos((radius_val**2 + d**2 - r**2) / (2 * radius_val * d))
        if a <= special_value:
            result = angle1 * r**2 + angle2 * radius_val**2 - d * r * np.sin(angle1)
        else:
            result = (
                angle1 * r**2
                + angle2 * radius_val**2
                - r**2 * np.sin(angle1) * np.cos(angle1)
                - radius_val**2 * np.sin(angle2) * np.cos(angle2)
            )
        return result - (a * np.pi * radius_val**2)

    def initial_guess(a, d, special_value):
        if a <= special_value:
            return d + 0.01
        elif a > 0.95:
            return 8
        else:
            return a * 10 - 1.5

    def calculate_radius(originaldf, radius_val):
        originaldf["radius"] = radius_val
        result_df = originaldf.copy()
        for index, row in originaldf.iterrows():
            for col in gene_list:
                a_value = row[col + "_original_ratio"]
                if a_value == 0:
                    cal_radius = 0
                    angle = 0
                else:
                    try:
                        guess = initial_guess(a_value, d, special_value)
                        cal_radius = fsolve(
                            equation,
                            guess,
                            args=(d, a_value, radius_val, special_value),
                        )[0]
                        angle = np.degrees(
                            np.arccos(
                                (cal_radius**2 + d**2 - radius_val**2)
                                / (2 * cal_radius * d)
                            )
                        )
                    except ValueError as e:
                        cal_radius = np.nan
                        angle = np.nan
                result_df.loc[index, f"{col}_radius"] = cal_radius
                result_df.loc[index, f"{col}_angle"] = angle

        return result_df

    def filter_and_merge(cell_ids, gene_names, sample_ids):
        results = {}
        for sample_id in sample_ids:
            if sample_id not in SAMPLES:
                raise ValueError("Sample not found.")
            else:
                adata = get_cached_adata(sample_id)

                if not cell_ids:
                    valid_cell_ids = adata.obs_names.tolist()
                else:
                    # Convert cell_ids to regular Python strings to avoid numpy string issues
                    if hasattr(cell_ids, '__iter__') and not isinstance(cell_ids, str):
                        cell_ids = [str(cell_id) for cell_id in cell_ids]
                    else:
                        cell_ids = [str(cell_ids)]
                    
                    valid_cell_ids = [
                        cell for cell in cell_ids if cell in adata.obs_names
                    ]

                valid_gene_names = [
                    gene for gene in gene_names if gene in adata.var_names
                ]

                if valid_gene_names:
                    filtered_adata = adata[valid_cell_ids, valid_gene_names].copy()
                    if issparse(filtered_adata.X):
                        expr_data = filtered_adata.X.toarray()
                    else:
                        expr_data = filtered_adata.X
                    expr_df = pd.DataFrame(
                        expr_data,
                        index=filtered_adata.obs_names,
                        columns=filtered_adata.var_names,
                    )
                else:
                    expr_df = pd.DataFrame(index=valid_cell_ids)

                expr_df = expr_df.reset_index().rename(columns={"index": "id"})

                missing_genes = [
                    gene for gene in gene_names if gene not in expr_df.columns
                ]
                for gene in missing_genes:
                    expr_df[gene] = 0

                expr_df = expr_df[["id"] + gene_names]

                coord_df = get_coordinates(sample_id).reset_index(drop=True)

                merged_df = pd.merge(expr_df, coord_df, on="id", how="inner")

                merged_df["total_expression"] = adata[valid_cell_ids, :].X.sum(axis=1)

                for gene in gene_names:
                    merged_df[f"{gene}_original_ratio"] = np.where(
                        merged_df["total_expression"] == 0,
                        0,
                        merged_df[gene] / merged_df["total_expression"],
                    )

                results[sample_id] = merged_df

        return results

    def get_coordinates(sample_id):
        if sample_id in SAMPLES:
            adata = get_cached_adata(sample_id)
            df = adata.obsm["spatial"].copy()
            df["cell_type"] = adata.obs["cell_type"]
            df["id"] = adata.obs.index
            return df

    position_cell_ratios_dict = filter_and_merge(cell_list, gene_list, sample_ids)
    results = {}

    for sample_id, merged_df in position_cell_ratios_dict.items():
        kosara_df = calculate_radius(merged_df, radius)
        formatted_results = []
        for _, row in kosara_df.iterrows():
            transformed_entry = {
                "id": row["id"],
                "cell_x": row["cell_x"],
                "cell_y": row["cell_y"],
                "cell_type": row["cell_type"],
                "total_expression": row["total_expression"],
                "angles": {},
                "radius": {},
                "ratios": {},
            }

            for gene in gene_list:
                transformed_entry["angles"][gene] = row.get(f"{gene}_angle", 0)
                transformed_entry["radius"][gene] = row.get(f"{gene}_radius", 0)
                transformed_entry["ratios"][gene] = row.get(f"{gene}_original_ratio", 0)

            formatted_results.append(transformed_entry)

        results[sample_id] = formatted_results

    return results


# get selected region's gene expression data
def get_selected_region_data(sample_id, cell_ids):
    return "not finished yet"


def get_umap_data(sample_id, cell_ids=None, n_neighbors=10, n_pcas=30, resolutions=1, adata_umap_title=None):
    """
    Generate UMAP data from gene expression data.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    sample_info = SAMPLES[sample_id]
    
    # Check if we have adata_path (new format)
    if "adata_path" in sample_info:
        # Get cached AnnData
        adata = get_cached_adata(sample_id)
        
        if cell_ids is not None:
            # Convert cell_ids to regular Python strings to avoid numpy string issues
            if hasattr(cell_ids, '__iter__') and not isinstance(cell_ids, str):
                # Handle lists, numpy arrays, pandas series, etc.
                cell_ids = [str(cell_id) for cell_id in cell_ids]
            else:
                cell_ids = str(cell_ids)
            adata = adata[cell_ids].copy()
        else:
            raise ValueError(f"No gene expression data available for sample {sample_id}")

        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)

        sc.tl.pca(adata, use_highly_variable=True)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcas)
        sc.tl.umap(adata)
        adata.obsm[f'X_umap_{adata_umap_title}'] = adata.obsm['X_umap'].copy()

        sc.tl.leiden(adata, resolution=resolutions, key_added=f'leiden_{adata_umap_title}')

        embedding = adata.obsm[f'X_umap_{adata_umap_title}'] # shape: (n_cells, 2)
        cluster_labels = adata.obs[f'leiden_{adata_umap_title}'].astype(int).values
        cell_ids_filtered = adata.obs_names.tolist()

        results = []
        for i in range(len(cell_ids_filtered)):
            results.append({
                'id': cell_ids_filtered[i],
                'x': float(embedding[i, 0]),
                'y': float(embedding[i, 1]),
                'cluster': f'Cluster {cluster_labels[i] + 1}'
            })

        # Add cluster information back to the original cached AnnData for GO analysis
        # but don't replace the entire object to preserve all cells for future UMAP generations
        original_adata = get_cached_adata(sample_id)
        
        # Add the leiden clustering results to the original cached AnnData
        # Only update cells that were part of this UMAP analysis
        leiden_col = f'leiden_{adata_umap_title}'
        if leiden_col not in original_adata.obs.columns:
            # Initialize the column with NaN for all cells
            original_adata.obs[leiden_col] = pd.NA
        
        # Update cluster assignments for the cells that were processed in this UMAP
        for cell_id, cluster_label in zip(cell_ids_filtered, cluster_labels):
            if cell_id in original_adata.obs_names:
                original_adata.obs.loc[cell_id, leiden_col] = cluster_label
        
        # Return the result
        return results
    else:
        raise ValueError(f"No gene expression data available for sample {sample_id}")


def perform_go_analysis(sample_id, cluster_id, adata_umap_title, top_n=5):
    """
    Perform GO analysis on the selected cluster cells and return top GO terms.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    sample_info = SAMPLES[sample_id]
    
    if "adata_path" in sample_info:
        # Get cached AnnData
        adata = get_cached_adata(sample_id)
        
        # Ensure the leiden column exists and is categorical for scanpy
        leiden_col = f'leiden_{adata_umap_title}'
        if leiden_col not in adata.obs.columns:
            raise ValueError(f"No clustering results found for {adata_umap_title}. Please generate UMAP first.")
        
        # Convert leiden column to categorical if it isn't already
        if not pd.api.types.is_categorical_dtype(adata.obs[leiden_col]):
            adata.obs[leiden_col] = adata.obs[leiden_col].astype('category')
        
        sc.tl.rank_genes_groups(adata, groupby=leiden_col, method='wilcoxon')
        cluster_name = str(cluster_id)
        top_genes = adata.uns['rank_genes_groups']['names'][cluster_name][:100].tolist()

        enr = gp.enrichr(gene_list=top_genes,
                     gene_sets='GO_Biological_Process_2023',
                     organism='Human',
                     cutoff=0.05)
        
        go_results = enr.results.sort_values(by='Combined Score', ascending=False)
        go_results.rename(columns={'Term': 'term', 'Adjusted P-value': 'adjusted_p_value', 'Combined Score': 'combined_score', 'P-value': 'p_value', "Genes": "genes", "Odds Ratio": "odds_ratio"}, inplace=True)
        top_df = go_results.head(top_n)

        return top_df.to_dict(orient='records')
    else:
        raise ValueError(f"No gene expression data available for sample {sample_id}")


def get_trajectory_gene_list(sample_id):
    """
    Get list of available genes from trajectory data for the given sample ID.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    # Check if sample has trajectory data
    if "horizontal_non_random_gene_trajectory_expression_path" not in SAMPLES[sample_id]:
        raise ValueError(f"No trajectory data available for sample {sample_id}")
    
    try:
        trajectory_path = SAMPLES[sample_id]["horizontal_non_random_gene_trajectory_expression_path"]
        df = pd.read_csv(trajectory_path, index_col=0)
        
        # Ensure variables column is string type
        df["variables"] = df["variables"].astype(str)
        
        # Return unique gene names
        gene_list = df["variables"].unique().tolist()
        return sorted(gene_list)  # Sort alphabetically for better UX
        
    except Exception as e:
        print(f"Error loading trajectory gene list for sample {sample_id}: {e}")
        raise ValueError(f"Error loading trajectory gene list: {str(e)}")


def get_trajectory_data(sample_id, selected_genes=None):
    """
    Get trajectory gene expression data for the given sample ID.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    # Check if sample has trajectory data
    if "horizontal_non_random_gene_trajectory_expression_path" not in SAMPLES[sample_id]:
        raise ValueError(f"No trajectory data available for sample {sample_id}")
    
    try:
        trajectory_path = SAMPLES[sample_id]["horizontal_non_random_gene_trajectory_expression_path"]
        df = pd.read_csv(trajectory_path, index_col=0)
        
        columns_to_drop = ["colour", "PANEL", "group", "fill", "linewidth", "linetype", "weight", "alpha", "flipped_aes"]
        df = df.drop(columns=[col for col in columns_to_drop if col in df.columns], errors='ignore')

        df["variables"] = df["variables"].astype(str)

        if selected_genes:
            df = df[df["variables"].isin(selected_genes)]
        
        # Group data by gene (variables)
        trajectory_data = {}
        
        for gene in df["variables"].unique():
            gene_data = df[df["variables"] == gene].copy()
            gene_data = gene_data.sort_values("x")  # Sort by x for proper line plotting
            
            trajectory_data[gene] = {
                "gene": gene,
                "data": gene_data[["x", "y", "ymin", "ymax", "se"]].to_dict(orient="records")
            }
        
        return trajectory_data
        
    except Exception as e:
        print(f"Error loading trajectory data for sample {sample_id}: {e}")
        raise ValueError(f"Error loading trajectory data: {str(e)}")


def get_pseudotime_data(sample_id, cell_ids, adata_umap_title, early_markers=None, n_neighbors=15, n_pcas=30, resolutions=1):
    """
    Get pseudotime data for the given sample ID and cluster ID.
    
    Parameters:
    - sample_id: ID of the sample
    - cell_ids: List of cell IDs to analyze
    - adata_umap_title: Title for the UMAP analysis
    - early_markers: Optional list of early marker genes for root identification
    - n_neighbors: Number of neighbors for neighbor graph construction (default: 15)
    - n_pcas: Number of principal components to use (default: 30)
    - resolutions: Resolution parameter for Leiden clustering (default: 1)
    """
    def identify_paga_roots(adata, early_markers=None, cluster_col='leiden'):
        # Calculate PAGA if not done
        if 'paga' not in adata.uns:
            sc.tl.paga(adata, groups=cluster_col)
        
        # Get connectivity metrics
        conn_matrix = adata.uns['paga']['connectivities'].toarray()
        out_degree = np.sum(conn_matrix, axis=1)
        in_degree = np.sum(conn_matrix, axis=0)
        connectivity_score = out_degree - in_degree
        
        # Biological validation if markers provided
        if early_markers and len(early_markers) > 0:
            bio_scores = []
            for i, cluster in enumerate(adata.obs[cluster_col].cat.categories):
                cluster_cells = adata.obs[cluster_col] == cluster
                if cluster_cells.sum() > 0:
                    # Check if markers exist in the dataset
                    valid_markers = [m for m in early_markers if m in adata.var_names]
                    if valid_markers:
                        marker_expr = adata[cluster_cells, valid_markers].X.mean()
                        bio_scores.append(marker_expr)
                    else:
                        bio_scores.append(0)
                else:
                    bio_scores.append(0)
            bio_scores = np.array(bio_scores)
            
            # Combine scores (normalize first) - avoid division by zero
            if connectivity_score.max() != connectivity_score.min():
                conn_norm = (connectivity_score - connectivity_score.min()) / (connectivity_score.max() - connectivity_score.min())
            else:
                conn_norm = np.zeros_like(connectivity_score)
                
            if bio_scores.max() != bio_scores.min():
                bio_norm = (bio_scores - bio_scores.min()) / (bio_scores.max() - bio_scores.min())
            else:
                bio_norm = np.zeros_like(bio_scores)
                
            combined_score = conn_norm + bio_norm
            root_cluster = np.argmax(combined_score)
        else:
            root_cluster = np.argmax(connectivity_score)
        
        return root_cluster, connectivity_score

    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    sample_info = SAMPLES[sample_id]
    
    if "adata_path" in sample_info:
        adata = get_cached_adata(sample_id)

        if cell_ids is not None:
            # Convert cell_ids to regular Python strings to avoid numpy string issues
            if hasattr(cell_ids, '__iter__') and not isinstance(cell_ids, str):
                # Handle lists, numpy arrays, pandas series, etc.
                cell_ids = [str(cell_id) for cell_id in cell_ids]
            else:
                cell_ids = str(cell_ids)
            adata = adata[cell_ids].copy()
        else:
            raise ValueError(f"No gene expression data available for sample {sample_id}")
            
        # Check if adata is empty after filtering
        if adata.n_obs == 0:
            raise ValueError("No cells remaining after filtering. Please check your cell_ids parameter.")

        # Apply the same preprocessing pipeline as get_umap_data to ensure proper data structure for DPT
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, use_highly_variable=True)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcas)
        sc.tl.umap(adata)
        
        # Store the UMAP for this subset
        adata.obsm[f'X_umap_{adata_umap_title}'] = adata.obsm['X_umap'].copy()
        
        # Perform clustering
        sc.tl.leiden(adata, resolution=resolutions, key_added=f'leiden_{adata_umap_title}')

        leiden_col = f'leiden_{adata_umap_title}'
        if not pd.api.types.is_categorical_dtype(adata.obs[leiden_col]):
            adata.obs[leiden_col] = adata.obs[leiden_col].astype('category')

        root_cluster, _ = identify_paga_roots(adata, early_markers=early_markers, cluster_col=leiden_col)

        # Convert root_cluster to match the category type
        categories = adata.obs[leiden_col].cat.categories
        if len(categories) > 0:
            # Try to match the type of the categories
            if isinstance(categories[0], str):
                root_cluster_key = str(root_cluster)
            else:
                root_cluster_key = int(root_cluster) if hasattr(root_cluster, 'item') else root_cluster
        else:
            root_cluster_key = str(root_cluster)
        
        root_cells = adata.obs[adata.obs[leiden_col] == root_cluster_key].index
        
        # Check if root cluster has any cells
        if len(root_cells) == 0:
            raise ValueError(f"Root cluster {root_cluster} has no cells. Available clusters: {list(adata.obs[leiden_col].cat.categories)}")
        
        # Find root cells in the current adata object
        root_indices = np.flatnonzero(adata.obs_names.isin(root_cells))
        if len(root_indices) == 0:
            # If no root cells found after filtering, use the first cell from any cluster
            adata.uns['iroot'] = 0
        else:
            adata.uns['iroot'] = root_indices[0]
        
        # Now that we have properly preprocessed data with a valid neighbor graph, run DPT
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)

        adata.obs[f'dpt_pseudotime_{adata_umap_title}'] = adata.obs['dpt_pseudotime'].copy()

        # Get pseudotime statistics for each cluster
        cluster_pseudotime = {}
        for cluster in adata.obs[leiden_col].cat.categories:
            cluster_cells = adata.obs[leiden_col] == cluster
            if cluster_cells.sum() > 0:
                cluster_pt = adata.obs.loc[cluster_cells, f'dpt_pseudotime_{adata_umap_title}']
                cluster_pseudotime[int(cluster)] = {
                    'mean_pseudotime': float(cluster_pt.mean()),
                    'median_pseudotime': float(cluster_pt.median()),
                    'std_pseudotime': float(cluster_pt.std()),
                    'min_pseudotime': float(cluster_pt.min()),
                    'max_pseudotime': float(cluster_pt.max()),
                    'n_cells': int(len(cluster_pt))
                }

        paga_connectivity = adata.uns['paga']['connectivities']
        
        # Convert PAGA connectivity to NetworkX graph - use integers for consistency
        G = nx.Graph()
        for i in range(paga_connectivity.shape[0]):
            for j in range(i+1, paga_connectivity.shape[1]):
                conn_value = paga_connectivity[i,j]
                if conn_value > 0.1:  # threshold for connectivity
                    # Avoid division by zero
                    weight = 1.0 / max(conn_value, 1e-6)
                    G.add_edge(i, j, weight=weight)

        # Find all possible paths from root cluster to leaf clusters
        leaf_clusters = [node for node in G.nodes() if G.degree(node) == 1 and node != root_cluster]

        # Store results as list of objects
        trajectory_objects = []

        for leaf in leaf_clusters:
            try:
                path = nx.shortest_path(G, root_cluster, leaf)
                path_pseudotimes = []
                
                for cluster in path:
                    if cluster in cluster_pseudotime:
                        path_pseudotimes.append(cluster_pseudotime[cluster]['mean_pseudotime'])
                    else:
                        path_pseudotimes.append(0.0)

                trajectory_obj = {
                    'path': [int(cluster) for cluster in path],
                    'pseudotimes': [f'{pt:.3f}' for pt in path_pseudotimes]
                }
                trajectory_objects.append(trajectory_obj)
            except nx.NetworkXNoPath:
                print(f"No path found to cluster {leaf}")
            except Exception as e:
                print(f"Error processing path to cluster {leaf}: {e}")

        return trajectory_objects
    else:
        raise ValueError(f"No gene expression data available for sample {sample_id}")