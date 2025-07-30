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
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
import umap

# Disable the PIL image limit entirely
Image.MAX_IMAGE_PIXELS = None


JSON_PATH = "./samples_list.json"
"""
    Load sample list from a JSON file.
"""
with open(JSON_PATH, "r") as f:
    SAMPLES = json.load(f)


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
                    adata = sc.read_h5ad(sample_info["adata_path"])
                    scalef = adata.uns['spatial']['skin_TXK6Z4X_A1']['scalefactors']['tissue_0.5_mpp_150_buffer_scalef']
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
                    adata = sc.read_h5ad(sample_info["adata_path"])
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
                adata_path = SAMPLES[sample_id]["adata"]
                adata = sc.read_h5ad(adata_path)

                if not cell_ids:
                    valid_cell_ids = adata.obs_names.tolist()
                else:
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
            adata = sc.read_h5ad(SAMPLES[sample_id]["adata"])
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
        adata = sc.read_h5ad(sample_info["adata_path"])
    else:
        raise ValueError(f"No gene expression data available for sample {sample_id}")
    
    if cell_ids is not None:
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

    sc.tl.leiden(adata, resolution=resolutions)
    adata.obs[f'leiden_{adata_umap_title}'] = adata.obs['leiden'].copy()

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

    # Return the result
    return results


def perform_go_analysis(sample_id, cluster_id, adata_umap_title, top_n=5):
    """
    Perform GO analysis on the selected cluster cells and return top GO terms.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    sample_info = SAMPLES[sample_id]
    
    if "adata_path" in sample_info:
        adata = sc.read_h5ad(sample_info["adata_path"])
    else:
        raise ValueError(f"No gene expression data available for sample {sample_id}")
    
    sc.tl.rank_genes_groups(adata, groupby=f'leiden_{adata_umap_title}', method='wilcoxon')
    cluster_name = str(cluster_id)
    top_genes = adata.uns['rank_genes_groups']['names'][cluster_name][:100].tolist()

    enr = gp.enrichr(gene_list=top_genes,
                 gene_sets='GO_Biological_Process_2023',
                 organism='Human',
                 cutoff=0.05)
    
    go_results = enr.results.sort_values(by='Combined Score', ascending=False)
    top_df = go_results.head(top_n)

    return top_df.to_dict(orient='records')


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