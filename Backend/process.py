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
            nuclei_df = pd.read_csv(SAMPLES[sample_id]["nuclei_path"], index_col=0)
            cell_nuclei_merged_df = nuclei_df[["id", "cell_x", "cell_y"]]

            cell_coordinate_result[sample_id] = cell_nuclei_merged_df.to_dict(
                orient="records"
            )

    return cell_coordinate_result


def get_gene_list(sample_ids):
    """
    Get the list of genes for the given sample IDs.
    """
    sample_gene_dict = {}

    for sample_id in sample_ids:
        gene_df = pd.read_csv(SAMPLES[sample_id]["gene_list_path"])
        gene_list = gene_df.columns.tolist()
        sample_gene_dict[sample_id] = gene_list

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


def get_umap_data(sample_id, n_neighbors=15, min_dist=0.1, n_components=2, n_clusters=5, random_state=42):
    """
    Generate UMAP data from gene expression data.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    # Load gene expression data
    gene_expression_path = SAMPLES[sample_id]["gene_expression_path"]
    df = pd.read_parquet(gene_expression_path)
    
    # Separate ID column from gene expression data
    cell_ids = df['id'].values
    gene_data = df.drop('id', axis=1).values
    
    # Filter out cells/genes with all zeros to improve UMAP performance
    # Remove cells with no expression
    cell_mask = gene_data.sum(axis=1) > 0
    gene_data_filtered = gene_data[cell_mask]
    cell_ids_filtered = cell_ids[cell_mask]
    
    # Remove genes with no expression across all cells
    gene_mask = gene_data_filtered.sum(axis=0) > 0
    gene_data_filtered = gene_data_filtered[:, gene_mask]
    
    if gene_data_filtered.shape[0] == 0:
        return []
    
    # Apply log transformation for better UMAP results
    # Add small constant to avoid log(0)
    gene_data_log = np.log1p(gene_data_filtered)
    
    # Add small random noise to help with spectral initialization
    noise_factor = 1e-6
    gene_data_log += np.random.normal(0, noise_factor, gene_data_log.shape)
    
    # Adjust n_neighbors if we have fewer cells than requested
    actual_n_neighbors = min(n_neighbors, gene_data_filtered.shape[0] - 1)
    
    # Perform UMAP dimensionality reduction
    reducer = umap.UMAP(
        n_neighbors=actual_n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=random_state,
        init='random'  # Use random initialization to avoid spectral warnings
    )
    
    embedding = reducer.fit_transform(gene_data_log)
    
    # Adjust n_clusters if we have fewer cells than requested clusters
    actual_n_clusters = min(n_clusters, gene_data_filtered.shape[0])
    
    # Perform clustering on the UMAP embedding
    kmeans = KMeans(n_clusters=actual_n_clusters, random_state=random_state, n_init=10)
    cluster_labels = kmeans.fit_predict(embedding)
    
    # Format results according to UMAP component structure
    results = []
    for i in range(len(cell_ids_filtered)):
        results.append({
            'id': cell_ids_filtered[i],
            'x': float(embedding[i, 0]),
            'y': float(embedding[i, 1]),
            'cluster': f'Cluster {cluster_labels[i] + 1}'
        })
    
    return results


def perform_go_analysis(sample_id, cell_ids, top_n=5):
    """
    Perform GO analysis on the selected cluster cells and return top GO terms.
    """
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample {sample_id} not found")
    
    try:
        gene_expression_path = SAMPLES[sample_id]["gene_expression_path"]
        df = pd.read_parquet(gene_expression_path)

        cluster_df = df[df['id'].isin(cell_ids)]
        
        if cluster_df.empty:
            return []
        
        # Calculate mean expression for each gene in the cluster
        gene_cols = [col for col in cluster_df.columns if col != 'id']
        mean_expression = cluster_df[gene_cols].mean()
        
        # Filter for significantly expressed genes (above threshold)
        expression_threshold = mean_expression.quantile(0.7)  # Top 30% of genes
        significant_genes = mean_expression[mean_expression > expression_threshold].index.tolist()
        
        if len(significant_genes) < 5:
            # If too few genes, use top 20 expressed genes
            significant_genes = mean_expression.nlargest(20).index.tolist()
        
        # Perform GO enrichment analysis using gseapy
        try:
            enr = gp.enrichr(
                gene_list=significant_genes,
                gene_sets=['GO_Biological_Process_2023'],
                organism='human',
                description='cluster_go_analysis',
                cutoff=0.05,
                no_plot=True
            )

            results = enr.results
            
            if results.empty:
                # If no significant results, try with less stringent criteria
                enr = gp.enrichr(
                    gene_list=significant_genes,
                    gene_sets=['GO_Biological_Process_2023'],
                    organism='human',
                    description='cluster_go_analysis',
                    cutoff=0.5,
                    no_plot=True
                )
                results = enr.results
            
            # Format and return top results
            go_results = []
            for i, (_, row) in enumerate(results.head(top_n).iterrows()):
                go_results.append({
                    'rank': i + 1,
                    'term': row['Term'],
                    'description': row['Term'].split('(')[0].strip(),
                    'p_value': float(row['P-value']),
                    'adjusted_p_value': float(row['Adjusted P-value']),
                    'odds_ratio': float(row['Odds Ratio']) if 'Odds Ratio' in row else 0,
                    'combined_score': float(row['Combined Score']) if 'Combined Score' in row else 0,
                    'genes': row['Genes'].split(';') if 'Genes' in row else []
                })
            
            return go_results
            
        except Exception as e:
            print(f"GO analysis error: {e}")
            return []
    
    except Exception as e:
        print(f"Error in GO analysis: {e}")
        return []