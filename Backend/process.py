import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import fsolve
import anndata as ad
import scanpy as sc
from PIL import Image
import tifffile as tifi
import squidpy as sq
from scipy.sparse import issparse

hirescalef = 0.10757315

SAMPLES = {
    "skin_TXK6Z4X_A1": {
        "id": "skin_TXK6Z4X_A1",
        "name": "skin_TXK6Z4X_A1",
        "adata": "../Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area|celltypist_cells_adata.h5",
        "wsi": "../Data/skin_TXK6Z4X_A1_processed/tmap/wsi.tif",
        "tiles": "../Data/skin_TXK6Z4X_A1_processed/skin_TXK6Z4X_A1_processed_tiles",
        "cells_layer": "../Data/skin_TXK6Z4X_A1_processed/cells_layer.png",
    },
    "skin_TXK6Z4X_D1": {
        "id": "skin_TXK6Z4X_D1",
        "name": "skin_TXK6Z4X_D1",
        "adata": "../Data/skin_TXK6Z4X_D1_processed/tmap/weighted_by_area|celltypist_cells_adata.h5",
        "wsi": "../Data/skin_TXK6Z4X_D1_processed/tmap/wsi.tif",
        "tiles": "../Data/skin_TXK6Z4X_D1_processed/skin_TXK6Z4X_D1_processed_tiles",
        "cells_layer": "../Data/skin_TXK6Z4X_D1_processed/cells_layer.png",
    },
}


# return sample list
def get_samples():
    return [{"value": sample["id"], "label": sample["name"]} for sample in SAMPLES.values()]


# return tissue width and height size
def get_hires_image_size(sample_ids):
    Image.MAX_IMAGE_PIXELS = None
    sizes = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            image = Image.open(SAMPLES[sample_id]["wsi"])
        sizes[sample_id] = image.size

    return sizes


# return unique cell types
def get_unique_cell_types(sample_ids):
    result = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            adata = sc.read_h5ad(SAMPLES[sample_id]["adata"])
            result[sample_id] = adata.obs["cell_type"].unique().tolist()

    return result


# return cell type, and cell coordinates
def get_cell_type_coordinates(sample_ids):
    result = {}

    for sample_id in sample_ids:
        if sample_id in SAMPLES:
            adata = sc.read_h5ad(SAMPLES[sample_id]["adata"])
            df = adata.obsm["spatial"].copy()
            df["cell_type"] = adata.obs["cell_type"]
            df["id"] = adata.obs.index
            result[sample_id] = df.to_dict(orient="records")

    return result


# return gene list
def get_gene_list(sample_names):
    sample_gene_dict = {}

    for sample_name in sample_names:
        sample_info = SAMPLES.get(sample_name)
        if not sample_info:
            continue

        h5ad_path = sample_info.get("adata")

        try:
            adata = sc.read_h5ad(h5ad_path)
        except Exception as e:
            print(f"Failed to read {h5ad_path}: {str(e)}")
            continue

        gene_sums = adata.X.sum(axis=0)
        gene_names = adata.var_names

        sample_gene_dict[sample_name] = {
            gene: float(gene_sums[i]) for i, gene in enumerate(gene_names)
        }

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
                a_value = row[col]
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
                    valid_cell_ids = [cell for cell in cell_ids if cell in adata.obs_names]

                valid_gene_names = [gene for gene in gene_names if gene in adata.var_names]

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

                missing_genes = [gene for gene in gene_names if gene not in expr_df.columns]
                for gene in missing_genes:
                    expr_df[gene] = 0

                expr_df = expr_df[['id'] + gene_names]

                coord_df = get_cell_type_coordinates(sample_id).reset_index(drop=True)

                merged_df = pd.merge(expr_df, coord_df, on="id", how="inner")

                merged_df["total_expression"] = adata.X.sum(axis=1)

                for gene in gene_names:
                    merged_df[f"{gene}_original_ratio"] = np.where(
                        merged_df["total_expression"] == 0,
                        0,
                        merged_df[gene] / merged_df["total_expression"],
                    )

                results[sample_id] = merged_df

        return results

    def filter_and_cumsum(merged_df, selected_columns):
        def conditional_cumsum(row, selected_columns):
            cumsum = 0
            for col in selected_columns:
                if row[col] != 0:
                    cumsum += row[col]
                    row[col] = cumsum
            return row

        filtered_df = merged_df[["id", "cell_x", "cell_y"] + selected_columns].copy()
        filtered_df[selected_columns] = filtered_df[selected_columns].apply(
            conditional_cumsum, axis=1, selected_columns=selected_columns
        )
        return filtered_df

    position_cell_ratios_dict = filter_and_merge(cell_list, gene_list, sample_ids)
    results = {}

    for sample_id, merged_df in position_cell_ratios_dict.items():
        # cumsum_df = filter_and_cumsum(merged_df, gene_list)

        # kosara_df = calculate_radius(cumsum_df, merged_df, radius)
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
                "ratios": {}
            }
            
            for gene in gene_list:
                transformed_entry["angles"][gene] = row.get(f"{gene}_angle", 0)
                transformed_entry["radius"][gene] = row.get(f"{gene}_radius", 0)
                transformed_entry["ratios"][gene] = row.get(f"{gene}_original_ratio", 0)
            
            formatted_results.append(transformed_entry)

        results[sample_id] = formatted_results

    return results
