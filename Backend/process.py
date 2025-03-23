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
    return [{"id": sample["id"], "name": sample["name"]} for sample in SAMPLES.values()]


# return tissue width and height size
def get_hires_image_size(sample_id):
    Image.MAX_IMAGE_PIXELS = None
    if sample_id == "skin_TXK6Z4X_A1":
        image = Image.open(SAMPLES["skin_TXK6Z4X_A1"]["wsi"])
    elif sample_id == "skin_TXK6Z4X_D1":
        image = Image.open(SAMPLES["skin_TXK6Z4X_D1"]["wsi"])

    return image.size


# return unique cell types
def get_unique_cell_types(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        adata = sc.read_h5ad(SAMPLES["skin_TXK6Z4X_A1"]["adata"])
    elif sample_id == "skin_TXK6Z4X_D1":
        adata = sc.read_h5ad(SAMPLES["skin_TXK6Z4X_D1"]["adata"])

    return adata.obs["cell_type"].unique().tolist()


# return cell type, and cell coordinates
def get_cell_type_coordinates(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        adata = sc.read_h5ad(SAMPLES["skin_TXK6Z4X_A1"]["adata"])
    elif sample_id == "skin_TXK6Z4X_D1":
        adata = sc.read_h5ad(SAMPLES["skin_TXK6Z4X_D1"]["adata"])

    df = adata.obsm["spatial"].copy()
    df["cell_type"] = adata.obs["cell_type"]
    df["id"] = adata.obs.index
    return df


# return gene list by sending a sample array
def get_gene_list(sample_names):
    grouped_options = []
    sample_gene_sets = {}

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

        sample_gene_sets[sample_name] = set(adata.var_names)

    if not sample_gene_sets:
        return grouped_options
    
    if len(sample_gene_sets) == 1:
        # Only one sample, return its gene list without common group
        sample_name = next(iter(sample_gene_sets))
        unique_gene_children = [
            {"value": f"{sample_name}-{gene}", "title": gene}
            for gene in sorted(sample_gene_sets[sample_name])
        ]
        grouped_options.append({"value": sample_name, "title": sample_name, "children": unique_gene_children})
    else:
        # Multiple samples, compute common and unique genes
        common_genes = set.intersection(*sample_gene_sets.values())
        unique_genes_per_sample = {
            sample: genes - common_genes for sample, genes in sample_gene_sets.items()
        }

        common_gene_children = [
            {"value": f"common-{gene}", "title": gene} for gene in sorted(common_genes)
        ]
        grouped_options.append(
            {
                "value": "common_genes",
                "title": "Common Genes",
                "children": common_gene_children,
            }
        )

        for sample_name, unique_genes in unique_genes_per_sample.items():
            unique_gene_children = [
                {"value": f"{sample_name}-{gene}", "title": gene}
                for gene in sorted(unique_genes)
            ]
            grouped_options.append(
                {
                    "value": sample_name,
                    "title": sample_name,
                    "children": unique_gene_children,
                }
            )
    
    return grouped_options


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

    def calculate_radius(df, originaldf, radius_val):
        originaldf = originaldf.copy()
        originaldf["radius"] = radius_val
        result_df = originaldf[["id", "cell_x", "cell_y", "radius"] + gene_list].copy()
        for index, row in df.iterrows():
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

                merged_df["total_expression"] = merged_df[gene_names].sum(axis=1)

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
        cumsum_df = filter_and_cumsum(merged_df, gene_list)

        kosara_df = calculate_radius(cumsum_df, merged_df, radius)

        results[sample_id] = kosara_df.to_dict(orient="records")

    return results
