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
def get_unique_columns(sample_names):
    columns_list = []

    for sample_name in sample_names:
        sample_info = SAMPLES.get(sample_name)

        h5ad_path = sample_info.get("adata")

        try:
            adata = sc.read_h5ad(h5ad_path)
        except Exception as e:
            print(f"Failed to read {h5ad_path}: {str(e)}")
            continue

        sample_columns = set(adata.var_names)
        columns_list.append(sample_columns)

    if not columns_list:
        return []

    if len(columns_list) == 1:
        unique_columns = columns_list[0]
    else:
        all_union = set().union(*columns_list)
        all_intersection = set(columns_list[0]).intersection(*columns_list[1:])
        unique_columns = all_union - all_intersection

    options = [
        {"value": col, "label": col} for i, col in enumerate(sorted(unique_columns))
    ]

    return options


def filter_and_merge(cell_ids, gene_names, sample_id):
    adata_path = SAMPLES[sample_id]["adata"]
    adata = sc.read_h5ad(adata_path)

    valid_cell_ids = [cell for cell in cell_ids if cell in adata.obs_names]
    valid_gene_names = [gene for gene in gene_names if gene in adata.var_names]

    if len(valid_cell_ids) == 0:
        raise ValueError("ID not found in AnnData.")
    if len(valid_gene_names) == 0:
        raise ValueError("Gene not found in AnnData.")

    filtered_adata = adata[valid_cell_ids, valid_gene_names].copy()

    if issparse(filtered_adata.X):
        expr_data = filtered_adata.X.toarray()
    else:
        expr_data = filtered_adata.X
    expr_df = pd.DataFrame(
        expr_data, index=filtered_adata.obs_names, columns=filtered_adata.var_names
    )

    expr_df = expr_df.reset_index().rename(columns={"index": "id"})

    coord_df = get_cell_type_coordinates(sample_id)

    merged_df = pd.merge(expr_df, coord_df, on="id", how="inner")

    merged_df["total_expression"] = merged_df[gene_names].sum(axis=1)

    for gene in gene_names:
        merged_df[f"{gene}_original_ratio"] = np.where(
            merged_df["total_expression"] == 0,
            0,
            merged_df[gene] / merged_df["total_expression"],
        )

    return merged_df


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

# get kosara data
def get_kosara_data(sample_id, gene_list, cell_list):
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

    position_cell_ratios_df = filter_and_merge(cell_list, gene_list, sample_id)
    kosara_df = calculate_radius(
        filter_and_cumsum(position_cell_ratios_df, gene_list),
        position_cell_ratios_df,
        radius,
    )

    return kosara_df