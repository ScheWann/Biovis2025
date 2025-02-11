import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import anndata as ad
import scanpy as sc
from PIL import Image
import tifffile as tifi
import squidpy as sq


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
    return [
        {"id": sample["id"], "name": sample["name"]}
        for sample in SAMPLES.values()
    ]

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
        {"value": col, "label": col}
        for i, col in enumerate(sorted(unique_columns))
    ]

    return options