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

skin_TXK6Z4X_A1_adata_path = "../Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area|celltypist_cells_adata.h5"
skin_TXK6Z4X_D1_adata_path = "../Data/skin_TXK6Z4X_D1_processed/tmap/weighted_by_area|celltypist_cells_adata.h5"

skin_TXK6Z4X_A1_wsi_path = "../Data/skin_TXK6Z4X_A1_processed/tmap/wsi.tif"
skin_TXK6Z4X_D1_wsi_path = "../Data/skin_TXK6Z4X_D1_processed/tmap/wsi.tif"

skin_TXK6Z4X_A1_wsi_tiles_dir = "../Data/skin_TXK6Z4X_A1_processed/skin_TXK6Z4X_A1_processed_tiles"
skin_TXK6Z4X_D1_wsi_tiles_dir = "../Data/skin_TXK6Z4X_D1_processed/skin_TXK6Z4X_D1_processed_tiles"

skin_TXK6Z4X_A1_cells_layer_image_path = "../Data/skin_TXK6Z4X_A1_processed/cells_layer.png"
skin_TXK6Z4X_D1_cells_layer_image_path = "../Data/skin_TXK6Z4X_D1_processed/cells_layer.png"


# return tissue width and height size
def get_hires_image_size(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        image = Image.open(skin_TXK6Z4X_A1_wsi_path)
    elif sample_id == "skin_TXK6Z4X_D1":
        image = Image.open(skin_TXK6Z4X_D1_wsi_path)

    return image.size

# return tissue image tiles
def get_tif_tiles(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        tiles = [f"/tiles/{file}" for file in os.listdir(skin_TXK6Z4X_A1_wsi_tiles_dir) if file.endswith('.tif')]
        return tiles
    elif sample_id == "skin_TXK6Z4X_D1":
        tiles = [f"/tiles/{file}" for file in os.listdir(skin_TXK6Z4X_D1_wsi_tiles_dir) if file.endswith('.tif')]
        return tiles


# return unique cell types
def get_unique_cell_types(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        adata = sc.read_h5ad(skin_TXK6Z4X_A1_adata_path)
    elif sample_id == "skin_TXK6Z4X_D1":
        adata = sc.read_h5ad(skin_TXK6Z4X_D1_adata_path)

    return adata.obs["cell_type"].unique().tolist()

# return cell type, and cell coordinates
def get_cell_type_coordinates(sample_id):
    if sample_id == "skin_TXK6Z4X_A1":
        adata = sc.read_h5ad(skin_TXK6Z4X_A1_adata_path)
    elif sample_id == "skin_TXK6Z4X_D1":
        adata = sc.read_h5ad(skin_TXK6Z4X_D1_adata_path)

    df = adata.obsm["spatial"].copy()
    df["cell_type"] = adata.obs["cell_type"]
    return df