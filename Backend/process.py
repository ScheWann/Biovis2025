import pandas as pd
from scipy.io import mmread
from scipy import sparse
import numpy as np
import scanpy as sc
from scipy.optimize import fsolve
from PIL import Image


position_df = pd.read_parquet("../Data/square_008um/spatial/tissue_positions.parquet")

um_008_adata = sc.read_10x_h5("../Data/square_008um/filtered_feature_bc_matrix.h5")


um_008_adata.var_names_make_unique()
um_008_matrix = um_008_adata.to_df()


# Hires image size
def get_hires_image_size():
    image = Image.open("../Data/spatial/tissue_hires_image.png")
    return image.size


# Positions and clusters for 008um
def get_um_008_positions_with_clusters(kmeans):
    kmeans_n_df = pd.read_csv(f"../Data/square_008um/analysis/clustering/gene_expression_kmeans_{kmeans}_clusters/clusters.csv")
    # GET ONLY CELLS IN TISSUE
    position_df_filtered = position_df[position_df["in_tissue"] == 1]

    merged_df = pd.merge(position_df_filtered, kmeans_n_df, left_on="barcode", right_on="Barcode")
    merged_df = merged_df.drop(columns=["in_tissue", "array_row", "array_col", "Barcode"])
    merged_df = merged_df.rename(columns={"pxl_row_in_fullres": "y", "pxl_col_in_fullres": "x", "Cluster": "cluster"})

    return merged_df