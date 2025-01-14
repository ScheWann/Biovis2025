import pandas as pd
from scipy.io import mmread
from scipy import sparse
import numpy as np
import scanpy as sc
from scipy.optimize import fsolve
from PIL import Image

# um_002_adata = sc.read_10x_h5(f"../Data/square_002um/filtered_feature_bc_matrix.h5")
um_008_adata = sc.read_10x_h5(f"../Data/square_008um/filtered_feature_bc_matrix.h5")
# um_016_adata = sc.read_10x_h5(f"../Data/square_016um/filtered_feature_bc_matrix.h5")

# um_002_adata.var_names_make_unique()
um_008_adata.var_names_make_unique()
# um_016_adata.var_names_make_unique()

um_008_feature_df = um_008_adata.to_df()


def current_kmeans(bin_size, kmeans):
    kmeans_n_df = pd.read_csv(f"../Data/square_{bin_size}um/analysis/clustering/gene_expression_kmeans_{kmeans}_clusters/clusters.csv")
    return kmeans_n_df


# Hires image size
def get_hires_image_size():
    image = Image.open("../Data/spatial/tissue_hires_image.png")
    return image.size


# Positions and clusters for 008um
def get_um_positions_with_clusters(bin_size, kmeans):
    position_df = pd.read_parquet(f"../Data/square_{bin_size}um/spatial/tissue_positions.parquet")

    kmeans_n_df = current_kmeans(bin_size, kmeans)
    scale_df = pd.read_json(f"../Data/square_{bin_size}um/spatial/scalefactors_json.json", typ="series")

    # GET ONLY CELLS IN TISSUE
    position_df_filtered = position_df[position_df["in_tissue"] == 1]

    merged_df = pd.merge(position_df_filtered, kmeans_n_df, left_on="barcode", right_on="Barcode")
    merged_df = merged_df.drop(columns=["in_tissue", "array_row", "array_col", "Barcode"])
    merged_df = merged_df.rename(columns={"pxl_row_in_fullres": "y", "pxl_col_in_fullres": "x", "Cluster": "cluster"})
    merged_df["x"] = merged_df["x"] * scale_df.tissue_hires_scalef
    merged_df["y"] = merged_df["y"] * scale_df.tissue_hires_scalef

    return merged_df


def get_umap_positions(bin_size, kmeans):
    umap_df = pd.read_csv(f"../Data/square_{bin_size}um/analysis/umap/gene_expression_2_components/projection.csv")
    kmeans_n_df = current_kmeans(bin_size, kmeans)
    umap_kmeans_merged_df = pd.merge(umap_df, kmeans_n_df, on="Barcode")
    return umap_kmeans_merged_df


def get_gene_list():
    # if bin_size == '008':
    #     gene_list = um_008_feature_df.columns.tolist()
    #     return gene_list
    gene_list = um_008_feature_df.columns.tolist()
    return gene_list
