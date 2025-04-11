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
import gseapy as gp
from scipy.sparse import issparse
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage, cophenet
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist

hirescalef = 0.10757315

SAMPLES = {
    "skin_TXK6Z4X_A1": {
        "id": "skin_TXK6Z4X_A1",
        "name": "skin_TXK6Z4X_A1",
        "adata": "../Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5",
        "wsi": "../Data/skin_TXK6Z4X_A1_processed/tmap/wsi.tif",
        "tiles": "../Data/skin_TXK6Z4X_A1_processed/skin_TXK6Z4X_A1_processed_tiles",
        "cells_layer": "../Data/skin_TXK6Z4X_A1_processed/cells_layer.png",
    },
    "skin_TXK6Z4X_D1": {
        "id": "skin_TXK6Z4X_D1",
        "name": "skin_TXK6Z4X_D1",
        "adata": "../Data/skin_TXK6Z4X_D1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5",
        "wsi": "../Data/skin_TXK6Z4X_D1_processed/tmap/wsi.tif",
        "tiles": "../Data/skin_TXK6Z4X_D1_processed/skin_TXK6Z4X_D1_processed_tiles",
        "cells_layer": "../Data/skin_TXK6Z4X_D1_processed/cells_layer.png",
    },
}


# return sample list
def get_samples():
    return [
        {"value": sample["id"], "label": sample["name"]} for sample in SAMPLES.values()
    ]


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


# return cell type
def get_cell_types(sample_id):
    sample_info = SAMPLES.get(sample_id)
    if not sample_info:
        return []

    adata = sc.read_h5ad(sample_info["adata"])
    cell_types = adata.obs["cell_type"].unique().tolist()

    return [{"value": ct, "label": ct} for ct in cell_types]


# return gene list
def get_gene_list_for_cell2cellinteraction(sample_id):
    sample_info_list = []

    sample_info = SAMPLES.get(sample_id)
    if not sample_info:
        return sample_info_list

    h5ad_path = sample_info.get("adata")

    try:
        adata = sc.read_h5ad(h5ad_path)

        for gene in adata.var_names:
            sample_info_list.append({
                'value': gene,
                'label': gene 
            })

    except Exception as e:
        print(f"Failed to read {h5ad_path}: {str(e)}")

    return sample_info_list


# return gene list(including gene numbers)
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
    if sample_id not in SAMPLES:
        raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES.")
    
    adata_path = SAMPLES[sample_id]["adata"]
    adata = sc.read_h5ad(adata_path)

    # filter cells based on cell_ids
    selected_cells_mask = adata.obs.index.isin(cell_ids)
    filtered_adata = adata[selected_cells_mask]
    
    all_genes = set()
    cell_expressions = {}

    for i, cell in enumerate(filtered_adata.obs.index):
        expression_values = filtered_adata[i].X.A[0] if hasattr(filtered_adata[i].X, 'A') else filtered_adata[i].X[0]
        nonzero_indices = np.where(expression_values > 0)[0]
        cell_expression = {filtered_adata.var.index[j]: float(expression_values[j]) for j in nonzero_indices}

        cell_expressions[cell] = cell_expression
        all_genes.update(cell_expression.keys())

    all_genes = sorted(all_genes)

    expression_data = [
        {
            "cell_id": cell,
            "expression": [cell_expressions[cell].get(gene, 0.0) for gene in all_genes]
        }
        for cell in filtered_adata.obs.index
    ]

    cell_type_annotations = filtered_adata.obs["cell_type"].to_dict()

    return {
        "metadata": {
            "cell_ids": list(filtered_adata.obs.index),
            "genes": all_genes,
            "cell_type_annotations": cell_type_annotations,
        },
        "expression_data": expression_data,
    }


def get_NMF_GO_data(sample_id, cell_list):
    # finding the best n_neighbors for leiden clustering
    def compute_silhouette_scores(adata, n_neighbors_list=[5, 10, 15, 20, 30]):
        silhouette_scores = {}

        for n_neighbors in n_neighbors_list:
            print(f"Trying n_neighbors = {n_neighbors}")
            # calculate neighbors
            sc.pp.neighbors(adata, use_rep='X_nmf', n_neighbors=n_neighbors)

            # leiden clustering
            sc.tl.leiden(adata, resolution=0.5)
            
            # calculate silhouette score
            labels = adata.obs['leiden'].astype(int)
            silhouette_avg = silhouette_score(adata.obsm['X_nmf'], labels)

            silhouette_scores[n_neighbors] = silhouette_avg
            print(f"Silhouette score for n_neighbors={n_neighbors}: {silhouette_avg:.3f}")
        
        return silhouette_scores

    # compute cophenetic correlation method
    def compute_cophenetic(W):
        try:
            dist = pdist(W.T)
            linkage_matrix = linkage(dist, method='average')
            coph_corr, _ = cophenet(linkage_matrix, dist)
            return coph_corr
        except Exception:
            return np.nan

    # finding the best k for NMF
    def auto_select_nmf_k_from_expr(expr_matrix, k_range=range(2, 21), n_repeats=5, random_state=42, coph_threshold=0.98):
        results = []

        for k in k_range:
            cophs = []
            errors = []
            for i in range(n_repeats):
                nmf = NMF(n_components=k, init='nndsvda', random_state=random_state+i, max_iter=1000)
                W = nmf.fit_transform(expr_matrix)
                H = nmf.components_
                recon = np.dot(W, H)
                error = np.linalg.norm(expr_matrix - recon)
                coph = compute_cophenetic(W)
                errors.append(error)
                cophs.append(coph)

            avg_coph = np.nanmean(cophs)
            avg_error = np.mean(errors)
            results.append((k, avg_coph, avg_error))

        # filtered cophenetic equals to nan
        valid_results = [(k, coph) for k, coph, _ in results if not np.isnan(coph)]

        # choose the min k with the first cophenetic >= 0.97
        for k, coph in valid_results:
            if coph >= coph_threshold:
                return k, results

        best_k = max(valid_results, key=lambda x: x[1])[0]
        return best_k, results

    # get top genes of each component(NMF)
    def get_top_genes(H, gene_names, top_n=10):
        top_genes = {}
        for i, comp in enumerate(H):
            # getting the indices of the top genes
            top_idx = np.argsort(comp)[::-1][:top_n]
            top_genes[f"Component_{i+1}"] = [gene_names[j] for j in top_idx]
        return top_genes

    if sample_id not in SAMPLES:
        raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES.")

    # ========== Load the data for the specified sample ID ========== 
    adata_path = SAMPLES[sample_id]["adata"]
    adata = sc.read_h5ad(adata_path)

    adata_region = adata[cell_list, :].copy()
    expr_matrix = adata_region.X
    if not isinstance(expr_matrix, np.ndarray):
        expr_matrix = expr_matrix.toarray()

    # ========== find the best component number for NMF ==========
    best_k, k_results = auto_select_nmf_k_from_expr(expr_matrix)

    for k, coph, err in k_results:
        print(f"k={k}, Cophenetic={coph:.3f}, Error={err:.2f}")

    # ========== NMF ==========
    n_components = best_k
    nmf_model = NMF(n_components=n_components, init='nndsvda', random_state=42)
    W = nmf_model.fit_transform(expr_matrix)
    H = nmf_model.components_ 

    # ========== clustering NMF result(M) ==========
    adata_region.obsm['X_nmf'] = W
    sil_scores = compute_silhouette_scores(adata_region, n_neighbors_list=[5, 10, 15, 20, 30])

    print("\nSilhouette scores for different n_neighbors:")
    for n, score in sil_scores.items():
        print(f"n_neighbors = {n}, silhouette score = {score:.3f}")

    # find the best n_neighbors based on silhouette score
    best_n_neighbors = max(sil_scores, key=sil_scores.get)
    print(f"\nBest n_neighbors based on silhouette score: {best_n_neighbors}")

    sc.pp.neighbors(adata_region, use_rep='X_nmf', n_neighbors=best_n_neighbors)
    sc.tl.leiden(adata_region, resolution=0.1)

    clusters = adata_region.obs['leiden']

    # convert W to a DataFrame for easier manipulation
    df_W = pd.DataFrame(W, index=adata_region.obs_names,
                        columns=[f"Component_{i+1}" for i in range(W.shape[1])])
    df_W["cluster"] = clusters.values

    # calculate the average activation of each NMF component for each cluster
    cluster_means = df_W.groupby("cluster").mean()

    # each cluster's cell ids
    cell_ids_by_cluster = {
        cluster: adata_region.obs.index[adata_region.obs['leiden'] == cluster].tolist()
        for cluster in adata_region.obs['leiden'].unique()
    }

    # ========== Go analysis based on the result of NMF(H) ==========
    gene_names = adata.var_names.tolist()
    top_genes = get_top_genes(H, gene_names, top_n=10)

    go_results = {}

    for comp, genes in top_genes.items():
        print(f"analyzing {comp} ...")
        enr = gp.enrich(
            gene_list=genes,
            gene_sets="../Data/c5.go.v2024.1.Hs.symbols.gmt",
            outdir=None,
            cutoff=0.5,
        )

        filtered = enr.results[enr.results["Adjusted P-value"] < 0.05]
        filtered = filtered.sort_values(by="Combined Score", ascending=False)
        filtered_top5 = filtered.head(5)
        if not filtered.empty:
            go_results[comp] = filtered_top5.to_dict(orient="records")
        else:
            print(f"{comp} no GO results found.")
    
    return {
        "NMF_matrix": W.tolist(),
        "GO_results": go_results,
        "cluster_means": cluster_means.to_dict(orient="records"),
        "cell_ids_by_cluster": cell_ids_by_cluster,
    }


def get_cell_cell_interaction_data(sample_id, receiver, sender, receiverGene, senderGene, cellIds):
    result = {}
    
    if sample_id not in SAMPLES:
        print(f"Error: Sample ID {sample_id} not found in SAMPLES.")
        return result
    
    adata_path = SAMPLES[sample_id]["adata"]
    adata = sc.read_h5ad(adata_path)
    
    filtered_adata = adata[adata.obs.index.isin(cellIds)]
    
    filtered_spatial = pd.DataFrame(
        filtered_adata.obsm["spatial"],
        columns=["cell_x", "cell_y"],
        index=filtered_adata.obs.index,
    )
    filtered_spatial["cell_type"] = filtered_adata.obs["cell_type"]
    filtered_spatial.rename(columns={"cell_x": "X", "cell_y": "Y"}, inplace=True)
    
    spatial_file = f"{sample_id}_spatial.txt"
    filtered_spatial.to_csv(spatial_file, sep="\t", index=True, index_label="")
    
    counts_data = filtered_adata.to_df()
    counts_file = f"{sample_id}_counts.txt"
    counts_data.to_csv(counts_file, sep="\t", index=True, index_label="")
    
    script_path = "../Spacia/spacia.py"
    output_path = "cell2cellinteractionOutput"

    if isinstance(receiverGene, list):
        receiverGene = "|".join(receiverGene)
    if isinstance(senderGene, list):
        senderGene = "|".join(senderGene)

    params = f'-rc {receiver} -sc {sender} -rf "{receiverGene}" -sf {senderGene} -d 30 -nc 20'
    
    cmd = f"python {script_path} {counts_file} {spatial_file} {params} -o {output_path}"
    print(f"Running command: {cmd}")
    
    os.system(cmd)
    
    interaction_file = os.path.join(output_path, "interaction.txt")
    if os.path.exists(interaction_file):
        interaction_df = pd.read_csv(interaction_file, sep="\t")
        result[sample_id] = interaction_df.to_dict(orient="records")
    else:
        print(f"Error: Interaction file for {sample_id} not found.")
    
    return result