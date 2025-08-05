from flask import Flask, request, jsonify, send_file
from process import SAMPLES
from flask_cors import CORS
import re
import os
from process import (
    get_samples_option,
    get_hires_image_size,
    get_coordinates,
    get_gene_list,
    get_kosara_data,
    get_selected_region_data,
    get_umap_data,
    perform_go_analysis,
    get_trajectory_data,
    get_trajectory_gene_list,
    load_adata_to_cache,
    clear_adata_cache,
    get_pseudotime_data,
    get_trajectory_gene_expression,
)


app = Flask(__name__)
CORS(app)

UPLOAD_FOLDER = "../Uploaded_Data"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)


@app.route("/api/get_samples_option", methods=["GET"])
def get_samples_option_route():
    """
    Get a list of available samples for the selector, grouped by cell scale(e.g., 2um, 8um)
    """
    return jsonify(get_samples_option())


@app.route("/api/load_adata_cache", methods=["POST"])
def load_adata_cache_route():
    """
    Load AnnData objects for the given sample IDs into the global cache.
    This should be called once when samples are confirmed.
    """
    sample_ids = request.json["sample_ids"]
    try:
        load_adata_to_cache(sample_ids)
        return jsonify({"status": "success", "message": f"Loaded AnnData for {len(sample_ids)} samples"})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/clear_adata_cache", methods=["POST"])
def clear_adata_cache_route():
    """
    Clear the global AnnData cache to free memory.
    """
    try:
        clear_adata_cache()
        return jsonify({"status": "success", "message": "AnnData cache cleared"})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_hires_image_size", methods=["POST"])
def get_hires_image_size_route():
    """
    Get high-resolution image size for the selected samples
    """
    sample_ids = request.json["sample_ids"]
    return jsonify(get_hires_image_size(sample_ids))


@app.route("/api/get_coordinates", methods=["POST"])
def get_coordinates_route():
    """
    Get coordinates for the selected samples
    """
    sample_ids = request.json["sample_ids"]
    return jsonify(get_coordinates(sample_ids))


@app.route("/api/get_gene_list", methods=["POST"])
def get_gene_list_route():
    """
    Get list of all genes for the selected samples
    """
    sample_ids = request.json["sample_ids"]
    return jsonify(get_gene_list(sample_ids))


@app.route("/api/get_kosara_data", methods=["POST"])
def get_kosara_data_route():
    """
    Get Kosara visualization format data
    """
    sample_ids = request.json["sample_ids"]
    gene_list = request.json["gene_list"]
    cell_list = request.json["cell_list"]
    return jsonify(get_kosara_data(sample_ids, gene_list, cell_list))


@app.route("/api/get_selected_region_data", methods=["POST"])
def get_selected_region_data_route():
    """
    Get gene expression data for the selected regions
    """
    sample_id = request.json["sample_id"]
    cell_list = request.json["cell_list"]
    return jsonify(get_selected_region_data(sample_id, cell_list))


@app.route("/api/get_umap_data", methods=["POST"])
def get_umap_data_route():
    """
    Generate UMAP data from gene expression data
    """
    sample_id = request.json["sample_id"]
    cell_ids = request.json.get("cell_ids", None)  # New parameter for specific cells
    n_neighbors = request.json.get("n_neighbors", 10)
    n_pcas = request.json.get("n_pcas", 30)
    resolutions = request.json.get("resolutions", 1)
    adata_umap_title = request.json.get("adata_umap_title", None)
    try:
        umap_data = get_umap_data(
            sample_id=sample_id,
            cell_ids=cell_ids,  # Pass cell_ids to the backend function
            n_neighbors=n_neighbors,
            n_pcas=n_pcas,
            resolutions=resolutions,
            adata_umap_title=adata_umap_title
        )
        return jsonify(umap_data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_go_analysis", methods=["POST"])
def get_go_analysis_route():
    """
    Perform GO analysis on selected cluster cells
    """
    sample_id = request.json["sample_id"]
    cluster_id = request.json["cluster_id"]
    adata_umap_title = request.json["adata_umap_title"]

    try:
        go_results = perform_go_analysis(
            sample_id=sample_id,
            cluster_id=cluster_id,
            adata_umap_title=adata_umap_title
        )
        return jsonify(go_results)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_trajectory_gene_list", methods=["POST"])
def get_trajectory_gene_list_route():
    """
    Get list of available genes from trajectory data
    """
    sample_id = request.json["sample_id"]
    is_vertical = request.json.get("is_vertical")

    try:
        gene_list = get_trajectory_gene_list(sample_id=sample_id, is_vertical=is_vertical)
        return jsonify(gene_list)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_trajectory_data", methods=["POST"])
def get_trajectory_data_route():
    """
    Get trajectory gene expression data for line chart visualization
    """
    sample_id = request.json["sample_id"]
    selected_genes = request.json.get("selected_genes", None)
    is_vertical = request.json.get("is_vertical")

    try:
        trajectory_data = get_trajectory_data(sample_id=sample_id, selected_genes=selected_genes, is_vertical=is_vertical)
        return jsonify(trajectory_data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_gene_name_search", methods=["POST"])
def get_gene_name_search():
    """
    Search for gene names based on a query string and return matching gene names.
    """
    sample_id = request.json["sample_id"]
    gene_name = request.json["gene_name"]

    if not gene_name or len(gene_name) < 2:
        return jsonify([])

    sample_gene_dict = get_gene_list([sample_id])

    if sample_id not in sample_gene_dict:
        return jsonify([])

    gene_list = sample_gene_dict[sample_id]

    # Filter genes that match the search query
    pattern = re.compile(re.escape(gene_name), re.IGNORECASE)
    matching_genes = [gene for gene in gene_list if pattern.search(gene)]

    return jsonify(matching_genes)


@app.route("/api/get_hires_image", methods=["POST"])
def get_hires_image_route():
    """
    Return the full high-resolution image for the given sample_id as JPEG.
    """
    sample_id = request.json["sample_id"]

    if sample_id not in SAMPLES:
        return jsonify({"error": f"Sample {sample_id} not found"}), 404

    image_path = SAMPLES[sample_id]["image_jpeg_path"]

    return send_file(image_path, mimetype="image/jpeg", as_attachment=False)


@app.route("/api/get_cell_boundary_image", methods=["POST"])
def get_cell_boundary_image_route():
    """
    Return the cell boundary image for the given sample_id as PNG.
    """
    sample_id = request.json["sample_id"]

    if sample_id not in SAMPLES:
        return jsonify({"error": f"Sample {sample_id} not found"}), 404

    cell_boundary_path = SAMPLES[sample_id]["cell_boundary_path"]

    return send_file(cell_boundary_path, mimetype="image/png", as_attachment=False)


@app.route("/api/upload_spaceranger", methods=["POST"])
def upload_spaceranger():
    """
    Upload Spaceranger output files and save them in the appropriate directory structure.
    """
    name = request.form.get("name")
    files = request.files.getlist("files")

    focus_patterns = [
        re.compile(r"binned_outputs/square_002um/filtered_feature_bc_matrix\.h5$"),
        re.compile(r"binned_outputs/square_008um/filtered_feature_bc_matrix\.h5$"),
        re.compile(r"binned_outputs/square_016um/filtered_feature_bc_matrix\.h5$"),
        re.compile(r"spatial/"),
    ]

    for file in files:
        rel_path = file.filename

        subdir = None
        for pattern in focus_patterns:
            if pattern.search(rel_path):
                if "binned_outputs/square_002um" in rel_path:
                    subdir = "binned_outputs/square_002um"
                elif "binned_outputs/square_008um" in rel_path:
                    subdir = "binned_outputs/square_008um"
                elif "binned_outputs/square_016um" in rel_path:
                    subdir = "binned_outputs/square_016um"
                elif "spatial/" in rel_path:

                    subdir = os.path.join(
                        "spatial", os.path.relpath(rel_path, start="spatial")
                    )
                break
        if subdir:
            if subdir.startswith("spatial"):
                save_path = os.path.join(UPLOAD_FOLDER, name, subdir)
            else:
                save_path = os.path.join(
                    UPLOAD_FOLDER, name, subdir, os.path.basename(rel_path)
                )
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            file.save(save_path)

    return jsonify({"status": "success"})


@app.route("/api/get_pseudotime_data", methods=["POST"])
def get_pseudotime_data_route():
    """
    Generate pseudotime analysis data using Slingshot trajectory inference
    """
    sample_id = request.json["sample_id"]
    cell_ids = request.json["cell_ids"]
    adata_umap_title = request.json["adata_umap_title"]
    early_markers = request.json.get("early_markers", None)
    n_neighbors = request.json.get("n_neighbors", 15)
    n_pcas = request.json.get("n_pcas", 30)
    resolutions = request.json.get("resolutions", 1)
    
    try:
        pseudotime_data = get_pseudotime_data(
            sample_id=sample_id,
            adata_umap_title=adata_umap_title,
            cell_ids=cell_ids,
            early_markers=early_markers,
            n_neighbors=n_neighbors,
            n_pcas=n_pcas,
            resolutions=resolutions
        )
        return jsonify(pseudotime_data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_trajectory_gene_expression", methods=["POST"])
def get_trajectory_gene_expression_route():
    """
    Get gene expression data along a specific trajectory path
    """
    sample_id = request.json["sample_id"]
    adata_umap_title = request.json["adata_umap_title"]
    gene_names = request.json["gene_names"]
    trajectory_path = request.json["trajectory_path"]
    
    try:
        gene_expression_data = get_trajectory_gene_expression(
            sample_id=sample_id,
            adata_umap_title=adata_umap_title,
            gene_names=gene_names,
            trajectory_path=trajectory_path
        )
        return jsonify(gene_expression_data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    app.run(debug=True, port=5003)
