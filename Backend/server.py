from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import re
import os
import io
from process import (
    get_samples_option,
    get_hires_image_size,
    get_coordinates,
    get_gene_list,
    get_kosara_data,
    get_selected_region_data,
    get_hires_image_crop,
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


@app.route("/api/get_hires_crop", methods=["POST"])
def get_hires_crop_route():
    """
    Get a high-resolution crop of the image for the specified region
    """
    data = request.json
    sample_id = data["sample_id"]
    center_x = data["center_x"]
    center_y = data["center_y"]
    crop_size = data.get("crop_size", 512)  # Default crop size
    
    try:
        image_bytes = get_hires_image_crop(sample_id, center_x, center_y, crop_size)
        return send_file(
            io.BytesIO(image_bytes),
            mimetype='image/jpeg',
            as_attachment=False
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/get_hires_image", methods=["POST"])
def get_hires_image_route():
    """
    Return the full high-resolution image for the given sample_id as JPEG.
    """
    data = request.json
    sample_id = data["sample_id"]
    from process import SAMPLES
    import PIL.Image
    if sample_id not in SAMPLES:
        return jsonify({"error": f"Sample {sample_id} not found"}), 404
    image_path = SAMPLES[sample_id]["image_tif_path"]
    try:
        image = PIL.Image.open(image_path)
        if image.mode != 'RGB':
            image = image.convert('RGB')
        img_byte_arr = io.BytesIO()
        image.save(img_byte_arr, format='JPEG', quality=85)
        img_byte_arr.seek(0)
        return send_file(
            img_byte_arr,
            mimetype='image/jpeg',
            as_attachment=False
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 500


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


if __name__ == "__main__":
    app.run(debug=True, port=5003)
