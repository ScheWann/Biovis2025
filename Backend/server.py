from flask import Flask, request, jsonify, send_from_directory, send_file
from flask_cors import CORS
import re
import os
import subprocess
import json
import sys
from process import (
    get_hires_image_size,
    get_unique_cell_types,
    get_cell_type_coordinates,
    get_samples,
    get_cell_types,
    get_gene_list,
    get_gene_list_for_cell2cellinteraction,
    get_kosara_data,
    get_selected_region_data,
    get_NMF_GO_data,
    get_cell_cell_interaction_data
)

app = Flask(__name__)
CORS(app)

UPLOAD_FOLDER = 'Data'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/', methods=['GET'])
def get_helloword():
    """Basic test endpoint"""
    return 'Hello World!'


@app.route('/get_available_samples', methods=['GET'])
def get_available_samples():
    """Get list of available samples"""
    return jsonify(get_samples())


@app.route('/get_hires_image_size', methods=['POST'])
def get_hires_image_size_route():
    """Get high-resolution image size for selected samples"""
    sample_ids = request.json['sample_ids']
    return jsonify(get_hires_image_size(sample_ids))


@app.route('/get_unique_cell_types', methods=['POST'])
def get_unique_cell_types_route():
    """Get unique cell types for selected samples"""
    sample_ids = request.json['sample_ids']
    return jsonify(get_unique_cell_types(sample_ids))


@app.route('/get_tile', methods=['GET'])
def serve_tile():
    """Serve image tiles for visualization"""
    sample_id = request.args.get('sample_id')
    x = request.args.get('x', type=int)
    y = request.args.get('y', type=int)

    if x is None or y is None:
        return jsonify({'error': 'Missing x or y parameters'}), 400

    filename = f"tile_{x*256}_{y*256}.tif"

    processed_sample_id = f"{sample_id}_processed"
    processed_sample_tiles = f"{sample_id}_processed_tiles"
    tile_dir = os.path.join("../Data", processed_sample_id, processed_sample_tiles)

    return send_from_directory(tile_dir, filename)


@app.route('/get_cell_type_coordinates', methods=['POST'])
def get_cell_type_coordinates_route():
    """Get cell type coordinates for selected samples"""
    sample_ids = request.json['sample_ids']
    return jsonify(get_cell_type_coordinates(sample_ids))


@app.route('/get_cell_types', methods=['POST'])
def get_cell_types_route():
    """Get cell types for selected samples"""
    sample_name = request.json['sample_name']
    return jsonify(get_cell_types(sample_name))


@app.route('/get_all_gene_list', methods=['POST'])
def get_all_gene_list():
    """Get list of all genes for selected samples"""
    sample_names = request.json['sample_names']
    return jsonify(get_gene_list(sample_names))


@app.route('/get_cell2cell_gene_list', methods=['POST'])
def get_cell2cell_gene_list_route():
    """Get list of all genes for selected samples"""
    sample_name = request.json['sample_name']
    return jsonify(get_gene_list_for_cell2cellinteraction(sample_name))


@app.route('/get_kosara_data', methods=['POST'])
def get_kosara_data_route():
    """Get Kosara visualization data"""
    sample_ids = request.json['sample_ids']
    gene_list = request.json['gene_list']
    cell_list = request.json['cell_list']
    return jsonify(get_kosara_data(sample_ids, gene_list, cell_list))


@app.route('/get_selected_region_data', methods=['POST'])
def get_selected_region_data_route():
    """Get gene expressiondata for selected regions"""
    sample_id = request.json['sample_id']
    cell_list = request.json['cell_list']
    return jsonify(get_selected_region_data(sample_id, cell_list))


@app.route('/get_NMF_GO_data', methods=['POST'])
def get_NMF_GO_data_route():
    """Get NMF GO data"""
    regions = request.json['regions']
    n_component = request.json['n_component']
    resolution = request.json['resolution']
    return jsonify(get_NMF_GO_data(regions, n_component, resolution))


@app.route('/get_cell_cell_interaction_data', methods=['POST'])
def get_cell_cell_interaction_data_route():
    """Get cell-cell interaction data"""
    regions = request.json['regions']
    receiver = request.json['receiver']
    sender = request.json['sender']
    receiverGene = request.json['receiverGene']
    senderGene = request.json['senderGene']
    return jsonify(get_cell_cell_interaction_data(regions, receiver, sender, receiverGene, senderGene))


@app.route('/get_gene_name_search')
def get_gene_name_search():
    """Search for genes by name"""
    query = request.args.get('q', '').strip().lower()
    gene_list = get_gene_list()
    
    # Return the full list if no query is provided
    if not query:
        return jsonify(gene_list)

    pattern = re.compile(re.escape(query), re.IGNORECASE)
    results = [item for item in gene_list if pattern.search(item)]
    
    return jsonify(results)


@app.route('/upload_spaceranger', methods=['POST'])
def upload_spaceranger():
    name = request.form.get('name')
    description = request.form.get('description')
    files = request.files.getlist('files')

    # Save files, preserving their relative paths
    for file in files:
        rel_path = file.filename  # This should be the relative path sent from frontend
        save_path = os.path.join(UPLOAD_FOLDER, name, rel_path)
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        file.save(save_path)

    # Optionally: validate required files exist, trigger downstream processing, etc.
    return jsonify({'status': 'success'})


if __name__ == "__main__":
    app.run(debug=True, port=5003)
