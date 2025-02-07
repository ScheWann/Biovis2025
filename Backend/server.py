from flask import Flask, request, jsonify, send_from_directory
import re
import os
from process import (
    # get_um_positions_with_clusters, 
    get_hires_image_size,
    get_unique_cell_types,
    get_cell_type_coordinates,
    get_samples,
    # get_umap_positions_with_clusters,
    # get_gene_list,
    # get_specific_gene_expression
)

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_helloword():
    return 'Hello World!'


@app.route('/get_available_samples', methods=['GET'])
def get_available_samples():
    return jsonify(get_samples())


@app.route('/get_hires_image_size', methods=['POST'])
def get_hires_image_size_route():
    sample_id = request.json['sample_id']
    return jsonify(get_hires_image_size(sample_id))


@app.route('/get_unique_cell_types', methods=['POST'])
def get_unique_cell_types_route():
    sample_id = request.json['sample_id']
    return jsonify(get_unique_cell_types(sample_id))


@app.route('/get_tile', methods=['GET'])
def serve_tile():
    sample_id = request.args.get('sample_id')
    x = request.args.get('x', type=int)
    y = request.args.get('y', type=int)

    if x is None or y is None:
        return jsonify({'error': 'Missing x or y parameters'}), 400

    filename = f"tile_{x*256}_{y*256}.tif"

    processed_sample_id = f"{sample_id}_processed"
    tile_dir = os.path.join("../Data", processed_sample_id, "skin_TXK6Z4X_A1_processed_tiles")

    return send_from_directory(tile_dir, filename)


@app.route('/get_cell_type_coordinates', methods=['POST'])
def get_cell_type_coordinates_route():
    sample_id = request.json['sample_id']
    return jsonify(get_cell_type_coordinates(sample_id).to_dict(orient='records'))


#################### OLD CODE ####################
@app.route('/get_um_positions_with_clusters', methods=['POST'])
def get_um_positions_with_clusters_route():
    bin_size = request.json['bin_size']
    kmeans = request.json['kmeans']
    return jsonify(get_um_positions_with_clusters(bin_size, kmeans).to_dict(orient='records'))


@app.route('/get_umap_positions', methods=['POST'])
def get_umap_positions_route():
    bin_size = request.json['bin_size']
    kmeans = request.json['kmeans']
    return jsonify(get_umap_positions_with_clusters(bin_size, kmeans).to_dict(orient='records'))


@app.route('/get_full_gene_list', methods=['GET'])
def get_full_gene_list():
    return jsonify(get_gene_list())


@app.route('/get_gene_name_search')
def get_gene_name_search():
    query = request.args.get('q', '').strip().lower()
    gene_list = get_gene_list()
    
    # Return the full list if no query is provided
    if not query:
        return jsonify(gene_list)

    pattern = re.compile(re.escape(query), re.IGNORECASE)
    results = [item for item in gene_list if pattern.search(item)]
    
    return jsonify(results)

@app.route('/get_specific_gene_expression', methods=['POST'])
def get_specific_gene_expression_route():
    bin_size = request.json['bin_size']
    gene_name = request.json['gene_name']
    return jsonify(get_specific_gene_expression(bin_size, gene_name).to_dict(orient='records'))

if __name__ == "__main__":
    app.run(debug=True, threaded=True)
