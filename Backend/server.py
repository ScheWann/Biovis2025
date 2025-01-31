from flask import Flask, request, jsonify
import re
from process import (
    get_um_positions_with_clusters, 
    get_hires_image_size,
    get_umap_positions_with_clusters,
    get_gene_list,
    get_specific_gene_expression
)

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_helloword():
    return 'Hello World!'


@app.route('/get_hires_image_size', methods=['GET'])
def get_hires_image_size_route():
    return jsonify(get_hires_image_size())


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
    app.run(debug=True)
