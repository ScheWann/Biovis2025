from flask import Flask, request, jsonify
from process import (
    get_um_positions_with_clusters, 
    get_hires_image_size,
    get_umap_positions
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
    return jsonify(get_umap_positions(bin_size, kmeans).to_dict(orient='records'))

if __name__ == "__main__":
    app.run(debug=True)
