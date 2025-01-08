from flask import Flask, request, jsonify
from process import (
    get_um_008_positions_with_clusters, 
    get_hires_image_size
)

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_helloword():
    return 'Hello World!'


@app.route('/get_hires_image_size', methods=['GET'])
def get_hires_image_size_route():
    return jsonify(get_hires_image_size())


@app.route('/get_um_008_positions_with_clusters', methods=['POST'])
def get_um_008_positions_with_clusters_route():
    kmeans = request.json['kmeans']
    return jsonify(get_um_008_positions_with_clusters(kmeans).to_dict(orient='records'))

if __name__ == "__main__":
    app.run(debug=True)
