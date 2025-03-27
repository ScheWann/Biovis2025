from flask import Flask, request, jsonify
from flask_cors import CORS
import scanpy as sc
import os
from DEAPLOG import run_deaplog_analysis

app = Flask(__name__)
CORS(app)

@app.route('/api/analyze', methods=['POST'])
def analyze():
    try:
        # Get parameters from request
        data = request.get_json()
        data_path = data.get('data_path')
        sample_percent = data.get('sample_percent')
        
        if not data_path:
            return jsonify({'error': 'data_path is required'}), 400
            
        # Load data
        adata = sc.read_h5ad(data_path)
        rdata = adata.copy()
        
        # Run analysis
        results = run_deaplog_analysis(rdata, adata, sample_percent)
        
        return jsonify(results)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000) 