from flask import Flask, request, jsonify, send_from_directory, send_file
from flask_cors import CORS
import re
import os
import subprocess
import json
from process import (
    # get_um_positions_with_clusters, 
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
    # get_umap_positions_with_clusters,
    # get_gene_list,
    # get_specific_gene_expression
)
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import sys
from functools import lru_cache
import time

# Add the Python directory to the system path for importing DEAPLOG module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Python'))

from DEAPLOG import run_deaplog_analysis

# Define workspace root for file paths
workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Initialize Flask app with static folder
app = Flask(__name__, static_folder='static')
CORS(app)  # Enable CORS for all routes

# Ensure static directories exist for storing figures
os.makedirs(os.path.join(app.static_folder, 'figures'), exist_ok=True)

# Cache for DEAPLOG results
@lru_cache(maxsize=10)
def get_cached_deaplog_results(sample_percent, step):
    """Cached version of DEAPLOG results"""
    try:
        # Path to the DEAPLOG script and data
        script_path = os.path.join(workspace_root, 'Python', 'DEAPLOG.py')
        data_path = os.path.join(workspace_root, 'Data', 'skin_TXK6Z4X_A1_processed', 'tmap', 'weighted_by_area_celltypist_cells_adata.h5')
        
        print(f"Debug - Parameters:")
        print(f"sample_percent: {sample_percent}")
        print(f"step: {step}")
        print(f"workspace_root: {workspace_root}")
        print(f"script_path: {script_path}")
        print(f"data_path: {data_path}")
        
        # Ensure the script exists
        if not os.path.exists(script_path):
            error_msg = f'DEAPLOG script not found at: {script_path}'
            print(f"Error: {error_msg}")
            return {'error': error_msg}, 500

        # Ensure the data file exists
        if not os.path.exists(data_path):
            error_msg = f'Data file not found at: {data_path}'
            print(f"Error: {error_msg}")
            return {'error': error_msg}, 500

        # Run the DEAPLOG script
        cmd = ['python', script_path, '--sample_percent', str(sample_percent), '--step', str(step), '--data_path', data_path]
        print(f"Debug - Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=workspace_root
        )
        
        print("Debug - Command output:")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        
        if result.returncode != 0:
            error_msg = f'DEAPLOG process failed with code {result.returncode}: {result.stderr}'
            print(f"Error: {error_msg}")
            return {'error': error_msg}, 500
            
        # Parse JSON output
        try:
            output_lines = result.stdout.strip().split('\n')
            json_data = None
            for line in reversed(output_lines):
                if line.strip().startswith('{'):
                    try:
                        json_data = json.loads(line)
                        break
                    except json.JSONDecodeError:
                        continue
                        
            if not json_data:
                error_msg = 'No valid JSON data found in DEAPLOG output'
                print(f"Error: {error_msg}")
                return {'error': error_msg}, 500
                
            print("Debug - Parsed JSON data:", json_data)
            return json_data
            
        except Exception as e:
            error_msg = f'Error parsing DEAPLOG output: {str(e)}'
            print(f"Error: {error_msg}")
            return {'error': error_msg}, 500
            
    except Exception as e:
        error_msg = f'Internal server error: {str(e)}'
        print(f"Error: {error_msg}")
        return {'error': error_msg}, 500

@app.route('/test-image')
def test_image():
    """Test endpoint to check available images in the figures directory"""
    try:
        # Get the absolute path to the figures directory
        figures_dir = os.path.join(workspace_root, 'Python', 'figures')
        print(f"Checking figures directory: {figures_dir}")
        
        # Check if the directory exists
        if not os.path.exists(figures_dir):
            print(f"Figures directory not found: {figures_dir}")
            return jsonify({'images': []})
            
        # List all PNG files in the directory
        images = [f for f in os.listdir(figures_dir) if f.endswith('.png')]
        print(f"Found images: {images}")
        
        return jsonify({'images': images})
    except Exception as e:
        print(f"Error listing images: {str(e)}")
        return jsonify({'images': []})

@app.route('/figures/<path:filename>')
def serve_figure(filename):
    """Serve figure files from the figures directory"""
    try:
        # Get the absolute path to the figures directory
        figures_dir = os.path.join(workspace_root, 'Python', 'figures')
        print(f"Serving figure from: {figures_dir}")
        
        # Construct the full path to the requested file
        file_path = os.path.join(figures_dir, filename)
        print(f"Requested file path: {file_path}")
        
        # Check if the file exists
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            return jsonify({'error': 'Image not found'}), 404
            
        # Serve the file with no caching
        response = send_file(file_path)
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '0'
        return response
    except Exception as e:
        print(f"Error serving figure: {str(e)}")
        return jsonify({'error': str(e)}), 500

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
    sample_id = request.json['sample_id']
    cell_list = request.json['cell_list']
    return jsonify(get_NMF_GO_data(sample_id, cell_list))

@app.route('/get_cell_cell_interaction_data', methods=['POST'])
def get_cell_cell_interaction_data_route():
    """Get cell-cell interaction data"""
    sample_id = request.json['sample_id']
    receiver = request.json['receiver']
    sender = request.json['sender']
    receiverGene = request.json['receiverGene']
    senderGene = request.json['senderGene']
    cellIds = request.json['cellIds']
    return jsonify(get_cell_cell_interaction_data(sample_id, receiver, sender, receiverGene, senderGene, cellIds))

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

@app.route('/get_specific_gene_expression', methods=['POST'])
def get_specific_gene_expression_route():
    bin_size = request.json['bin_size']
    gene_name = request.json['gene_name']
    return jsonify(get_specific_gene_expression(bin_size, gene_name).to_dict(orient='records'))

@app.route('/get_deaplog_results', methods=['GET'])
def get_deaplog_results():
    try:
        sample_percent = request.args.get('sample_percent', default=0.01, type=float)
        step = request.args.get('step', default=0, type=int)
        
        # Get cached results
        results = get_cached_deaplog_results(sample_percent, step)
        
        # If results is a tuple (error case), return it directly
        if isinstance(results, tuple):
            return jsonify(results[0]), results[1]
            
        return jsonify(results)
        
    except Exception as e:
        error_msg = f'Internal server error: {str(e)}'
        print(f"Error: {error_msg}")
        return jsonify({'error': error_msg}), 500

@app.route('/run_deaplog', methods=['POST'])
def run_deaplog():
    try:
        # Get parameters from request
        data = request.get_json()
        sample_percent = float(data.get('sample_percent', 0.1))
        
        # Run DEAPLOG analysis
        results = run_deaplog_analysis(adata, sample_percent)
        
        # Save results to files
        with open('static/deaplog_results.json', 'w') as f:
            json.dump(results, f)
        
        # Save UMAP plot
        plt.figure(figsize=(10, 10))
        sc.pl.umap(adata, color='leiden', save='_deaplog.png')
        
        # Save PAGA plot
        plt.figure(figsize=(10, 10))
        sc.pl.paga(adata, save='_deaplog.png')
        
        return jsonify({
            'status': 'success',
            'message': 'DEAPLOG analysis completed successfully',
            'results': results
        })
        
    except Exception as e:
        print(f"Error in DEAPLOG analysis: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

@app.route('/api/analyze', methods=['POST'])
def analyze():
    """Run DEAPLOG analysis on the provided data"""
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

if __name__ == "__main__":
    app.run(debug=True, port=5003)
