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
    get_gene_list,
    get_kosara_data,
    # get_umap_positions_with_clusters,
    # get_gene_list,
    # get_specific_gene_expression
)
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

# Define workspace root
workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

app = Flask(__name__, static_folder='static')
CORS(app)  # Enable CORS for all routes

# Ensure static directories exist
os.makedirs(os.path.join(app.static_folder, 'figures'), exist_ok=True)

@app.route('/test-image')
def test_image():
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
            
        # Serve the file
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
    processed_sample_tiles = f"{sample_id}_processed_tiles"
    tile_dir = os.path.join("../Data", processed_sample_id, processed_sample_tiles)

    return send_from_directory(tile_dir, filename)


@app.route('/get_cell_type_coordinates', methods=['POST'])
def get_cell_type_coordinates_route():
    sample_id = request.json['sample_id']
    return jsonify(get_cell_type_coordinates(sample_id).to_dict(orient='records'))


@app.route('/get_all_gene_list', methods=['POST'])
def get_all_gene_list():
    sample_names = request.json['sample_names']
    return jsonify(get_gene_list(sample_names))


@app.route('/get_kosara_data', methods=['POST'])
def get_kosara_data_route():
    sample_ids = request.json['sample_ids']
    gene_list = request.json['gene_list']
    cell_list = request.json['cell_list']
    return jsonify(get_kosara_data(sample_ids, gene_list, cell_list))

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
        sample_percent = request.args.get('sample_percent', default=1.0, type=float)
        step = request.args.get('step', default=0, type=int)
        
        # Path to the DEAPLOG script and data
        workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        script_path = os.path.join(workspace_root, 'Python', 'DEAPLOG.py')
        data_path = 'Data/skin_TXK6Z4X_A1_processed/tmap/weighted_by_area_celltypist_cells_adata.h5'  # Using relative path
        
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
            return jsonify({'error': error_msg}), 500

        # Ensure the data file exists
        abs_data_path = os.path.join(workspace_root, data_path)
        if not os.path.exists(abs_data_path):
            error_msg = f'Data file not found at: {abs_data_path}'
            print(f"Error: {error_msg}")
            return jsonify({'error': error_msg}), 500

        # Run the DEAPLOG script from workspace root
        cmd = ['python', script_path, '--sample_percent', str(sample_percent), '--step', str(step), '--data_path', data_path]
        print(f"Debug - Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=workspace_root  # Set working directory to workspace root
        )
        
        print("Debug - Command output:")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        
        if result.returncode != 0:
            error_msg = f'DEAPLOG process failed with code {result.returncode}: {result.stderr}'
            print(f"Error: {error_msg}")
            return jsonify({'error': error_msg}), 500
        
        # Try to parse the JSON output from the script
        try:
            # Split the output into lines and find the last line that looks like JSON
            lines = result.stdout.strip().split('\n')
            response_data = {
                'sample_percent': sample_percent,
                'step': step,
                'n_cells': 0,
                'n_clusters': 0,
                'marker_genes': [],
            }
            
            for line in reversed(lines):
                line = line.strip()
                if line.startswith('{') and line.endswith('}'):
                    try:
                        parsed_data = json.loads(line)
                        print("Debug - Parsed JSON data:", parsed_data)
                        response_data.update(parsed_data)
                        break
                    except json.JSONDecodeError as e:
                        print(f"Warning - Failed to parse JSON line: {e}")
                        continue
            
            print("Debug - Final response data:", response_data)
            return jsonify(response_data)
            
        except Exception as e:
            error_msg = f'Error parsing DEAPLOG output: {str(e)}'
            print(f"Error: {error_msg}")
            return jsonify({'error': error_msg}), 500
            
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

if __name__ == "__main__":
    app.run(debug=True, port=5000)
