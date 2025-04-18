import json
from collections import Counter
import os

# --- Configuration ---
input_filename = 'go_results2.json' # Corrected input filename based on user path
output_filename = 'sorted_gene2_counts.json'
# ---------------------

def process_gene_data(input_path, output_path):
    """Reads GO results, counts gene frequencies per component, sorts them, 
       and saves the sorted counts to a JSON file.
    """
    # Check if input file exists
    if not os.path.exists(input_path):
        print(f"Error: Input file not found at '{input_path}'") # Keep this check
        return

    # Add a print statement right before attempting to open
    print(f"Attempting to open verified input file at: {input_path}") 

    # Load the JSON data
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{input_path}'. Please ensure it's valid JSON.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{input_path}': {e}")
        return

    # Dictionary to hold final sorted gene counts per component
    sorted_component_gene_counts = {}

    # Loop through components 1 to 10
    for i in range(1, 11):
        component_key = f'Component_{i}'
        
        if component_key in data and isinstance(data[component_key], list):
            gene_counter = Counter()
            
            # Iterate through entries in the component list
            for entry in data[component_key]:
                if isinstance(entry, dict):
                    genes_str = entry.get("Genes", "")
                    if genes_str: # Ensure the string is not empty
                        # Split potentially semicolon-separated genes and remove empty strings
                        individual_genes = [gene.strip() for gene in genes_str.split(';') if gene.strip()]
                        gene_counter.update(individual_genes)
            
            # Sort genes by count (most frequent first)
            sorted_genes = dict(gene_counter.most_common())
            
            # Store in the result dictionary with lowercase component key
            sorted_component_gene_counts[component_key.lower()] = sorted_genes
        else:
            print(f"Warning: Component '{component_key}' not found or not a list in the input file.")

    # Save the sorted results to the output JSON file
    try:
        with open(output_path, 'w', encoding='utf-8') as out_file:
            json.dump(sorted_component_gene_counts, out_file, indent=2, ensure_ascii=False)
        print(f"Successfully processed data. Sorted gene counts saved to '{output_path}'")
    except Exception as e:
        print(f"An error occurred while writing to '{output_path}': {e}")

# --- Run the processing ---
if __name__ == "__main__":
    process_gene_data(input_filename, output_filename)
# ------------------------- 