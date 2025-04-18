import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import os
from collections import Counter

# --- Configuration ---
NMF_FILENAME = 'region1_NMF.csv'
GENE_COUNTS_FILENAME = 'sorted_gene1_counts.json'
OUTPUT_FILENAME = 'cluster_component_bubble_chart.png'
NUM_CLUSTERS = 10
NUM_COMPONENTS = 10
NUM_TOP_GENES_OVERALL = 10 # Number of top genes to assign unique colors for the legend
NUM_TOP_GENES_PER_COMPONENT = 5 # Number of top genes considered per component
CLUSTER_COLOR = 'blue'
DEFAULT_COMPONENT_COLOR = 'lightgrey'
BUBBLE_SIZE_SCALE = 1500  # Adjust to control component bubble sizes
CLUSTER_POSITIONS = np.array([(i % 5, -(i // 5)) for i in range(NUM_CLUSTERS)]) * 3 # Simple grid layout
COMPONENT_RADIUS = 0.8 # Radius for placing components around clusters
# ---------------------

def load_data(nmf_path, gene_path):
    """Loads NMF data from CSV and gene counts from JSON."""
    if not os.path.exists(nmf_path):
        print(f"Error: NMF file not found at 
{nmf_path}
")
        return None, None
    if not os.path.exists(gene_path):
        print(f"Error: Gene counts file not found at 
{gene_path}
")
        return None, None
    
    try:
        nmf_df = pd.read_csv(nmf_path, index_col=0) # Assuming first column is cluster index
        # Ensure index/columns match expected format (e.g., Cluster_1, Component_1)
        # Optional: Add validation here if needed
    except Exception as e:
        print(f"Error reading NMF file 
{nmf_path}
: {e}")
        return None, None

    try:
        with open(gene_path, 'r', encoding='utf-8') as f:
            gene_data = json.load(f)
    except Exception as e:
        print(f"Error reading or parsing JSON file 
{gene_path}
: {e}")
        return None, None
        
    return nmf_df, gene_data

def get_overall_top_genes(gene_data, num_top_overall, num_top_component):
    """Identifies the most frequent top genes across all components."""
    overall_gene_counts = Counter()
    for i in range(1, NUM_COMPONENTS + 1):
        comp_key = f'component_{i}' # Assumes lowercase keys in JSON
        if comp_key in gene_data:
            # Take top N genes as per the pre-sorted JSON structure
            component_top_genes = list(gene_data[comp_key].items())[:num_top_component]
            for gene, count in component_top_genes:
                 overall_gene_counts[gene] += count # Use actual counts for overall ranking
                 
    # Get the most common unique genes overall
    top_overall_genes = [gene for gene, count in overall_gene_counts.most_common(num_top_overall)]
    return top_overall_genes

def create_visualization(nmf_df, gene_data):
    """Creates the bubble chart visualization."""
    
    # --- Data Preparation ---
    # Normalize NMF data for bubble size (0 to 1 range scaling)
    max_nmf_val = nmf_df.values.max()
    min_nmf_val = nmf_df.values.min()
    if max_nmf_val == min_nmf_val: # Avoid division by zero
        normalized_nmf = nmf_df.applymap(lambda x: 0.5) # Assign a default size
    else:
        normalized_nmf = (nmf_df - min_nmf_val) / (max_nmf_val - min_nmf_val)
    
    # Get top overall genes and create a color map
    top_genes_list = get_overall_top_genes(gene_data, NUM_TOP_GENES_OVERALL, NUM_TOP_GENES_PER_COMPONENT)
    # Use a visually distinct colormap
    colors = plt.cm.get_cmap('tab10', len(top_genes_list)) 
    gene_color_map = {gene: colors(i) for i, gene in enumerate(top_genes_list)}

    # --- Plotting Setup ---
    fig, ax = plt.subplots(figsize=(18, 12))
    ax.set_aspect('equal', adjustable='box')
    ax.margins(0.1)

    # --- Plotting Clusters and Components ---
    component_legend_handles = {} # To store handles for gene color legend

    for cluster_idx in range(NUM_CLUSTERS):
        cluster_name = f'Cluster_{cluster_idx + 1}'
        cluster_name_json = f'component_{cluster_idx + 1}' # JSON keys are lowercase
        cluster_pos = CLUSTER_POSITIONS[cluster_idx]
        
        # Plot Cluster Center
        ax.scatter(cluster_pos[0], cluster_pos[1], s=400, color=CLUSTER_COLOR, 
                   label='Cluster' if cluster_idx == 0 else "", zorder=2, alpha=0.7)
        ax.text(cluster_pos[0], cluster_pos[1] + 0.3, f'C{cluster_idx + 1}', 
                ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Plot Components around the cluster
        for comp_idx in range(NUM_COMPONENTS):
            component_name = f'Component_{comp_idx + 1}'
            component_name_json = f'component_{comp_idx + 1}' # JSON keys are lowercase
            
            # Get component value for this cluster
            try:
                nmf_value = normalized_nmf.loc[cluster_name, component_name]
            except KeyError:
                print(f"Warning: Could not find NMF value for {cluster_name}, {component_name}. Skipping bubble.")
                continue
                
            bubble_size = nmf_value * BUBBLE_SIZE_SCALE
            if bubble_size < 1: continue # Don't plot tiny bubbles
            
            # Determine component color based on top gene
            bubble_color = DEFAULT_COMPONENT_COLOR
            top_gene_for_comp = None
            if component_name_json in gene_data and gene_data[component_name_json]:
                # Get the #1 gene for this component (JSON assumed sorted)
                top_gene_for_comp = list(gene_data[component_name_json].keys())[0]
                if top_gene_for_comp in gene_color_map:
                    bubble_color = gene_color_map[top_gene_for_comp]
                    # Store legend handle if this is the first time seeing this gene color
                    if top_gene_for_comp not in component_legend_handles:
                         component_legend_handles[top_gene_for_comp] = Line2D([0], [0], marker='o', color='w',
                                                                             label=f'{top_gene_for_comp}',
                                                                             markerfacecolor=bubble_color, markersize=10)
                else:
                     bubble_color = DEFAULT_COMPONENT_COLOR # Use default if top gene not in overall top list
            
            # Calculate component position
            angle = 2 * np.pi * comp_idx / NUM_COMPONENTS
            comp_pos_x = cluster_pos[0] + COMPONENT_RADIUS * np.cos(angle)
            comp_pos_y = cluster_pos[1] + COMPONENT_RADIUS * np.sin(angle)
            
            # Plot Component Bubble
            ax.scatter(comp_pos_x, comp_pos_y, s=bubble_size, color=bubble_color, alpha=0.6, zorder=3)
            # Optional: Add component number text (can get cluttered)
            # ax.text(comp_pos_x, comp_pos_y, f'{comp_idx+1}', ha='center', va='center', fontsize=6, color='white' if np.mean(bubble_color[:3]) < 0.5 else 'black')

    # --- Create Legends ---
    # Cluster Legend
    cluster_legend = Line2D([0], [0], marker='o', color='w', label='Cluster Center', 
                             markerfacecolor=CLUSTER_COLOR, markersize=10)
    
    # Component Size Legend (Proxy Bubbles)
    size_legend_handles = []
    size_values = np.linspace(normalized_nmf.values.min(), normalized_nmf.values.max(), 4)
    for val in size_values:
         size = val * BUBBLE_SIZE_SCALE
         if size < 1: continue
         size_legend_handles.append(Line2D([0], [0], marker='o', color='w', 
                                           label=f'{val:.2f} (norm. value)',
                                           markerfacecolor=DEFAULT_COMPONENT_COLOR, 
                                           markersize=np.sqrt(size / np.pi)*0.5 )) # Approximate marker size
    
    # Combine Legends
    handles = [cluster_legend] + list(component_legend_handles.values()) 
    labels = [h.get_label() for h in handles]
    leg1 = ax.legend(handles=handles, labels=labels, title="Top Genes & Cluster", loc='upper left', bbox_to_anchor=(1.02, 1))
    ax.add_artist(leg1)
    
    # Add size legend separately
    if size_legend_handles:
        leg2 = ax.legend(handles=size_legend_handles, title="Component Size", loc='lower left', bbox_to_anchor=(1.02, 0)) 
        ax.add_artist(leg2)
        
    # --- Final Touches ---
    ax.set_title('Cluster-Component Relationship with Top Genes', fontsize=16)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for legends

    # Save the plot
    try:
        plt.savefig(OUTPUT_FILENAME, dpi=300, bbox_inches='tight')
        print(f"Chart saved successfully to 
{OUTPUT_FILENAME}
")
    except Exception as e:
        print(f"Error saving chart: {e}")
    plt.show()

# --- Main Execution ---
if __name__ == "__main__":
    nmf_data, gene_counts = load_data(NMF_FILENAME, GENE_COUNTS_FILENAME)
    if nmf_data is not None and gene_counts is not None:
        create_visualization(nmf_data, gene_counts)
# --------------------- 