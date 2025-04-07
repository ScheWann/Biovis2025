import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as path_effects

# Set style
plt.style.use('seaborn')
sns.set_context("notebook", font_scale=1.2)

# Read the JSON file
with open('output/umap_data.json', 'r') as f:
    data = json.load(f)

# Extract data
umap_coords = np.array(data['umap_coordinates'])
cell_types = np.array(data['cell_types'])
cell_clusters = np.array(data['cell_clusters'])
pseudotime = np.array(data['pseudotime'])
gene_locations = data['gene_locations']

# Create figure with subplots
fig = plt.figure(figsize=(20, 6))

# Plot 1: Cell Types
ax1 = fig.add_subplot(131)
unique_types = np.unique(cell_types)
colors = sns.color_palette("husl", n_colors=len(unique_types))
for i, cell_type in enumerate(unique_types):
    mask = cell_types == cell_type
    ax1.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                c=[colors[i]], label=cell_type, alpha=0.6, s=10)
ax1.set_title('Cell Types Distribution', pad=20)
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.set_xlabel('UMAP1')
ax1.set_ylabel('UMAP2')

# Plot 2: Clusters
ax2 = fig.add_subplot(132)
unique_clusters = np.unique(cell_clusters)
colors = sns.color_palette("husl", n_colors=len(unique_clusters))
for i, cluster in enumerate(unique_clusters):
    mask = cell_clusters == cluster
    ax2.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                c=[colors[i]], label=f'Cluster {cluster}', alpha=0.6, s=10)
ax2.set_title('Cell Clusters Distribution', pad=20)
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.set_xlabel('UMAP1')
ax2.set_ylabel('UMAP2')

# Plot 3: Pseudotime
ax3 = fig.add_subplot(133)
scatter = ax3.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                     c=pseudotime, cmap='viridis', alpha=0.6, s=10)
plt.colorbar(scatter, ax=ax3, label='Pseudotime')
ax3.set_title('Pseudotime Distribution', pad=20)
ax3.set_xlabel('UMAP1')
ax3.set_ylabel('UMAP2')

# Adjust layout and save
plt.tight_layout()
plt.savefig('output/deaplog_visualization.png', dpi=300, bbox_inches='tight')
plt.close()

# Create gene locations plot
plt.figure(figsize=(12, 8))
scatter = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                     c=pseudotime, cmap='viridis', alpha=0.1, s=10)
plt.colorbar(scatter, label='Pseudotime')

# Plot gene locations
for gene, info in gene_locations.items():
    plt.scatter(info['x'], info['y'], c='red', s=100, alpha=0.6)
    text = plt.text(info['x'], info['y'], gene, fontsize=8)
    text.set_path_effects([path_effects.withStroke(linewidth=3, foreground='white')])

plt.title('Gene Locations on UMAP')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.tight_layout()
plt.savefig('output/gene_locations.png', dpi=300, bbox_inches='tight')
plt.close()

print("Visualization completed. Check the output directory for the generated plots.") 