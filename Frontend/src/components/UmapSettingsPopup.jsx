import React, { useState, useEffect } from 'react';
import '../styles/UmapSettingsPopup.css';

export const UmapSettingsPopup = ({
  visible,
  setVisible,
  position,
  currentSettings,
  onUpdateSettings,
  onPseudotimeAnalysis,
  onLoadingStart,
  sampleId,
  cellIds,
  adata_umap_title
}) => {
  const [settings, setSettings] = useState({
    n_neighbors: 10,
    n_pcas: 30,
    resolutions: 1
  });
  const [loading, setLoading] = useState(false);
  const [pseudotimeLoading, setPseudotimeLoading] = useState(false);

  // Initialize settings from current adata_umap_title when popup opens
  useEffect(() => {
    if (visible && adata_umap_title) {
      try {
        const parts = adata_umap_title.split('_');
        if (parts.length >= 5) {
          // Extract the last 3 parts as the parameters
          const n_neighbors = parseInt(parts[parts.length - 3]) || 10;
          const n_pcas = parseInt(parts[parts.length - 2]) || 30;
          const resolutions = parseFloat(parts[parts.length - 1]) || 1;
          
          setSettings({
            n_neighbors,
            n_pcas,
            resolutions
          });
        }
      } catch (error) {
        console.warn("Could not parse parameters from adata_umap_title, using defaults", error);
      }
    }
  }, [visible, adata_umap_title]);

  const handleInputChange = (field, value) => {
    setSettings(prev => ({
      ...prev,
      [field]: value
    }));
  };

  const handleUpdateUMAP = async () => {
    if (!sampleId || !cellIds || cellIds.length === 0) {
      alert('No cells selected for UMAP analysis');
      return;
    }

    // Validate settings
    if (settings.n_neighbors < 1 || settings.n_neighbors > 50) {
      alert('Number of neighbors must be between 1 and 50');
      return;
    }
    if (settings.n_pcas < 1 || settings.n_pcas > 100) {
      alert('Number of PCs must be between 1 and 100');
      return;
    }
    if (settings.resolutions < 0.1 || settings.resolutions > 5) {
      alert('Resolution must be between 0.1 and 5.0');
      return;
    }

    setLoading(true);
    
    // Start the loading animation in the parent component
    if (onLoadingStart) {
      onLoadingStart();
    }
    
    try {
      // Generate new adata_umap_title with updated parameters
      const parts = adata_umap_title.split('_');
      const baseName = parts.slice(0, -3).join('_'); // Remove the last 3 parts (parameters)
      const newAdataUmapTitle = `${baseName}_${settings.n_neighbors}_${settings.n_pcas}_${settings.resolutions}`;

      // Call the API to get updated UMAP data
      const response = await fetch('/api/get_umap_data', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sample_id: sampleId,
          cell_ids: cellIds,
          n_neighbors: settings.n_neighbors,
          n_pcas: settings.n_pcas,
          resolutions: settings.resolutions,
          adata_umap_title: newAdataUmapTitle
        })
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const data = await response.json();
      
      if (data.error) {
        throw new Error(data.error);
      }

      // Call the parent callback with new data and title
      onUpdateSettings(data, newAdataUmapTitle, settings);
      
      // Close the popup
      setVisible(false);
      
    } catch (error) {
      console.error('Error updating UMAP:', error);
      alert(`Failed to update UMAP: ${error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handlePseudotimeAnalysis = async () => {
    if (!sampleId || !cellIds || cellIds.length === 0) {
      alert('No cells selected for pseudotime analysis');
      return;
    }

    setPseudotimeLoading(true);
    
    try {
      await onPseudotimeAnalysis(sampleId, cellIds);
    } catch (error) {
      console.error('Error in pseudotime analysis:', error);
      alert(`Failed to perform pseudotime analysis: ${error.message}`);
    } finally {
      setPseudotimeLoading(false);
    }
  };

  if (!visible) return null;

  return (
    <div 
      className="umap-settings-popup-overlay"
      onClick={(e) => {
        if (e.target === e.currentTarget) {
          setVisible(false);
        }
      }}
    >
      <div 
        className="umap-settings-popup"
        style={{
          left: Math.min(position.x, window.innerWidth - 420),
          top: Math.min(position.y, window.innerHeight - 400)
        }}
      >
        <div className="umap-settings-header">
          <h3>UMAP Settings</h3>
          <div className="umap-settings-current">
            <small>Current: {settings.n_neighbors} neighbors, {settings.n_pcas} PCs, {settings.resolutions} resolution</small>
          </div>
          <button 
            className="umap-settings-close"
            onClick={() => setVisible(false)}
          >
            ×
          </button>
        </div>
        
        <div className="umap-settings-content">
          <div className="umap-setting-group">
            <label htmlFor="n_neighbors">Number of Neighbors:</label>
            <input
              id="n_neighbors"
              type="number"
              min="1"
              max="50"
              value={settings.n_neighbors}
              onChange={(e) => handleInputChange('n_neighbors', parseInt(e.target.value) || 10)}
            />
            <small>Controls the number of neighbors for neighbor graph construction (1-50)</small>
          </div>

          <div className="umap-setting-group">
            <label htmlFor="n_pcas">Number of PCs:</label>
            <input
              id="n_pcas"
              type="number"
              min="1"
              max="100"
              value={settings.n_pcas}
              onChange={(e) => handleInputChange('n_pcas', parseInt(e.target.value) || 30)}
            />
            <small>Number of principal components to use (1-100)</small>
          </div>

          <div className="umap-setting-group">
            <label htmlFor="resolutions">Resolution:</label>
            <input
              id="resolutions"
              type="number"
              min="0.1"
              max="5"
              step="0.1"
              value={settings.resolutions}
              onChange={(e) => handleInputChange('resolutions', parseFloat(e.target.value) || 1)}
            />
            <small>Resolution parameter for Leiden clustering (0.1-5.0)</small>
          </div>

          <div className="umap-settings-actions">
            <button
              className="umap-settings-btn umap-settings-btn-primary"
              onClick={handleUpdateUMAP}
              disabled={loading}
            >
              {loading ? 'Updating...' : 'Update UMAP'}
            </button>
            
            <button
              className="umap-settings-btn umap-settings-btn-secondary"
              onClick={handlePseudotimeAnalysis}
              disabled={pseudotimeLoading}
            >
              {pseudotimeLoading ? 'Analyzing...' : 'Pseudotime Analysis'}
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};
