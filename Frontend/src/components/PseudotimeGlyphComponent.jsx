import { PseudotimeGlyph } from './PseudotimeGlyph';
import { Empty, Spin, Select, Button } from 'antd';
import { useState, useEffect, useRef, useMemo } from 'react';

export const PseudotimeGlyphComponent = ({
    adata_umap_title,
    pseudotimeDataSets,
    pseudotimeLoadingStates,
    relatedSampleIds,
    clusterColorMappings,
    hoveredTrajectory,
    setHoveredTrajectory,
    umapDataSets,
}) => {
    // State for tracking selected glyphs
    const [selectedGlyphs, setSelectedGlyphs] = useState(new Set());
    
    // State for selected genes (multiple selection)
    const [selectedGenes, setSelectedGenes] = useState([]);
    
    // State for highly variable genes
    const [highVariableGenes, setHighVariableGenes] = useState([]);
    
    // State for gene expression analysis
    const [geneExpressionData, setGeneExpressionData] = useState([]);
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [componentError, setComponentError] = useState(null);
    
    // Ref to track if genes have been fetched for current sample IDs
    const fetchedSampleIdsRef = useRef(new Set());
    
    // Memoize the sample IDs to prevent unnecessary re-fetching
    const memoizedSampleIds = useMemo(() => {
        return [...relatedSampleIds].sort(); // Sort for consistent comparison
    }, [relatedSampleIds]);
    
    // Load highly variable genes only when sample IDs actually change
    useEffect(() => {
        const fetchHighVariableGenes = async () => {
            // Create a key from current sample IDs to check if we need to fetch
            const currentSampleKey = memoizedSampleIds.join(',');
            
            // Only fetch if sample IDs are different from what we've already fetched
            if (memoizedSampleIds.length === 0) {
                setHighVariableGenes([]);
                fetchedSampleIdsRef.current.clear();
                return;
            }
            
            // Check if we've already fetched for these exact sample IDs
            if (fetchedSampleIdsRef.current.has(currentSampleKey)) {
                return; // Skip fetching, we already have the data
            }
            
            try {
                const response = await fetch("/api/get_highly_variable_genes", {
                    method: "POST",
                    headers: {
                        "Content-Type": "application/json",
                    },
                    body: JSON.stringify({
                        sample_ids: memoizedSampleIds,
                        top_n: 20
                    }),
                });
                
                if (response.ok) {
                    const data = await response.json();
                    
                    // Convert the response to the format expected by the select component
                    const geneList = [];
                    Object.entries(data).forEach(([sampleId, genes]) => {
                        genes.forEach(gene => {
                            geneList.push({
                                value: `${sampleId}_${gene}`,
                                label: gene,
                                sampleId: sampleId
                            });
                        });
                    });
                    
                    setHighVariableGenes(geneList);
                    
                    // Mark these sample IDs as fetched
                    fetchedSampleIdsRef.current.add(currentSampleKey);
                } else {
                    console.error("Failed to fetch highly variable genes");
                    setHighVariableGenes([]);
                }
            } catch (error) {
                console.error("Error fetching highly variable genes:", error);
                setHighVariableGenes([]);
            }
        };
        
        fetchHighVariableGenes();
    }, [memoizedSampleIds]);

    // Group genes by sample ID for the select options
    const geneOptions = highVariableGenes.reduce((groups, gene) => {
        if (!groups[gene.sampleId]) {
            groups[gene.sampleId] = [];
        }
        groups[gene.sampleId].push({
            label: gene.label,
            value: gene.value
        });
        return groups;
    }, {});
    
    // Convert to antd Select option format with groups
    const selectOptions = Object.entries(geneOptions).map(([sampleId, genes]) => ({
        label: sampleId,
        options: genes
    }));
    
    // Handle glyph selection
    const handleGlyphSelection = (glyphIndex, isSelected) => {
        const newSelected = new Set(selectedGlyphs);
        if (isSelected) {
            newSelected.add(glyphIndex);
        } else {
            newSelected.delete(glyphIndex);
        }
        setSelectedGlyphs(newSelected);
    };
    
    // Handle confirmation button click
    const handleAnalyzeGeneExpression = async () => {
        if (selectedGlyphs.size === 0) {
            console.warn('No glyphs selected');
            return;
        }
        
        if (selectedGenes.length === 0) {
            console.warn('No genes selected');
            return;
        }
        
        setIsAnalyzing(true);
        
        try {
            const analysisRequests = [];
            
            // Extract clean gene names from the selected gene values (format: sampleId_geneName)
            const geneNames = selectedGenes.map(geneValue => {
                // The gene name is the label part, which should be just the gene name without sample ID
                const geneInfo = highVariableGenes.find(g => g.value === geneValue);
                return geneInfo ? geneInfo.label : geneValue; // Use label (clean gene name) or fallback to value
            });

            // Create analysis requests for each selected glyph
            Array.from(selectedGlyphs).forEach(glyphIndex => {
                const trajectoryData = allPseudotimeData[glyphIndex];
                if (trajectoryData && !trajectoryData.isPlaceholder && trajectoryData.path) {
                    // Use the related sample IDs directly - typically there should be one main sample
                    // For pseudotime analysis, we should use the same sample that generated the trajectory
                    const sampleId = relatedSampleIds.length > 0 ? relatedSampleIds[0] : 'unknown';
                    
                    // The trajectory ID should match the glyph index or trajectory data structure
                    const trajectoryId = glyphIndex;

                    analysisRequests.push({
                        sample_id: sampleId,
                        genes: geneNames,
                        adata_umap_title: adata_umap_title,
                        trajectory_id: trajectoryId,
                        trajectory_path: trajectoryData.path
                    });
                }
            });

            // Make API calls to get_trajectory_gene_expression for each request
            const analysisResults = [];
            
            for (const request of analysisRequests) {
                try {
                    const response = await fetch("/api/get_trajectory_gene_expression", {
                        method: "POST",
                        headers: {
                            "Content-Type": "application/json",
                        },
                        body: JSON.stringify({
                            sample_id: request.sample_id,
                            adata_umap_title: request.adata_umap_title,
                            gene_names: request.genes,
                            trajectory_path: request.trajectory_path
                        }),
                    });
                    
                    if (response.ok) {
                        const data = await response.json();
                        // Validate and sanitize the response data to prevent memory issues
                        if (data && typeof data === 'object') {
                            analysisResults.push({
                                ...request,
                                gene_expression_data: data
                            });
                        } else {
                            console.warn(`Invalid gene expression data format for sample ${request.sample_id}`);
                        }
                    } else {
                        const errorText = await response.text();
                        console.error(`Failed to get gene expression data for sample ${request.sample_id}: ${response.status} ${errorText}`);
                    }
                } catch (error) {
                    console.error(`Error fetching gene expression data for sample ${request.sample_id}:`, error);
                }
            }
            
            // Only update state if we have successful results
            if (analysisResults.length > 0) {
                setGeneExpressionData(analysisResults);
            } else {
                console.warn('No gene expression data was successfully retrieved');
            }
        } catch (error) {
            console.error('Error during gene expression analysis:', error);
            // Don't let errors break the entire component - just log them and continue
        } finally {
            setIsAnalyzing(false);
        }
    };
    
    // Convert pseudotimeDataSets object to array of all trajectory data
    const allPseudotimeData = [];
    
    try {
        Object.entries(pseudotimeDataSets).forEach(([title, dataArray]) => {
            if (Array.isArray(dataArray)) {
                dataArray.forEach((trajectoryData, index) => {
                    // Find the corresponding UMAP dataset title
                    let displayTitle = title;
                    
                    // Look for matching UMAP dataset by adata_umap_title
                    if (umapDataSets && Array.isArray(umapDataSets)) {
                        const matchingUmapDataset = umapDataSets.find(dataset => 
                            dataset.adata_umap_title === title || dataset.title === title
                        );
                        
                        if (matchingUmapDataset) {
                            displayTitle = matchingUmapDataset.title;
                        }
                    }
                    
                    // Ensure sampleId is set - try multiple fallback strategies
                    let sampleId = trajectoryData.sampleId;
                    
                    // First fallback: try relatedSampleIds
                    if (!sampleId && relatedSampleIds.length > 0) {
                        sampleId = relatedSampleIds[0];
                    }
                    
                    // Second fallback: try to extract from title (format: prefix_sampleId_suffix)
                    if (!sampleId) {
                        const titleToCheck = title || adata_umap_title;
                        if (titleToCheck && typeof titleToCheck === 'string') {
                            // Look for patterns like skin_TXK6Z4X_A1 or similar
                            const match = titleToCheck.match(/(skin_[A-Z0-9]+_[A-Z0-9]+)/);
                            if (match) {
                                sampleId = match[1];
                            }
                        }
                    }
                    
                    allPseudotimeData.push({
                        ...trajectoryData,
                        source_title: title,
                        display_title: `${displayTitle} - ${index + 1}`,
                        isLoading: pseudotimeLoadingStates[title] || false
                    });
                });
            }
        });

        // Add loading placeholders for datasets that are currently being loaded
        Object.entries(pseudotimeLoadingStates).forEach(([title, isLoading]) => {
            if (isLoading && !pseudotimeDataSets[title]) {
                // Find the corresponding UMAP dataset title for loading placeholder
                let displayTitle = title;
                if (umapDataSets && Array.isArray(umapDataSets)) {
                    const matchingUmapDataset = umapDataSets.find(dataset => 
                        dataset.adata_umap_title === title || dataset.title === title
                    );
                    
                    if (matchingUmapDataset) {
                        displayTitle = matchingUmapDataset.title;
                    }
                }
                
                allPseudotimeData.push({
                    source_title: title,
                    display_title: `${displayTitle} - Loading...`,
                    isLoading: true,
                    isPlaceholder: true
                });
            }
        });
    } catch (error) {
        console.error('Error processing pseudotime data:', error);
        // Continue with empty data rather than breaking the component
    }

    // Check if there's any global loading happening and no data exists yet
    const anyLoading = Object.values(pseudotimeLoadingStates).some(loading => loading);
    const hasNoData = Object.keys(pseudotimeDataSets).length === 0;

    // If there's a component error, show error state
    if (componentError) {
        return (
            <div style={{ 
                textAlign: 'center', 
                marginTop: '15%', 
                color: '#ff4d4f',
                padding: '20px'
            }}>
                <div style={{ fontSize: '16px', marginBottom: '10px' }}>Component Error</div>
                <div style={{ fontSize: '12px', marginBottom: '10px' }}>{componentError.message}</div>
                <Button 
                    size="small" 
                    onClick={() => setComponentError(null)}
                    type="primary"
                >
                    Retry
                </Button>
            </div>
        );
    }

    // If allPseudotimeData is not an array or is empty, show loading or empty state
    if (anyLoading && hasNoData) {
        return (
            <div style={{ textAlign: 'center', marginTop: '15%' }}>
                <Spin size="large" />
            </div>  
        );
    }

    if (!allPseudotimeData || allPseudotimeData.length === 0) {
        return (
            <Empty
                description="No pseudotime data available"
                image={Empty.PRESENTED_IMAGE_SIMPLE}
            />
        );
    }

    // Calculate responsive dimensions
    const numGlyphs = allPseudotimeData.length;
    const maxPerRow = 3;
    const glyphsPerRow = Math.min(numGlyphs, maxPerRow);
    const numRows = Math.ceil(numGlyphs / maxPerRow);

    const gapSize = 10; // gap in px

    // Wrap the entire render in error handling to prevent white screen
    try {
        return (
            <div style={{
                padding: '10px',
                width: '100%',
                height: '100%',
                boxSizing: 'border-box',
                position: 'relative'
            }}>
                {/* Gene Selection Dropdown and Confirmation Button */}
                <div style={{
                    position: 'absolute',
                    top: '10px',
                    right: '10px',
                    zIndex: 1000,
                    display: 'flex',
                    gap: '8px',
                    alignItems: 'center'
                }}>
                    <Select
                        placeholder="Select genes"
                        value={selectedGenes}
                        onChange={setSelectedGenes}
                        style={{ width: '200px' }}
                        size="small"
                        options={selectOptions}
                        mode="multiple"
                        showSearch
                        filterOption={(input, option) =>
                            (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
                        }
                        maxTagCount="responsive"
                    />
                    <Button
                        onClick={handleAnalyzeGeneExpression}
                        disabled={selectedGlyphs.size === 0 || selectedGenes.length === 0 || isAnalyzing}
                        loading={isAnalyzing}
                        type="primary"
                        size="small"
                    >
                        Analyze
                    </Button>
                </div>

                <div style={{
                    display: 'grid',
                    gridTemplateColumns: `repeat(${glyphsPerRow}, 1fr)`,
                    gridTemplateRows: `repeat(${numRows}, 1fr)`,
                    gap: `${gapSize}px`,
                    width: '100%',
                    height: `calc(100% - 30px)`,
                    justifyItems: 'center',
                    alignItems: 'center',
                    marginTop: '35px'
                }}>
                    {allPseudotimeData.map((trajectoryData, index) => (
                        <div
                            key={index}
                            style={{
                                width: '100%',
                                height: '100%',
                                minWidth: '200px',
                                minHeight: '200px',
                                textAlign: 'center'
                            }}
                        >
                            {trajectoryData.isPlaceholder ? (
                                <div style={{
                                    width: '100%',
                                    height: '100%',
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    border: '2px dashed #ccc',
                                    borderRadius: '8px',
                                    color: '#666',
                                    fontSize: '8px'
                                }}>
                                    Loading {trajectoryData.source_title}...
                                </div>
                            ) : (
                                (() => {
                                    try {
                                        const geneDataForGlyph = geneExpressionData.find(data => data.trajectory_id === index)?.gene_expression_data || null;
                                        
                                        // Get cluster color mapping for this trajectory data
                                        // Try different key formats to find the correct mapping
                                        let clusterColors = null;
                                        if (clusterColorMappings) {
                                            // Try the source_title first
                                            clusterColors = clusterColorMappings[trajectoryData.source_title]?.clusters;
                                            
                                            // If not found, try the adata_umap_title
                                            if (!clusterColors && adata_umap_title) {
                                                clusterColors = clusterColorMappings[adata_umap_title]?.clusters;
                                            }
                                            
                                            // If still not found, try with sample_id prefix
                                            if (!clusterColors && trajectoryData.sampleId) {
                                                const keyWithSample = `${trajectoryData.sampleId}_${trajectoryData.source_title || adata_umap_title}`;
                                                clusterColors = clusterColorMappings[keyWithSample]?.clusters;
                                            }
                                        }

                                        return (
                                            <PseudotimeGlyph
                                                adata_umap_title={trajectoryData.display_title}
                                                pseudotimeData={[trajectoryData]}
                                                pseudotimeLoading={trajectoryData.isLoading}
                                                isSelected={selectedGlyphs.has(index)}
                                                onSelectionChange={(isSelected) => handleGlyphSelection(index, isSelected)}
                                                geneExpressionData={geneDataForGlyph}
                                                clusterColors={clusterColors}
                                                hoveredTrajectory={hoveredTrajectory}
                                                setHoveredTrajectory={setHoveredTrajectory}
                                                trajectoryIndex={index}
                                                source_title={trajectoryData.source_title || adata_umap_title}
                                            />
                                        );
                                    } catch (error) {
                                        console.error(`Error rendering glyph ${index}:`, error);
                                        return (
                                            <div style={{
                                                width: '100%',
                                                height: '100%',
                                                display: 'flex',
                                                alignItems: 'center',
                                                justifyContent: 'center',
                                                border: '2px solid #ff4d4f',
                                                borderRadius: '8px',
                                                color: '#ff4d4f',
                                                fontSize: '12px'
                                            }}>
                                                Error rendering glyph
                                            </div>
                                        );
                                    }
                                })()
                            )}
                        </div>
                    ))}
                </div>
            </div>
        );
    } catch (error) {
        console.error('Critical component render error:', error);
        return (
            <div style={{ 
                textAlign: 'center', 
                marginTop: '15%', 
                color: '#ff4d4f',
                padding: '20px'
            }}>
                <div style={{ fontSize: '16px', marginBottom: '10px' }}>Component Render Error</div>
                <div style={{ fontSize: '12px', marginBottom: '10px' }}>{error.message}</div>
                <Button 
                    size="small" 
                    onClick={() => window.location.reload()}
                    type="primary"
                >
                    Reload Page
                </Button>
            </div>
        );
    }
};