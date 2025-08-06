import PseudotimeGlyph from './PseudotimeGlyph';
import { Empty, Spin, Select, Button } from 'antd';
import { useState, useEffect } from 'react';

export const PseudotimeGlyphComponent = ({
    adata_umap_title,
    pseudotimeDataSets,
    pseudotimeLoadingStates,
    relatedSampleIds,
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
    
    // Load highly variable genes on component mount
    useEffect(() => {
        const fetchHighVariableGenes = async () => {
            try {
                if (relatedSampleIds.length > 0) {
                    const response = await fetch("/api/get_highly_variable_genes", {
                        method: "POST",
                        headers: {
                            "Content-Type": "application/json",
                        },
                        body: JSON.stringify({
                            sample_ids: relatedSampleIds,
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
                    } else {
                        console.error("Failed to fetch highly variable genes");
                        // Fallback to mock data
                        setHighVariableGenes([]);
                    }
                } else {
                    setHighVariableGenes([]);
                }
            } catch (error) {
                console.error("Error fetching highly variable genes:", error);
                setHighVariableGenes([]);
            }
        };
        
        if (Object.keys(pseudotimeDataSets).length > 0) {
            fetchHighVariableGenes();
        }
    }, [pseudotimeDataSets, relatedSampleIds]);

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
        console.log('Selected glyphs:', Array.from(newSelected));
    };
    
    // Debug selected genes changes
    useEffect(() => {
        if (selectedGenes.length > 0) {
            console.log('Selected genes:', selectedGenes);
        }
    }, [selectedGenes]);
    
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
                        analysisResults.push({
                            ...request,
                            gene_expression_data: data
                        });
                    } else {
                        console.error(`Failed to get gene expression data for sample ${request.sample_id}`);
                    }
                } catch (error) {
                    console.error(`Error fetching gene expression data for sample ${request.sample_id}:`, error);
                }
            }
            
            setGeneExpressionData(analysisResults);

            // Show success message with the results
            if (analysisResults.length > 0) {
                console.log(`Successfully analyzed ${analysisResults.length} trajectory(ies) for genes ${geneNames.join(', ')}`);
                // The results are now stored in geneExpressionData state and logged to console
                // Each result has the structure: { sample_id, genes, adata_umap_title, trajectory_id, trajectory_path, gene_expression_data }
                // where gene_expression_data contains: [{ gene: "SOX2", timePoints: [0.0, 0.4, 0.7, 1.0], expressions: [0.8, 0.6, 0.3, 0.2] }, ...]
            }
            
        } catch (error) {
            console.error('Error during gene expression analysis:', error);
        } finally {
            setIsAnalyzing(false);
        }
    };
    
    // Convert pseudotimeDataSets object to array of all trajectory data
    const allPseudotimeData = [];
    Object.entries(pseudotimeDataSets).forEach(([title, dataArray]) => {
        if (Array.isArray(dataArray)) {
            dataArray.forEach((trajectoryData, index) => {
                allPseudotimeData.push({
                    ...trajectoryData,
                    source_title: title,
                    display_title: `${title} - Trajectory ${index + 1}`,
                    isLoading: pseudotimeLoadingStates[title] || false
                });
            });
        }
    });

    // Add loading placeholders for datasets that are currently being loaded
    Object.entries(pseudotimeLoadingStates).forEach(([title, isLoading]) => {
        if (isLoading && !pseudotimeDataSets[title]) {
            allPseudotimeData.push({
                source_title: title,
                display_title: `${title} - Loading...`,
                isLoading: true,
                isPlaceholder: true
            });
        }
    });

    // Check if there's any global loading happening and no data exists yet
    const anyLoading = Object.values(pseudotimeLoadingStates).some(loading => loading);
    const hasNoData = Object.keys(pseudotimeDataSets).length === 0;

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
                            <PseudotimeGlyph
                                adata_umap_title={trajectoryData.display_title || `${adata_umap_title} - Trajectory ${index + 1}`}
                                pseudotimeData={[trajectoryData]}
                                pseudotimeLoading={trajectoryData.isLoading}
                                isSelected={selectedGlyphs.has(index)}
                                onSelectionChange={(isSelected) => handleGlyphSelection(index, isSelected)}
                            />
                        )}
                    </div>
                ))}
            </div>
        </div>
    );
};