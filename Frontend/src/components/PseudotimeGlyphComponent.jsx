import PseudotimeGlyph from './PseudotimeGlyph';
import { Empty, Spin, Select } from 'antd';
import { useState, useEffect } from 'react';

export const PseudotimeGlyphComponent = ({
    adata_umap_title,
    pseudotimeDataSets,
    pseudotimeLoadingStates,
    relatedSampleIds,
}) => {
    // State for tracking selected glyphs
    const [selectedGlyphs, setSelectedGlyphs] = useState(new Set());
    
    // State for selected gene
    const [selectedGene, setSelectedGene] = useState(null);
    
    // State for highly variable genes
    const [highVariableGenes, setHighVariableGenes] = useState([]);
    const [genesLoading, setGenesLoading] = useState(false);
    
    // Load highly variable genes on component mount
    useEffect(() => {
        const fetchHighVariableGenes = async () => {
            setGenesLoading(true);
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
            } finally {
                setGenesLoading(false);
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
    
    // Debug selected gene changes
    useEffect(() => {
        if (selectedGene) {
            console.log('Selected gene:', selectedGene);
        }
    }, [selectedGene]);
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
            {/* Gene Selection Dropdown in upper right corner */}
            <div style={{
                position: 'absolute',
                top: '10px',
                right: '10px',
                zIndex: 1000,
                width: '200px'
            }}>
                <Select
                    placeholder={genesLoading ? "Loading genes..." : "Select a gene"}
                    value={selectedGene}
                    onChange={setSelectedGene}
                    style={{ width: '100%' }}
                    size="small"
                    options={selectOptions}
                    loading={genesLoading}
                    disabled={genesLoading}
                    showSearch
                    filterOption={(input, option) =>
                        (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
                    }
                />
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