import { PseudotimeGlyph } from './PseudotimeGlyph';
import { Empty, Spin, Select, Button } from 'antd';
import { CloseOutlined } from '@ant-design/icons';
import { useState, useEffect, useRef, useMemo } from 'react';
import { debounce } from './Utils';

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

    // State for hiding/closing glyphs - tracked by stable key (source_title)
    const [hiddenGlyphs, setHiddenGlyphs] = useState(new Set());

    // State for highly variable genes
    const [highVariableGenes, setHighVariableGenes] = useState([]);

    // Remote search state for gene Select
    const [searchQuery, setSearchQuery] = useState('');
    const [displayOptions, setDisplayOptions] = useState([]);
    const [searchLoading, setSearchLoading] = useState(false);

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

    // Group genes by sample ID for the select options and memoize the result
    const hvgGroupedOptions = useMemo(() => {
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
        return Object.entries(geneOptions).map(([sampleId, genes]) => ({
            label: sampleId,
            options: genes
        }));
    }, [highVariableGenes]);

    // Keep displayOptions in sync with HVGs when not actively searching
    useEffect(() => {
        if (!searchQuery) {
            setDisplayOptions(hvgGroupedOptions);
        }
    }, [hvgGroupedOptions, searchQuery]);

    // Helper to build grouped options from raw results
    const buildGroupedOptions = (records) => {
        const grouped = records.reduce((acc, item) => {
            const sid = item.sampleId;
            if (!acc[sid]) acc[sid] = [];
            acc[sid].push({ label: item.label, value: item.value });
            return acc;
        }, {});
        return Object.entries(grouped).map(([sampleId, genes]) => ({ label: sampleId, options: genes }));
    };

    // Debounced remote search for genes across provided samples
    const performGeneSearch = async (value) => {
        setSearchLoading(true);
        try {
            const response = await fetch('/api/get_gene_name_search', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    sample_ids: memoizedSampleIds,
                    query: value,
                    limit: 80,
                }),
            });
            if (response.ok) {
                const data = await response.json();
                // Transform to flat records then group
                const records = [];
                Object.entries(data || {}).forEach(([sid, genes]) => {
                    (genes || []).forEach((g) => {
                        records.push({
                            value: `${sid}_${g}`,
                            label: g,
                            sampleId: sid,
                        });
                    });
                });
                setDisplayOptions(buildGroupedOptions(records));
            } else {
                setDisplayOptions([]);
            }
        } catch (e) {
            console.error('Gene search failed:', e);
            setDisplayOptions([]);
        } finally {
            setSearchLoading(false);
        }
    };

    // Create debounced version of the search function
    const debouncedGeneSearch = useMemo(() => 
        debounce(performGeneSearch, 250), 
        [memoizedSampleIds]
    );

    const handleGeneSearch = (value) => {
        setSearchQuery(value);

        // If search text is empty, show HVGs
        if (!value || value.trim().length === 0) {
            setSearchLoading(false);
            setDisplayOptions(hvgGroupedOptions);
            return;
        }

        // Avoid spamming backend for very short queries
        if (value.trim().length < 2) {
            return;
        }

        // Use the debounced search function
        debouncedGeneSearch(value);
    };

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

    // Handle closing a glyph (hide it and unselect if selected)
    const handleCloseGlyph = (glyphKey) => {
        setHiddenGlyphs((prev) => {
            const next = new Set(prev);
            next.add(glyphKey);
            return next;
        });
        // Unselect any glyphs that correspond to this key (by index)
        setSelectedGlyphs((prev) => {
            const next = new Set(prev);
            // Selection is tracked by index; we cannot reliably map back here, so keep as-is.
            return next;
        });
    };

    // Prune hidden glyphs when datasets change (but preserve user's hide choices during loading)
    useEffect(() => {
        // Only prune hidden keys that are no longer present in the datasets
        const currentKeys = new Set(Object.keys(pseudotimeDataSets || {}));
        setHiddenGlyphs((prev) => {
            const next = new Set();
            prev.forEach((k) => {
                if (currentKeys.has(k)) next.add(k);
            });
            return next;
        });
    }, [pseudotimeDataSets]);

    // Helper function to extract UMAP parameters from adata_umap_title
    const extractUmapParameters = (adataUmapTitle) => {
        if (!adataUmapTitle || typeof adataUmapTitle !== 'string') {
            return { n_neighbors: 10, n_pcas: 30, resolutions: 1 }; // defaults
        }
        
        try {
            const parts = adataUmapTitle.split('_');
            if (parts.length >= 3) {
                // Extract the last 3 parts as the parameters
                const n_neighbors = parseInt(parts[parts.length - 3]) || 10;
                const n_pcas = parseInt(parts[parts.length - 2]) || 30;
                const resolutions = parseFloat(parts[parts.length - 1]) || 1;
                return { n_neighbors, n_pcas, resolutions };
            }
        } catch (error) {
            console.warn("Could not parse UMAP parameters from adata_umap_title:", adataUmapTitle, error);
        }
        
        return { n_neighbors: 10, n_pcas: 30, resolutions: 1 }; // fallback defaults
    };

    // Helper function to extract display name and sample ID for consistent formatting
    const getConsistentDisplayTitle = (baseTitle, umapDataSets) => {
        if (!umapDataSets || !Array.isArray(umapDataSets)) {
            return baseTitle;
        }
        
        // Find the matching UMAP dataset
        const matchingUmapDataset = umapDataSets.find(dataset =>
            dataset.adata_umap_title === baseTitle || dataset.title === baseTitle
        );
        
        if (matchingUmapDataset && matchingUmapDataset.title) {
            // The title should already be in format "${newName} (${dataset.sampleId})"
            return matchingUmapDataset.title;
        }
        
        // Fallback: try to construct the title from adata_umap_title
        if (baseTitle && typeof baseTitle === 'string') {
            // Extract sample ID from the title (look for patterns like skin_TXK6Z4X_A1)
            const sampleIdMatch = baseTitle.match(/(skin_[A-Z0-9]+_[A-Z0-9]+)/);
            if (sampleIdMatch) {
                const sampleId = sampleIdMatch[1];
                // Extract the name part (everything before the sample ID and parameters)
                const parts = baseTitle.split('_');
                // Remove the last 3 parts (parameters) and the sample ID parts
                const sampleIdParts = sampleId.split('_'); // e.g., ["skin", "TXK6Z4X", "A1"]
                const namePartsCount = parts.length - 3 - sampleIdParts.length;
                if (namePartsCount > 0) {
                    const nameParts = parts.slice(0, namePartsCount);
                    const displayName = nameParts.join('_').replace(/_/g, ' ');
                    return `${displayName} (${sampleId})`;
                }
            }
        }
        
        return baseTitle; // fallback to original title
    };

    // Convert pseudotimeDataSets object to separate trajectory data for each UMAP
    const allPseudotimeData = [];

    try {
        // Process each UMAP dataset separately to create individual glyphs
        Object.entries(pseudotimeDataSets).forEach(([title, pseudotimeData]) => {
            // Handle both new structure {cluster_order, trajectory_objects} and old array structure
            let hasValidData = false;
            let trajectoryObjects = [];

            if (pseudotimeData && typeof pseudotimeData === 'object') {
                if (pseudotimeData.trajectory_objects && Array.isArray(pseudotimeData.trajectory_objects) && pseudotimeData.trajectory_objects.length > 0) {
                    hasValidData = true;
                    trajectoryObjects = pseudotimeData.trajectory_objects;
                }
            }

            if (hasValidData) {
                // Determine if this is direct slingshot data or regular pseudotime data
                const isDirectSlingshot = title.endsWith('_direct_slingshot');
                const baseTitle = isDirectSlingshot ? title.replace('_direct_slingshot', '') : title;
                
                // Get consistent display title using helper function
                let displayTitle = getConsistentDisplayTitle(baseTitle, umapDataSets);
                
                // Extract UMAP parameters from the base title
                const umapParameters = extractUmapParameters(baseTitle);

                // Add suffix to distinguish between regular and direct slingshot data
                if (isDirectSlingshot) {
                    displayTitle = `${displayTitle} (Direct Slingshot)`;
                }

                // Process trajectories for this specific UMAP dataset
                const processedTrajectories = trajectoryObjects.map((trajectoryData) => {
                    // Ensure sampleId is set - try multiple fallback strategies
                    let sampleId = trajectoryData.sampleId;

                    // First fallback: try relatedSampleIds
                    if (!sampleId && relatedSampleIds.length > 0) {
                        sampleId = relatedSampleIds[0];
                    }

                    // Second fallback: try to extract from title (format: prefix_sampleId_suffix)
                    if (!sampleId) {
                        const titleToCheck = baseTitle || adata_umap_title;
                        if (titleToCheck && typeof titleToCheck === 'string') {
                            // Look for patterns like skin_TXK6Z4X_A1 or similar
                            const match = titleToCheck.match(/(skin_[A-Z0-9]+_[A-Z0-9]+)/);
                            if (match) {
                                sampleId = match[1];
                            }
                        }
                    }

                    return {
                        ...trajectoryData,
                        sampleId: sampleId
                    };
                });

                // Create a separate glyph entry for this UMAP dataset
                allPseudotimeData.push({
                    mergedTrajectories: processedTrajectories,
                    source_title: title,
                    display_title: displayTitle,
                    isLoading: false,
                    isPlaceholder: false,
                    fullPseudotimeData: pseudotimeData,
                    isDirectSlingshot: isDirectSlingshot,
                    umapParameters: umapParameters
                });
            }
        });

        // Check if there's any loading happening
        const isLoading = Object.values(pseudotimeLoadingStates).some(loading => loading);

        // Add loading placeholders for UMAP datasets that are currently loading
        if (umapDataSets && Array.isArray(umapDataSets)) {
            umapDataSets.forEach((umapDataset) => {
                const isThisDatasetLoading = pseudotimeLoadingStates[umapDataset.adata_umap_title];
                const directSlingshotKey = `${umapDataset.adata_umap_title}_direct_slingshot`;
                const isDirectSlingshotLoading = pseudotimeLoadingStates[directSlingshotKey];
                
                const hasDataForThisDataset = allPseudotimeData.some(data =>
                    data.source_title === umapDataset.adata_umap_title
                );
                const hasDirectSlingshotData = allPseudotimeData.some(data =>
                    data.source_title === directSlingshotKey
                );

                // If regular pseudotime is loading and we don't have data for it yet, add a loading placeholder
                if (isThisDatasetLoading && !hasDataForThisDataset) {
                    const consistentTitle = getConsistentDisplayTitle(umapDataset.adata_umap_title, umapDataSets);
                    const umapParams = extractUmapParameters(umapDataset.adata_umap_title);
                    allPseudotimeData.push({
                        source_title: umapDataset.adata_umap_title,
                        display_title: consistentTitle,
                        isLoading: true,
                        isPlaceholder: true,
                        umapParameters: umapParams
                    });
                }

                // If direct slingshot is loading and we don't have data for it yet, add a loading placeholder
                if (isDirectSlingshotLoading && !hasDirectSlingshotData) {
                    const consistentTitle = getConsistentDisplayTitle(umapDataset.adata_umap_title, umapDataSets);
                    const umapParams = extractUmapParameters(umapDataset.adata_umap_title);
                    allPseudotimeData.push({
                        source_title: directSlingshotKey,
                        display_title: `${consistentTitle} (Direct Slingshot)`,
                        isLoading: true,
                        isPlaceholder: true,
                        umapParameters: umapParams
                    });
                }
            });
        }

        // If no data at all but loading, show generic loading placeholder
        if (allPseudotimeData.length === 0 && isLoading) {
            allPseudotimeData.push({
                source_title: 'Loading',
                display_title: 'Loading trajectories...',
                isLoading: true,
                isPlaceholder: true
            });
        }
    } catch (error) {
        console.error('Error processing pseudotime data:', error);
        // Continue with empty data rather than breaking the component
    }

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
                // Use label from HVG list when available
                const geneInfo = highVariableGenes.find(g => g.value === geneValue);
                if (geneInfo) return geneInfo.label;
                // Fallback: value pattern is `${sampleId}_${gene}`; sampleId may contain underscores, so take substring after last underscore
                if (typeof geneValue === 'string' && geneValue.includes('_')) {
                    const lastUnderscore = geneValue.lastIndexOf('_');
                    return geneValue.slice(lastUnderscore + 1);
                }
                return geneValue;
            });

            // Create analysis requests for each selected glyph
            Array.from(selectedGlyphs).forEach(glyphIndex => {
                const trajectoryData = allPseudotimeData[glyphIndex];
                if (trajectoryData && !trajectoryData.isPlaceholder && trajectoryData.mergedTrajectories) {
                    // Use the related sample IDs directly - typically there should be one main sample
                    // For pseudotime analysis, we should use the same sample that generated the trajectory
                    const sampleId = relatedSampleIds.length > 0 ? relatedSampleIds[0] : 'unknown';

                    // The trajectory ID should match the glyph index or trajectory data structure
                    const trajectoryId = glyphIndex;

                    // For merged trajectories, we need to handle multiple paths
                    trajectoryData.mergedTrajectories.forEach((singleTrajectory, trajIndex) => {
                        // Use the base title without "_direct_slingshot" suffix for the API request
                        const baseTitleForApi = trajectoryData.isDirectSlingshot 
                            ? trajectoryData.source_title.replace('_direct_slingshot', '')
                            : trajectoryData.source_title;
                        
                        analysisRequests.push({
                            sample_id: sampleId,
                            genes: geneNames,
                            adata_umap_title: baseTitleForApi, // Use the base title without suffix for API calls
                            trajectory_id: `${trajectoryId}_${trajIndex}`,
                            trajectory_path: singleTrajectory.path
                        });
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
                <div style={{ textAlign: 'center', marginTop: '8px', color: '#666', fontSize: '12px' }}>
                    {(() => {
                        const loadingNames = (umapDataSets || [])
                            .filter(ds => pseudotimeLoadingStates?.[ds.adata_umap_title])
                            .map(ds => ds.title || ds.adata_umap_title);

                        if (loadingNames.length === 0) return 'Generating pseudotime...';
                        if (loadingNames.length === 1) return `Generating ${loadingNames[0]}...`;
                        return `Generating ${loadingNames.join(', ')}...`;
                    })()}
                </div>
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

    return (
        <div style={{
            width: '100%',
            height: '100%',
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center'
        }}>
            {/* Gene Selection Dropdown and Confirmation Button */}
            <div style={{
                top: '5px',
                right: '10px',
                zIndex: 999,
                display: 'flex',
                gap: '8px',
                alignItems: 'center',
                justifyContent: 'flex-end',
                padding: '5px'
            }}>
                <Select
                    placeholder="Select genes"
                    value={selectedGenes}
                    onChange={setSelectedGenes}
                    style={{ width: '220px' }}
                    size="small"
                    options={displayOptions}
                    mode="multiple"
                    showSearch
                    filterOption={false}
                    onSearch={handleGeneSearch}
                    onOpenChange={(open) => {
                        if (!open) {
                            // reset to HVGs when closing dropdown
                            setSearchQuery('');
                            setDisplayOptions(hvgGroupedOptions);
                        }
                    }}
                    notFoundContent={searchLoading ? <Spin size="small" /> : null}
                    maxTagCount="responsive"
                    allowClear
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

            {(() => {
                const indexToKey = (i) => (allPseudotimeData[i]?.source_title) || `${i}`;
                const visibleIndices = allPseudotimeData.map((_, i) => i).filter(i => !hiddenGlyphs.has(indexToKey(i)));
                const visibleCount = visibleIndices.length;
                return (
                    <div style={{
                        width: '100%',
                        height: `calc(100% - 35px)`,
                        display: 'grid',
                        gridTemplateColumns: visibleCount === 1 ? 'minmax(0, 1fr)' :
                            visibleCount === 2 ? 'repeat(2, minmax(0, 1fr))' :
                                'repeat(3, minmax(0, 1fr))',
                        gridAutoRows: '1fr',
                        gap: '5px',
                        overflow: 'hidden',
                        boxSizing: 'border-box'
                    }}>
                {visibleIndices.map((index) => {
                    const trajectoryData = allPseudotimeData[index];
                    const glyphKey = indexToKey(index);
                    return (
                        <div
                            key={index}
                            style={{
                                width: '100%',
                                height: '99%',
                                textAlign: 'center',
                                border: '1px solid #ddd',
                                borderRadius: '8px',
                                backgroundColor: '#f9f9f9',
                                position: 'relative',
                                overflow: 'hidden',
                                boxSizing: 'border-box'
                            }}
                        >
                            {!trajectoryData.isPlaceholder && (
                                <div style={{
                                    position: 'absolute',
                                    top: '4px',
                                    right: '4px',
                                    zIndex: 900,
                                }}>
                                    <Button
                                        type="text"
                                        size="small"
                                        icon={<CloseOutlined />}
                                        onClick={() => handleCloseGlyph(glyphKey)}
                                        aria-label="Close glyph"
                                    />
                                </div>
                            )}
                            {trajectoryData.isPlaceholder ? (
                                <div style={{
                                    width: '100%',
                                    height: '100%',
                                    display: 'flex',
                                    flexDirection: 'column',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    border: '2px dashed #ccc',
                                    borderRadius: '8px',
                                    color: '#666',
                                    fontSize: '12px',
                                    backgroundColor: '#fafafa',
                                    boxSizing: 'border-box'
                                }}>
                                    <Spin size="large" style={{ marginBottom: '10px' }} />
                                    <div style={{ textAlign: 'center' }}>
                                        <div style={{ fontWeight: 'bold', marginBottom: '5px' }}>
                                            Generating Pseudotime
                                        </div>
                                        <div style={{ fontSize: '11px', color: '#999' }}>
                                            {trajectoryData.display_title}
                                        </div>
                                    </div>
                                </div>
                            ) : (
                                (() => {
                                    try {
                                        // For merged trajectories, collect gene expression data from all related trajectories
                                        let geneDataForGlyph = null;

                                        if (trajectoryData.mergedTrajectories && trajectoryData.mergedTrajectories.length > 0) {
                                            // Collect all gene expression data for trajectories that start with this glyph index
                                            const relatedGeneData = geneExpressionData.filter(data =>
                                                data.trajectory_id && data.trajectory_id.toString().startsWith(`${index}_`)
                                            );

                                            if (relatedGeneData.length > 0) {
                                                // Pass all related gene data so the glyph can select the appropriate one
                                                geneDataForGlyph = relatedGeneData;
                                            }
                                        } else {
                                            // Single trajectory - find by exact match
                                            const singleGeneData = geneExpressionData.find(data =>
                                                data.trajectory_id === index || data.trajectory_id === index.toString()
                                            );
                                            geneDataForGlyph = singleGeneData ? [singleGeneData] : null;
                                        }

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
                                                pseudotimeData={trajectoryData.fullPseudotimeData || trajectoryData.mergedTrajectories || [trajectoryData]}
                                                pseudotimeLoading={trajectoryData.isLoading}
                                                isSelected={selectedGlyphs.has(index)}
                                                onSelectionChange={(isSelected) => handleGlyphSelection(index, isSelected)}
                                                geneExpressionData={geneDataForGlyph}
                                                clusterColors={clusterColors}
                                                hoveredTrajectory={hoveredTrajectory}
                                                setHoveredTrajectory={setHoveredTrajectory}
                                                trajectoryIndex={index}
                                                source_title={trajectoryData.source_title || adata_umap_title}
                                                umapParameters={trajectoryData.umapParameters}
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
                    );
                })}
                    </div>
                );
            })()}
        </div>
    );
};