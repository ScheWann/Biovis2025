import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox } from "antd";
import { CloseOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, PathLayer } from '@deck.gl/layers';
import { booleanPointInPolygon } from '@turf/turf';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { debounce } from 'lodash';
import { convertHSLtoRGB, convertHEXToRGB, hashStringToHue } from './ColorUtils';
import "../styles/MultiSampleViewer.css";

// Generate random color for new regions
const generateRandomColor = () => {
    const randomHue = Math.floor(Math.random() * 360);
    return convertHSLtoRGB(randomHue, 100, 50);
};

export const MultiSampleViewerNew = ({
    setSampleDataLoading,
    selectedSamples,
    coordinatesData,
    interestedRegions,
    setInterestedRegions,
    analyzedRegion,
    setAnalyzedRegion,
}) => {
    const containerRef = useRef(null);
    const layersRef = useRef([]);
    const [mainViewState, setMainViewState] = useState(null);
    const [containerSize, setContainerSize] = useState({ width: 0, height: 0 });
    const [imageSizes, setImageSizes] = useState({});

    // Drawing and region states
    const [features, setFeatures] = useState(
        selectedSamples.reduce((acc, sample) => ({
            ...acc,
            [sample.id]: { type: 'FeatureCollection', features: [] }
        }), {})
    );
    const [tempRegions, setTempRegions] = useState({});
    const [sampleOffsets, setSampleOffsets] = useState({});
    const [regionName, setRegionName] = useState('');
    const [regionColor, setRegionColor] = useState(generateRandomColor());
    const [isDrawingActive, setIsDrawingActive] = useState(false);
    const [activeDrawingSample, setActiveDrawingSample] = useState(null);
    const [currentZoom, setCurrentZoom] = useState(-3);
    const [isZooming, setIsZooming] = useState(false);

    // Gene-related states
    const [geneList, setGeneList] = useState({});
    const [selectedGenes, setSelectedGenes] = useState([]);
    const [geneExpressionData, setGeneExpressionData] = useState({});
    const [hoveredCell, setHoveredCell] = useState(null);

    const debouncedSetZoom = useCallback(
        debounce((zoom) => {
            setCurrentZoom(zoom);
            setIsZooming(false);
        }, 50),
        []
    );

    // Views
    const mainView = new OrthographicView({
        id: 'main',
        controller: true
    });

    const minimapView = new OrthographicView({
        id: 'minimap',
        x: containerSize.width - 210,
        y: containerSize.height - 210,
        width: 200,
        height: 200,
        controller: false
    });

    // Set container size
    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        const resizeObserver = new ResizeObserver((entries) => {
            for (const entry of entries) {
                const { width, height } = entry.contentRect;
                setContainerSize({ width, height });
            }
        });

        resizeObserver.observe(container);
        return () => resizeObserver.disconnect();
    }, []);

    // Fetch gene list for selected samples
    useEffect(() => {
        const sampleIds = selectedSamples.map(item => item.id);

        fetch('/api/get_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_ids: sampleIds })
        })
            .then(res => res.json())
            .then(data => {
                setGeneList(prev => ({ ...prev, ...data }));
            })
            .catch(error => console.error('Error fetching gene list:', error));
    }, [selectedSamples]);

    // Get image sizes for all samples
    useEffect(() => {
        const fetchImageSizes = () => {
            const sampleIds = selectedSamples.map(sample => sample.id);

            fetch('/api/get_hires_image_size', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ sample_ids: sampleIds })
            })
                .then(res => res.json())
                .then(data => {
                    setImageSizes(data);
                })
        };

        if (selectedSamples.length) {
            fetchImageSizes();
        }
    }, [selectedSamples]);

    // Calculate sample offsets and initial view state
    useEffect(() => {
        if (Object.keys(imageSizes).length === selectedSamples.length) {
            const offsets = {};
            let currentX = 0;

            selectedSamples.forEach((sample, index) => {
                offsets[sample.id] = [currentX, 0];
                if (imageSizes[sample.id]) {
                    currentX += imageSizes[sample.id][0] + 500; // 500px gap between samples
                }
            });

            setSampleOffsets(offsets);
        }
    }, [imageSizes, selectedSamples]);

    // Set initial view state
    useEffect(() => {
        if (!selectedSamples.length || !sampleOffsets || !imageSizes) return;

        const firstSample = selectedSamples[0];
        if (!firstSample) return;

        const offset = sampleOffsets[firstSample.id] ?? [0, 0];
        const size = imageSizes[firstSample.id] ?? [0, 0];

        setMainViewState({
            target: [
                offset[0] + size[0] / 2,
                offset[1] + size[1] / 2,
                0
            ],
            zoom: -3,
            maxZoom: 2.5,
            minZoom: -5
        });
    }, [selectedSamples, sampleOffsets, imageSizes]);

    // Calculate overall bounds for minimap
    const overallBounds = useMemo(() => {
        if (!selectedSamples.length) return null;

        let xMin = Infinity, yMin = Infinity, xMax = -Infinity, yMax = -Infinity;

        selectedSamples.forEach(s => {
            const offset = sampleOffsets[s.id] || [0, 0];
            const size = imageSizes[s.id] || [0, 0];

            xMin = Math.min(xMin, offset[0]);
            yMin = Math.min(yMin, offset[1]);
            xMax = Math.max(xMax, offset[0] + size[0]);
            yMax = Math.max(yMax, offset[1] + size[1]);
        });

        return { xMin, yMin, xMax, yMax };
    }, [selectedSamples, sampleOffsets, imageSizes]);

    // Filter cell data for scatter plots
    const filteredCellData = useMemo(() => {
        return selectedSamples.reduce((acc, sample) => {
            const cellData = coordinatesData[sample.id] || [];
            const offset = sampleOffsets[sample.id] || [0, 0];

            acc[sample.id] = cellData.map(cell => ({
                ...cell,
                x: cell.cell_x + offset[0],
                y: cell.cell_y + offset[1]
            }));

            return acc;
        }, {});
    }, [selectedSamples, coordinatesData, sampleOffsets]);

    // Generate tissue image layers
    const generateImageLayers = useCallback(() => {
        return selectedSamples.map(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];

            if (!imageSize || imageSize.length < 2) return null;

            return new BitmapLayer({
                id: `tissue-image-${sample.id}`,
                image: `/${sample.id}_full.jpg`,
                bounds: [
                    offset[0],
                    offset[1] + imageSize[1],
                    offset[0] + imageSize[0],
                    offset[1]
                ],
                opacity: 0.8,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [selectedSamples, imageSizes, sampleOffsets]);

    // Generate cell scatter layers
    const generateCellLayers = useCallback(() => {
        return selectedSamples.map(sample => {
            const cellData = filteredCellData[sample.id] || [];

            return new ScatterplotLayer({
                id: `cells-${sample.id}`,
                data: cellData,
                getPosition: d => [d.x, d.y],
                getRadius: 3,
                getColor: d => {
                    // Color cells based on selected genes or default
                    if (selectedGenes.length > 0 && geneExpressionData[sample.id]) {
                        const cellExpression = geneExpressionData[sample.id][d.id];
                        if (cellExpression) {
                            const intensity = Math.min(cellExpression.total_expression || 0, 10) / 10;
                            return [255 * intensity, 100, 100, 200];
                        }
                    }
                    // Default gray color
                    return [150, 150, 150, 100];
                },
                pickable: true,
                radiusUnits: 'pixels',
                radiusScale: 1,
                radiusMinPixels: 1,
                radiusMaxPixels: 10,
                updateTriggers: {
                    getColor: [selectedGenes, geneExpressionData]
                },
                transitions: {
                    getColor: 0  // Disable color transitions to prevent flashing
                }
            });
        }).filter(Boolean);
    }, [selectedSamples, filteredCellData, selectedGenes, geneExpressionData]);

    // Generate editable layers for drawing
    const generateEditLayers = useCallback(() => {
        return selectedSamples.map(sample => {
            if (!isDrawingActive || activeDrawingSample !== sample.id) {
                return null;
            }

            return new EditableGeoJsonLayer({
                id: `editable-${sample.id}`,
                data: features[sample.id] || { type: 'FeatureCollection', features: [] },
                mode: DrawPolygonMode,
                selectedFeatureIndexes: [],
                onEdit: ({ updatedData }) => {
                    setFeatures(prev => ({
                        ...prev,
                        [sample.id]: updatedData
                    }));

                    // Update temp regions
                    handleRegionUpdate(sample.id, updatedData);
                },
                getLineColor: regionColor,
                getFillColor: [...regionColor, 50],
                getLineWidth: 2,
                pickable: true
            });
        }).filter(Boolean);
    }, [selectedSamples, features, isDrawingActive, activeDrawingSample, regionColor]);

    // Handle region updates during drawing
    const handleRegionUpdate = (sampleId, updatedData) => {
        if (!updatedData.features.length) return;

        const polygon = updatedData.features[0];
        if (!polygon || !polygon.geometry) return;

        const cellData = filteredCellData[sampleId] || [];
        const cellsInRegion = cellData.filter(cell => {
            const point = [cell.x, cell.y];
            return booleanPointInPolygon(point, polygon);
        });

        setTempRegions(prev => ({
            ...prev,
            [sampleId]: {
                polygon,
                cellIds: cellsInRegion.map(cell => cell.id),
                cellCount: cellsInRegion.length
            }
        }));
    };

    // Save drawn region
    const handleSaveRegion = () => {
        if (!regionName.trim() || !activeDrawingSample || !tempRegions[activeDrawingSample]) {
            return;
        }

        const newRegion = {
            id: Date.now().toString(),
            name: regionName.trim(),
            sampleId: activeDrawingSample,
            cellIds: tempRegions[activeDrawingSample].cellIds,
            color: `rgb(${regionColor.join(',')})`,
            polygon: tempRegions[activeDrawingSample].polygon
        };

        setInterestedRegions(prev => [...prev, newRegion]);

        // Reset drawing state
        setIsDrawingActive(false);
        setActiveDrawingSample(null);
        setRegionName('');
        setRegionColor(generateRandomColor());
        setFeatures(prev => ({
            ...prev,
            [activeDrawingSample]: { type: 'FeatureCollection', features: [] }
        }));
        setTempRegions({});
    };

    // Delete region
    const handleDeleteRegion = (regionId) => {
        setInterestedRegions(prev => prev.filter(r => r.id !== regionId));
        if (analyzedRegion === regionId) {
            setAnalyzedRegion(null);
        }
    };

    // Gene selection handlers
    const onGeneSelect = (gene) => {
        if (!selectedGenes.includes(gene)) {
            setSelectedGenes(prev => [...prev, gene]);
        }
    };

    const clearGeneSelection = () => {
        setSelectedGenes([]);
        setGeneExpressionData({});
    };

    const confirmGeneSelection = (sampleId) => {
        if (selectedGenes.length === 0) return;

        setSampleDataLoading(true);

        fetch('/api/get_gene_expression', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                sample_id: sampleId,
                genes: selectedGenes
            })
        })
            .then(res => res.json())
            .then(data => {
                setGeneExpressionData(prev => ({
                    ...prev,
                    [sampleId]: data
                }));
                setSampleDataLoading(false);
            })
            .catch(error => {
                console.error('Error fetching gene expression:', error);
                setSampleDataLoading(false);
            });
    };

    // Minimap layers
    const minimapSampleOutlines = useMemo(() => selectedSamples.map(s => {
        const offset = sampleOffsets[s.id] || [0, 0];
        const size = imageSizes[s.id] || [0, 0];

        return {
            path: [
                [offset[0], offset[1]],
                [offset[0] + size[0], offset[1]],
                [offset[0] + size[0], offset[1] + size[1]],
                [offset[0], offset[1] + size[1]],
                [offset[0], offset[1]]
            ]
        };
    }), [selectedSamples, sampleOffsets, imageSizes]);

    const minimapRegionOutlines = useMemo(() => interestedRegions.map(r => ({
        path: r.polygon.geometry.coordinates[0],
        color: r.color
    })), [interestedRegions]);

    const minimapLayers = [
        new PathLayer({
            id: 'minimap-samples',
            data: minimapSampleOutlines,
            getPath: d => d.path,
            getColor: [200, 200, 200, 100],
            getWidth: 1,
            viewId: 'minimap'
        }),
        new PathLayer({
            id: 'minimap-regions',
            data: minimapRegionOutlines,
            getPath: d => d.path,
            getColor: d => [...convertHEXToRGB(d.color), 200],
            getWidth: 2,
            viewId: 'minimap'
        })
    ];

    // Combine all layers with stable references during zoom
    const layers = useMemo(() => {
        const newLayers = [
            ...generateImageLayers(),
            ...generateCellLayers(),
            ...generateEditLayers(),
            ...minimapLayers
        ];
        
        // During zoom operations, use cached layers to prevent flashing
        if (isZooming && layersRef.current.length > 0) {
            // Only update non-interactive layers during zoom
            const interactiveLayers = newLayers.filter(layer => 
                layer.id.includes('editable') || layer.id.includes('minimap')
            );
            const staticLayers = layersRef.current.filter(layer => 
                !layer.id.includes('editable') && !layer.id.includes('minimap')
            );
            
            layersRef.current = [...staticLayers, ...interactiveLayers];
            return layersRef.current;
        }
        
        layersRef.current = newLayers;
        return newLayers;
    }, [
        generateImageLayers,
        generateCellLayers,
        generateEditLayers,
        minimapLayers,
        isZooming
    ]);

    const handleViewStateChange = useCallback(({ viewState, viewId }) => {
        if (viewId === 'main') {
            setMainViewState(viewState);
            if (!isZooming) {
                setIsZooming(true);
            }
            debouncedSetZoom(viewState.zoom);
        }
    }, [debouncedSetZoom, isZooming]);

    // Gene settings component
    const GeneSettings = ({ geneList, sampleId }) => {
        const [searchText, setSearchText] = useState('');

        const filteredGenes = useMemo(() =>
            Object.entries(geneList)
                .sort(([, countA], [, countB]) => countB - countA)
                .filter(([gene]) => gene.toLowerCase().includes(searchText.toLowerCase()))
            || []
            , [geneList, searchText]);

        return (
            <div style={{ maxHeight: 400 }}>
                <Input.Search
                    size="small"
                    placeholder="Search genes"
                    value={searchText}
                    onChange={e => setSearchText(e.target.value)}
                    style={{ marginBottom: 8 }}
                />
                <div style={{ maxHeight: 340, marginBottom: 10, overflowY: 'auto' }}>
                    {filteredGenes.map(([gene, count]) => (
                        <div key={gene} style={{
                            display: 'flex',
                            alignItems: 'center',
                            padding: '8px 0',
                            borderBottom: '1px solid #f0f0f0',
                        }}>
                            <Checkbox
                                style={{ marginRight: 8 }}
                                checked={selectedGenes.includes(gene)}
                                onChange={e => {
                                    if (e.target.checked) {
                                        onGeneSelect(gene);
                                    } else {
                                        setSelectedGenes(prev => prev.filter(g => g !== gene));
                                    }
                                }}
                            />
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', width: '100%' }}>
                                <span style={{ fontSize: 12 }}>{gene}</span>
                                <span style={{ fontSize: 12, color: '#666' }}>{count}</span>
                            </div>
                        </div>
                    ))}
                </div>
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 5 }}>
                    <Button size='small' style={{ width: '50%' }} onClick={clearGeneSelection}>Clear</Button>
                    <Button size='small' style={{ width: '50%' }} onClick={() => confirmGeneSelection(sampleId)}>Confirm</Button>
                </div>
            </div>
        );
    };

    const collapseItems = selectedSamples.map((sample) => ({
        key: sample.id,
        label: sample.name,
        children: (
            <GeneSettings
                geneList={geneList[sample.id] || {}}
                sampleId={sample.id}
            />
        )
    }));

    return (
        <div style={{ width: '100%', height: '100%', display: 'flex', borderRight: '2px solid #e8e8e8' }}>
            <div
                ref={containerRef}
                style={{
                    flex: 1,
                    position: 'relative',
                    cursor: isDrawingActive ? 'crosshair' : (hoveredCell ? 'pointer' : 'grab')
                }}
            >
                <DeckGL
                    layers={layers}
                    views={[mainView, minimapView]}
                    viewState={{
                        main: mainViewState,
                        minimap: overallBounds ? {
                            target: [
                                (overallBounds.xMin + overallBounds.xMax) / 2,
                                (overallBounds.yMin + overallBounds.yMax) / 2,
                                0
                            ],
                            zoom: Math.log2(Math.min(
                                200 / (overallBounds.xMax - overallBounds.xMin),
                                200 / (overallBounds.yMax - overallBounds.yMin)
                            ))
                        } : undefined
                    }}
                    onViewStateChange={handleViewStateChange}
                    onClick={info => {
                        if (info.picked && info.layer?.id?.includes('cells-')) {
                            setHoveredCell({
                                ...info.object,
                                sampleId: info.layer.id.replace('cells-', ''),
                                x: info.x,
                                y: info.y
                            });
                        }
                    }}
                    onHover={info => {
                        if (info.picked && info.layer?.id?.includes('cells-')) {
                            setHoveredCell({
                                ...info.object,
                                sampleId: info.layer.id.replace('cells-', ''),
                                x: info.x,
                                y: info.y
                            });
                        } else {
                            setHoveredCell(null);
                        }
                    }}
                    controller={{
                        inertia: false,  // Disable inertia to prevent momentum-based flashing
                        scrollZoom: {
                            speed: 0.005,  // Much slower zoom speed for better control
                            smooth: false  // Disable smooth zoom to prevent intermediate frames
                        },
                        dragPan: !isDrawingActive,
                        dragRotate: false,
                        touchAction: 'pan-y'  // Better touch handling
                    }}
                />

                {/* Hover tooltip */}
                {hoveredCell && (
                    <div style={{
                        position: 'absolute',
                        left: hoveredCell.x,
                        top: hoveredCell.y - 40,
                        pointerEvents: 'none',
                        backgroundColor: 'rgba(255, 255, 255, 0.9)',
                        padding: 8,
                        borderRadius: 4,
                        boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
                        transform: 'translateX(-50%)',
                        fontSize: 12,
                        zIndex: 1000,
                        textAlign: 'left'
                    }}>
                        <div><strong>Sample:</strong> {hoveredCell.sampleId}</div>
                        <div><strong>Cell Type:</strong> {hoveredCell.cell_type}</div>
                        {hoveredCell.id && (
                            <div><strong>Cell ID:</strong> {hoveredCell.id}</div>
                        )}
                    </div>
                )}

                {/* Gene list panel */}
                <div style={{ position: 'absolute', top: 10, left: 10, zIndex: 10 }}>
                    <Collapse
                        items={collapseItems}
                        defaultActiveKey={['0']}
                        style={{ background: '#ffffff', width: 300, opacity: 0.9 }}
                    />
                </div>

                {/* Drawing controls */}
                <div style={{
                    position: 'absolute',
                    top: 10,
                    right: 10,
                    zIndex: 10,
                    width: 250
                }}>
                    <div style={{
                        background: 'rgba(255,255,255,0.9)',
                        padding: '10px',
                        borderRadius: '10px',
                        display: 'flex',
                        flexDirection: 'column',
                        marginBottom: 10
                    }}>
                        <Input
                            size='small'
                            allowClear
                            placeholder="New region name"
                            value={regionName}
                            onChange={e => setRegionName(e.target.value)}
                            style={{ marginBottom: 8 }}
                        />
                        <div style={{ display: 'flex', gap: 5, marginBottom: 8 }}>
                            <ColorPicker
                                size='small'
                                value={`rgb(${regionColor.join(',')})`}
                                onChange={color => {
                                    const rgb = color.toRgb();
                                    setRegionColor([rgb.r, rgb.g, rgb.b]);
                                }}
                            />
                            <Button
                                size='small'
                                onClick={() => {
                                    if (isDrawingActive && regionName) {
                                        handleSaveRegion();
                                    } else {
                                        setIsDrawingActive(true);
                                        setActiveDrawingSample(selectedSamples[0]?.id);
                                    }
                                }}
                                block
                            >
                                {isDrawingActive && regionName ? 'Save' : 'Draw'}
                            </Button>
                        </div>
                        {tempRegions[activeDrawingSample] && (
                            <div style={{ fontSize: 12, color: '#666' }}>
                                Cells in region: {tempRegions[activeDrawingSample].cellCount}
                            </div>
                        )}
                    </div>

                    {/* Saved regions list */}
                    {interestedRegions.length > 0 && (
                        <div style={{
                            background: 'rgba(255,255,255,0.9)',
                            padding: '10px',
                            borderRadius: '10px',
                            maxHeight: 300,
                            overflowY: 'auto'
                        }}>
                            <div style={{ marginBottom: 8, fontWeight: 'bold' }}>Saved Regions</div>
                            {interestedRegions.map(region => (
                                <div key={region.id} style={{
                                    display: 'flex',
                                    justifyContent: 'space-between',
                                    alignItems: 'center',
                                    padding: '5px 0',
                                    borderBottom: '1px solid #f0f0f0'
                                }}>
                                    <div style={{ flex: 1 }}>
                                        <div style={{ fontSize: 12, fontWeight: 'bold' }}>{region.name}</div>
                                        <div style={{ fontSize: 10, color: '#666' }}>
                                            {region.cellIds.length} cells
                                        </div>
                                    </div>
                                    <Button
                                        size="small"
                                        type="text"
                                        icon={<CloseOutlined />}
                                        onClick={() => handleDeleteRegion(region.id)}
                                    />
                                </div>
                            ))}
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
};

export default MultiSampleViewerNew;
