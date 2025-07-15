import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { GeneSettings } from './GeneList';
import { Collapse, Radio, Button } from "antd";
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, PolygonLayer, LineLayer } from '@deck.gl/layers';


export const SampleViewer = ({
    selectedSamples,
    coordinatesData,
}) => {
    const containerRef = useRef(null);
    const [mainViewState, setMainViewState] = useState(null);
    const [containerSize, setContainerSize] = useState({ width: 0, height: 0 });
    const [imageSizes, setImageSizes] = useState({});
    const [availableGenes, setAvailableGenes] = useState([]); // All genes that have been added to the list
    const [selectedGenes, setSelectedGenes] = useState([]); // Currently selected (checked) genes
    const [radioCellGeneModes, setRadioCellGeneModes] = useState(
        selectedSamples.reduce((acc, sample) => ({ ...acc, [sample.id]: 'cellTypes' }), {})
    );

    // Drawing state
    const [isDrawing, setIsDrawing] = useState(false);
    const [drawingPoints, setDrawingPoints] = useState([]);
    const [customAreas, setCustomAreas] = useState([]);
    const [currentDrawingSample, setCurrentDrawingSample] = useState(null);
    const [mousePosition, setMousePosition] = useState(null);

    const radioOptions = [
        {
            label: 'Cell Type',
            value: 'cellTypes',
        },
        {
            label: 'Genes',
            value: 'genes',
        },
    ];

    // Calculate sample offsets based on image sizes
    const sampleOffsets = useMemo(() => {
        if (selectedSamples.length <= 1) return {};

        const offsets = {};
        let currentX = 0;

        selectedSamples.forEach((sample) => {
            offsets[sample.id] = [currentX, 0];
            if (imageSizes[sample.id]) {
                currentX += imageSizes[sample.id][0] + 500; // 500px gap between samples
            }
        });

        return offsets;
    }, [selectedSamples, imageSizes]);

    // Main view
    const mainView = new OrthographicView({
        id: 'main',
        controller: true
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

    // Consolidated effect for managing samples and view state
    useEffect(() => {
        if (selectedSamples.length === 0) return;

        const sampleIds = selectedSamples.map(sample => sample.id);

        // Fetch image sizes
        fetch('/api/get_hires_image_size', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_ids: sampleIds })
        })
            .then(res => res.json())
            .then(data => {
                setImageSizes(data);
            });
    }, [selectedSamples.length]);

    // Initialize radio modes for new samples
    useEffect(() => {
        setRadioCellGeneModes(prev => {
            const newModes = { ...prev };
            selectedSamples.forEach(sample => {
                if (!newModes[sample.id]) {
                    newModes[sample.id] = 'cellTypes';
                }
            });
            return newModes;
        });
    }, [selectedSamples]);

    // Set initial view state when image sizes or offsets change
    useEffect(() => {
        if (!selectedSamples.length || !imageSizes[selectedSamples[0]?.id]) return;

        const firstSample = selectedSamples[0];
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
    }, [selectedSamples, imageSizes, sampleOffsets]);

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

    const changeCellGeneMode = (sampleId, e) => {
        setRadioCellGeneModes(prev => ({
            ...prev,
            [sampleId]: e.target.value
        }));
    };

    // Drawing methods
    const startDrawing = (sampleId) => {
        setIsDrawing(true);
        setCurrentDrawingSample(sampleId);
        setDrawingPoints([]);
    };

    const finishDrawing = () => {
        if (drawingPoints.length >= 3 && currentDrawingSample) {
            const newArea = {
                id: `area-${Date.now()}`,
                sampleId: currentDrawingSample,
                points: [...drawingPoints],
                name: `Custom Area ${customAreas.length + 1}`
            };
            setCustomAreas(prev => [...prev, newArea]);
        }
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
    };

    const cancelDrawing = () => {
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
    };

    const deleteArea = (areaId) => {
        setCustomAreas(prev => prev.filter(area => area.id !== areaId));
    };

    // Undo last point
    const undoLastPoint = () => {
        if (drawingPoints.length > 0) {
            setDrawingPoints(prev => prev.slice(0, -1));
        }
    };

    // Check if point should snap to first point (auto-close)
    const shouldSnapToFirst = useCallback((currentPoint) => {
        if (drawingPoints.length < 3 || !mainViewState) return false;

        const firstPoint = drawingPoints[0];

        // Calculate distance in world coordinates
        const worldDistance = Math.sqrt(
            Math.pow(currentPoint[0] - firstPoint[0], 2) +
            Math.pow(currentPoint[1] - firstPoint[1], 2)
        );

        // Dynamic snap distance based on zoom level
        // At zoom 0, snap distance is ~50 world units
        // At zoom -3, snap distance is ~400 world units
        // At zoom 2, snap distance is ~6 world units
        const baseSnapDistance = 50;
        const zoomFactor = Math.pow(2, -mainViewState.zoom);
        const dynamicSnapDistance = baseSnapDistance * zoomFactor;

        return worldDistance < dynamicSnapDistance;
    }, [drawingPoints, mainViewState]);

    const handleMapClick = useCallback((info) => {
        if (!isDrawing || !info.coordinate) return;

        const [x, y] = info.coordinate;
        const currentPoint = [x, y];

        // Check for auto-close (snap to first point)
        if (shouldSnapToFirst(currentPoint)) {
            finishDrawing();
            return;
        }

        setDrawingPoints(prev => [...prev, currentPoint]);
    }, [isDrawing, shouldSnapToFirst]);

    // Track mouse movement for preview
    const handleMouseMove = useCallback((info) => {
        if (!isDrawing) {
            setMousePosition(null);
            return;
        }

        if (info.coordinate) {
            setMousePosition(info.coordinate);
        }
    }, [isDrawing]);

    const handleKeyPress = useCallback((event) => {
        if (!isDrawing) return;

        if (event.key === 'Enter') {
            finishDrawing();
        } else if (event.key === 'Escape') {
            cancelDrawing();
        } else if (event.key === 'Backspace' || event.key === 'Delete') {
            event.preventDefault();
            undoLastPoint();
        }
    }, [isDrawing, drawingPoints, currentDrawingSample]);

    // Add keyboard event listener
    useEffect(() => {
        document.addEventListener('keydown', handleKeyPress);
        return () => document.removeEventListener('keydown', handleKeyPress);
    }, [handleKeyPress]);

    // Create collapse items for each sample
    const collapseItems = selectedSamples.map((sample, index) => ({
        key: sample.id,
        label: sample.name,
        children: (
            <>
                <Radio.Group
                    block
                    options={radioOptions}
                    size='small'
                    value={radioCellGeneModes[sample.id]}
                    optionType="button"
                    style={{ marginBottom: 10 }}
                    onChange={(e) => changeCellGeneMode(sample.id, e)}
                />

                {/* Drawing controls */}
                <div style={{ marginBottom: 10 }}>
                    <Button
                        size="small"
                        type={isDrawing && currentDrawingSample === sample.id ? "primary" : "default"}
                        onClick={() => {
                            if (isDrawing && currentDrawingSample === sample.id) {
                                finishDrawing();
                            } else {
                                startDrawing(sample.id);
                            }
                        }}
                        style={{ marginRight: 5 }}
                    >
                        {isDrawing && currentDrawingSample === sample.id ? "Finish Area" : "Draw Area"}
                    </Button>
                    {isDrawing && currentDrawingSample === sample.id && (
                        <>
                            <Button size="small" onClick={undoLastPoint} style={{ marginRight: 5 }}>
                                Undo
                            </Button>
                            <Button size="small" onClick={cancelDrawing}>
                                Cancel
                            </Button>
                        </>
                    )}
                </div>

                {/* Show custom areas for this sample */}
                {customAreas.filter(area => area.sampleId === sample.id).length > 0 && (
                    <div style={{ marginBottom: 10 }}>
                        <div style={{ fontSize: '12px', fontWeight: 'bold', marginBottom: 5 }}>
                            Custom Areas:
                        </div>
                        {customAreas.filter(area => area.sampleId === sample.id).map(area => (
                            <div key={area.id} style={{
                                display: 'flex',
                                justifyContent: 'space-between',
                                alignItems: 'center',
                                fontSize: '11px',
                                marginBottom: 3
                            }}>
                                <span>{area.name}</span>
                                <Button
                                    size="small"
                                    type="text"
                                    danger
                                    onClick={() => deleteArea(area.id)}
                                    style={{ fontSize: '10px', padding: '0 4px' }}
                                >
                                    Delete
                                </Button>
                            </div>
                        ))}
                    </div>
                )}

                {radioCellGeneModes[sample.id] === 'cellTypes' ? (
                    <div style={{ padding: '10px 0' }}>
                        <div>Cell Type view selected</div>
                    </div>
                ) : (
                    <GeneSettings
                        sampleId={sample.id}
                        availableGenes={availableGenes}
                        setAvailableGenes={setAvailableGenes}
                        selectedGenes={selectedGenes}
                        setSelectedGenes={setSelectedGenes}
                    />
                )}
            </>
        )
    }));

    const handleViewStateChange = useCallback(({ viewState, viewId }) => {
        if (viewId === 'main') {
            setMainViewState(viewState);
        }
    }, []);

    // Generate tissue image layers
    const generateImageLayers = useCallback(() => {
        return selectedSamples.map(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];

            if (!imageSize) return null;

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
                getFillColor: [150, 150, 150, 100],
                pickable: true,
                radiusUnits: 'pixels',
            });
        }).filter(Boolean);
    }, [selectedSamples, filteredCellData]);

    // Generate custom area layers
    const generateCustomAreaLayers = useCallback(() => {
        const layers = [];

        // Completed custom areas
        customAreas.forEach(area => {
            // const offset = sampleOffsets[area.sampleId] || [0, 0];

            layers.push(new PolygonLayer({
                id: `custom-area-${area.id}`,
                data: [{ polygon: area.points }],
                getPolygon: d => d.polygon,
                getFillColor: [255, 0, 0, 50],
                getLineColor: [255, 0, 0, 200],
                getLineWidth: 2,
                lineWidthUnits: 'pixels',
                pickable: true,
            }));
        });

        // Current drawing visualization
        if (isDrawing && currentDrawingSample) {
            // Draw completed points
            if (drawingPoints.length > 0) {
                layers.push(new ScatterplotLayer({
                    id: 'drawing-points',
                    data: drawingPoints.map((point, index) => ({
                        position: point,
                        index,
                        isFirst: index === 0
                    })),
                    getPosition: d => d.position,
                    getRadius: d => {
                        if (d.isFirst && drawingPoints.length >= 3) {
                            // Make first point larger and pulsing when it can be snapped to
                            const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                            return shouldSnap ? 12 : 8;
                        }
                        return 5;
                    },
                    getFillColor: d => {
                        if (d.isFirst && drawingPoints.length >= 3) {
                            const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                            return shouldSnap ? [255, 255, 0, 255] : [255, 255, 0, 200];
                        }
                        return [0, 255, 0, 200];
                    },
                    getLineColor: d => {
                        if (d.isFirst && drawingPoints.length >= 3) {
                            const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                            return shouldSnap ? [255, 165, 0, 255] : [200, 200, 0, 255];
                        }
                        return [0, 150, 0, 255];
                    },
                    getLineWidth: d => d.isFirst && drawingPoints.length >= 3 ? 3 : 2,
                    radiusUnits: 'pixels',
                    lineWidthUnits: 'pixels',
                    pickable: false,
                }));
            }

            // Draw mouse cursor when hovering
            if (mousePosition) {
                const shouldSnap = shouldSnapToFirst(mousePosition);
                layers.push(new ScatterplotLayer({
                    id: 'drawing-cursor',
                    data: [{ position: mousePosition }],
                    getPosition: d => d.position,
                    getRadius: shouldSnap ? 10 : 4,
                    getFillColor: shouldSnap ? [255, 165, 0, 200] : [0, 255, 0, 150],
                    getLineColor: shouldSnap ? [255, 140, 0, 255] : [0, 255, 0, 200],
                    getLineWidth: shouldSnap ? 3 : 1,
                    radiusUnits: 'pixels',
                    lineWidthUnits: 'pixels',
                    pickable: false,
                }));

                // Draw preview line from last point to mouse
                if (drawingPoints.length > 0) {
                    const lastPoint = drawingPoints[drawingPoints.length - 1];
                    const targetPoint = shouldSnap ? drawingPoints[0] : mousePosition;

                    layers.push(new LineLayer({
                        id: 'drawing-preview-line',
                        data: [{
                            sourcePosition: lastPoint,
                            targetPosition: targetPoint
                        }],
                        getSourcePosition: d => d.sourcePosition,
                        getTargetPosition: d => d.targetPosition,
                        getColor: shouldSnap ? [255, 165, 0, 200] : [0, 255, 0, 150],
                        getWidth: shouldSnap ? 4 : 2,
                        widthUnits: 'pixels',
                        pickable: false,
                    }));
                }

                // Show snap zone indicator when close to first point
                if (shouldSnap && drawingPoints.length >= 3) {
                    const firstPoint = drawingPoints[0];
                    layers.push(new ScatterplotLayer({
                        id: 'snap-zone-indicator',
                        data: [{ position: firstPoint }],
                        getPosition: d => d.position,
                        getRadius: 15,
                        getFillColor: [255, 165, 0, 50],
                        getLineColor: [255, 165, 0, 150],
                        getLineWidth: 2,
                        radiusUnits: 'pixels',
                        lineWidthUnits: 'pixels',
                        pickable: false,
                    }));
                }
            }

            // Draw lines connecting the points
            if (drawingPoints.length > 1) {
                const lineSegments = [];
                for (let i = 0; i < drawingPoints.length - 1; i++) {
                    lineSegments.push({
                        sourcePosition: drawingPoints[i],
                        targetPosition: drawingPoints[i + 1]
                    });
                }

                layers.push(new LineLayer({
                    id: 'drawing-lines',
                    data: lineSegments,
                    getSourcePosition: d => d.sourcePosition,
                    getTargetPosition: d => d.targetPosition,
                    getColor: [0, 255, 0, 200],
                    getWidth: 2,
                    widthUnits: 'pixels',
                    pickable: false,
                }));
            }

            // Show polygon preview when we have at least 3 points
            if (drawingPoints.length >= 3) {
                layers.push(new PolygonLayer({
                    id: 'drawing-preview',
                    data: [{ polygon: drawingPoints }],
                    getPolygon: d => d.polygon,
                    getFillColor: [0, 255, 0, 30],
                    getLineColor: [0, 255, 0, 0],
                    getLineWidth: 0,
                    pickable: false,
                }));
            }
        }

        return layers;
    }, [customAreas, isDrawing, drawingPoints, currentDrawingSample, sampleOffsets, mousePosition, shouldSnapToFirst]);

    // Combine all layers
    const layers = useMemo(() => [
        ...generateImageLayers(),
        ...generateCellLayers(),
        ...generateCustomAreaLayers()
    ], [generateImageLayers, generateCellLayers, generateCustomAreaLayers]);

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div
                ref={containerRef}
                style={{
                    width: '100%',
                    height: '100%',
                    position: 'relative'
                }}
            >
                <DeckGL
                    layers={layers}
                    views={[mainView]}
                    viewState={{
                        main: mainViewState
                    }}
                    onViewStateChange={handleViewStateChange}
                    onClick={handleMapClick}
                    onHover={handleMouseMove}
                    controller={!isDrawing ? true : { dragPan: false, dragRotate: false, doubleClickZoom: false }}
                    getCursor={({ isHovering, isDragging }) => {
                        if (isDrawing) {
                            if (mousePosition && shouldSnapToFirst(mousePosition)) {
                                return 'pointer';
                            }
                            return 'crosshair';
                        }
                        return isHovering ? 'pointer' : 'grab';
                    }}
                />

                {/* Sample controls */}
                <div style={{ position: 'absolute', top: 10, left: 10, zIndex: 10 }}>
                    <Collapse
                        items={collapseItems}
                        defaultActiveKey={[selectedSamples[0]?.id]}
                        style={{ background: '#ffffff', width: 300, opacity: 0.9 }}
                    />
                </div>

                {/* Drawing instructions */}
                {isDrawing && (
                    <div style={{
                        position: 'absolute',
                        top: 10,
                        right: 10,
                        zIndex: 10,
                        background: '#ffffff',
                        padding: '10px',
                        borderRadius: '4px',
                        opacity: 0.9,
                        fontSize: '12px'
                    }}>
                        <div style={{ fontWeight: 'bold', marginBottom: '5px' }}>
                            Drawing Mode Active
                        </div>
                        <div>• Click to add points</div>
                        <div>• Press Enter to finish</div>
                        <div>• Press Escape to cancel</div>
                        <div>• Backspace/Delete to undo</div>
                        <div>• Click near first point to auto-close</div>
                        <div>• Points: {drawingPoints.length}</div>
                        {drawingPoints.length >= 3 && mousePosition && shouldSnapToFirst(mousePosition) && (
                            <div style={{ color: 'orange', fontWeight: 'bold' }}>
                                Click to close polygon!
                            </div>
                        )}
                        {drawingPoints.length >= 3 && (!mousePosition || !shouldSnapToFirst(mousePosition)) && (
                            <div style={{ color: 'green', fontWeight: 'bold' }}>
                                Ready to finish!
                            </div>
                        )}
                    </div>
                )}
            </div>
        </div>
    );
};