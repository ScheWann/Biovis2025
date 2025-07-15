import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { GeneSettings } from './GeneList';
import { Collapse, Radio, Button, Input, ColorPicker } from "antd";
import { CloseOutlined } from '@ant-design/icons';
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
    
    // Area customization tooltip state
    const [isAreaTooltipVisible, setIsAreaTooltipVisible] = useState(false);
    const [pendingArea, setPendingArea] = useState(null);
    const [areaName, setAreaName] = useState('');
    const [areaColor, setAreaColor] = useState('#ff0000');
    const [tooltipPosition, setTooltipPosition] = useState({ x: 0, y: 0 });

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
            // Create pending area and show tooltip for customization
            const newPendingArea = {
                id: `area-${Date.now()}`,
                sampleId: currentDrawingSample,
                points: [...drawingPoints],
                name: `Custom Area ${customAreas.length + 1}`,
                color: '#ff0000'
            };
            
            // Find the rightmost point of the drawn area for tooltip positioning
            const rightmostPoint = findRightmostPoint(drawingPoints);
            // Find the vertical center of the drawn area
            const verticalCenter = findVerticalCenter(drawingPoints);
            
            setPendingArea(newPendingArea);
            setAreaName(newPendingArea.name);
            setAreaColor(newPendingArea.color);
            // Use rightmost x-coordinate but vertical center for y-coordinate
            setTooltipPosition({ x: rightmostPoint.x, y: verticalCenter.y });
            setIsAreaTooltipVisible(true);
            
            // Keep drawing state until tooltip is closed so area remains visible
            // Don't clear drawing state here - it will be cleared when tooltip is closed
        } else {
            // If not enough points, just clear the drawing state
            setIsDrawing(false);
            setDrawingPoints([]);
            setCurrentDrawingSample(null);
        }
    };

    const cancelDrawing = () => {
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
    };

    const deleteArea = (areaId) => {
        setCustomAreas(prev => prev.filter(area => area.id !== areaId));
    };

    // Handle area tooltip actions
    const handleAreaTooltipSave = () => {
        if (pendingArea) {
            const finalArea = {
                ...pendingArea,
                name: areaName || pendingArea.name,
                color: areaColor
            };
            setCustomAreas(prev => [...prev, finalArea]);
        }
        
        // Clear all drawing and tooltip state
        setIsAreaTooltipVisible(false);
        setPendingArea(null);
        setAreaName('');
        setAreaColor('#ff0000');
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
        setMousePosition(null);
    };

    const handleAreaTooltipCancel = () => {
        // Clear all drawing and tooltip state without saving
        setIsAreaTooltipVisible(false);
        setPendingArea(null);
        setAreaName('');
        setAreaColor('#ff0000');
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
        setMousePosition(null);
    };

    // Calculate centroid of a polygon for tooltip positioning
    const calculateCentroid = (points) => {
        if (points.length === 0) return { x: 0, y: 0 };
        
        const sum = points.reduce((acc, point) => ({
            x: acc.x + point[0],
            y: acc.y + point[1]
        }), { x: 0, y: 0 });
        
        return {
            x: sum.x / points.length,
            y: sum.y / points.length
        };
    };

    // Find the rightmost point of a polygon for tooltip positioning
    const findRightmostPoint = (points) => {
        if (points.length === 0) return { x: 0, y: 0 };
        
        const rightmost = points.reduce((max, point) => {
            return point[0] > max[0] ? point : max;
        }, points[0]);
        
        return {
            x: rightmost[0],
            y: rightmost[1]
        };
    };

    // Find the vertical center of a polygon for tooltip positioning
    const findVerticalCenter = (points) => {
        if (points.length === 0) return { x: 0, y: 0 };
        
        const yCoordinates = points.map(point => point[1]);
        const minY = Math.min(...yCoordinates);
        const maxY = Math.max(...yCoordinates);
        const centerY = (minY + maxY) / 2;
        
        return {
            x: 0, // x is not used, only y matters for vertical center
            y: centerY
        };
    };

    // Convert world coordinates to screen coordinates for tooltip positioning
    const worldToScreen = (worldX, worldY) => {
        if (!mainViewState || !containerRef.current) return { x: 0, y: 0 };
        
        const container = containerRef.current;
        const rect = container.getBoundingClientRect();
        
        // Get the center of the view
        const centerX = rect.width / 2;
        const centerY = rect.height / 2;
        
        // Calculate the scale based on zoom
        const scale = Math.pow(2, mainViewState.zoom);
        
        // Transform world coordinates to screen coordinates
        const screenX = centerX + (worldX - mainViewState.target[0]) * scale;
        const screenY = centerY + (worldY - mainViewState.target[1]) * scale; // Removed inversion
        
        return { x: screenX, y: screenY };
    };

    // Calculate tooltip position with real-time updates based on current view state
    const getTooltipPosition = useCallback(() => {
        if (!isAreaTooltipVisible || !pendingArea || !containerRef.current || !mainViewState) {
            return { left: 0, top: 0 };
        }

        // Recalculate screen position based on current view state
        const screenPos = worldToScreen(tooltipPosition.x, tooltipPosition.y);
        const container = containerRef.current;
        const rect = container.getBoundingClientRect();
        
        const tooltipWidth = 280;
        const tooltipHeight = 180;
        
        // Calculate position relative to viewport (for fixed positioning)
        const left = rect.left + screenPos.x + 20; // 20px to the right of the rightmost point
        let top = rect.top + screenPos.y - tooltipHeight / 2; // Center vertically on the area's vertical center
        
        // Only constrain the top position to stay within the viewport bounds
        if (top < 10) {
            top = 10; // Force to top edge of viewport
        }
        
        if (top + tooltipHeight > window.innerHeight - 10) {
            top = window.innerHeight - tooltipHeight - 10; // Force to bottom edge of viewport
        }
        
        return { left, top };
    }, [isAreaTooltipVisible, pendingArea, tooltipPosition, mainViewState]);

    // Convert hex color to RGB array
    const hexToRgb = (hex) => {
        const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? [
            parseInt(result[1], 16),
            parseInt(result[2], 16),
            parseInt(result[3], 16)
        ] : [255, 0, 0];
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
                                marginBottom: 3,
                                padding: '4px 6px',
                                backgroundColor: '#f8f9fa',
                                borderRadius: 4,
                                border: `1px solid ${area.color || '#ff0000'}`,
                                borderLeft: `4px solid ${area.color || '#ff0000'}`
                            }}>
                                <span style={{ fontWeight: 500 }}>{area.name}</span>
                                <Button
                                    size="small"
                                    type="text"
                                    danger
                                    onClick={() => deleteArea(area.id)}
                                    style={{ fontSize: '10px', padding: '0 4px', height: 20 }}
                                >
                                    ×
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
            const areaColor = hexToRgb(area.color || '#ff0000');
            
            layers.push(new PolygonLayer({
                id: `custom-area-${area.id}`,
                data: [{ polygon: area.points }],
                getPolygon: d => d.polygon,
                getFillColor: [...areaColor, 50],
                getLineColor: [...areaColor, 200],
                getLineWidth: 2,
                lineWidthUnits: 'pixels',
                pickable: true,
            }));
        });

        // Current drawing visualization
        if ((isDrawing || isAreaTooltipVisible) && (drawingPoints.length > 0 || pendingArea)) {
            const pointsToShow = pendingArea ? pendingArea.points : drawingPoints;
            
            // When tooltip is visible, show the finalized area without animation effects
            if (isAreaTooltipVisible && pendingArea) {
                // Show finalized polygon with pending area color
                const pendingAreaColor = hexToRgb(areaColor || pendingArea.color);
                layers.push(new PolygonLayer({
                    id: 'pending-area-preview',
                    data: [{ polygon: pendingArea.points }],
                    getPolygon: d => d.polygon,
                    getFillColor: [...pendingAreaColor, 50],
                    getLineColor: [...pendingAreaColor, 200],
                    getLineWidth: 2,
                    lineWidthUnits: 'pixels',
                    pickable: false,
                }));
                
                // Show clean points without animation
                layers.push(new ScatterplotLayer({
                    id: 'pending-area-points',
                    data: pendingArea.points.map((point, index) => ({
                        position: point,
                        index
                    })),
                    getPosition: d => d.position,
                    getRadius: 4,
                    getFillColor: pendingAreaColor,
                    getLineColor: [0, 0, 0, 150],
                    getLineWidth: 1,
                    radiusUnits: 'pixels',
                    lineWidthUnits: 'pixels',
                    pickable: false,
                }));
            } else {
                // Normal drawing mode with animation effects
                // Draw completed points
                if (pointsToShow.length > 0) {
                    layers.push(new ScatterplotLayer({
                        id: 'drawing-points',
                        data: pointsToShow.map((point, index) => ({
                            position: point,
                            index,
                            isFirst: index === 0
                        })),
                        getPosition: d => d.position,
                        getRadius: d => {
                            if (d.isFirst && pointsToShow.length >= 3) {
                                // Make first point larger and pulsing when it can be snapped to
                                const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                                return shouldSnap ? 12 : 8;
                            }
                            return 5;
                        },
                        getFillColor: d => {
                            if (d.isFirst && pointsToShow.length >= 3) {
                                const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                                return shouldSnap ? [255, 255, 0, 255] : [255, 255, 0, 200];
                            }
                            return [0, 255, 0, 200];
                        },
                        getLineColor: d => {
                            if (d.isFirst && pointsToShow.length >= 3) {
                                const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                                return shouldSnap ? [255, 165, 0, 255] : [200, 200, 0, 255];
                            }
                            return [0, 150, 0, 255];
                        },
                        getLineWidth: d => d.isFirst && pointsToShow.length >= 3 ? 3 : 2,
                        radiusUnits: 'pixels',
                        lineWidthUnits: 'pixels',
                        pickable: false,
                    }));
                }

                // Draw lines connecting the points
                if (pointsToShow.length > 1) {
                    const lineSegments = [];
                    for (let i = 0; i < pointsToShow.length - 1; i++) {
                        lineSegments.push({
                            sourcePosition: pointsToShow[i],
                            targetPosition: pointsToShow[i + 1]
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
                if (pointsToShow.length >= 3) {
                    layers.push(new PolygonLayer({
                        id: 'drawing-preview',
                        data: [{ polygon: pointsToShow }],
                        getPolygon: d => d.polygon,
                        getFillColor: [0, 255, 0, 30],
                        getLineColor: [0, 255, 0, 0],
                        getLineWidth: 0,
                        pickable: false,
                    }));
                }
            }

            // Only show mouse interactions if still actively drawing (not just showing tooltip)
            if (isDrawing && !isAreaTooltipVisible) {
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
            }
        }

        return layers;
    }, [customAreas, isDrawing, isAreaTooltipVisible, drawingPoints, pendingArea, areaColor, currentDrawingSample, sampleOffsets, mousePosition, shouldSnapToFirst]);

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
                    controller={
                        isAreaTooltipVisible ? false : // Disable all interactions when tooltip is visible
                        !isDrawing ? true : { dragPan: false, dragRotate: false, doubleClickZoom: false }
                    }
                    getCursor={({ isHovering, isDragging }) => {
                        if (isAreaTooltipVisible) {
                            return 'default'; // Normal cursor when tooltip is visible
                        }
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
                {/* {isDrawing && (
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
                )} */}

                {/* Area Customization Tooltip */}
                {isAreaTooltipVisible && pendingArea && (
                    <>
                        {/* Overlay to prevent interactions with the map */}
                        <div 
                            style={{
                                position: 'absolute',
                                top: 0,
                                left: 0,
                                right: 0,
                                bottom: 0,
                                zIndex: 9000, // Lower than tooltip but higher than deck.gl
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                cursor: 'default',
                                pointerEvents: 'auto'
                            }}
                            onClick={(e) => e.stopPropagation()}
                        />
                        
                        <div 
                            style={{
                                position: 'fixed', // Changed from absolute to fixed to allow positioning outside container
                                left: getTooltipPosition().left,
                                top: getTooltipPosition().top,
                                zIndex: 10000, // Increased z-index to ensure it's above everything
                                background: '#ffffff',
                                border: '1px solid #d9d9d9',
                                borderRadius: 8,
                                boxShadow: '0 6px 16px rgba(0, 0, 0, 0.12)',
                                padding: 12,
                                minWidth: 240,
                                maxWidth: 280,
                                pointerEvents: 'auto' // Ensure tooltip is interactive
                            }}
                        >
                        {/* Close button */}
                        <div style={{ 
                            position: 'absolute', 
                            top: 8, 
                            right: 8,
                            cursor: 'pointer',
                            padding: 4,
                            borderRadius: 4,
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center'
                        }}
                        onClick={handleAreaTooltipCancel}
                        onMouseEnter={(e) => e.target.style.backgroundColor = '#f5f5f5'}
                        onMouseLeave={(e) => e.target.style.backgroundColor = 'transparent'}
                        >
                            <CloseOutlined style={{ fontSize: 12, color: '#666' }} />
                        </div>

                        {/* Title */}
                        <div style={{ 
                            fontWeight: 'bold', 
                            marginBottom: 12, 
                            fontSize: 14,
                            color: '#262626',
                            paddingRight: 20
                        }}>
                            Customize Area
                        </div>

                        {/* Area Name Input */}
                        <div style={{ marginBottom: 12 }}>
                            <label style={{ 
                                display: 'block', 
                                marginBottom: 6, 
                                fontSize: 12,
                                fontWeight: 500,
                                color: '#595959'
                            }}>
                                Area Name:
                            </label>
                            <Input
                                value={areaName}
                                onChange={(e) => setAreaName(e.target.value)}
                                placeholder="Enter area name"
                                maxLength={50}
                                size="small"
                            />
                        </div>
                        
                        {/* Color Picker */}
                        <div style={{ marginBottom: 12 }}>
                            <label style={{ 
                                display: 'block', 
                                marginBottom: 6,
                                fontSize: 12,
                                fontWeight: 500,
                                color: '#595959'
                            }}>
                                Area Color:
                            </label>
                            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                <ColorPicker
                                    value={areaColor}
                                    onChange={(color) => setAreaColor(color.toHexString())}
                                    size="small"
                                />
                                <div 
                                    style={{ 
                                        width: 24, 
                                        height: 24, 
                                        backgroundColor: areaColor,
                                        border: '1px solid #d9d9d9',
                                        borderRadius: 4
                                    }}
                                />
                            </div>
                        </div>

                        {/* Action Buttons */}
                        <div style={{ display: 'flex', gap: 8, justifyContent: 'flex-end' }}>
                            <Button 
                                size="small" 
                                onClick={handleAreaTooltipCancel}
                            >
                                Cancel
                            </Button>
                            <Button 
                                size="small" 
                                type="primary"
                                onClick={handleAreaTooltipSave}
                            >
                                Save Area
                            </Button>
                        </div>
                        </div>
                    </>
                )}
            </div>
        </div>
    );
};