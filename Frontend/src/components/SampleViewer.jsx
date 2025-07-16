import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { GeneSettings } from './GeneList';
import { Collapse, Radio, Button, Input, ColorPicker, AutoComplete } from "antd";
import { CloseOutlined, EditOutlined, RedoOutlined, BorderOutlined } from '@ant-design/icons';
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
    const [areaColor, setAreaColor] = useState('#f72585');
    const [tooltipPosition, setTooltipPosition] = useState({ x: 0, y: 0 });

    // Area edit/delete popup state
    const [isAreaEditPopupVisible, setIsAreaEditPopupVisible] = useState(false);
    const [selectedAreaForEdit, setSelectedAreaForEdit] = useState(null);
    const [editAreaName, setEditAreaName] = useState('');
    const [editAreaColor, setEditAreaColor] = useState('#f72585');
    const [editPopupPosition, setEditPopupPosition] = useState({ x: 0, y: 0 });
    const [editNeighbors, setEditNeighbors] = useState(10);
    const [editNPcas, setEditNPcas] = useState(50);
    const [editResolution, setEditResolution] = useState(0.5);

    // Minimap state
    const [minimapVisible, setMinimapVisible] = useState(true);
    const [minimapAnimating, setMinimapAnimating] = useState(false);
    const minimapRef = useRef(null);

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

    // Main view
    const mainView = new OrthographicView({
        id: 'main',
        controller: true
    });

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

    // Change cell/gene mode for a sample
    const changeCellGeneMode = (sampleId, e) => {
        setRadioCellGeneModes(prev => ({
            ...prev,
            [sampleId]: e.target.value
        }));
    };

    // Reset view to initial position and zoom
    const resetView = () => {
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
            x: 0,
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

    // Calculate minimap viewport bounds based on current view state
    const getMinimapViewportBounds = useCallback(() => {
        if (!mainViewState || !containerSize.width || !selectedSamples.length) return null;

        const firstSample = selectedSamples[0];
        const imageSize = imageSizes[firstSample.id];
        const offset = sampleOffsets[firstSample.id] || [0, 0];

        if (!imageSize) return null;

        // Calculate the world bounds of the current viewport
        const scale = Math.pow(2, mainViewState.zoom);
        const halfWidth = containerSize.width / (2 * scale);
        const halfHeight = containerSize.height / (2 * scale);

        const viewportBounds = {
            left: mainViewState.target[0] - halfWidth,
            right: mainViewState.target[0] + halfWidth,
            top: mainViewState.target[1] - halfHeight,
            bottom: mainViewState.target[1] + halfHeight
        };

        // Convert to relative coordinates within the first sample's image
        const relativeBounds = {
            left: Math.max(0, (viewportBounds.left - offset[0]) / imageSize[0]),
            right: Math.min(1, (viewportBounds.right - offset[0]) / imageSize[0]),
            top: Math.max(0, (viewportBounds.top - offset[1]) / imageSize[1]),
            bottom: Math.min(1, (viewportBounds.bottom - offset[1]) / imageSize[1])
        };

        return relativeBounds;
    }, [mainViewState, containerSize, selectedSamples, imageSizes, sampleOffsets]);

    // Handle minimap click to navigate
    const handleMinimapClick = useCallback((event) => {
        if (!minimapRef.current || !selectedSamples.length) return;

        const firstSample = selectedSamples[0];
        const imageSize = imageSizes[firstSample.id];
        const offset = sampleOffsets[firstSample.id] || [0, 0];

        if (!imageSize) return;

        const rect = minimapRef.current.getBoundingClientRect();
        const x = (event.clientX - rect.left) / rect.width;
        const y = (event.clientY - rect.top) / rect.height;

        // Convert relative coordinates to world coordinates
        const worldX = offset[0] + x * imageSize[0];
        const worldY = offset[1] + y * imageSize[1];

        // Update main view to center on clicked position
        setMainViewState(prev => ({
            ...prev,
            target: [worldX, worldY, 0]
        }));
    }, [selectedSamples, imageSizes, sampleOffsets]);

    // Toggle minimap with fade animation
    const toggleMinimapVisible = useCallback(() => {
        if (minimapVisible) {
            setMinimapVisible(false);
            setMinimapAnimating(true);

            // Stop animating
            setTimeout(() => {
                setMinimapAnimating(false);
            }, 300);
        } else {
            setMinimapAnimating(true);

            setTimeout(() => {
                setMinimapVisible(true);
            }, 10);

            // Reset animating state after transition completes
            setTimeout(() => {
                setMinimapAnimating(false);
            }, 310);
        }
    }, [minimapVisible]);

    // Toggle drawing mode
    const toggleDrawingMode = () => {
        if (isDrawing) {
            // If currently drawing, finish or cancel
            if (drawingPoints.length >= 3) {
                finishDrawing();
            } else {
                cancelDrawing();
            }
        } else {
            // Start drawing mode - sample will be determined on first click
            setIsDrawing(true);
            setCurrentDrawingSample(null);
            setDrawingPoints([]);
        }
    };

    const finishDrawing = () => {
        if (drawingPoints.length >= 3 && currentDrawingSample) {
            // Create pending area and show tooltip for customization
            const newPendingArea = {
                id: `area-${Date.now()}`,
                sampleId: currentDrawingSample,
                points: [...drawingPoints],
                name: `Custom Area ${customAreas.length + 1}`,
                color: '#f72585'
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
        } else {
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

        setIsAreaTooltipVisible(false);
        setPendingArea(null);
        setAreaName('');
        setAreaColor('#f72585');
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
        setMousePosition(null);
    };

    // Clear all drawing and tooltip state without saving
    const handleAreaTooltipCancel = () => {
        setIsAreaTooltipVisible(false);
        setPendingArea(null);
        setAreaName('');
        setAreaColor('#f72585');
        setIsDrawing(false);
        setDrawingPoints([]);
        setCurrentDrawingSample(null);
        setMousePosition(null);
    };

    // Calculate tooltip position with real-time updates based on current view state
    const getTempAreaCompleteTooltipPosition = useCallback(() => {
        if (!isAreaTooltipVisible || !pendingArea || !containerRef.current || !mainViewState) {
            return { left: 0, top: 0 };
        }

        // Recalculate screen position based on current view state
        const screenPos = worldToScreen(tooltipPosition.x, tooltipPosition.y);
        const container = containerRef.current;
        const rect = container.getBoundingClientRect();

        const tooltipWidth = 280;
        const tooltipHeight = 180;

        // 20px to the right of the rightmost point
        const left = rect.left + screenPos.x + 20;
        let top = rect.top + screenPos.y - tooltipHeight / 2;

        // Only constrain the top position to stay within the viewport bounds
        if (top < 10) {
            top = 10;
        }

        if (top + tooltipHeight > window.innerHeight - 10) {
            // Force to bottom edge of viewport
            top = window.innerHeight - tooltipHeight - 10;
        }

        return { left, top };
    }, [isAreaTooltipVisible, pendingArea, tooltipPosition, mainViewState]);

    // Undo last point
    const undoLastPoint = () => {
        if (drawingPoints.length > 0) {
            setDrawingPoints(prev => prev.slice(0, -1));
        }
    };

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

    // Handle area click for editing
    const handleAreaClick = useCallback((info) => {
        if (isDrawing || isAreaTooltipVisible) return;

        // Find which area was clicked
        const clickedObject = info.object;
        if (clickedObject) {
            // Extract area ID from layer ID
            const layerId = info.layer?.id;
            if (layerId && layerId.startsWith('custom-area-')) {
                const areaId = layerId.replace('custom-area-', '');
                const area = customAreas.find(a => a.id === areaId);

                if (area) {
                    setSelectedAreaForEdit(area);
                    setEditAreaName(area.name);
                    setEditAreaColor(area.color);
                    setEditNeighbors(area.neighbors || 10);
                    setEditNPcas(area.n_pcas || 50);
                    setEditResolution(area.resolution || 0.5);

                    // Find the rightmost point and vertical center of the area
                    const rightmostPoint = findRightmostPoint(area.points);
                    const verticalCenter = findVerticalCenter(area.points);
                    const areaPosition = { x: rightmostPoint.x, y: verticalCenter.y };

                    const screenPos = worldToScreen(areaPosition.x, areaPosition.y);
                    const container = containerRef.current;
                    const rect = container.getBoundingClientRect();

                    // Popup dimensions
                    const popupWidth = 280;
                    const popupHeight = 300;

                    let left = rect.left + screenPos.x + 20;
                    let top = rect.top + screenPos.y - popupHeight / 2;

                    // Check if popup would go off the right edge of the screen
                    if (left + popupWidth > window.innerWidth - 10) {
                        left = rect.left + screenPos.x - popupWidth - 20;
                    }

                    // Ensure popup doesn't go off the left edge of the screen
                    if (left < 10) {
                        left = 10;
                    }

                    // Constrain the top position to stay within the viewport bounds
                    if (top < 10) {
                        top = 10;
                    }

                    if (top + popupHeight > window.innerHeight - 10) {
                        // Force to bottom edge of viewport
                        top = window.innerHeight - popupHeight - 10;
                    }

                    setEditPopupPosition({
                        x: left,
                        y: top
                    });

                    setIsAreaEditPopupVisible(true);
                }
            }
        }
    }, [isDrawing, isAreaTooltipVisible, customAreas, worldToScreen]);

    // Handle area edit save
    const handleAreaEditSave = () => {
        if (selectedAreaForEdit) {
            setCustomAreas(prev => prev.map(area =>
                area.id === selectedAreaForEdit.id
                    ? {
                        ...area,
                        name: editAreaName,
                        color: editAreaColor,
                        neighbors: editNeighbors,
                        n_pcas: editNPcas,
                        resolution: editResolution
                    }
                    : area
            ));
        }
        handleAreaEditCancel();
    };

    // Handle area deletion
    const handleAreaDelete = () => {
        if (selectedAreaForEdit) {
            setCustomAreas(prev => prev.filter(area => area.id !== selectedAreaForEdit.id));
        }
        handleAreaEditCancel();
    };

    // Cancel area edit popup
    const handleAreaEditCancel = () => {
        setIsAreaEditPopupVisible(false);
        setSelectedAreaForEdit(null);
        setEditAreaName('');
        setEditAreaColor('#f72585');
        setEditNeighbors(10);
        setEditNPcas(50);
        setEditResolution(0.5);
    };

    // Convert hex color to RGB array
    const hexToRgb = (hex) => {
        const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? [
            parseInt(result[1], 16),
            parseInt(result[2], 16),
            parseInt(result[3], 16)
        ] : [255, 0, 0];
    };

    // Determine which sample contains the given coordinate
    const getSampleAtCoordinate = useCallback((x, y) => {
        for (const sample of selectedSamples) {
            const offset = sampleOffsets[sample.id] || [0, 0];
            const imageSize = imageSizes[sample.id];

            if (imageSize) {
                const [offsetX, offsetY] = offset;
                const [width, height] = imageSize;

                if (x >= offsetX && x <= offsetX + width &&
                    y >= offsetY && y <= offsetY + height) {
                    return sample.id;
                }
            }
        }
        return null;
    }, [selectedSamples, sampleOffsets, imageSizes]);

    const handleMapClick = useCallback((info) => {
        // First check if we clicked on a custom area (for editing)
        if (!isDrawing && !isAreaTooltipVisible) {
            handleAreaClick(info);
            return;
        }

        if (!isDrawing || !info.coordinate) return;

        const [x, y] = info.coordinate;
        const currentPoint = [x, y];

        // If this is the first point and no sample is selected, determine the sample
        if (!currentDrawingSample && drawingPoints.length === 0) {
            const sampleId = getSampleAtCoordinate(x, y);
            if (!sampleId) return; // Click outside any sample
            setCurrentDrawingSample(sampleId);
        }

        // Check for auto-close (snap to first point)
        if (shouldSnapToFirst(currentPoint)) {
            finishDrawing();
            return;
        }

        setDrawingPoints(prev => [...prev, currentPoint]);
    }, [isDrawing, isAreaTooltipVisible, shouldSnapToFirst, currentDrawingSample, drawingPoints.length, getSampleAtCoordinate, handleAreaClick]);

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

            // Calculate dynamic radius based on zoom level
            // At zoom 0, radius = 5 pixels
            // At zoom -3, radius = 1 pixel (zoomed out)
            // At zoom 2, radius = 20 pixels (zoomed in)
            const baseRadius = 5;
            const zoomFactor = Math.pow(2, mainViewState?.zoom || 0);
            const dynamicRadius = Math.max(1, Math.min(20, baseRadius * zoomFactor));

            return new ScatterplotLayer({
                id: `cells-${sample.id}`,
                data: cellData,
                getPosition: d => [d.x, d.y],
                getRadius: dynamicRadius,
                getFillColor: [0, 0, 0, 0],
                getLineColor: [150, 150, 150, 255], // Outline color
                getLineWidth: 1,
                lineWidthUnits: 'pixels',
                pickable: true,
                radiusUnits: 'pixels',
                stroked: true,
                filled: false,
            });
        }).filter(Boolean);
    }, [selectedSamples, filteredCellData, mainViewState]);

    // Generate custom area layers
    const generateCustomAreaLayers = useCallback(() => {
        const layers = [];

        // Completed custom areas
        customAreas.forEach(area => {
            const areaColor = hexToRgb(area.color || '#ff0000');

            layers.push(new PolygonLayer({
                id: `custom-area-${area.id}`,
                data: [{ polygon: area.points, areaId: area.id }],
                getPolygon: d => d.polygon,
                getFillColor: [...areaColor, 50],
                getLineColor: [...areaColor, 200],
                getLineWidth: 2,
                lineWidthUnits: 'pixels',
                pickable: true,
                onClick: (info) => handleAreaClick(info),
            }));
        });

        // Current drawing visualization
        if ((isDrawing || isAreaTooltipVisible) && (drawingPoints.length > 0 || pendingArea)) {
            const pointsToShow = pendingArea ? pendingArea.points : drawingPoints;

            // When tooltip is visible, show the finalized area without animation effects
            if (isAreaTooltipVisible && pendingArea) {
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
                                return shouldSnap ? 10 : 4;
                            }
                            return 5;
                        },
                        getFillColor: d => {
                            if (d.isFirst && pointsToShow.length >= 3) {
                                const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                                return shouldSnap ? [255, 202, 58, 255] : [255, 202, 58, 200];
                            }
                            return [0, 255, 150, 200];
                        },
                        getLineColor: d => {
                            if (d.isFirst && pointsToShow.length >= 3) {
                                const shouldSnap = mousePosition && shouldSnapToFirst(mousePosition);
                                return shouldSnap ? [255, 165, 0, 255] : [200, 200, 0, 255];
                            }
                            return [0, 255, 150, 255];
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
                        getColor: [0, 255, 150, 200],
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
                        getFillColor: [0, 255, 150, 30],
                        getLineColor: [0, 255, 150, 0],
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
                        getRadius: shouldSnap ? 8 : 5,
                        getFillColor: shouldSnap ? [255, 202, 58, 200] : [0, 255, 150, 150],
                        getLineColor: shouldSnap ? [255, 140, 0, 255] : [0, 255, 150, 200],
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
                            getColor: shouldSnap ? [255, 202, 58, 200] : [0, 255, 150, 150],
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

    // Get image sizes for selected samples
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

    // Initialize modes(cell type or genes) for samples
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

    // Add keyboard event listener
    useEffect(() => {
        document.addEventListener('keydown', handleKeyPress);
        return () => document.removeEventListener('keydown', handleKeyPress);
    }, [handleKeyPress]);

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
                        isAreaTooltipVisible || isAreaEditPopupVisible ? false :
                            !isDrawing ? true : { dragPan: false, dragRotate: false, doubleClickZoom: false }
                    }
                    getCursor={({ isHovering, isDragging }) => {
                        if (isAreaTooltipVisible || isAreaEditPopupVisible) {
                            return 'default';
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
                <div style={{ position: 'absolute', top: 10, left: 10, zIndex: 20 }}>
                    <Collapse
                        items={collapseItems}
                        defaultActiveKey={[selectedSamples[0]?.id]}
                        style={{ background: '#ffffff', width: 300, opacity: 0.9 }}
                    />
                </div>

                {/* Global Drawing Control */}
                <div style={{ position: 'absolute', top: 10, right: 10, zIndex: 20, display: 'flex', flexDirection: 'column', gap: 8 }}>
                    <div style={{ display: 'flex', justifyContent: 'flex-end', gap: 8 }}>
                        {/* Reset View Button */}
                        <Button
                            size="big"
                            onClick={resetView}
                            style={{
                                backgroundColor: '#ffffff',
                                borderColor: '#d9d9d9',
                                color: '#000000',
                                boxShadow: '0 2px 8px rgba(0, 0, 0, 0.15)',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center'
                            }}
                            icon={<RedoOutlined style={{ fontSize: '18px' }} />}
                            title="Reset view to initial position and zoom"
                        />

                        {/* Minimap Toggle Button */}
                        <Button
                            size="big"
                            onClick={toggleMinimapVisible}
                            style={{
                                backgroundColor: minimapVisible ? '#1890ff' : '#ffffff',
                                borderColor: minimapVisible ? '#1890ff' : '#d9d9d9',
                                color: minimapVisible ? '#ffffff' : '#000000',
                                boxShadow: '0 2px 8px rgba(0, 0, 0, 0.15)',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center'
                            }}
                            icon={<BorderOutlined style={{ fontSize: '18px' }} />}
                            title={minimapVisible ? 'Hide minimap' : 'Show minimap'}
                        />

                        {/* Drawing Toggle Button */}
                        <Button
                            size="big"
                            onClick={toggleDrawingMode}
                            style={{
                                backgroundColor: isDrawing ? '#1890ff' : '#ffffff',
                                borderColor: isDrawing ? '#1890ff' : '#ffffff',
                                color: isDrawing ? '#ffffff' : '#000000',
                                boxShadow: '0 2px 8px rgba(0, 0, 0, 0.15)',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center'
                            }}
                            icon={<EditOutlined style={{ fontSize: '18px' }} />}
                            title={isDrawing ? 'Click to finish/cancel drawing' : 'Click to start drawing areas'}
                        />
                    </div>

                    {/* Keyboard Shortcuts Panel */}
                    <div style={{
                        opacity: isDrawing ? 1 : 0,
                        visibility: isDrawing ? 'visible' : 'hidden',
                        transform: isDrawing ? 'translateY(0)' : 'translateY(-10px)',
                        transition: 'opacity 0.3s ease-in-out, visibility 0.3s ease-in-out, transform 0.3s ease-in-out',
                        backgroundColor: 'rgba(255, 255, 255, 0.8)',
                        color: '#000000',
                        padding: '8px 12px',
                        borderRadius: 6,
                        fontSize: '12px',
                        lineHeight: '1.4',
                        boxShadow: '0 2px 8px rgba(0, 0, 0, 0.2)',
                        minWidth: '200px',
                        textAlign: 'left'
                    }}>
                        <div style={{ fontWeight: 'bold', marginBottom: 4 }}>Keyboard Shortcuts:</div>
                        <div style={{ marginBottom: 2, display: 'flex', alignItems: 'center', gap: 6 }}>
                            <kbd style={{
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                padding: '2px 4px',
                                borderRadius: 3,
                                fontSize: '11px',
                            }}>Enter</kbd>
                            <span>Finish drawing</span>
                        </div>
                        <div style={{ marginBottom: 2, display: 'flex', alignItems: 'center', gap: 6 }}>
                            <kbd style={{
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                padding: '2px 4px',
                                borderRadius: 3,
                                fontSize: '11px'
                            }}>Esc</kbd>
                            <span>Cancel drawing</span>
                        </div>
                        <div style={{ marginBottom: 2, display: 'flex', alignItems: 'center', gap: 6 }}>
                            <kbd style={{
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                padding: '2px 4px',
                                borderRadius: 3,
                                fontSize: '11px'
                            }}>Backspace</kbd>
                            <span>Undo last point</span>
                        </div>
                    </div>
                </div>

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
                                zIndex: 999,
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                cursor: 'default',
                                pointerEvents: 'auto'
                            }}
                            onClick={(e) => e.stopPropagation()}
                        />

                        <div
                            style={{
                                position: 'fixed',
                                left: getTempAreaCompleteTooltipPosition().left,
                                top: getTempAreaCompleteTooltipPosition().top,
                                zIndex: 1000,
                                background: '#ffffff',
                                border: '1px solid #d9d9d9',
                                borderRadius: 8,
                                boxShadow: '0 6px 16px rgba(0, 0, 0, 0.12)',
                                padding: 12,
                                minWidth: 240,
                                maxWidth: 280,
                                pointerEvents: 'auto'
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
                                marginBottom: 5,
                                fontSize: 14,
                                color: '#262626',
                                paddingRight: 20,
                                textAlign: 'left'
                            }}>
                                Customize Area
                            </div>

                            {/* Area Name Input */}
                            <div style={{ marginBottom: 8 }}>
                                <label style={{
                                    display: 'block',
                                    marginBottom: 6,
                                    fontSize: 12,
                                    fontWeight: 500,
                                    color: '#595959',
                                    textAlign: 'left'
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
                            <div style={{ marginBottom: 8 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                    <label style={{
                                        fontSize: 12,
                                        fontWeight: 500,
                                        color: '#595959',
                                        minWidth: 'fit-content'
                                    }}>
                                        Area Color:
                                    </label>
                                    <ColorPicker
                                        value={areaColor}
                                        onChange={(color) => setAreaColor(color.toHexString())}
                                        size="small"
                                    />
                                </div>
                            </div>

                            {/* Action Buttons */}
                            <div style={{ display: 'flex', gap: 5, justifyContent: 'flex-end' }}>
                                <Button
                                    size="small"
                                    color="pink"
                                    variant="solid"
                                    onClick={handleAreaTooltipSave}
                                >
                                    Save Area
                                </Button>
                            </div>
                        </div>
                    </>
                )}

                {/* Minimap */}
                {(minimapVisible || minimapAnimating) && selectedSamples.length > 0 && imageSizes[selectedSamples[0]?.id] && (
                    <div
                        style={{
                            position: 'absolute',
                            bottom: 10,
                            right: 10,
                            width: 200,
                            height: 150,
                            zIndex: 15,
                            backgroundColor: 'rgba(255, 255, 255, 0.9)',
                            border: '2px solid #d9d9d9',
                            borderRadius: 8,
                            boxShadow: '0 4px 12px rgba(0, 0, 0, 0.15)',
                            overflow: 'hidden',
                            cursor: 'pointer',
                            opacity: minimapVisible ? 1 : 0,
                            transition: 'opacity 0.3s ease-in-out, transform 0.3s ease-in-out',
                            transform: minimapVisible ? 'translateY(0) scale(1)' : 'translateY(20px) scale(0.95)',
                            pointerEvents: minimapVisible ? 'auto' : 'none'
                        }}
                        ref={minimapRef}
                        onClick={handleMinimapClick}
                    >
                        {/* Minimap background image */}
                        <img
                            src={`/${selectedSamples[0].id}_full.jpg`}
                            alt="Minimap"
                            style={{
                                width: '100%',
                                height: '100%',
                                objectFit: 'cover',
                                display: 'block'
                            }}
                            draggable={false}
                        />

                        {/* Viewport indicator */}
                        {(() => {
                            const viewportBounds = getMinimapViewportBounds();
                            if (!viewportBounds) return null;

                            const left = viewportBounds.left * 100;
                            const top = viewportBounds.top * 100;
                            const width = (viewportBounds.right - viewportBounds.left) * 100;
                            const height = (viewportBounds.bottom - viewportBounds.top) * 100;

                            return (
                                <div
                                    style={{
                                        position: 'absolute',
                                        left: `${Math.max(0, Math.min(100, left))}%`,
                                        top: `${Math.max(0, Math.min(100, top))}%`,
                                        width: `${Math.max(0, Math.min(100 - left, width))}%`,
                                        height: `${Math.max(0, Math.min(100 - top, height))}%`,
                                        border: '2px solid #1890ff',
                                        backgroundColor: 'rgba(24, 144, 255, 0.2)',
                                        pointerEvents: 'none',
                                        boxSizing: 'border-box'
                                    }}
                                />
                            );
                        })()}

                        {/* Close button */}
                        <div
                            style={{
                                position: 'absolute',
                                top: 4,
                                right: 4,
                                width: 18,
                                height: 18,
                                backgroundColor: 'rgba(255, 255, 255, 0.9)',
                                border: '1px solid #d9d9d9',
                                borderRadius: 3,
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center',
                                cursor: 'pointer',
                                fontSize: 12,
                                color: '#666'
                            }}
                            onClick={(e) => {
                                e.stopPropagation();
                                toggleMinimapVisible();
                            }}
                            onMouseEnter={(e) => e.target.style.backgroundColor = 'rgba(245, 245, 245, 0.9)'}
                            onMouseLeave={(e) => e.target.style.backgroundColor = 'rgba(255, 255, 255, 0.9)'}
                        >
                            <CloseOutlined style={{ fontSize: 8, color: '#666' }} />
                        </div>
                    </div>
                )}

                {/* Area Edit/Delete Popup */}
                {isAreaEditPopupVisible && selectedAreaForEdit && (
                    <>
                        {/* Overlay to prevent interactions with the map */}
                        <div
                            style={{
                                position: 'absolute',
                                top: 0,
                                left: 0,
                                right: 0,
                                bottom: 0,
                                zIndex: 999,
                                backgroundColor: 'rgba(0, 0, 0, 0.1)',
                                cursor: 'default',
                                pointerEvents: 'auto'
                            }}
                            onClick={(e) => {
                                e.stopPropagation();
                                handleAreaEditCancel();
                            }}
                        />

                        <div
                            style={{
                                position: 'fixed',
                                left: editPopupPosition.x,
                                top: editPopupPosition.y,
                                zIndex: 1000,
                                background: '#ffffff',
                                border: '1px solid #d9d9d9',
                                borderRadius: 8,
                                boxShadow: '0 6px 16px rgba(0, 0, 0, 0.12)',
                                padding: 12,
                                minWidth: 240,
                                maxWidth: 280,
                                pointerEvents: 'auto'
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
                                onClick={handleAreaEditCancel}
                                onMouseEnter={(e) => e.target.style.backgroundColor = '#f5f5f5'}
                                onMouseLeave={(e) => e.target.style.backgroundColor = 'transparent'}
                            >
                                <CloseOutlined style={{ fontSize: 12, color: '#666' }} />
                            </div>

                            {/* Title */}
                            <div style={{
                                fontWeight: 'bold',
                                marginBottom: 5,
                                fontSize: 14,
                                color: '#262626',
                                paddingRight: 20,
                                textAlign: 'left'
                            }}>
                                Edit Area
                            </div>

                            {/* Area Name Input */}
                            <div style={{ marginBottom: 8 }}>
                                <label style={{
                                    display: 'block',
                                    marginBottom: 6,
                                    fontSize: 12,
                                    fontWeight: 500,
                                    color: '#595959',
                                    textAlign: 'left'
                                }}>
                                    Area Name:
                                </label>
                                <Input
                                    value={editAreaName}
                                    onChange={(e) => setEditAreaName(e.target.value)}
                                    placeholder="Enter area name"
                                    maxLength={50}
                                    size="small"
                                />
                            </div>

                            {/* Color Picker */}
                            <div style={{ marginBottom: 8 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                    <label style={{
                                        fontSize: 12,
                                        fontWeight: 500,
                                        minWidth: '70px',
                                        textAlign: 'left',
                                        color: '#595959',
                                    }}>
                                        Area Color:
                                    </label>
                                    <ColorPicker
                                        value={editAreaColor}
                                        onChange={(color) => setEditAreaColor(color.toHexString())}
                                        size="small"
                                    />
                                </div>
                            </div>

                            {/* Neighbors Input */}
                            <div style={{ marginBottom: 8 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                    <label style={{
                                        fontSize: 12,
                                        fontWeight: 500,
                                        color: '#595959',
                                        minWidth: '70px',
                                        textAlign: 'left'
                                    }}>
                                        Neighbors:
                                    </label>
                                    <AutoComplete
                                        value={editNeighbors.toString()}
                                        onChange={(value) => setEditNeighbors(parseInt(value) || 10)}
                                        options={[
                                            { value: '5' },
                                            { value: '10' },
                                            { value: '15' },
                                            { value: '20' },
                                            { value: '25' },
                                            { value: '30' }
                                        ]}
                                        size="small"
                                        style={{ flex: 1 }}
                                        placeholder="10"
                                    />
                                </div>
                            </div>

                            {/* N PCAs Input */}
                            <div style={{ marginBottom: 8 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                    <label style={{
                                        fontSize: 12,
                                        fontWeight: 500,
                                        color: '#595959',
                                        minWidth: '70px',
                                        textAlign: 'left'
                                    }}>
                                        N PCAs:
                                    </label>
                                    <AutoComplete
                                        value={editNPcas.toString()}
                                        onChange={(value) => setEditNPcas(parseInt(value) || 50)}
                                        options={[
                                            { value: '10' },
                                            { value: '20' },
                                            { value: '30' },
                                            { value: '40' },
                                            { value: '50' },
                                            { value: '75' },
                                            { value: '100' }
                                        ]}
                                        size="small"
                                        style={{ flex: 1 }}
                                        placeholder="50"
                                    />
                                </div>
                            </div>

                            {/* Resolution Input */}
                            <div style={{ marginBottom: 8 }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                                    <label style={{
                                        fontSize: 12,
                                        fontWeight: 500,
                                        color: '#595959',
                                        minWidth: '70px',
                                        textAlign: 'left'
                                    }}>
                                        Resolution:
                                    </label>
                                    <AutoComplete
                                        value={editResolution.toString()}
                                        onChange={(value) => setEditResolution(parseFloat(value) || 0.5)}
                                        options={[
                                            { value: '0.1' },
                                            { value: '0.2' },
                                            { value: '0.3' },
                                            { value: '0.4' },
                                            { value: '0.5' },
                                            { value: '0.6' },
                                            { value: '0.7' },
                                            { value: '0.8' },
                                            { value: '0.9' },
                                            { value: '1.0' }
                                        ]}
                                        size="small"
                                        style={{ flex: 1 }}
                                        placeholder="0.5"
                                    />
                                </div>
                            </div>

                            {/* Action Buttons */}
                            <div style={{ display: 'flex', gap: 5, justifyContent: 'space-between' }}>
                                <Button
                                    size="small"
                                    danger
                                    onClick={handleAreaDelete}
                                    style={{ flex: 1 }}
                                >
                                    Delete
                                </Button>
                                <Button
                                    size="small"
                                    type="primary"
                                    onClick={handleAreaEditSave}
                                    style={{ flex: 1 }}
                                >
                                    Save
                                </Button>
                            </div>
                        </div>
                    </>
                )}
            </div>
        </div>
    );
};