import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { GeneSettings } from './GeneList';
import { CellSettings } from './CellList';
import { Collapse, Radio, Button, Input, ColorPicker, AutoComplete, Spin, Switch, message } from "antd";
import { CloseOutlined, EditOutlined, RedoOutlined, BorderOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, PolygonLayer, LineLayer } from '@deck.gl/layers';
import { convertHEXToRGB, COLOR_PALETTE, getSequentialColor } from './Utils';


export const SampleViewer = ({
    selectedSamples,
    coordinatesData,
    cellTypesData,
    selectedCellTypes,
    setSelectedCellTypes,
    cellTypeColors,
    setCellTypeColors,
    umapDataSets,
    setUmapDataSets,
    umapLoading,
    setUmapLoading,
    hoveredCluster,
    onImagesLoaded,
    kosaraDisplayEnabled = true,
    trajectoryGenes = [],
    trajectoryGenesSample = null
}) => {
    const containerRef = useRef(null);
    const lastLoadedTrajectoryRef = useRef(null); // Track the last loaded trajectory gene combination to prevent redundant API calls
    const viewStatePendingRef = useRef(null);
    const viewStateRafRef = useRef(null);
    const kosaraLoadingSamplesRef = useRef({});
    const fetchingImages = useRef(new Set()); // Track which images are currently being fetched
    const imagesLoadedCallbackCalled = useRef(false); // Track if callback has been called for current samples

    const [mainViewState, setMainViewState] = useState(null);
    const [containerSize, setContainerSize] = useState({ width: 0, height: 0 });
    const [imageSizes, setImageSizes] = useState({});
    const [availableGenes, setAvailableGenes] = useState([]); // All genes that have been added to the list
    const [selectedGenes, setSelectedGenes] = useState([]); // Currently selected (checked) genes
    const [geneColorMap, setGeneColorMap] = useState({}); // { geneName: '#RRGGBB' }

    const [radioCellGeneModes, setRadioCellGeneModes] = useState(
        selectedSamples.reduce((acc, sample) => ({ ...acc, [sample.id]: 'cellTypes' }), {})
    );

    // Previous modes for each sample when Kosara display is toggled off
    const [previousModes, setPreviousModes] = useState(
        selectedSamples.reduce((acc, sample) => ({ ...acc, [sample.id]: 'cellTypes' }), {})
    );

    // Previous gene selections when Kosara display is toggled off
    const [previousGeneSelections, setPreviousGeneSelections] = useState({});

    const [previousKosaraData, setPreviousKosaraData] = useState({});

    // Kosara gene expression data per sample returned from backend
    const [kosaraDataBySample, setKosaraDataBySample] = useState({}); // { sampleId: [ { id, cell_x, cell_y, cell_type, total_expression, angles:{}, radius:{}, ratios:{} }, ... ] }

    // Loading state for gene mode switching
    const [isKosaraLoading, setIsKosaraLoading] = useState(false);

    // Track which samples currently have an in-flight kosara request
    const [kosaraLoadingSamples, setKosaraLoadingSamples] = useState({}); // { sampleId: true }

    // Single gene expression data per sample for sequential coloring
    const [singleGeneDataBySample, setSingleGeneDataBySample] = useState({}); // { sampleId: { cells: [...], min_expression: ..., max_expression: ... } }

    // Drawing state
    const [isDrawing, setIsDrawing] = useState(false);
    const [drawingPoints, setDrawingPoints] = useState([]);
    const [customAreas, setCustomAreas] = useState([]);
    const [currentDrawingSample, setCurrentDrawingSample] = useState(null);
    const [mousePosition, setMousePosition] = useState(null);
    const [hoveredCell, setHoveredCell] = useState(null);

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
    const [editNPcas, setEditNPcas] = useState(30);
    const [editResolutions, setEditResolutions] = useState(1);

    // Minimap state
    const [minimapVisible, setMinimapVisible] = useState(true);
    const [minimapAnimating, setMinimapAnimating] = useState(false);
    const minimapRef = useRef(null);

    // High-definition magnifying glass state
    const [magnifierVisible, setMagnifierVisible] = useState(false);
    const [magnifierData, setMagnifierData] = useState(null);
    const [magnifierViewport, setMagnifierViewport] = useState({ x: 0.5, y: 0.5, size: 200 });
    const [magnifierMousePos, setMagnifierMousePos] = useState({ x: 0, y: 0 });
    const [keyPressed, setKeyPressed] = useState(false);

    const magnifierRef = useRef(null);

    // Track previous image URL to only reset loading state when it truly changes
    const prevMagnifierUrlRef = useRef(null);
    // Track if current magnifier image has finished loading
    const [magnifierImageLoaded, setMagnifierImageLoaded] = useState(false);
    const [magnifierImageVersion, setMagnifierImageVersion] = useState(0);

    // Add state for preloaded high-res images
    const [hiresImages, setHiresImages] = useState({}); // { sampleId: imageUrl }


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
    const mainView = useMemo(() => new OrthographicView({
        id: 'main',
        controller: true
    }), []);

    // Stable controller instance to avoid re-initialization flashes
    const deckController = useMemo(() => {
        if (isAreaTooltipVisible || isAreaEditPopupVisible) return false;
        if (!isDrawing) {
            // Use DeckGL defaults to avoid controller re-inits that may cause flashes
            return true;
        }
        return {
            dragPan: false,
            dragRotate: false,
            doubleClickZoom: false,
            scrollZoom: true
        };
    }, [isAreaTooltipVisible, isAreaEditPopupVisible, isDrawing]);

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
            const cellData = coordinatesData && coordinatesData[sample.id] ? coordinatesData[sample.id] : [];
            const offset = sampleOffsets && sampleOffsets[sample.id] ? sampleOffsets[sample.id] : [0, 0];

            acc[sample.id] = cellData.map(cell => ({
                ...cell,
                x: cell.cell_x + offset[0],
                y: cell.cell_y + offset[1]
            }));

            return acc;
        }, {});
    }, [selectedSamples, coordinatesData, sampleOffsets]);

    // Function to load Kosara data for a specific sample
    const loadKosaraDataForSample = async (sampleId, genes) => {
        if (!genes || genes.length === 0) return;

        // Check if this sample is already loading to prevent duplicate requests
        if (kosaraLoadingSamples[sampleId]) {
            return;
        }

        try {
            setIsKosaraLoading(true);
            setKosaraLoadingSamples(prev => ({ ...prev, [sampleId]: true }));

            const response = await fetch('/api/get_kosara_data', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    sample_ids: [sampleId],
                    gene_list: genes
                })
            });

            if (!response.ok) {
                console.error('Failed to fetch Kosara data:', response.status, response.statusText);
                return;
            }

            const data = await response.json();
            if (data && data[sampleId]) {
                setKosaraDataBySample(prev => ({
                    ...prev,
                    [sampleId]: Array.isArray(data[sampleId]) ? data[sampleId] : []
                }));

                // Set the sample to gene mode
                setRadioCellGeneModes(prev => ({ ...prev, [sampleId]: 'genes' }));
            }
        } catch (err) {
            console.error('Error fetching Kosara data:', err);
        } finally {
            setKosaraLoadingSamples(prev => {
                const next = { ...prev };
                delete next[sampleId];
                if (Object.keys(next).length === 0) {
                    setIsKosaraLoading(false);
                }
                return next;
            });
        }
    };

    // Function to load single gene expression data for sequential coloring
    const loadSingleGeneDataForSample = async (sampleId, geneName) => {
        if (!geneName) return;

        try {
            setIsKosaraLoading(true);
            setKosaraLoadingSamples(prev => ({ ...prev, [sampleId]: true }));

            const response = await fetch('/api/get_single_gene_expression', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    sample_ids: [sampleId],
                    gene_name: geneName
                })
            });

            if (!response.ok) {
                console.error('Failed to fetch single gene expression data:', response.status, response.statusText);
                return;
            }

            const data = await response.json();
            if (data && data[sampleId]) {
                setSingleGeneDataBySample(prev => ({
                    ...prev,
                    [sampleId]: data[sampleId]
                }));

                // Set the sample to gene mode
                setRadioCellGeneModes(prev => ({ ...prev, [sampleId]: 'genes' }));
            }
        } catch (err) {
            console.error('Error fetching single gene expression data:', err);
        } finally {
            setKosaraLoadingSamples(prev => {
                const next = { ...prev };
                delete next[sampleId];
                if (Object.keys(next).length === 0) {
                    setIsKosaraLoading(false);
                }
                return next;
            });
        }
    };

    // Change cell/gene mode for a sample
    const changeCellGeneMode = (sampleId, e) => {
        const newMode = e.target.value;

        // Update previous modes when manually changing mode
        setPreviousModes(prev => ({ ...prev, [sampleId]: radioCellGeneModes[sampleId] }));

        // If switching TO genes with existing kosara data, show spinner first, then defer the expensive mode switch
        if (newMode === 'genes' && kosaraDataBySample[sampleId]?.length) {
            // Start spinner immediately so it can paint before heavy polygon generation
            setIsKosaraLoading(true);
            // Provide a short fallback auto-hide ONLY if no real fetch starts (no sample in kosaraLoadingSamples)
            if (spinnerFallbackTimeoutRef.current) clearTimeout(spinnerFallbackTimeoutRef.current);
            spinnerFallbackTimeoutRef.current = setTimeout(() => {
                if (Object.keys(kosaraLoadingSamplesRef.current).length === 0) {
                    setIsKosaraLoading(false);
                }
            }, 400);
            // Defer the mode change to next animation frame (after spinner paints)
            const schedule = typeof requestAnimationFrame === 'function' ? requestAnimationFrame : (fn) => setTimeout(fn, 0);
            schedule(() => {
                setRadioCellGeneModes(prev => ({ ...prev, [sampleId]: newMode }));
            });
        } else {
            // Normal immediate mode change (cellTypes or genes without cached data)
            setRadioCellGeneModes(prev => ({ ...prev, [sampleId]: newMode }));
        }
    };

    // Receive data from GeneList and store & switch to gene mode
    const handleKosaraData = useCallback((sampleId, dataArray, dataType = 'kosara') => {
        if (dataType === 'single_gene') {
            // Handle single gene expression data
            setSingleGeneDataBySample(prev => ({
                ...prev,
                [sampleId]: dataArray
            }));
            // Clear kosara data for this sample since we're using single gene mode
            setKosaraDataBySample(prev => {
                const updated = { ...prev };
                delete updated[sampleId];
                return updated;
            });
        } else {
            // Handle kosara data
            setKosaraDataBySample(prev => ({
                ...prev,
                [sampleId]: Array.isArray(dataArray) ? dataArray : []
            }));
            // Clear single gene data for this sample since we're using kosara mode
            setSingleGeneDataBySample(prev => {
                const updated = { ...prev };
                delete updated[sampleId];
                return updated;
            });
        }

        setRadioCellGeneModes(prev => ({ ...prev, [sampleId]: 'genes' }));
        // Mark this sample's request complete
        setKosaraLoadingSamples(prev => {
            const next = { ...prev };
            delete next[sampleId];
            // If no more in-flight requests, stop spinner
            if (Object.keys(next).length === 0) {
                setIsKosaraLoading(false);
            }
            return next;
        });
    }, []);

    // Start loading when confirm button is clicked (before fetching data)
    const handleKosaraLoadingStart = useCallback((sampleId) => {
        // Clear any fallback hide for cached-mode switches
        if (spinnerFallbackTimeoutRef.current) clearTimeout(spinnerFallbackTimeoutRef.current);
        setKosaraLoadingSamples(prev => ({ ...prev, [sampleId]: true }));
        setIsKosaraLoading(true);
    }, []);

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

    // Kosara path generation utilities
    const generateCirclePoints = useCallback((cx, cy, r, steps = 50) => {
        const points = [];
        const angleStep = (2 * Math.PI) / steps;
        for (let i = 0; i < steps; i++) {
            const theta = i * angleStep;
            points.push([cx + r * Math.cos(theta), cy + r * Math.sin(theta)]);
        }
        points.push(points[0]);
        return points;
    }, []);

    const generateSingleArcPoints = useCallback((startX, startY, endX, endY, r, largeArcFlag, sweepFlag) => {
        const dx = endX - startX;
        const dy = endY - startY;
        const d = Math.hypot(dx, dy);
        if (d > 2 * r) {
            return [];
        }

        const midX = (startX + endX) / 2;
        const midY = (startY + endY) / 2;
        const h = Math.sqrt(r * r - (d / 2) * (d / 2));

        const ux = -dy / d;
        const uy = dx / d;

        const cx1 = midX + h * ux;
        const cy1 = midY + h * uy;
        const cx2 = midX - h * ux;
        const cy2 = midY - h * uy;

        const computeAngles = (cx, cy) => {
            const startAngle = Math.atan2(startY - cy, startX - cx);
            const endAngle = Math.atan2(endY - cy, endX - cx);
            let delta = endAngle - startAngle;
            delta = ((delta % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
            return { startAngle, endAngle, delta };
        };

        const cand1 = computeAngles(cx1, cy1);
        const cand2 = computeAngles(cx2, cy2);

        const effectiveDelta = candidate => (sweepFlag === 1 ? candidate.delta : (2 * Math.PI - candidate.delta));
        const eff1 = effectiveDelta(cand1);
        const eff2 = effectiveDelta(cand2);

        let chosen, cx, cy;
        if (largeArcFlag === 0) {
            if (eff1 <= Math.PI && eff2 > Math.PI) {
                chosen = cand1; cx = cx1; cy = cy1;
            } else if (eff2 <= Math.PI && eff1 > Math.PI) {
                chosen = cand2; cx = cx2; cy = cy2;
            } else {
                chosen = cand1; cx = cx1; cy = cy1;
            }
        } else {
            if (eff1 >= Math.PI && eff2 < Math.PI) {
                chosen = cand1; cx = cx1; cy = cy1;
            } else if (eff2 >= Math.PI && eff1 < Math.PI) {
                chosen = cand2; cx = cx2; cy = cy2;
            } else {
                chosen = cand1; cx = cx1; cy = cy1;
            }
        }

        const deltaEffective = sweepFlag === 1 ? effectiveDelta(chosen) : -effectiveDelta(chosen);

        const steps = 20;
        const points = [];
        for (let i = 0; i <= steps; i++) {
            const t = i / steps;
            const theta = chosen.startAngle + t * deltaEffective;
            const x = cx + r * Math.cos(theta);
            const y = cy + r * Math.sin(theta);
            points.push([x, y]);
        }
        return points;
    }, []);

    const generateComplexArcPoints = useCallback((startX, startY, endX, endY, outerRadius, innerRadius, outerFlags = { large: 0, sweep: 1 }, innerFlags = { large: 0, sweep: 0 }) => {
        const outerArc = generateSingleArcPoints(startX, startY, endX, endY, outerRadius, outerFlags.large, outerFlags.sweep);
        const innerArc = generateSingleArcPoints(endX, endY, startX, startY, innerRadius, innerFlags.large, innerFlags.sweep);
        return [...outerArc, ...innerArc];
    }, [generateSingleArcPoints]);

    const generateKosaraPath = useCallback((pointX, pointY, angles, ratios, cal_radius) => {
        const baseRadius = 5;
        const paths = [];
        const cellTypes = selectedGenes;

        let startpointX, startpointY, endpointX, endpointY;
        let lastStartPointX, lastStartPointY, lastEndPointX, lastEndPointY, lastCircleRadius = 0;
        const originalPointX = pointX - baseRadius * Math.cos(45 * Math.PI / 180);
        const originalPointY = pointY + baseRadius * Math.sin(45 * Math.PI / 180);

        const cellIndices = ratios
            .filter(item => item[1] !== 0 && cellTypes.includes(item[0]))
            .sort((a, b) => cellTypes.indexOf(a[0]) - cellTypes.indexOf(b[0]))
            .slice(0, 9)
            .map(item => item[0]);

        // we no longer compute colors here; defer to layer accessor using gene ids
        let cellAngles = cellIndices.map(index => angles.find(item => item[0] === index));
        let cellRadius = cellIndices.map(index => cal_radius.find(item => item[0] === index));

        const ratioSum = ratios.reduce((acc, item) => acc + item[1], 0);

        if (cellAngles.length === 0) {
            const circlePoints = generateCirclePoints(pointX, pointY, baseRadius, 50);
            paths.push({ path: circlePoints, color: '#FFFFFF' });
        } else {
            cellAngles = cellAngles.map(angle => [angle[0], angle[1]]);
            cellRadius = cellRadius.map(rad => [rad[0], rad[1]]);

            cellAngles.forEach((angle, index) => {
                const cal_cell_radius = cellRadius[index][1];
                let points = [];

                startpointX = originalPointX + Math.abs(cal_cell_radius * Math.cos((angle[1] + 45) * Math.PI / 180));
                startpointY = originalPointY - Math.abs(cal_cell_radius * Math.sin((angle[1] + 45) * Math.PI / 180));
                endpointX = originalPointX + Math.abs(cal_cell_radius * Math.cos((angle[1] - 45) * Math.PI / 180));
                endpointY = originalPointY - Math.abs(cal_cell_radius * Math.sin((angle[1] - 45) * Math.PI / 180));

                if (index === 0) {
                    const isLargeArcInner = cal_cell_radius > Math.sqrt(3) * baseRadius;
                    points = generateComplexArcPoints(
                        startpointX,
                        startpointY,
                        endpointX,
                        endpointY,
                        cal_cell_radius,
                        baseRadius,
                        { large: 0, sweep: 1 },
                        { large: isLargeArcInner ? 1 : 0, sweep: 1 }
                    );
                } else if (index === cellAngles.length - 1 && ratioSum === 1) {
                    const isLargeArcInner = lastCircleRadius <= Math.sqrt(3) * baseRadius;
                    points = generateComplexArcPoints(
                        lastStartPointX,
                        lastStartPointY,
                        lastEndPointX,
                        lastEndPointY,
                        lastCircleRadius,
                        baseRadius,
                        { large: 0, sweep: 1 },
                        { large: isLargeArcInner ? 1 : 0, sweep: 0 }
                    );
                } else {
                    const pointsSegment1 = generateSingleArcPoints(
                        lastStartPointX, lastStartPointY,
                        lastEndPointX, lastEndPointY,
                        lastCircleRadius,
                        0,
                        1
                    );

                    const pointsSegment2 = generateSingleArcPoints(
                        lastEndPointX, lastEndPointY,
                        endpointX, endpointY,
                        baseRadius,
                        0,
                        0
                    );

                    const pointsSegment3 = generateSingleArcPoints(
                        endpointX, endpointY,
                        startpointX, startpointY,
                        cal_cell_radius,
                        0,
                        0
                    );

                    const pointsSegment4 = generateSingleArcPoints(
                        startpointX, startpointY,
                        lastStartPointX, lastStartPointY,
                        baseRadius,
                        0,
                        0
                    );
                    points = [...pointsSegment1, ...pointsSegment2, ...pointsSegment3, ...pointsSegment4];
                }

                paths.push({
                    path: points,
                    gene: cellIndices[index]
                });

                lastCircleRadius = cal_cell_radius;
                lastStartPointX = startpointX;
                lastStartPointY = startpointY;
                lastEndPointX = endpointX;
                lastEndPointY = endpointY;
            });

            if (ratioSum < 1) {
                const isLargeArcInner = lastCircleRadius <= Math.sqrt(3) * baseRadius;
                const points = generateComplexArcPoints(
                    lastStartPointX,
                    lastStartPointY,
                    lastEndPointX,
                    lastEndPointY,
                    lastCircleRadius,
                    baseRadius,
                    { large: 0, sweep: 1 },
                    { large: isLargeArcInner ? 1 : 0, sweep: 0 }
                );
                paths.push({ path: points, color: '#333333' });
            }
        }
        return paths;
    }, [generateCirclePoints, generateComplexArcPoints, generateSingleArcPoints, selectedGenes]);

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

        // Calculate total world bounds for all samples
        let totalLeft = Infinity;
        let totalRight = -Infinity;
        let totalTop = Infinity;
        let totalBottom = -Infinity;

        selectedSamples.forEach(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            if (imageSize) {
                totalLeft = Math.min(totalLeft, offset[0]);
                totalRight = Math.max(totalRight, offset[0] + imageSize[0]);
                totalTop = Math.min(totalTop, offset[1]);
                totalBottom = Math.max(totalBottom, offset[1] + imageSize[1]);
            }
        });

        if (totalLeft === Infinity) return null;

        const totalWidth = totalRight - totalLeft;
        const totalHeight = totalBottom - totalTop;

        // Convert to relative coordinates within the total combined area
        const relativeBounds = {
            left: Math.max(0, (viewportBounds.left - totalLeft) / totalWidth),
            right: Math.min(1, (viewportBounds.right - totalLeft) / totalWidth),
            top: Math.max(0, (viewportBounds.top - totalTop) / totalHeight),
            bottom: Math.min(1, (viewportBounds.bottom - totalTop) / totalHeight)
        };

        return relativeBounds;
    }, [mainViewState, containerSize, selectedSamples, imageSizes, sampleOffsets]);

    // Function to detect which sample was clicked in the minimap
    const getSampleFromMinimapClick = useCallback((clickX, clickY) => {
        if (!selectedSamples.length) return null;

        // Calculate total world bounds for all samples
        let totalLeft = Infinity;
        let totalRight = -Infinity;
        let totalTop = Infinity;
        let totalBottom = -Infinity;

        selectedSamples.forEach(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            if (imageSize) {
                totalLeft = Math.min(totalLeft, offset[0]);
                totalRight = Math.max(totalRight, offset[0] + imageSize[0]);
                totalTop = Math.min(totalTop, offset[1]);
                totalBottom = Math.max(totalBottom, offset[1] + imageSize[1]);
            }
        });

        if (totalLeft === Infinity) return null;

        const totalWidth = totalRight - totalLeft;
        const totalHeight = totalBottom - totalTop;

        // Check which sample contains the click coordinates
        for (const sample of selectedSamples) {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];

            if (!imageSize) continue;

            // Calculate relative position and size within the total bounds
            const relativeLeft = (offset[0] - totalLeft) / totalWidth;
            const relativeTop = (offset[1] - totalTop) / totalHeight;
            const relativeWidth = imageSize[0] / totalWidth;
            const relativeHeight = imageSize[1] / totalHeight;

            // Check if click is within this sample's bounds
            if (clickX >= relativeLeft && clickX <= relativeLeft + relativeWidth &&
                clickY >= relativeTop && clickY <= relativeTop + relativeHeight) {
                return sample;
            }
        }

        return null;
    }, [selectedSamples, imageSizes, sampleOffsets]);

    // Handle minimap click to navigate
    const handleMinimapClick = useCallback((event) => {
        if (!minimapRef.current || !selectedSamples.length) return;

        const rect = minimapRef.current.getBoundingClientRect();
        const x = (event.clientX - rect.left) / rect.width;
        const y = (event.clientY - rect.top) / rect.height;

        // First, check if we clicked on a specific sample
        const clickedSample = getSampleFromMinimapClick(x, y);

        if (clickedSample) {
            // Center the view on the clicked sample
            const imageSize = imageSizes[clickedSample.id];
            const offset = sampleOffsets[clickedSample.id] || [0, 0];

            if (imageSize) {
                const centerX = offset[0] + imageSize[0] / 2;
                const centerY = offset[1] + imageSize[1] / 2;

                setMainViewState(prev => ({
                    ...prev,
                    target: [centerX, centerY, 0]
                }));
                return;
            }
        }

        // Fallback to original behavior: pan to clicked position
        // Calculate total world bounds for all samples
        let totalLeft = Infinity;
        let totalRight = -Infinity;
        let totalTop = Infinity;
        let totalBottom = -Infinity;

        selectedSamples.forEach(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            if (imageSize) {
                totalLeft = Math.min(totalLeft, offset[0]);
                totalRight = Math.max(totalRight, offset[0] + imageSize[0]);
                totalTop = Math.min(totalTop, offset[1]);
                totalBottom = Math.max(totalBottom, offset[1] + imageSize[1]);
            }
        });

        if (totalLeft === Infinity) return;

        const totalWidth = totalRight - totalLeft;
        const totalHeight = totalBottom - totalTop;

        // Convert relative coordinates to world coordinates
        const worldX = totalLeft + x * totalWidth;
        const worldY = totalTop + y * totalHeight;

        // Update main view to center on clicked position
        setMainViewState(prev => ({
            ...prev,
            target: [worldX, worldY, 0]
        }));
    }, [selectedSamples, imageSizes, sampleOffsets, getSampleFromMinimapClick]);

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

    // Update magnifier viewport based on mouse position
    const updateMagnifierViewport = useCallback((worldX, worldY, sampleId) => {
        if (!sampleId) return;

        const offset = sampleOffsets[sampleId] || [0, 0];
        const imageSize = imageSizes[sampleId];

        if (!imageSize) return;

        // Convert world coordinates to image coordinates
        const imageX = worldX - offset[0];
        const imageY = worldY - offset[1];

        // Calculate viewport position (as percentage of image)
        const viewportX = imageX / imageSize[0];
        const viewportY = imageY / imageSize[1];

        const clampedViewport = {
            x: Math.max(0, Math.min(1, viewportX)),
            y: Math.max(0, Math.min(1, viewportY)),
            size: 200 // Fixed viewport size in magnifier pixels
        };

        // Only update if values actually changed to avoid render → hover → setState loops
        setMagnifierViewport(prev => {
            if (
                prev &&
                Math.abs(prev.x - clampedViewport.x) < 1e-4 &&
                Math.abs(prev.y - clampedViewport.y) < 1e-4 &&
                prev.size === clampedViewport.size
            ) {
                return prev;
            }
            return clampedViewport;
        });

        setMagnifierMousePos(prev => {
            if (prev && prev.x === worldX && prev.y === worldY) return prev;
            return { x: worldX, y: worldY };
        });
    }, [sampleOffsets, imageSizes]);

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
                    setEditNPcas(area.n_pcas || 30);
                    setEditResolutions(area.resolutions || 1);

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
                        resolutions: editResolutions
                    }
                    : area
            ));

            // Update corresponding UMAP datasets with new area color and name
            setUmapDataSets(prev => prev.map(dataset => {
                // Check if this dataset corresponds to the edited area by matching the original area name and sample
                if (dataset.areaName === selectedAreaForEdit.name && 
                    dataset.sampleId === selectedAreaForEdit.sampleId) {
                    
                    // If the area name changed, we need to update the adata_umap_title as well
                    let newAdataUmapTitle = dataset.adata_umap_title;
                    if (editAreaName !== selectedAreaForEdit.name) {
                        // Extract the current adata_umap_title and replace the name part
                        const oldFormattedName = selectedAreaForEdit.name.split(' ').join('_');
                        const newFormattedName = editAreaName.split(' ').join('_');
                        newAdataUmapTitle = dataset.adata_umap_title.replace(oldFormattedName, newFormattedName);
                    }
                    
                    return {
                        ...dataset,
                        areaColor: editAreaColor,
                        areaName: editAreaName,
                        adata_umap_title: newAdataUmapTitle,
                        title: `${editAreaName} (${dataset.sampleId})` // Update display title too
                    };
                }
                return dataset;
            }));
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
        setEditNPcas(30);
        setEditResolutions(1);
    };

    // Memoize getSampleAtCoordinate to prevent infinite effect loops
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
        if (isDrawing) {
            if (info.coordinate) {
                setMousePosition(info.coordinate);
            }
        } else {
            // Always track mouse position for magnifier, but only set drawing mousePosition when drawing
            if (info.coordinate) {
                // Update mouse position for magnifier tracking
                setMagnifierMousePos(prev => {
                    const [worldX, worldY] = info.coordinate;
                    if (prev && prev.x === worldX && prev.y === worldY) return prev;
                    return { x: worldX, y: worldY };
                });

                // Handle magnifier when key is pressed
                if (magnifierVisible && !isAreaTooltipVisible && !isAreaEditPopupVisible) {
                    const [worldX, worldY] = info.coordinate;

                    // Determine which sample the mouse is over
                    const hoveredSample = getSampleAtCoordinate(worldX, worldY);

                    if (hoveredSample) {
                        // Update magnifier viewport and mouse position
                        updateMagnifierViewport(worldX, worldY, hoveredSample);
                    }
                }
            }

            // Clear drawing mousePosition when not drawing
            if (mousePosition !== null) setMousePosition(null);
        }
    }, [isDrawing, mousePosition, magnifierVisible, isAreaTooltipVisible, isAreaEditPopupVisible, getSampleAtCoordinate, updateMagnifierViewport]);

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
                    value={radioCellGeneModes && radioCellGeneModes[sample.id] ? radioCellGeneModes[sample.id] : 'cellTypes'}
                    optionType="button"
                    style={{ marginBottom: 10 }}
                    onChange={(e) => changeCellGeneMode(sample.id, e)}
                />

                {(radioCellGeneModes && radioCellGeneModes[sample.id] ? radioCellGeneModes[sample.id] : 'cellTypes') === 'cellTypes' ? (
                    <CellSettings
                        cellTypesData={cellTypesData && cellTypesData[sample.id] ? cellTypesData[sample.id] : []}
                        selectedCellTypes={selectedCellTypes && selectedCellTypes[sample.id] ? selectedCellTypes[sample.id] : []}
                        setSelectedCellTypes={(newSelectedTypes) => {
                            setSelectedCellTypes(prev => ({
                                ...prev,
                                [sample.id]: newSelectedTypes
                            }));
                        }}
                        cellTypeColors={cellTypeColors}
                        setCellTypeColors={setCellTypeColors}
                    />
                ) : (
                    <GeneSettings
                        sampleId={sample.id}
                        availableGenes={availableGenes}
                        setAvailableGenes={setAvailableGenes}
                        selectedGenes={selectedGenes}
                        setSelectedGenes={setSelectedGenes}
                        geneColorMap={geneColorMap}
                        setGeneColorMap={setGeneColorMap}
                        onKosaraData={handleKosaraData}
                        onKosaraLoadingStart={handleKosaraLoadingStart}
                    />
                )}
            </>
        )
    }));

    const handleViewStateChange = useCallback(({ viewState, viewId }) => {
        if (viewId !== 'main') return;
        // Throttle setState to animation frames to avoid flicker during zoom
        const nextState = (prev => ({
            ...viewState,
            maxZoom: prev?.maxZoom ?? 2.5,
            minZoom: prev?.minZoom ?? -5
        }))(mainViewState);
        viewStatePendingRef.current = nextState;
        if (viewStateRafRef.current == null) {
            viewStateRafRef.current = requestAnimationFrame(() => {
                viewStateRafRef.current = null;
                const pending = viewStatePendingRef.current;
                viewStatePendingRef.current = null;
                if (!pending) return;
                setMainViewState(prev => {
                    if (
                        prev &&
                        prev.zoom === pending.zoom &&
                        prev.maxZoom === pending.maxZoom &&
                        prev.minZoom === pending.minZoom &&
                        prev.target && pending.target &&
                        prev.target[0] === pending.target[0] &&
                        prev.target[1] === pending.target[1] &&
                        prev.target[2] === pending.target[2]
                    ) {
                        return prev;
                    }
                    return pending;
                });
            });
        }
    }, [mainViewState]);

    // Helper function to check if two areas are the same
    const arePointsSimilar = (points1, points2, tolerance = 1) => {
        if (!points1 || !points2 || points1.length !== points2.length) {
            return false;
        }
        
        // Check if all points are within tolerance distance
        return points1.every((point1, index) => {
            const point2 = points2[index];
            const dx = point1[0] - point2[0];
            const dy = point1[1] - point2[1];
            const distance = Math.sqrt(dx * dx + dy * dy);
            return distance <= tolerance;
        });
    };

    // Helper function to check for duplicate UMAP
    const findDuplicateUmap = (sampleId, areaPoints, neighbors, nPcas, resolutions) => {
        if (!umapDataSets || umapDataSets.length === 0) return null;
        
        return umapDataSets.find(dataset => {
            // Check if sample ID matches
            if (dataset.sampleId !== sampleId) return false;
            
            // Check if parameters match
            const datasetTitle = dataset.adata_umap_title || '';
            const titleParts = datasetTitle.split('_');
            if (titleParts.length >= 3) {
                const datasetNeighbors = parseInt(titleParts[titleParts.length - 3]) || 0;
                const datasetNPcas = parseInt(titleParts[titleParts.length - 2]) || 0;
                const datasetResolutions = parseFloat(titleParts[titleParts.length - 1]) || 0;
                
                if (datasetNeighbors !== neighbors || 
                    datasetNPcas !== nPcas || 
                    Math.abs(datasetResolutions - resolutions) > 0.001) {
                    return false;
                }
            }
            
            // Check if area points are similar
            return arePointsSimilar(dataset.areaPoints, areaPoints, 5); // 5 pixel tolerance
        });
    };

    const generateUmap = () => {
        if (!selectedAreaForEdit) return;

        // Get cells that are within the selected area BEFORE making the API call
        const sampleCells = coordinatesData && coordinatesData[selectedAreaForEdit.sampleId] ? coordinatesData[selectedAreaForEdit.sampleId] : [];
        const offset = sampleOffsets && sampleOffsets[selectedAreaForEdit.sampleId] ? sampleOffsets[selectedAreaForEdit.sampleId] : [0, 0];

        // Filter cells that are within the drawn polygon
        const cellsInArea = sampleCells.filter(cell => {
            const localX = cell.cell_x;
            const localY = cell.cell_y;

            // Simple point-in-polygon check
            return isPointInPolygon([localX, localY], selectedAreaForEdit.points.map(p => [p[0] - offset[0], p[1] - offset[1]]));
        });

        const cellIdsInArea = cellsInArea.map(cell => cell.id);

        // Check if we have any cells in the selected area
        if (cellIdsInArea.length === 0) {
            message.warning('No cells found in the selected area');
            return;
        }

        // Check for duplicate UMAP before generating
        const duplicateUmap = findDuplicateUmap(
            selectedAreaForEdit.sampleId,
            selectedAreaForEdit.points,
            editNeighbors,
            editNPcas,
            editResolutions
        );

        if (duplicateUmap) {
            message.warning({
                content: `A UMAP with the same parameters and area already exists. Please modify the parameters or select a different area.`,
                duration: 6, // Show for 6 seconds since it's a longer message
            });
            return;
        }

        // Generate a unique ID for this UMAP dataset
        const umapId = `${selectedAreaForEdit.sampleId}_${selectedAreaForEdit.name}_${Date.now()}`;
        const umapTitle = `${selectedAreaForEdit.name} (${selectedAreaForEdit.sampleId})`;

        const name = selectedAreaForEdit.name;
        const formattedName = name.split(' ').join('_');
        const adata_umap_title = `${formattedName}_${selectedAreaForEdit.sampleId}_${editNeighbors}_${editNPcas}_${editResolutions}`;

        // Add a new loading dataset entry
        setUmapDataSets(prev => [
            ...prev,
            {
                id: umapId,
                title: umapTitle,
                adata_umap_title: adata_umap_title,
                data: [],
                loading: true,
                sampleId: selectedAreaForEdit.sampleId,
                areaPoints: selectedAreaForEdit.points,
                areaColor: selectedAreaForEdit.color,
                areaName: selectedAreaForEdit.name
            }
        ]);

        setUmapLoading(true);

        fetch('/api/get_umap_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                sample_id: selectedAreaForEdit.sampleId,
                cell_ids: cellIdsInArea,  // Pass the cell IDs to the backend
                n_neighbors: editNeighbors,
                n_pcas: editNPcas,
                resolutions: editResolutions,
                adata_umap_title: adata_umap_title
            })
        })
            .then(res => res.json())
            .then(data => {
                // No need to filter data since backend only returns data for specified cells
                setUmapDataSets(prev =>
                    prev.map(dataset =>
                        dataset.id === umapId
                            ? { ...dataset, data: data, loading: false }
                            : dataset
                    )
                );
                setUmapLoading(false);
            })
            .catch(error => {
                console.error('Error generating UMAP:', error);

                // Remove the failed dataset entry
                setUmapDataSets(prev => prev.filter(dataset => dataset.id !== umapId));
                setUmapLoading(false);
            });
    }

    // Helper function for point-in-polygon check
    const isPointInPolygon = (point, polygon) => {
        const [x, y] = point;
        let inside = false;

        for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
            const [xi, yi] = polygon[i];
            const [xj, yj] = polygon[j];

            if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
                inside = !inside;
            }
        }

        return inside;
    };

    // Generate tissue image layers
    const generateImageLayers = useCallback(() => {
        const layers = selectedSamples.map(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            const hasImage = !!hiresImages[sample.id];

            if (!imageSize) return null;

            if (!hasImage) {
                return null;
            }

            const layer = new BitmapLayer({
                id: `tissue-image-${sample.id}`,
                image: hiresImages[sample.id],
                bounds: [
                    offset[0],
                    offset[1] + imageSize[1],
                    offset[0] + imageSize[0],
                    offset[1]
                ],
                opacity: 0.8,
                parameters: { depthTest: false },
                // Keep layer stable; no custom updateTriggers/transitions to avoid flicker on zoom
            });

            return layer;
        }).filter(Boolean);

        return layers;
    }, [selectedSamples, imageSizes, sampleOffsets, hiresImages]);

    // Generate cell scatter layers and kosara polygons (gene mode)
    const kosaraPolygonsBySample = useMemo(() => {
        const result = {};
        selectedSamples.forEach(sample => {
            const sampleId = sample.id;
            const mode = radioCellGeneModes[sampleId];
            // Only generate kosara polygons if Kosara display is enabled
            if (kosaraDisplayEnabled && mode === 'genes' && selectedGenes.length > 0 && (kosaraDataBySample[sampleId]?.length > 0)) {
                const offset = sampleOffsets[sampleId] || [0, 0];
                const optimizedPathData = kosaraDataBySample[sampleId].flatMap(d => {
                    const angles = Object.entries(d.angles || {});
                    const ratios = Object.entries(d.ratios || {});
                    const radius = Object.entries(d.radius || {});
                    return generateKosaraPath(
                        (d.cell_x || 0) + offset[0],
                        (d.cell_y || 0) + offset[1],
                        angles,
                        ratios,
                        radius
                    ).map(path => ({
                        id: d.id,
                        cell_type: d.cell_type,
                        points: path.path,
                        color: path.color,
                        gene: path.gene,
                        total_expression: d.total_expression,
                        ratios: d.ratios,
                    }));
                });
                result[sampleId] = optimizedPathData;
            }
        });
        return result;
    }, [selectedSamples, radioCellGeneModes, selectedGenes, geneColorMap, kosaraDataBySample, sampleOffsets, generateKosaraPath, kosaraDisplayEnabled]);

    // Precompute hovered ID sets per sample for efficient matching
    const hoveredIdsSetBySample = useMemo(() => {
        if (!hoveredCluster || !hoveredCluster.sampleId || !Array.isArray(hoveredCluster.cellIds)) return {};
        return { [hoveredCluster.sampleId]: new Set(hoveredCluster.cellIds.map(String)) };
    }, [hoveredCluster]);

    const generateCellLayers = useCallback(() => {
        return selectedSamples.flatMap(sample => {
            const sampleId = sample.id;
            const mode = radioCellGeneModes[sampleId];

            // If in gene mode and single gene data available, draw single gene expression visualization
            if (kosaraDisplayEnabled && mode === 'genes' && singleGeneDataBySample[sampleId]?.cells?.length > 0) {
                const singleGeneData = singleGeneDataBySample[sampleId];
                const offset = sampleOffsets[sampleId] || [0, 0];
                const baseRadius = 5;
                const zoomFactor = mainViewState ? Math.pow(2, mainViewState.zoom * 0.8) : 1;
                const dynamicRadius = baseRadius * zoomFactor;

                const expressionData = singleGeneData.cells.map(cell => ({
                    ...cell,
                    x: (cell.cell_x || 0) + offset[0],
                    y: (cell.cell_y || 0) + offset[1]
                }));

                const layers = [new ScatterplotLayer({
                    id: `single-gene-expression-${sampleId}`,
                    data: expressionData,
                    getPosition: d => [d.x, d.y],
                    getRadius: dynamicRadius,
                    getFillColor: d => {
                        // Get the selected gene and its color from the color map
                        const selectedGene = selectedGenes[0]; // For single gene mode
                        const baseColor = geneColorMap[selectedGene] || "#d73027";

                        const color = getSequentialColor(
                            d.expression,
                            singleGeneData.min_expression,
                            singleGeneData.max_expression,
                            baseColor
                        );
                        // Binary opacity: low for zero expression, high for any expression
                        const opacity = d.expression === 0 ? 0 : 255;
                        return [...color, opacity];
                    },
                    pickable: true,
                    stroked: false,
                    radiusUnits: 'pixels',
                    parameters: { depthTest: false, blend: true },
                    updateTriggers: {
                        data: [singleGeneData, sampleId],
                        getFillColor: [singleGeneData.min_expression, singleGeneData.max_expression, geneColorMap],
                        getRadius: [sampleId, mainViewState?.zoom]
                    },
                    transitions: {
                        getPosition: 0,
                        getFillColor: 0,
                        getRadius: 0
                    }
                })];

                // Add a highlight overlay for hovered cells
                const hoveredSet = hoveredIdsSetBySample[sampleId] || null;
                if (hoveredSet) {
                    const highlightData = expressionData.filter(d => hoveredSet.has(String(d.id)));

                    layers.push(new ScatterplotLayer({
                        id: `single-gene-highlight-${sampleId}`,
                        data: highlightData,
                        getPosition: d => [d.x, d.y],
                        getRadius: dynamicRadius * 1.6,
                        getFillColor: [255, 215, 0, 220],
                        getLineColor: [255, 140, 0, 255],
                        getLineWidth: 2,
                        lineWidthUnits: 'pixels',
                        radiusUnits: 'pixels',
                        pickable: false,
                        stroked: true,
                        updateTriggers: {
                            getRadius: [sampleId, mainViewState?.zoom, hoveredCluster],
                        },
                        parameters: { depthTest: false }
                    }));
                }

                return layers;
            }
            // If in gene mode and kosara data available and Kosara display is enabled, draw kosara polygons + optional highlight overlay
            else if (kosaraDisplayEnabled && mode === 'genes' && kosaraPolygonsBySample[sampleId]?.length > 0) {
                const optimizedPathData = kosaraPolygonsBySample[sampleId];
                const layers = [new PolygonLayer({
                    id: `kosara-polygons-${sampleId}`,
                    data: optimizedPathData,
                    getPolygon: d => d.points,
                    getFillColor: d => {
                        // if this polygon corresponds to a gene slice, color via map/palette; otherwise use fixed color
                        if (d.gene) {
                            const hex = geneColorMap[d.gene] || (() => {
                                const pos = selectedGenes.indexOf(d.gene);
                                const fallback = COLOR_PALETTE[(pos >= 0 ? pos : 0) % COLOR_PALETTE.length];
                                return fallback;
                            })();
                            const rgbColor = convertHEXToRGB(hex);
                            return [...rgbColor, 255];
                        }
                        const rgbColor = convertHEXToRGB(d.color || '#333333');
                        return [...rgbColor, 255];
                    },
                    pickable: true,
                    stroked: false,
                    parameters: { depthTest: false, blend: true },
                    updateTriggers: { data: [kosaraPolygonsBySample[sampleId], sampleId], getFillColor: [geneColorMap, selectedGenes] },
                    transitions: {
                        getPolygon: 0,
                        getFillColor: 0
                    }
                })];

                // Add a highlight overlay for hovered cells to ensure cross-highlighting works in gene mode
                const hoveredSet = hoveredIdsSetBySample[sampleId] || null;
                if (hoveredSet) {
                    const highlightData = (filteredCellData[sampleId] || []).filter(d => hoveredSet.has(String(d.id ?? d.cell_id)));
                    const baseRadius = 5;
                    const zoomFactor = mainViewState ? Math.pow(2, mainViewState.zoom * 0.8) : 1;
                    const dynamicRadius = baseRadius * zoomFactor;

                    layers.push(new ScatterplotLayer({
                        id: `cells-highlight-${sampleId}`,
                        data: highlightData,
                        getPosition: d => [d.x, d.y],
                        getRadius: dynamicRadius * 1.6,
                        getFillColor: [255, 215, 0, 220],
                        getLineColor: [255, 140, 0, 255],
                        getLineWidth: 2,
                        lineWidthUnits: 'pixels',
                        radiusUnits: 'pixels',
                        pickable: false,
                        stroked: true,
                        updateTriggers: {
                            getRadius: [sampleId, mainViewState?.zoom, hoveredCluster],
                        },
                        parameters: { depthTest: false }
                    }));
                }

                return layers;
            }

            // Otherwise, default cell scatter for cell type highlighting
            const cellData = filteredCellData[sampleId] || [];
            const baseRadius = 5;
            // Make radius responsive to zoom: smaller when zoomed out, larger when zoomed in
            const zoomFactor = mainViewState ? Math.pow(2, mainViewState.zoom * 0.8) : 1;
            const dynamicRadius = baseRadius * zoomFactor;

            const hoveredSet = hoveredIdsSetBySample[sampleId] || null;

            const hasConfirmedGeneData = (singleGeneDataBySample[sampleId]?.cells?.length > 0) || (kosaraDataBySample[sampleId]?.length > 0);
            const isGeneModeWithoutSelection = mode === 'genes' && (!selectedGenes || selectedGenes.length === 0);
            const isGeneModeWithoutConfirmedData = mode === 'genes' && selectedGenes.length > 0 && !hasConfirmedGeneData;
            const shouldShowCellTypes = mode === 'cellTypes' || isGeneModeWithoutSelection || isGeneModeWithoutConfirmedData;

            return [new ScatterplotLayer({
                id: `cells-${sampleId}`,
                data: cellData,
                getPosition: d => [d.x, d.y],
                getRadius: d => {
                    const localId = d.id ?? d.cell_id;
                    const isHoveredSample = !!hoveredSet;
                    if (isHoveredSample && hoveredSet.has(String(localId))) {
                        return dynamicRadius * 1.5;
                    }
                    return dynamicRadius;
                },
                getFillColor: d => {
                    const cellType = d.cell_type;

                    if (shouldShowCellTypes) {
                        const sampleSelectedCellTypes = selectedCellTypes && selectedCellTypes[sampleId] ? selectedCellTypes[sampleId] : [];
                        if (cellType && sampleSelectedCellTypes.includes(cellType)) {
                            const color = cellTypeColors[cellType];
                            if (color) {
                                const rgb = color.match(/\w\w/g)?.map(x => parseInt(x, 16)) || [100, 100, 100];
                                return [...rgb, 200];
                            }
                        }
                    }
                    
                    const localId = d.id ?? d.cell_id;
                    if (hoveredSet && hoveredSet.has(String(localId))) {
                        return [255, 215, 0, 200];
                    }
                    if (hoveredSet) {
                        return [150, 150, 150, 50];
                    }
                    return [0, 0, 0, 0];
                },
                getLineColor: d => {
                    const localId = d.id ?? d.cell_id;
                    if (hoveredSet && hoveredSet.has(String(localId))) {
                        return [255, 140, 0, 255];
                    }
                    if (hoveredSet) {
                        return [150, 150, 150, 100];
                    }
                    return [0, 0, 0, 0];
                },
                getLineWidth: d => {
                    const localId = d.id ?? d.cell_id;
                    if (hoveredSet && hoveredSet.has(String(localId))) {
                        return 2;
                    }
                    if (hoveredSet) {
                        return 1;
                    }
                    return 0;
                },
                lineWidthUnits: 'pixels',
                pickable: true,
                radiusUnits: 'pixels',
                stroked: true,
                filled: (!!hoveredSet) || (shouldShowCellTypes && selectedCellTypes && selectedCellTypes[sampleId] && selectedCellTypes[sampleId].length > 0),
                updateTriggers: {
                    getFillColor: [hoveredCluster, selectedCellTypes && selectedCellTypes[sampleId] ? selectedCellTypes[sampleId] : [], cellTypeColors, sampleId, shouldShowCellTypes],
                    getLineColor: [hoveredCluster, sampleId],
                    getRadius: [sampleId, mainViewState?.zoom, hoveredCluster],
                    getLineWidth: [hoveredCluster, sampleId],
                },
                transitions: {
                    getPosition: 0,
                    getRadius: 0
                }
            })];
        }).filter(Boolean);
    }, [selectedSamples, filteredCellData, hoveredCluster, hoveredIdsSetBySample, selectedCellTypes, cellTypeColors, radioCellGeneModes, kosaraPolygonsBySample, singleGeneDataBySample, mainViewState, kosaraDisplayEnabled]);

    // Generate custom area layers
    const generateCustomAreaLayers = useCallback(() => {
        const layers = [];

        // Completed custom areas
        customAreas.forEach(area => {
            const areaColor = convertHEXToRGB(area.color || '#ff0000');

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
                const pendingAreaColor = convertHEXToRGB(areaColor || pendingArea.color);
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
    const layers = useMemo(() => {
        // Build layers in stable order to prevent reordering artifacts
        const imgLayers = generateImageLayers();
        const cellLayers = generateCellLayers();
        const areaLayers = generateCustomAreaLayers();
        return [...imgLayers, ...cellLayers, ...areaLayers];
    }, [generateImageLayers, generateCellLayers, generateCustomAreaLayers]);

    // Stable viewState wrapper to avoid creating a new object on every render
    const deckViewState = useMemo(() => ({ main: mainViewState }), [mainViewState]);

    useEffect(() => { kosaraLoadingSamplesRef.current = kosaraLoadingSamples; }, [kosaraLoadingSamples]);

    // Update radioCellGeneModes when selectedSamples changes
    useEffect(() => {
        setRadioCellGeneModes(prev => {
            const newModes = { ...prev };
            selectedSamples.forEach(sample => {
                if (!(sample.id in newModes)) {
                    newModes[sample.id] = 'cellTypes';
                }
            });
            return newModes;
        });

        // Initialize previous modes for new samples
        setPreviousModes(prev => {
            const newPreviousModes = { ...prev };
            selectedSamples.forEach(sample => {
                if (!(sample.id in newPreviousModes)) {
                    newPreviousModes[sample.id] = 'cellTypes';
                }
            });
            return newPreviousModes;
        });

        // Reset trajectory tracking when samples change to allow loading for new samples
        lastLoadedTrajectoryRef.current = null;
    }, [selectedSamples]);

    // Clear gene expression data when no genes are selected
    useEffect(() => {
        if (selectedGenes.length === 0) {
            // Clear both single gene and kosara data when no genes are selected
            setSingleGeneDataBySample({});
            setKosaraDataBySample({});
        }
    }, [selectedGenes]);

    // Handle Kosara display toggle changes from parent
    useEffect(() => {
        if (!kosaraDisplayEnabled) {
            // Reset trajectory tracking when kosara is disabled
            lastLoadedTrajectoryRef.current = null;

            // When Kosara display is turned OFF, save current state and restore previous state

            // Save current gene selections
            setPreviousGeneSelections(prev => ({ ...prev, current: selectedGenes }));
            setPreviousKosaraData(prev => ({ ...prev, current: kosaraDataBySample }));

            // Save current modes as previous modes for next time
            setPreviousModes(radioCellGeneModes);

            // Restore previous state or default to cellTypes
            const restoredModes = {};
            selectedSamples.forEach(sample => {
                restoredModes[sample.id] = previousModes[sample.id] || 'cellTypes';
            });
            setRadioCellGeneModes(restoredModes);

            // Clear current gene displays when kosara is disabled
            setSelectedGenes([]);
            setGeneColorMap({});
            setKosaraDataBySample({});
            setSingleGeneDataBySample({});
        } else {
            // When Kosara display is turned ON, restore saved current state

            // Save current state as previous
            setPreviousModes(radioCellGeneModes);

            // Restore the "current" state that was saved when we turned off Kosara
            if (previousGeneSelections.current) {
                setSelectedGenes(previousGeneSelections.current);
            }

            if (previousKosaraData.current) {
                setKosaraDataBySample(previousKosaraData.current);

                // Set modes back to genes if there was kosara data
                const restoredModes = {};
                selectedSamples.forEach(sample => {
                    if (previousKosaraData.current[sample.id]?.length > 0) {
                        restoredModes[sample.id] = 'genes';
                    } else {
                        restoredModes[sample.id] = radioCellGeneModes[sample.id];
                    }
                });
                setRadioCellGeneModes(restoredModes);
            }
        }
    }, [kosaraDisplayEnabled]);

    // Handle trajectory gene selection changes
    useEffect(() => {
        if (kosaraDisplayEnabled && trajectoryGenes.length > 0 && trajectoryGenesSample) {
            // Create a key to track the current trajectory selection
            const currentTrajectoryKey = `${trajectoryGenesSample}:${trajectoryGenes.sort().join(',')}`;

            // Only proceed if this is a new trajectory selection
            if (lastLoadedTrajectoryRef.current === currentTrajectoryKey) {
                return; // Skip if we've already loaded this exact combination
            }

            // Check if the trajectory genes sample matches any of our selected samples
            const matchingSample = selectedSamples.find(sample => sample.id === trajectoryGenesSample);
            if (matchingSample) {
                // Update the tracking reference
                lastLoadedTrajectoryRef.current = currentTrajectoryKey;

                // Update available genes and selected genes for the trajectory sample
                setAvailableGenes(trajectoryGenes);
                setSelectedGenes(trajectoryGenes);



                // Load data for the trajectory sample - use single gene mode if only one gene
                if (trajectoryGenes.length === 1) {
                    loadSingleGeneDataForSample(trajectoryGenesSample, trajectoryGenes[0]);
                } else {
                    loadKosaraDataForSample(trajectoryGenesSample, trajectoryGenes);
                }
            }
        }
    }, [trajectoryGenes, trajectoryGenesSample, kosaraDisplayEnabled, selectedSamples]);

    // Preload high-res images for all selected samples
    useEffect(() => {
        let isMounted = true;

        // Reset callback flag when samples change
        imagesLoadedCallbackCalled.current = false;

        // Clean up images that are no longer needed
        setHiresImages(prev => {
            const currentSampleIds = new Set(selectedSamples.map(s => s.id));
            const filteredImages = {};
            Object.keys(prev).forEach(sampleId => {
                if (currentSampleIds.has(sampleId)) {
                    filteredImages[sampleId] = prev[sampleId];
                } else {
                    // Revoke the object URL to free memory
                    try { URL.revokeObjectURL(prev[sampleId]); } catch (e) { }
                }
            });
            return filteredImages;
        });

        selectedSamples.forEach(sample => {
            // Check if image is already loaded or currently being fetched
            setHiresImages(prev => {
                // If image is already loaded, don't fetch again
                if (prev[sample.id]) {
                    return prev;
                }

                // If already fetching this image, don't start another request
                if (fetchingImages.current.has(sample.id)) {
                    return prev;
                }

                // Mark as being fetched
                fetchingImages.current.add(sample.id);

                // Start fetching asynchronously
                fetch('/api/get_hires_image', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ sample_id: sample.id })
                })
                    .then(response => {
                        return response.ok ? response.blob() : null;
                    })
                    .then(blob => {
                        if (blob && isMounted) {
                            const imageUrl = URL.createObjectURL(blob);
                            setHiresImages(currentState => {
                                const newState = {
                                    ...currentState,
                                    [sample.id]: imageUrl
                                };
                                return newState;
                            });
                        }
                    })
                    .catch(error => {
                        console.error(`Error fetching image for ${sample.id}:`, error);
                    })
                    .finally(() => {
                        // Remove from fetching set when done
                        fetchingImages.current.delete(sample.id);
                    });

                // Return current state immediately (fetch is async)
                return prev;
            });
        });

        return () => {
            isMounted = false;
        };
    }, [selectedSamples]);

    // Check if all images are loaded and call callback
    useEffect(() => {
        if (!onImagesLoaded || imagesLoadedCallbackCalled.current) return;

        const timeoutId = setTimeout(() => {
            const allImagesLoaded = selectedSamples.length > 0 && selectedSamples.every(sample => hiresImages[sample.id]);

            if (allImagesLoaded) {
                imagesLoadedCallbackCalled.current = true;
                onImagesLoaded();
            }
        }, 2000);

        return () => clearTimeout(timeoutId);
    }, [hiresImages, selectedSamples, onImagesLoaded]);

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

    // Add keyboard event listener for both drawing and magnifier
    useEffect(() => {
        const handleKeyDown = (event) => {
            // Skip if user is typing in an input field
            if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
                return;
            }

            // Handle drawing keys
            handleKeyPress(event);

            // Handle magnifier keys
            if ((event.code === 'Space') && !keyPressed && !isDrawing) {
                event.preventDefault();
                event.stopPropagation();
                setKeyPressed(true);
                setMagnifierVisible(true);
            }
        };

        const handleKeyUp = (event) => {
            // Skip if user is typing in an input field
            if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
                return;
            }

            if ((event.code === 'Space') && keyPressed) {
                event.preventDefault();
                event.stopPropagation();
                setKeyPressed(false);
                setMagnifierVisible(false);
            }
        };

        document.addEventListener('keydown', handleKeyDown, true); // Use capture phase
        document.addEventListener('keyup', handleKeyUp, true); // Use capture phase

        return () => {
            document.removeEventListener('keydown', handleKeyDown, true);
            document.removeEventListener('keyup', handleKeyUp, true);
        };
    }, [handleKeyPress, keyPressed, isDrawing]);

    // Add native mouse event listener to ensure mouse position is always tracked
    useEffect(() => {
        const handleNativeMouseMove = (event) => {
            if (!containerRef.current || !mainViewState) return;

            const rect = containerRef.current.getBoundingClientRect();
            const x = event.clientX - rect.left;
            const y = event.clientY - rect.top;

            // Convert screen coordinates to world coordinates using DeckGL's viewport
            try {
                // Create a viewport from the current view state
                const viewport = new OrthographicView({ id: 'main' }).makeViewport({
                    width: rect.width,
                    height: rect.height,
                    viewState: mainViewState
                });

                const worldCoords = viewport.unproject([x, y]);

                if (worldCoords && worldCoords.length >= 2) {
                    setMagnifierMousePos(prev => {
                        if (prev && prev.x === worldCoords[0] && prev.y === worldCoords[1]) return prev;
                        return { x: worldCoords[0], y: worldCoords[1] };
                    });
                }
            } catch (error) {
                // Silently ignore projection errors
            }
        };

        // Only add listener when magnifier could be used
        if (containerRef.current) {
            containerRef.current.addEventListener('mousemove', handleNativeMouseMove);
            return () => {
                if (containerRef.current) {
                    containerRef.current.removeEventListener('mousemove', handleNativeMouseMove);
                }
            };
        }
    }, [mainViewState]);

    // Initialize magnifier position when it becomes visible
    useEffect(() => {
        if (magnifierVisible && !isDrawing && mainViewState && selectedSamples.length > 0) {
            // Try to use current mouse position first
            if (magnifierMousePos && magnifierMousePos.x !== undefined && magnifierMousePos.y !== undefined) {
                const hoveredSample = getSampleAtCoordinate(magnifierMousePos.x, magnifierMousePos.y);
                if (hoveredSample) {
                    // Use current mouse position
                    updateMagnifierViewport(magnifierMousePos.x, magnifierMousePos.y, hoveredSample);
                    return;
                }
            }

            // Fallback: Initialize magnifier at the center of the first sample only if no mouse position
            const firstSample = selectedSamples[0];
            const offset = sampleOffsets[firstSample.id] ?? [0, 0];
            const size = imageSizes[firstSample.id] ?? [0, 0];

            if (size[0] > 0 && size[1] > 0) {
                const centerX = offset[0] + size[0] / 2;
                const centerY = offset[1] + size[1] / 2;

                // Update magnifier position
                updateMagnifierViewport(centerX, centerY, firstSample.id);
            }
        }
    }, [magnifierVisible, isDrawing, mainViewState, selectedSamples, sampleOffsets, imageSizes, updateMagnifierViewport, magnifierMousePos, getSampleAtCoordinate]);

    // Cleanup effect for magnifier and images
    useEffect(() => {
        return () => {
            Object.values(hiresImages).forEach(url => {
                try { URL.revokeObjectURL(url); } catch (e) { }
            });

            if (viewStateRafRef.current) {
                cancelAnimationFrame(viewStateRafRef.current);
                viewStateRafRef.current = null;
            }
            viewStatePendingRef.current = null;
        };
    }, [hiresImages]);

    // In handleMouseMove, update magnifier logic to use hiresImages
    useEffect(() => {
        if (!magnifierVisible || isDrawing || isAreaTooltipVisible || isAreaEditPopupVisible) {
            setMagnifierData(prev => (prev === null ? prev : null));
            return;
        }

        if (!magnifierMousePos || !selectedSamples.length) {
            setMagnifierData(prev => (prev === null ? prev : null));
            return;
        }

        const { x: worldX, y: worldY } = magnifierMousePos;
        const hoveredSample = getSampleAtCoordinate(worldX, worldY);

        if (hoveredSample && imageSizes[hoveredSample]) {
            const imageUrl = hiresImages[hoveredSample];
            const size = imageSizes[hoveredSample];
            if (!imageUrl) {
                setMagnifierData(prev => (prev === null ? prev : null));
            } else {
                setMagnifierData(prev => {
                    if (
                        prev &&
                        prev.sampleId === hoveredSample &&
                        prev.imageUrl === imageUrl &&
                        prev.imageSize && size &&
                        prev.imageSize[0] === size[0] &&
                        prev.imageSize[1] === size[1]
                    ) {
                        return prev;
                    }
                    return {
                        imageUrl,
                        sampleId: hoveredSample,
                        imageSize: size
                    };
                });
            }
        } else {
            setMagnifierData(prev => (prev === null ? prev : null));
        }
    }, [magnifierVisible, magnifierMousePos, selectedSamples, hiresImages, imageSizes, isDrawing, isAreaTooltipVisible, isAreaEditPopupVisible, getSampleAtCoordinate]);

    useEffect(() => {
        const currentUrl = magnifierData?.imageUrl || null;
        if (!magnifierVisible) {
            // When hidden, clear loaded flag (spinner hidden by visibility anyway)
            if (magnifierImageLoaded) setMagnifierImageLoaded(false);
            prevMagnifierUrlRef.current = currentUrl; // store for next open
            return;
        }
        // Visible: only reset if URL actually changed
        if (prevMagnifierUrlRef.current !== currentUrl) {
            prevMagnifierUrlRef.current = currentUrl;
            if (currentUrl) {
                setMagnifierImageLoaded(false);
                setMagnifierImageVersion(v => v + 1); // trigger remount so onLoad always fires
            }
        }
    }, [magnifierVisible, magnifierData?.imageUrl, magnifierImageLoaded]);

    // Fallback initialize magnifierData if visible but not yet set (ensures spinner/image sequence shows)
    useEffect(() => {
        if (!magnifierVisible) return;
        if (magnifierData) return;
        if (!selectedSamples.length) return;
        const firstSample = selectedSamples[0];
        const imgUrl = hiresImages[firstSample.id];
        const size = imageSizes[firstSample.id];
        if (imgUrl && size) {
            setMagnifierData({ imageUrl: imgUrl, sampleId: firstSample.id, imageSize: size });
        }
    }, [magnifierVisible, magnifierData, selectedSamples, hiresImages, imageSizes]);

    return (
        <div style={{ width: '100%', height: '100%', position: 'relative' }}>
            {/* Main content */}
            <div
                ref={containerRef}
                style={{
                    width: '100%',
                    height: '100%',
                    position: 'relative',
                    overflow: 'hidden'
                }}
            >
                <DeckGL
                    layers={layers}
                    views={[mainView]}
                    viewState={deckViewState}
                    onViewStateChange={handleViewStateChange}
                    onClick={handleMapClick}
                    onHover={(info) => {
                        // Preserve existing mouse-move behavior
                        handleMouseMove(info);

                        // Show tooltip for cells and kosara polygons
                        if (info && info.object && info.layer && info.layer.id) {
                            const layerId = info.layer.id;
                            if (layerId.startsWith('kosara-polygons-')) {
                                const sampleId = layerId.replace('kosara-polygons-', '');
                                const { id, cell_type, ratios, total_expression } = info.object || {};
                                setHoveredCell({
                                    id,
                                    sampleId,
                                    cell_type,
                                    ratios,
                                    total_expression,
                                    x: info.x,
                                    y: info.y
                                });
                            } else if (layerId.startsWith('single-gene-expression-')) {
                                const sampleId = layerId.replace('single-gene-expression-', '');
                                const { id, cell_type, expression } = info.object || {};
                                setHoveredCell({
                                    id,
                                    sampleId,
                                    cell_type,
                                    expression,
                                    x: info.x,
                                    y: info.y
                                });
                            } else if (layerId.startsWith('cells-')) {
                                const sampleId = layerId.split('-')[1];
                                const { id, cell_type } = info.object || {};
                                setHoveredCell({
                                    id,
                                    sampleId,
                                    cell_type,
                                    x: info.x,
                                    y: info.y
                                });
                            } else {
                                setHoveredCell(null);
                            }
                        } else {
                            setHoveredCell(null);
                        }
                    }}
                    controller={deckController}
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

                {hoveredCell && (
                    <div style={{
                        position: 'absolute',
                        left: hoveredCell.x + 12,
                        top: hoveredCell.y - 40,
                        pointerEvents: 'none',
                        backgroundColor: 'rgba(255, 255, 255, 0.9)',
                        padding: 8,
                        borderRadius: 4,
                        boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
                        transform: 'none',
                        whiteSpace: 'nowrap',
                        willChange: 'left, top',
                        fontSize: 12,
                        zIndex: 1000,
                        textAlign: 'left'
                    }}>
                        {hoveredCell.id ? (
                            <>
                                <div><strong>Sample:</strong> {hoveredCell.sampleId}</div>
                                <div><strong>Cell Type:</strong> {hoveredCell.cell_type}</div>
                                {hoveredCell.total_expression !== undefined && (
                                    <div><strong>Total Expression:</strong> {Number(hoveredCell.total_expression).toFixed(5)}</div>
                                )}
                                {hoveredCell.ratios && Object.entries(hoveredCell.ratios).map(([gene, expression]) => (
                                    <div key={gene}><strong>{gene}:</strong> {Number(expression).toFixed(5) * 100}%</div>
                                ))}
                            </>
                        ) : (
                            <>
                                <div>Sample: {hoveredCell.sampleId}</div>
                                <div>Cell Type: {hoveredCell.cell_type}</div>
                            </>
                        )}
                    </div>
                )}

                {/* Sample controls */}
                <div style={{
                    position: 'absolute',
                    top: minimapVisible ? 170 : 10,
                    left: 10,
                    zIndex: 20,
                    transition: 'top 0.3s ease-in-out'
                }}>
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
                            type={minimapVisible ? "primary" : "default"}
                            onClick={toggleMinimapVisible}
                            style={{
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
                            type={isDrawing ? "primary" : "default"}
                            onClick={toggleDrawingMode}
                            style={{
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
                        <div style={{ fontWeight: 'bold', marginBottom: 4 }}>Drawing Shortcuts:</div>
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
                {(minimapVisible || minimapAnimating) && selectedSamples.length > 0 && (
                    <div
                        style={{
                            position: 'absolute',
                            top: 10,
                            left: 10,
                            width: 296,
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
                        {/* Minimap background - composite view for multiple samples */}
                        <div style={{
                            width: '100%',
                            height: '100%',
                            position: 'relative',
                            backgroundColor: '#f0f0f0'
                        }}>
                            {(() => {
                                // Calculate total world bounds for positioning
                                let totalLeft = Infinity;
                                let totalRight = -Infinity;
                                let totalTop = Infinity;
                                let totalBottom = -Infinity;

                                selectedSamples.forEach(sample => {
                                    const imageSize = imageSizes[sample.id];
                                    const offset = sampleOffsets[sample.id] || [0, 0];
                                    if (imageSize) {
                                        totalLeft = Math.min(totalLeft, offset[0]);
                                        totalRight = Math.max(totalRight, offset[0] + imageSize[0]);
                                        totalTop = Math.min(totalTop, offset[1]);
                                        totalBottom = Math.max(totalBottom, offset[1] + imageSize[1]);
                                    }
                                });

                                if (totalLeft === Infinity) return null;

                                const totalWidth = totalRight - totalLeft;
                                const totalHeight = totalBottom - totalTop;

                                // Render each sample image positioned within the composite minimap
                                return selectedSamples.map(sample => {
                                    const imageSize = imageSizes[sample.id];
                                    const offset = sampleOffsets[sample.id] || [0, 0];

                                    if (!imageSize) return null;

                                    // Calculate relative position and size within the total bounds
                                    const relativeLeft = ((offset[0] - totalLeft) / totalWidth) * 100;
                                    const relativeTop = ((offset[1] - totalTop) / totalHeight) * 100;
                                    const relativeWidth = (imageSize[0] / totalWidth) * 100;
                                    const relativeHeight = (imageSize[1] / totalHeight) * 100;

                                    return (
                                        <div
                                            key={sample.id}
                                            style={{
                                                position: 'absolute',
                                                left: `${relativeLeft}%`,
                                                top: `${relativeTop}%`,
                                                width: `${relativeWidth}%`,
                                                height: `${relativeHeight}%`,
                                                cursor: 'pointer',
                                                border: '1px solid rgba(24, 144, 255, 0.3)',
                                                borderRadius: '2px',
                                                overflow: 'hidden',
                                                transition: 'all 0.2s ease',
                                                boxSizing: 'border-box'
                                            }}
                                            onMouseEnter={(e) => {
                                                e.target.style.border = '2px solid #1890ff';
                                                e.target.style.boxShadow = '0 2px 8px rgba(24, 144, 255, 0.4)';
                                                e.target.style.transform = 'scale(1.02)';
                                                e.target.style.zIndex = '3';
                                            }}
                                            onMouseLeave={(e) => {
                                                e.target.style.border = '1px solid rgba(24, 144, 255, 0.3)';
                                                e.target.style.boxShadow = 'none';
                                                e.target.style.transform = 'scale(1)';
                                                e.target.style.zIndex = '1';
                                            }}
                                            title={`Click to center on ${sample.name}`}
                                        >
                                            <img
                                                src={`/${sample.id}_full.jpg`}
                                                alt={`Minimap ${sample.name}`}
                                                style={{
                                                    width: '100%',
                                                    height: '100%',
                                                    objectFit: 'cover',
                                                    display: 'block',
                                                    pointerEvents: 'none'
                                                }}
                                                draggable={false}
                                            />
                                            {/* Sample label overlay */}
                                            <div
                                                style={{
                                                    position: 'absolute',
                                                    bottom: '2px',
                                                    left: '2px',
                                                    right: '2px',
                                                    fontSize: '8px',
                                                    fontWeight: 'bold',
                                                    color: 'white',
                                                    backgroundColor: 'rgba(0, 0, 0, 0.7)',
                                                    padding: '1px 3px',
                                                    borderRadius: '2px',
                                                    textAlign: 'center',
                                                    textShadow: '1px 1px 1px rgba(0, 0, 0, 0.8)',
                                                    whiteSpace: 'nowrap',
                                                    overflow: 'hidden',
                                                    textOverflow: 'ellipsis',
                                                    pointerEvents: 'none'
                                                }}
                                            >
                                                {sample.name || sample.id}
                                            </div>
                                        </div>
                                    );
                                });
                            })()}
                        </div>

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
                                        boxSizing: 'border-box',
                                        zIndex: 10
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

                {/* Magnifying Glass */}
                {magnifierVisible && (
                    <div
                        ref={magnifierRef}
                        style={{
                            position: 'absolute',
                            bottom: 10,
                            left: 10,
                            width: 300,
                            height: 300,
                            zIndex: 15,
                            backgroundColor: 'rgba(255, 255, 255, 0.95)',
                            border: '2px solid #1890ff',
                            borderRadius: 8,
                            boxShadow: '0 4px 12px rgba(0, 0, 0, 0.2)',
                            overflow: 'hidden',
                            opacity: magnifierVisible ? 1 : 0,
                            transition: 'opacity 0.2s ease-in-out',
                            pointerEvents: 'none' // allow underlying map to keep receiving hover events
                        }}
                    >
                        {/* Header */}
                        <div style={{
                            padding: '6px 12px',
                            borderBottom: '1px solid #e8e8e8',
                            fontSize: '11px',
                            fontWeight: 'bold',
                            color: '#262626',
                            backgroundColor: '#f0f8ff',
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                            pointerEvents: 'auto' // header still interactive if needed
                        }}>
                            <span>Magnifier - {magnifierData?.sampleId || ''}</span>
                            <span style={{ fontSize: '9px', color: '#666' }}>
                                Hold Space
                            </span>
                        </div>

                        {/* Magnifier View */}
                        <div style={{
                            position: 'relative',
                            width: '100%',
                            height: 280,
                            overflow: 'hidden',
                            backgroundColor: '#ffffff',
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center'
                        }}>
                            {/* Image (hidden until loaded) */}
                            {magnifierData && (
                                <img
                                    key={`magnifier-img-${magnifierImageVersion}`}
                                    src={magnifierData.imageUrl}
                                    alt="HD Magnifier"
                                    style={{
                                        position: 'absolute',
                                        width: magnifierData.imageSize[0] * 2,
                                        height: magnifierData.imageSize[1] * 2,
                                        left: -(magnifierViewport.x * magnifierData.imageSize[0] * 2) + 150,
                                        top: -(magnifierViewport.y * magnifierData.imageSize[1] * 2) + 140,
                                        imageRendering: 'crisp-edges',
                                        transition: 'left 0.1s ease-out, top 0.1s ease-out, opacity 0.2s',
                                        opacity: magnifierImageLoaded ? 1 : 0
                                    }}
                                    draggable={false}
                                    onLoad={() => setMagnifierImageLoaded(true)}
                                    onError={() => setMagnifierImageLoaded(true)}
                                />
                            )}

                            {/* Spinner overlay while loading */}
                            {magnifierVisible && !magnifierImageLoaded && (
                                <div style={{
                                    position: 'absolute',
                                    inset: 0,
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    background: 'rgba(255,255,255,0.8)',
                                    zIndex: 9999,
                                    pointerEvents: 'none'
                                }}>
                                    <Spin />
                                </div>
                            )}

                            {/* Crosshairs always visible */}
                            <div style={{
                                position: 'absolute',
                                left: 150,
                                top: 0,
                                width: 1,
                                height: '100%',
                                backgroundColor: '#ff4d4f',
                                pointerEvents: 'none',
                                boxShadow: '0 0 2px rgba(0,0,0,0.5)',
                                zIndex: 4
                            }} />
                            <div style={{
                                position: 'absolute',
                                left: 0,
                                top: 140,
                                width: '100%',
                                height: 1,
                                backgroundColor: '#ff4d4f',
                                pointerEvents: 'none',
                                boxShadow: '0 0 2px rgba(0,0,0,0.5)',
                                zIndex: 4
                            }} />
                            <div style={{
                                position: 'absolute',
                                left: 147,
                                top: 137,
                                width: 6,
                                height: 6,
                                borderRadius: '50%',
                                backgroundColor: '#ff4d4f',
                                border: '1px solid white',
                                pointerEvents: 'none',
                                boxShadow: '0 0 4px rgba(0,0,0,0.5)',
                                zIndex: 5
                            }} />
                            {/* Coordinates indicator */}
                            <div style={{
                                position: 'absolute',
                                bottom: 15,
                                right: 5,
                                padding: '2px 6px',
                                backgroundColor: 'rgba(0, 0, 0, 0.7)',
                                color: 'white',
                                fontSize: '9px',
                                borderRadius: 3,
                                fontFamily: 'monospace'
                            }}>
                                X: {Math.round(magnifierMousePos?.x || 0)} Y: {Math.round(magnifierMousePos?.y || 0)}
                            </div>
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
                            <div style={{ marginBottom: 5 }}>
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
                                        value={editResolutions.toString()}
                                        onChange={(value) => setEditResolutions(parseFloat(value) || 0.5)}
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
                            <Button
                                size="small"
                                style={{ marginBottom: 5, width: '100%' }}
                                color="pink"
                                variant="outlined"
                                onClick={generateUmap}
                                loading={umapLoading}
                            >
                                {umapLoading ? 'Generating...' : 'Generate UMAP'}
                            </Button>
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

            {/* Loading overlay */}
            {isKosaraLoading && (
                <div style={{
                    position: 'absolute',
                    top: 0,
                    left: 0,
                    right: 0,
                    bottom: 0,
                    backgroundColor: 'rgba(255, 255, 255, 0.5)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    zIndex: 1000,
                    flexDirection: 'column',
                    gap: '16px'
                }}>
                    <Spin size="large" />
                    <div style={{ fontSize: '16px', color: '#666' }}>
                        Loading Kosara visualization...
                    </div>
                </div>
            )}
        </div>
    );
};