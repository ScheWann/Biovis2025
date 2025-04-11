import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox, TreeSelect, message, Switch, Radio } from "antd";
import { CloseOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, TextLayer, GeoJsonLayer, PolygonLayer, PathLayer } from '@deck.gl/layers';
import { booleanPointInPolygon, centroid, sample } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';
import "../styles/MultiSampleViewer.css";

const CELL_COLORS = [
    '#FF6B6B', '#4ECDC4', '#45B7D1',
    '#96CEB4', '#FFEEAD', '#D4A5A5',
    '#88D8B0', '#FF9999', '#99CCFF'
];

const stringToHash = (str) => {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
        const char = str.charCodeAt(i);
        hash = (hash << 5) - hash + char;
        hash |= 0;
    }
    return Math.abs(hash);
};

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

// HSL to RGB
const hslToRgb = (h, s, l) => {
    h /= 360;
    s /= 100;
    l /= 100;
    let r, g, b;
    if (s === 0) {
        r = g = b = l;
    } else {
        const hue2rgb = (p, q, t) => {
            if (t < 0) t += 1;
            if (t > 1) t -= 1;
            if (t < 1 / 6) return p + (q - p) * 6 * t;
            if (t < 1 / 2) return q;
            if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
            return p;
        };
        const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        const p = 2 * l - q;
        r = hue2rgb(p, q, h + 1 / 3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1 / 3);
    }
    return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
};

const hexToRgbs = (hex) => {
    const r = parseInt(hex.slice(1, 3), 16);
    const g = parseInt(hex.slice(3, 5), 16);
    const b = parseInt(hex.slice(5, 7), 16);
    return [r, g, b];
};

// generate random color, call this function after saving region to update the color for next drawing
const generateRandomColor = () => {
    return hslToRgb(Math.floor(Math.random() * 360), 100, 50);
};

const debounce = (fn, delay) => {
    let timer;
    return (...args) => {
        clearTimeout(timer);
        timer = setTimeout(() => {
            fn(...args);
        }, delay);
    };
};

export const MultiSampleViewer = ({
    setLoading,
    samples,
    cellTypeCoordinatesData,
    cellTypeDir,
    regions,
    setRegions,
    setNMFGOData,
    setNMFGODataLoading,
    analyzedRegion,
    setAnalyzedRegion,
    NMFclusterCells,
    setSelectedRegionGeneExpressionData
}) => {
    const [imageSizes, setImageSizes] = useState({});
    const [tileSize] = useState(256);
    const [features, setFeatures] = useState(
        samples.reduce((acc, sample) => ({
            ...acc,
            [sample.id]: { type: 'FeatureCollection', features: [] }
        }), {})
    );
    const [tempRegions, setTempRegions] = useState({});
    const [visibleCellTypes, setVisibleCellTypes] = useState({});
    const [colorMaps, setColorMaps] = useState({});
    const [hoveredCell, setHoveredCell] = useState(null);
    const [sampleOffsets, setSampleOffsets] = useState({});
    const [activeSample, setActiveSample] = useState(samples[0]?.id);
    const [visibleSamples, setVisibleSamples] = useState(
        samples.reduce((acc, sample) => ({ ...acc, [sample.id]: true }), {})
    );
    const [geneList, setGeneList] = useState({});
    const [selectedGenes, setSelectedGenes] = useState([]);
    const [regionName, setRegionName] = useState('');
    const [regionColor, setRegionColor] = useState(generateRandomColor());
    const [isDrawingActive, setIsDrawingActive] = useState(false);
    const [activeDrawingSample, setActiveDrawingSample] = useState(null);
    const [currentZoom, setCurrentZoom] = useState(-3);
    const [radioCellGeneModes, setRadioCellGeneModes] = useState(samples.reduce((acc, sample) => ({ ...acc, [sample.id]: 'cellTypes' }), {}));
    const [partWholeMode, setPartWholeMode] = useState(true); // defaultly showing gene expression value in the whole regions
    const [geneExpressionData, setGeneExpressionData] = useState([]);
    const TILE_LOAD_ZOOM_THRESHOLD = 0;

    // fetch gene list
    useEffect(() => {
        const sampleNames = samples.map(item => item.id);
        const missingSampleNames = sampleNames.filter(id => !geneList[id]);

        if (missingSampleNames.length === 0) return;

        fetch('/get_all_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_names: missingSampleNames })
        })
            .then(res => res.json())
            .then(data => {
                setGeneList(prev => ({
                    ...prev,
                    ...data
                }));
            });
    }, [samples]);

    // Initialize color map and visible cell types for each sample
    useEffect(() => {
        const initialStates = samples.reduce((acc, sample) => {
            const initialColorMap = {};
            const initialVisibility = {};
            const sampleCellTypes = cellTypeDir?.[sample.id] || [];

            sampleCellTypes.forEach((cellType) => {
                const hash = stringToHash(cellType);
                const hue = hash % 360;
                initialColorMap[cellType] = hslToRgb(hue, 100, 50);
                initialVisibility[cellType] = true;
            });

            return {
                colorMaps: { ...acc.colorMaps, [sample.id]: initialColorMap },
                visibleCellTypes: { ...acc.visibleCellTypes, [sample.id]: initialVisibility }
            };
        }, { colorMaps: {}, visibleCellTypes: {} });

        setColorMaps(initialStates.colorMaps);
        setVisibleCellTypes(initialStates.visibleCellTypes);
    }, [samples, cellTypeDir]);

    // get image sizes for all samples
    useEffect(() => {
        const fetchImageSizes = () => {
            const sampleIds = samples.map(sample => sample.id);

            fetch('/get_hires_image_size', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ sample_ids: sampleIds })
            })
                .then(res => res.json())
                .then(data => {
                    setImageSizes(data);
                })
        };

        if (samples.length) {
            fetchImageSizes();
        }
    }, [samples]);

    useEffect(() => {
        if (Object.keys(imageSizes).length === samples.length) {
            let currentX = 0;
            const margin = 100;
            const newOffsets = {};
            samples.forEach(sample => {
                const size = imageSizes[sample.id];
                newOffsets[sample.id] = [currentX, 0];
                currentX += size[0] + margin;
            });
            setSampleOffsets(newOffsets);
        }
        setVisibleSamples(
            samples.reduce((acc, sample) => ({ ...acc, [sample.id]: true }), {})
        );
    }, [imageSizes, samples]);

    const analyzedRegionData = useMemo(() => {
        if (!analyzedRegion) return null;
        const region = regions.find(r => r.name === analyzedRegion);
        if (!region) return null;
        return {
            sampleId: region.sampleId,
            cellIds: region.cellIds,
        };
    }, [analyzedRegion, regions]);

    const fetchGeneExpressionData = (sampleId, cell_ids) => {
        fetch('/get_selected_region_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: sampleId, cell_list: cell_ids })
        })
            .then(res => res.json())
            .then(data => {
                setSelectedRegionGeneExpressionData(data);
            });
    }

    const fetchNMFGOExpressionData = (sampleId, cell_ids) => {
        setNMFGODataLoading(true);
        fetch('/get_NMF_GO_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: sampleId, cell_list: cell_ids })
        })
            .then(res => res.json())
            .then(data => {
                setNMFGODataLoading(false);
                setNMFGOData(data);
            });
    }

    const fetchCell2CellInteractionData = (regionName, sampleId, cell_ids) => {
        const regions = { [regionName]: { sampleId, cellIds: cell_ids } };
        fetch('/get_cell_cell_interaction_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ regions: regions })
        })
            .then(res => res.json())
            .then(data => {
                setSelectedRegionGeneExpressionData(data);
            });
    };


    // save filtered cell data, only recalculate when dependencies change
    const filteredCellData = useMemo(() => {
        return samples.reduce((acc, sample) => {
            const data = cellTypeCoordinatesData[sample.id] || [];
            const filtered = data.filter(cell => visibleCellTypes[sample.id]?.[cell.cell_type] ?? true);
            return { ...acc, [sample.id]: filtered };
        }, {});
    }, [samples, cellTypeCoordinatesData, visibleCellTypes]);

    // debounce zoom level update
    const debouncedSetZoom = useMemo(() => debounce((zoom) => {
        setCurrentZoom(zoom);
    }, 100), []);

    const changeCellGeneMode = (sampleId, e) => {
        setRadioCellGeneModes(prev => ({
            ...prev,
            [sampleId]: e.target.value
        }));
    };

    const onVisibilityGeneChange = (gene) => {
        setSelectedGenes([...selectedGenes, gene]);
    }

    const cleanGeneSelection = () => {
        setSelectedGenes([]);
    }

    const confirmKosaraPlot = (sampleId) => {
        const sample = [sampleId];
        // const sampleList = samples.map(sample => sample.id);
        // const cleanedGenes = selectedGenes.map(gene => gene.split('-')[1] || gene);
        setLoading(true);
        if (partWholeMode) {
            regions.forEach(region => {
                fetch('/get_kosara_data', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ sample_ids: sample, gene_list: selectedGenes, cell_list: region.cellIds })
                })
                    .then(res => res.json())
                    .then(data => {
                        setGeneExpressionData(data[sample]);
                        setLoading(false);
                    });
            });
        } else {
            fetch('/get_kosara_data', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ sample_ids: sample, gene_list: selectedGenes, cell_list: [] })
            })
                .then(res => res.json())
                .then(data => {
                    setLoading(false);
                    setGeneExpressionData(data[sample]);
                });
        }
    }

    const generateCirclePoints = (cx, cy, r, steps = 20) => {
        const points = [];
        const angleStep = (2 * Math.PI) / steps;
        for (let i = 0; i < steps; i++) {
            const theta = i * angleStep;
            points.push([cx + r * Math.cos(theta), cy + r * Math.sin(theta)]);
        }

        points.push(points[0]);
        return points;
    }

    const generateKosaraPath = (pointX, pointY, angles, ratios, cal_radius) => {
        const baseRadius = 5;
        let paths = [];
        let cellTypes = selectedGenes;
        let startpointX, startpointY, endpointX, endpointY = 0;
        let lastStartPointX, lastStartPointY, lastEndPointX, lastEndPointY, lastCircleRadius = 0;
        let originalPointX = pointX - baseRadius * Math.cos(45 * Math.PI / 180);
        let originalPointY = pointY + baseRadius * Math.sin(45 * Math.PI / 180);

        let cellIndices = ratios
            .filter(item => item[1] !== 0 && cellTypes.includes(item[0]))
            .sort((a, b) => cellTypes.indexOf(a[0]) - cellTypes.indexOf(b[0]))
            .slice(0, 9)
            .map(item => item[0]);

        let cellColors = cellIndices.map(index => {
            const positionInSelectedGenes = selectedGenes.indexOf(index);
            return CELL_COLORS[positionInSelectedGenes % CELL_COLORS.length];
        });
        let cellAngles = cellIndices.map(index => angles.find(item => item[0] === index));
        let cellRadius = cellIndices.map(index => cal_radius.find(item => item[0] === index));

        const ratioSum = ratios.reduce((acc, item) => acc + item[1], 0);

        // If no selected cells are shown, draw an empty circle
        if (cellAngles.length === 0) {
            const circlePoints = generateCirclePoints(originalPointX, originalPointY, baseRadius, 50);
            paths.push({ path: circlePoints, color: '#FFFFFF' });
        } else {
            cellAngles = cellAngles.map(angle => [angle[0], angle[1]]);
            cellRadius = cellRadius.map(rad => [rad[0], rad[1]]);

            cellAngles.forEach((angle, index) => {
                let cal_cell_radius = cellRadius[index][1];
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
                }
                else if (index === cellAngles.length - 1 && ratioSum == 1) {
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
                }
                else {
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
                    color: cellColors[index]
                });

                lastCircleRadius = cal_cell_radius;
                lastStartPointX = startpointX;
                lastStartPointY = startpointY;
                lastEndPointX = endpointX;
                lastEndPointY = endpointY;
            });

            if (ratioSum < 1) {
                let points = [];
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

                paths.push({
                    path: points,
                    color: '#333333'
                });
            }
        }
        return paths;
    }

    const generateComplexArcPoints = (
        startX,
        startY,
        endX,
        endY,
        outerRadius,
        innerRadius,
        outerFlags = { large: 0, sweep: 1 },
        innerFlags = { large: 0, sweep: 0 }
    ) => {
        // outerRadius
        const outerArc = generateSingleArcPoints(
            startX, startY,
            endX, endY,
            outerRadius,
            outerFlags.large,
            outerFlags.sweep
        );

        // innerRadius
        const innerArc = generateSingleArcPoints(
            endX, endY,
            startX, startY,
            innerRadius,
            innerFlags.large,
            innerFlags.sweep
        );

        return [...outerArc, ...innerArc];
    };

    const generateSingleArcPoints = (startX, startY, endX, endY, r, largeArcFlag, sweepFlag) => {
        const dx = endX - startX;
        const dy = endY - startY;
        const d = Math.hypot(dx, dy);
        if (d > 2 * r) {
            console.error("Chord length is greater than the diameter");
            return [];
        }

        // middle point and h value
        const midX = (startX + endX) / 2;
        const midY = (startY + endY) / 2;
        const h = Math.sqrt(r * r - (d / 2) * (d / 2));

        // unit vector perpendicular to the chord (two candidates)
        const ux = -dy / d;
        const uy = dx / d;

        // two candidate centers
        const cx1 = midX + h * ux;
        const cy1 = midY + h * uy;
        const cx2 = midX - h * ux;
        const cy2 = midY - h * uy;

        // given the center of the circle, calculate the start point, end point angle and arc length (without sweepFlag correction, the angle difference is in [0, 2π))
        const computeAngles = (cx, cy) => {
            const startAngle = Math.atan2(startY - cy, startX - cx);
            const endAngle = Math.atan2(endY - cy, endX - cx);
            let delta = endAngle - startAngle;
            // nomoralize to [0, 2π)
            delta = ((delta % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
            return { startAngle, endAngle, delta };
        }

        const cand1 = computeAngles(cx1, cy1);
        const cand2 = computeAngles(cx2, cy2);

        // get effective delta, if sweepFlag is 1, then effective delta is delta; if sweepFlag is 0, then effective delta is 2π - delta
        const effectiveDelta = (candidate) => {
            return sweepFlag === 1 ? candidate.delta : (2 * Math.PI - candidate.delta);
        }

        const eff1 = effectiveDelta(cand1);
        const eff2 = effectiveDelta(cand2);

        // when largeArcFlag === 0, if eff1 <= Math.PI && eff2 > Math.PI, choose cand1; if eff2 <= Math.PI && eff1 > Math.PI, choose cand2
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

        // generate points by sampling
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
    }

    const collapseItems = samples.map((sample) => ({
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
                    <CellTypeSettings
                        cellTypes={cellTypeDir[sample.id] || []}
                        cellData={cellTypeCoordinatesData[sample.id]}
                        colorMap={colorMaps[sample.id] || {}}
                        visibleMap={visibleCellTypes[sample.id] || {}}
                        onColorChange={(type, color) => {
                            setColorMaps(prev => ({
                                ...prev,
                                [sample.id]: { ...prev[sample.id], [type]: color }
                            }));
                        }}
                        onVisibilityCellTypeChange={(type, visible) => {
                            setVisibleCellTypes(prev => ({
                                ...prev,
                                [sample.id]: { ...prev[sample.id], [type]: visible }
                            }));
                        }}
                    />
                ) : (
                    <GeneSettings
                        geneList={geneList[sample.id] || {}}
                        sampleId={sample.id}
                        onVisibilityGeneChange={onVisibilityGeneChange}
                        cleanGeneSelection={cleanGeneSelection}
                        confirmKosaraPlot={confirmKosaraPlot}
                    />
                )}
            </>
        )
    }));

    // TileLayer
    const generateTileLayers = useCallback(() => {
        if (currentZoom < TILE_LOAD_ZOOM_THRESHOLD) {
            return [];
        }
        return samples.map(sample => {
            const imageSize = imageSizes[sample.id];
            if (!imageSize || imageSize.length < 2) return null;
            const offset = sampleOffsets[sample.id] || [0, 0];
            return new TileLayer({
                id: `tif-tiles-${sample.id}`,
                visible: visibleSamples[sample.id],
                tileSize,
                extent: [
                    offset[0],
                    offset[1],
                    offset[0] + imageSize[0],
                    offset[1] + imageSize[1]
                ],
                minZoom: 0,
                maxZoom: 0,
                getTileData: async ({ index }) => {
                    try {
                        const { x, y } = index;
                        const tileOriginX = Math.floor(offset[0] / tileSize);
                        const tileOriginY = Math.floor(offset[1] / tileSize);
                        const relativeX = x - tileOriginX;
                        const relativeY = y - tileOriginY;
                        const response = await fetch(`/get_tile?sample_id=${sample.id}&x=${relativeX}&y=${relativeY}`);
                        const blob = await response.blob();
                        const tiff = await fromBlob(blob);
                        const image = await tiff.getImage();
                        const rgba = await image.readRasters({ interleave: true });
                        const canvas = document.createElement('canvas');
                        [canvas.width, canvas.height] = [image.getWidth(), image.getHeight()];
                        const ctx = canvas.getContext('2d');
                        ctx.translate(canvas.width, 0);
                        ctx.scale(-1, 1);
                        ctx.putImageData(new ImageData(new Uint8ClampedArray(rgba), canvas.width, canvas.height), 0, 0);
                        return {
                            image: canvas,
                            bounds: [
                                [relativeX * tileSize + offset[0], (relativeY + 1) * tileSize + offset[1]],
                                [relativeX * tileSize + offset[0], relativeY * tileSize + offset[1]],
                                [(relativeX + 1) * tileSize + offset[0], relativeY * tileSize + offset[1]],
                                [(relativeX + 1) * tileSize + offset[0], (relativeY + 1) * tileSize + offset[1]]
                            ]
                        };
                    } catch (error) {
                        console.error('Loading tile failed:', error);
                        return null;
                    }
                },
                renderSubLayers: props => props.data && new BitmapLayer(props, {
                    image: props.data.image,
                    bounds: props.data.bounds
                })
            });
        }).filter(Boolean);
    }, [samples, imageSizes, tileSize, visibleSamples, sampleOffsets, currentZoom]);

    const generateWholePngLayers = useCallback(() => {
        return samples.map(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            if (!imageSize || imageSize.length < 2) return null;
            return new BitmapLayer({
                id: `png-replacement-${sample.id}`,
                visible: visibleSamples[sample.id] && currentZoom < TILE_LOAD_ZOOM_THRESHOLD,
                image: `/${sample.id}_full.jpg`,
                bounds: [offset[0], offset[1] + imageSize[1], offset[0] + imageSize[0], offset[1]],
                opacity: 1,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [imageSizes, visibleSamples, sampleOffsets, currentZoom]);

    // cell boundaries image layers
    const generateMarkerImageLayers = useCallback(() => {
        return samples.map(sample => {
            const imageSize = imageSizes[sample.id];
            const offset = sampleOffsets[sample.id] || [0, 0];
            if (!imageSize || imageSize.length < 2) return null;
            return new BitmapLayer({
                id: `marker-image-${sample.id}`,
                visible: visibleSamples[sample.id],
                image: `/${sample.id}_cells_layer.png`,
                bounds: [offset[0], offset[1] + imageSize[1], offset[0] + imageSize[0], offset[1]],
                opacity: 1,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [samples, imageSizes, visibleSamples]);

    // nucleus coordinates scatterplot layers
    const generateCellLayers = useCallback(() => {
        return samples.map(sample => {
            const sampleId = sample.id;
            const useGeneData = sampleId === samples[0]?.id &&
                selectedGenes.length > 0 &&
                geneExpressionData.length > 0;

            const analyzedRegionCells = analyzedRegionData?.cellIds || [];
            const isPartMode = partWholeMode && analyzedRegionCells.length > 0;

            let data;
            let getFillColor;

            if (useGeneData) {
                let regionLayers = [];
                if (isPartMode) {
                    const inRegionData = geneExpressionData.filter(d => analyzedRegionCells.includes(d.id));
                    const outRegionData = filteredCellData[sampleId].filter(d => !analyzedRegionCells.includes(d.id));

                    // kosara path layer for inRegionData
                    if (inRegionData.length > 0) {
                        const optimizedPathData = inRegionData.flatMap(d => {
                            const angles = Object.entries(d.angles);
                            const ratios = Object.entries(d.ratios);
                            const radius = Object.entries(d.radius);
                            const offset = sampleOffsets[d.sampleId] || [0, 0];

                            return generateKosaraPath(
                                d.cell_x + offset[0],
                                d.cell_y + offset[1],
                                angles,
                                ratios,
                                radius
                            ).map(path => ({
                                id: d.id,
                                cell_type: d.cell_type,
                                points: path.path,
                                color: path.color,
                                total_expression: d.total_expression,
                                ratios: d.ratios,
                            }));
                        });

                        regionLayers.push(new PolygonLayer({
                            id: `Scatters-${sampleId}`,
                            data: optimizedPathData,
                            getPolygon: d => d.points,
                            getFillColor: d => {
                                const rgbColor = hexToRgbs(d.color);
                                return [...rgbColor, 255];
                            },
                            pickable: true,
                            stroked: false,
                            parameters: { depthTest: false, blend: true },
                            updateTriggers: {
                                data: [inRegionData, sampleOffsets]
                            }
                        }));
                    }

                    // scatterplot layer for outRegionData
                    if (outRegionData.length > 0) {
                        const offset = sampleOffsets[sampleId] || [0, 0];
                        regionLayers.push(new ScatterplotLayer({
                            id: `Scatters-${sampleId}-out`,
                            data: outRegionData,
                            getPosition: d => [d.cell_x + offset[0], d.cell_y + offset[1]],
                            getFillColor: d => {
                                const defaultColor = colorMaps[sampleId]?.[d.cell_type] || [0, 0, 0];
                                return [...defaultColor, 255];
                            },
                            getRadius: 5,
                            pickable: true,
                            parameters: { depthTest: false },
                            updateTriggers: {
                                getFillColor: [colorMaps[sampleId], visibleCellTypes[sampleId]],
                                data: [outRegionData]
                            }
                        }));
                    }
                } else {
                    if (geneExpressionData.length > 0) {
                        const optimizedPathData = geneExpressionData.flatMap(d => {
                            const angles = Object.entries(d.angles);
                            const ratios = Object.entries(d.ratios);
                            const radius = Object.entries(d.radius);
                            const offset = sampleOffsets[d.sampleId] || [0, 0];

                            return generateKosaraPath(
                                d.cell_x + offset[0],
                                d.cell_y + offset[1],
                                angles,
                                ratios,
                                radius
                            ).map(path => ({
                                id: d.id,
                                cell_type: d.cell_type,
                                points: path.path,
                                color: path.color,
                                total_expression: d.total_expression,
                                ratios: d.ratios,
                            }));
                        });

                        regionLayers.push(new PolygonLayer({
                            id: `Scatters-${sampleId}`,
                            data: optimizedPathData,
                            getPolygon: d => d.points,
                            getFillColor: d => {
                                const rgbColor = hexToRgbs(d.color);
                                return [...rgbColor, 255];
                            },
                            pickable: true,
                            stroked: false,
                            parameters: { depthTest: false, blend: true },
                            updateTriggers: {
                                data: [geneExpressionData, sampleOffsets]
                            }
                        }));
                    }
                }
                return regionLayers;
            } else {
                data = filteredCellData[sampleId] || [];

                getFillColor = d => {
                    const defaultColor = colorMaps[sampleId]?.[d.cell_type] || [0, 0, 0];

                    if (analyzedRegionData && sampleId === analyzedRegionData.sampleId) {
                        const isInAnalyzedRegion = analyzedRegionData.cellIds.includes(d.id);

                        if (isInAnalyzedRegion) {
                            if (NMFclusterCells && NMFclusterCells.length > 0) {
                                const isHighlighted = NMFclusterCells.includes(d.id);
                                return isHighlighted
                                    ? [...defaultColor, 255]
                                    : [...defaultColor, 50];
                            }

                            return [...defaultColor, 255];
                        }
                    }

                    return [...defaultColor, 255];
                };

                const offset = sampleOffsets[sampleId] || [0, 0];
                return new ScatterplotLayer({
                    id: `Scatters-${sampleId}`,
                    visible: visibleSamples[sampleId],
                    data: data,
                    getPosition: d => [d.cell_x + offset[0], d.cell_y + offset[1]],
                    getRadius: 5,
                    getFillColor: getFillColor,
                    pickable: true,
                    parameters: { depthTest: false },
                    updateTriggers: {
                        getFillColor: [colorMaps[sampleId], visibleCellTypes[sampleId], NMFclusterCells, analyzedRegionData],
                        data: [data]
                    }
                });
            }
        }).filter(Boolean);
    }, [samples, filteredCellData, colorMaps, visibleSamples, sampleOffsets, visibleCellTypes, geneExpressionData, selectedGenes, analyzedRegion, NMFclusterCells, partWholeMode, regions]);

    // update current region data
    const handleRegionUpdate = (sampleId, updatedData) => {
        if (!updatedData || !Array.isArray(updatedData.features)) {
            console.error('Invalid data:', updatedData);
            return;
        }
        setFeatures(prev => ({
            ...prev,
            [sampleId]: updatedData
        }));
        const lastFeature = updatedData.features[updatedData.features.length - 1];
        if (lastFeature) {
            const offset = sampleOffsets[sampleId] || [0, 0];
            const localFeature = {
                ...lastFeature,
                geometry: {
                    ...lastFeature.geometry,
                    coordinates: lastFeature.geometry.coordinates.map(ring =>
                        ring.map(([x, y]) => [x - offset[0], y - offset[1]])
                    )
                }
            };
            setTempRegions(prev => ({
                ...prev,
                [sampleId]: {
                    ...prev[sampleId],
                    data: cellTypeCoordinatesData?.[sampleId]?.filter(cell =>
                        booleanPointInPolygon([cell.cell_x, cell.cell_y], localFeature.geometry)
                    )
                }
            }));
            setActiveDrawingSample(sampleId);
        }
    };

    // Showing genes in the whole region or the selected region
    const changeGeneShowRange = () => {
        setPartWholeMode(!partWholeMode);
    };

    // save region, save region name and region color, and update region color to new random color
    const handleSaveRegion = () => {
        if (!regionName) {
            message.error('Please enter a region name');
            return;
        }

        let drawnFeature = null;
        for (const sample of samples) {
            if (features[sample.id]?.features?.length > 0) {
                drawnFeature = features[sample.id].features[features[sample.id].features.length - 1];
                break;
            }
        }

        if (!drawnFeature) {
            message.error('No region available to save');
            return;
        }

        const centroidFeature = centroid(drawnFeature);
        const [cx, cy] = centroidFeature.geometry.coordinates;

        let targetSampleId = null;
        for (const sample of samples) {
            const offset = sampleOffsets[sample.id] || [0, 0];
            const size = imageSizes[sample.id];
            if (size && size.length === 2) {
                const [w, h] = size;
                if (cx >= offset[0] && cx <= offset[0] + w && cy >= offset[1] && cy <= offset[1] + h) {
                    targetSampleId = sample.id;
                    break;
                }
            }
        }

        if (!targetSampleId) {
            message.error('Drawn region does not fall within any sample');
            return;
        }

        const targetOffset = sampleOffsets[targetSampleId] || [0, 0];
        const localFeature = {
            ...drawnFeature,
            geometry: {
                ...drawnFeature.geometry,
                coordinates: drawnFeature.geometry.coordinates.map(ring =>
                    ring.map(([x, y]) => [x - targetOffset[0], y - targetOffset[1]])
                )
            }
        };

        const cellsInRegion = (cellTypeCoordinatesData[targetSampleId] || []).filter(cell =>
            booleanPointInPolygon([cell.cell_x, cell.cell_y], localFeature.geometry)
        );
        const cellIds = cellsInRegion.map(cell => cell.id);

        const newRegion = {
            id: `${targetSampleId}-${regionName}-${Date.now()}`,
            name: regionName,
            color: regionColor,
            sampleId: targetSampleId,
            feature: {
                type: 'FeatureCollection',
                features: [localFeature]
            },
            cellIds: cellIds
        };

        setRegions(prev => [...prev, newRegion]);
        const newFeatures = {};
        samples.forEach(sample => {
            newFeatures[sample.id] = { type: 'FeatureCollection', features: [] };
        });
        setFeatures(newFeatures);
        setIsDrawingActive(false);
        message.success('Region saved successfully');
        setRegionName('');
        setRegionColor(generateRandomColor());
    };

    // generate EditableGeoJsonLayer: only make the layer of the current activeDrawingSample visible
    const generateEditLayers = useCallback(() => {
        return samples.map(sample => {
            const tempRegion = tempRegions[sample.id] || {};
            return new EditableGeoJsonLayer({
                id: `edit-layer-${sample.id}`,
                data: features[sample.id] || { type: 'FeatureCollection', features: [] },
                mode: isDrawingActive ? DrawPolygonMode : undefined,
                _subLayerProps: {
                    'polygons-stroke': {
                        getCursor: () => 'crosshair'
                    }
                },
                selectedFeatureIndexes: [],
                onEdit: ({ updatedData }) => {
                    handleRegionUpdate(sample.id, updatedData);
                },
                visible: isDrawingActive && sample.id === activeDrawingSample,
                getLineColor: tempRegion.color ? [...tempRegion.color, 200] : [255, 0, 0, 200],
                getFillColor: tempRegion.color ? [...tempRegion.color, 50] : [255, 255, 0, 50],
                lineWidthMinPixels: 2
            });
        });
    }, [samples, features, tempRegions, isDrawingActive, activeDrawingSample]);

    // Calculate the number of cells in the region
    const getCellCount = (region) => {
        return region.cellIds ? region.cellIds.length : 0;
    };

    const handleDeleteRegion = (regionId) => {
        setRegions(prev => prev.filter(region => region.id !== regionId));
        message.success('Region deleted!');
    };

    const layers = useMemo(() => {
        const regionLabelLayer = new TextLayer({
            id: 'region-labels',
            data: regions
                .map(region => {
                    if (!region.feature?.features) return null;
                    const offset = sampleOffsets[region.sampleId] || [0, 0];
                    const globalFeature = {
                        ...region.feature.features[0],
                        geometry: {
                            ...region.feature.features[0].geometry,
                            coordinates: region.feature.features[0].geometry.coordinates.map(ring =>
                                ring.map(([x, y]) => [x + offset[0], y + offset[1]])
                            )
                        }
                    };
                    const center = centroid(globalFeature);
                    return {
                        position: center.geometry.coordinates,
                        text: region.name
                    };
                })
                .filter(Boolean),
            getPosition: d => d.position,
            getText: d => d.text,
            getColor: () => [0, 0, 0, 255],
            getSize: 16,
            sizeUnits: 'pixels'
        });
        return [
            ...generateWholePngLayers(),
            ...generateTileLayers(),
            ...generateMarkerImageLayers(),
            ...generateEditLayers(),
            regionLabelLayer,
            ...regions.map(region => {
                const offset = sampleOffsets[region.sampleId] || [0, 0];
                const globalFeature = {
                    ...region.feature,
                    features: region.feature.features.map(f => ({
                        ...f,
                        geometry: {
                            ...f.geometry,
                            coordinates: f.geometry.coordinates.map(ring =>
                                ring.map(([x, y]) => [x + offset[0], y + offset[1]])
                            )
                        },
                        properties: {
                            ...(f.properties || {}),
                            __regionMeta: {
                                id: region.id,
                                name: region.name,
                                sampleId: region.sampleId,
                                cell_ids: region.cellIds,
                            }
                        },
                    }))
                };
                return new GeoJsonLayer({
                    id: `Selected-region-${region.id}`,
                    data: globalFeature,
                    pickable: true,
                    stroked: true,
                    filled: true,
                    lineWidthMinPixels: 2,
                    getLineColor: () => [...region.color, 200],
                    getFillColor: () => [...region.color, 30],
                });
            }),
            ...generateCellLayers(),
        ].filter(Boolean);
    }, [
        generateWholePngLayers,
        generateTileLayers,
        generateMarkerImageLayers,
        generateCellLayers,
        generateEditLayers,
        regions,
        sampleOffsets
    ]);

    return (
        <div style={{ width: '100%', height: '100%', display: 'flex', borderRight: '2px solid #e8e8e8' }}>
            <div
                style={{
                    flex: 1,
                    position: 'relative',
                    cursor: isDrawingActive ? 'crosshair' : (hoveredCell ? 'pointer' : 'grab')
                }}
                onMouseEnter={() => {
                    const canvas = document.querySelector('.deckgl-wrapper canvas');
                    if (canvas) {
                        canvas.style.cursor = isDrawingActive ? 'crosshair' : (hoveredCell ? 'pointer' : 'grab');
                    }
                }}
            >
                <DeckGL
                    layers={layers}
                    views={new OrthographicView({ controller: true })}
                    initialViewState={{
                        target: (() => {
                            const firstSample = samples[0];
                            if (!firstSample) return [0, 0, 0];
                            const sampleId = firstSample.id;
                            const offset = sampleOffsets[sampleId] ?? [0, 0];
                            const size = imageSizes[sampleId] ?? [0, 0];
                            return [
                                offset[0] + size[0] / 2,
                                offset[1] + size[1] / 2,
                                0
                            ];
                        })(),
                        zoom: -3,
                        maxZoom: 2.5,
                        minZoom: -5
                    }}
                    onViewStateChange={({ viewState }) => {
                        debouncedSetZoom(viewState.zoom);
                    }}
                    onClick={info => {
                        if (info.object && info.layer.id.startsWith('Selected-')) {
                            // fetchGeneExpressionData(info.object.properties.__regionMeta.sampleId, info.object.properties.__regionMeta.cell_ids);
                            setAnalyzedRegion(info.object.properties.__regionMeta.name);
                            fetchNMFGOExpressionData(info.object.properties.__regionMeta.sampleId, info.object.properties.__regionMeta.cell_ids);
                            // fetchCell2CellInteractionData(info.object.properties.__regionMeta.name, info.object.properties.__regionMeta.sampleId, info.object.properties.__regionMeta.cell_ids);
                        }
                    }}
                    onHover={info => {
                        if (info.object && info.layer.id.startsWith('Scatters-')) {
                            const sampleId = info.layer.id.split('-')[1];
                            const isGeneDataMode = selectedGenes.length > 0 && geneExpressionData.length > 0;

                            if (isGeneDataMode) {
                                const id = info.object.id;
                                const ratios = info.object["ratios"]; // {gene name1: expression,...}
                                const cell_type = info.object["cell_type"];
                                const total_expression = info.object["total_expression"];
                                setHoveredCell({
                                    id,
                                    sampleId,
                                    cell_type,
                                    ratios,
                                    total_expression,
                                    x: info.x,
                                    y: info.y
                                });
                            } else {
                                setHoveredCell({
                                    ...info.object,
                                    sampleId,
                                    x: info.x,
                                    y: info.y
                                });
                            }
                        } else {
                            setHoveredCell(null);
                        }
                    }}
                    controller={{
                        inertia: true,
                        scrollZoom: {
                            speed: 0.05,
                            smooth: true
                        },
                        dragPan: !isDrawingActive,
                        dragRotate: false
                    }}
                />

                {/* tooltip */}
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
                        {hoveredCell.id ? (
                            <>
                                <div><strong>Sample:</strong> {hoveredCell.sampleId}</div>
                                <div><strong>Cell Type:</strong> {hoveredCell.cell_type}</div>
                                {hoveredCell.total_expression && (
                                    <div><strong>Total Expression:</strong> {hoveredCell.total_expression}</div>
                                )}
                                {hoveredCell.ratios &&
                                    Object.entries(hoveredCell.ratios).map(([gene, expression]) => (
                                        <div key={gene}>
                                            <strong>{gene}:</strong> {expression?.toFixed(2)}
                                        </div>
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

                {/* Sample list Collapse */}
                <div style={{ position: 'absolute', top: 10, left: 10, zIndex: 10 }}>
                    <Collapse items={collapseItems} defaultActiveKey={['0']} style={{ background: '#ffffff', width: 300, opacity: 0.8 }} />
                </div>

                {/* rightCornerControls: Saved regions and drawing controls */}
                <div style={{
                    position: 'absolute',
                    top: 10,
                    right: 10,
                    zIndex: 10,
                    width: 250
                }}>
                    {/* region drawing controls */}
                    <div style={{ background: 'rgba(255,255,255,0.8)', padding: '10px', borderRadius: '10px', display: 'flex', flexDirection: 'column' }}>
                        <Input
                            size='small'
                            allowClear
                            placeholder="New region name"
                            value={regionName}
                            onChange={e => setRegionName(e.target.value)}
                            style={{ marginBottom: 8 }}
                        />
                        <div style={{ display: 'flex', gap: 5 }}>
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
                                variant="outlined"
                                onClick={() => {
                                    if (!isDrawingActive) {
                                        setActiveDrawingSample(activeSample);
                                        setFeatures(prev => ({
                                            ...prev,
                                            [activeSample]: { type: 'FeatureCollection', features: [] }
                                        }));
                                        setIsDrawingActive(true);
                                    } else {
                                        handleSaveRegion();
                                    }
                                }}
                                block
                            >
                                {isDrawingActive && regionName ? 'Save' : 'Draw'}
                            </Button>
                        </div>
                    </div>

                    {/* saved regions */}
                    <Collapse style={{ marginTop: 10, background: 'rgba(255,255,255,0.8)' }}>
                        <Collapse.Panel header={`Selected Region (${regions.length})`} key="selected-region">
                            {/* <Switch
                                defaultChecked
                                onChange={changeGeneShowRange}
                                checked={partWholeMode}
                                checkedChildren="Part"
                                unCheckedChildren="Whole"
                                style={{
                                    backgroundColor: partWholeMode ? '#ED9121' : '#74C365',
                                    width: '100%',
                                    marginBottom: 10
                                }}
                            />
                            <TreeSelect
                                showSearch
                                style={{
                                    width: '100%',
                                    marginBottom: 10
                                }}
                                value={selectedGenes}
                                dropdownStyle={{
                                    maxHeight: 400,
                                    overflow: 'auto',
                                }}
                                size='small'
                                placeholder="Select genes"
                                allowClear
                                multiple
                                onChange={setSelectedGenes}
                                treeData={geneList}
                            />
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 5, marginBottom: 5 }}>
                                <Button size='small' style={{ width: '50%' }} onClick={cleanGeneSelection}>Clear</Button>
                                <Button size='small' style={{ width: '50%' }} onClick={confirmKosaraPlot}>Confirm</Button>
                            </div> */}
                            {regions.length > 0 && (
                                <>
                                    {regions.map(region => {
                                        const sampleId = region.sampleId;
                                        const cellCount = getCellCount(region);
                                        return (
                                            <div key={region.id} style={{
                                                display: 'flex',
                                                justifyContent: 'space-between',
                                                alignItems: 'center',
                                                marginBottom: 4,
                                                fontSize: '12px',
                                                color: '#555'
                                            }}>
                                                <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'flex-start', justifyContent: 'center' }}>
                                                    <span style={{ fontSize: 14, marginBottom: 10 }}>{region.name}</span>
                                                    <span style={{ fontSize: 10, color: "#666" }}>Sample: {sampleId}</span>
                                                    <span style={{ fontSize: 10, color: "#666" }}>Cells: {cellCount}</span>
                                                </div>
                                                <CloseOutlined
                                                    onClick={(e) => { e.stopPropagation(); handleDeleteRegion(region.id); }}
                                                    style={{ fontSize: '12px' }}
                                                />
                                            </div>
                                        );
                                    })}
                                </>
                            )}
                            {regions.length === 0 && (
                                <div style={{ fontSize: '12px', color: '#999' }}>No regions selected</div>
                            )}
                        </Collapse.Panel>
                    </Collapse>
                </div>
            </div>
        </div>
    );
};

// cell types and related counts display
const CellTypeSettings = ({
    cellTypes,
    cellData,
    colorMap,
    visibleMap,
    onColorChange,
    onVisibilityCellTypeChange,
}) => {
    const [searchText, setSearchText] = useState('');

    const filteredTypes = useMemo(() => {
        return (cellTypes || [])
            .map(type => ({
                type,
                count: (cellData || []).filter(cell => cell.cell_type === type).length
            }))
            .filter(item => item.type.toLowerCase().includes(searchText.toLowerCase()))
            .sort((a, b) => b.count - a.count);
    }, [cellTypes, cellData, searchText]);

    return (
        <div style={{ maxHeight: 400, overflowY: 'auto' }}>
            <Input.Search
                size='small'
                placeholder="Search cell types"
                value={searchText}
                onChange={e => setSearchText(e.target.value)}
                style={{ marginBottom: 5 }}
            />

            {filteredTypes.map(({ type, count }) => (
                <div key={type} style={{
                    display: 'flex',
                    alignItems: 'center',
                    padding: '8px 0',
                    borderBottom: '1px solid #f0f0f0'
                }}>
                    <Checkbox
                        checked={visibleMap[type] ?? true}
                        onChange={e => onVisibilityCellTypeChange(type, e.target.checked)}
                        style={{ marginRight: 8 }}
                    />
                    <ColorPicker
                        size="small"
                        value={`rgb(${(colorMap[type] || [0, 0, 0]).join(',')})`}
                        onChange={color => {
                            const rgb = color.toRgb();
                            onColorChange(type, [rgb.r, rgb.g, rgb.b]);
                        }}
                        style={{ marginRight: 12 }}
                    />
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', width: '100%' }}>
                        <span style={{ fontSize: 14 }}>{type}</span>
                        <span style={{ fontSize: 12, color: '#666' }}>{count}</span>
                    </div>
                </div>
            ))}
        </div>
    );
};

// gene names display
const GeneSettings = ({ geneList, sampleId, onVisibilityGeneChange, cleanGeneSelection, confirmKosaraPlot }) => {
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
                        <Checkbox style={{ marginRight: 8 }}
                            onChange={e => onVisibilityGeneChange(gene)}
                        />
                        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', width: '100%' }}>
                            <span style={{ fontSize: 12 }}>{gene}</span>
                            <span style={{ fontSize: 12, color: '#666' }}>{count}</span>
                        </div>
                    </div>
                ))}
            </div>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 5 }}>
                <Button size='small' style={{ width: '50%' }} onClick={cleanGeneSelection}>Clear</Button>
                <Button size='small' style={{ width: '50%' }} onClick={() => confirmKosaraPlot(sampleId)}>Confirm</Button>
            </div>
        </div>
    );
};
