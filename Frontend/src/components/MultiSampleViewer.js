import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox, Select, message } from "antd";
import { CloseOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, TextLayer, GeoJsonLayer } from '@deck.gl/layers';
import { booleanPointInPolygon, centroid } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';
import "../styles/MultiSampleViewer.css";


const stringToHash = (str) => {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
        const char = str.charCodeAt(i);
        hash = (hash << 5) - hash + char;
        hash |= 0;
    }
    return Math.abs(hash);
};

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
    samples,
    cellTypeCoordinatesData,
    cellTypeDir,
    regions,
    setRegions
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
    const [geneList, setGeneList] = useState([]);
    const [selectedGenes, setSelectedGenes] = useState([]);
    const [regionName, setRegionName] = useState('');
    const [regionColor, setRegionColor] = useState(generateRandomColor());
    const [isDrawingActive, setIsDrawingActive] = useState(false);
    const [activeDrawingSample, setActiveDrawingSample] = useState(null);
    const [currentZoom, setCurrentZoom] = useState(-3);
    const TILE_LOAD_ZOOM_THRESHOLD = -2;

    // fetch gene list
    useEffect(() => {
        const sampleNames = samples.map(item => item.id);

        fetch('/get_all_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_names: sampleNames })
        })
            .then(res => res.json())
            .then(data => {
                setGeneList(data);
            });
    }, [setGeneList]);

    // intialize color map and visible cell types for each sample
    useEffect(() => {
        const initialStates = samples.reduce((acc, sample) => {
            const initialColorMap = {};
            const initialVisibility = {};
            cellTypeDir?.forEach((cellType, index) => {
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
        const fetchImageSizes = async () => {
            const samplesToFetch = samples.filter(sample => !imageSizes.hasOwnProperty(sample.id));
            if (samplesToFetch.length === 0) return;

            const newSizes = {};
            await Promise.all(samplesToFetch.map(async (sample) => {
                const response = await fetch("/get_hires_image_size", {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ sample_id: sample.id })
                });
                const data = await response.json();
                newSizes[sample.id] = data;
            }));

            setImageSizes(prev => ({
                ...prev,
                ...newSizes
            }));
        };

        if (samples.length) {
            fetchImageSizes();
        }
    }, [samples, imageSizes]);

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
    }, [samples, imageSizes, visibleSamples, sampleOffsets, currentZoom]);

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
                opacity: 0.05,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [samples, imageSizes, visibleSamples]);

    // nucleus coordinates scatterplot layers
    const generateCellLayers = useCallback(() => {
        return samples.map(sample => {
            const data = filteredCellData[sample.id];
            if (!data || data.length === 0) return null;
            const offset = sampleOffsets[sample.id] || [0, 0];
            return new ScatterplotLayer({
                id: `cells-${sample.id}`,
                visible: visibleSamples[sample.id],
                data: data,
                getPosition: d => [d.cell_x + offset[0], d.cell_y + offset[1]],
                getRadius: 5,
                getFillColor: d => colorMaps[sample.id]?.[d.cell_type] || [0, 0, 0],
                pickable: true,
                parameters: { depthTest: false },
                updateTriggers: {
                    getFillColor: [colorMaps[sample.id], visibleCellTypes[sample.id]],
                    data: [data]
                }
            });
        }).filter(Boolean);
    }, [samples, filteredCellData, colorMaps, visibleSamples, sampleOffsets, visibleCellTypes]);

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

    const confirmKosaraPlot = () => {
        console.log(selectedGenes, regions, '//////')
        regions.forEach(region => {
            fetch('/get_kosara_data', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ sample_id: region.sampleId, gene_list: selectedGenes, cell_list: region.cellIds})
            })
                .then(res => res.json())
                .then(data => {
                    console.log(data, 'kosara data')
                });
        })
    }

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
            ...generateCellLayers(),
            ...generateEditLayers(),
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
                        }
                    }))
                };
                return new GeoJsonLayer({
                    id: `region-${region.id}`,
                    data: globalFeature,
                    stroked: true,
                    filled: true,
                    lineWidthMinPixels: 2,
                    getLineColor: () => [...region.color, 200],
                    getFillColor: () => [...region.color, 50],
                    pickable: false
                });
            }),
            regionLabelLayer
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
        <div style={{ height: '100%', display: 'flex' }}>
            <div className="controls" style={{ position: 'absolute', top: 10, left: 10, zIndex: 10 }}>
                <Collapse style={{ background: '#ffffff', width: 300, opacity: 0.8 }}>
                    {samples.map(sample => (
                        <Collapse.Panel key={sample.id} header={sample.name}>
                            <Button
                                block
                                onClick={() => setVisibleSamples(prev => ({
                                    ...prev,
                                    [sample.id]: !prev[sample.id]
                                }))}
                            >
                                {visibleSamples[sample.id] ? 'Hide Sample' : 'Show Sample'}
                            </Button>
                            <Collapse style={{ marginTop: 16 }}>
                                <Collapse.Panel key={`cell-types-${sample.id}`} header="Cell Types">
                                    <CellTypeSettings
                                        cellTypes={cellTypeDir}
                                        cellData={cellTypeCoordinatesData[sample.id]}
                                        colorMap={colorMaps[sample.id] || {}}
                                        visibleMap={visibleCellTypes[sample.id] || {}}
                                        onColorChange={(type, color) => {
                                            setColorMaps(prev => ({
                                                ...prev,
                                                [sample.id]: { ...prev[sample.id], [type]: color }
                                            }));
                                        }}
                                        onVisibilityChange={(type, visible) => {
                                            setVisibleCellTypes(prev => ({
                                                ...prev,
                                                [sample.id]: { ...prev[sample.id], [type]: visible }
                                            }));
                                        }}
                                    />
                                </Collapse.Panel>
                            </Collapse>
                        </Collapse.Panel>
                    ))}
                </Collapse>
            </div>

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
                    onHover={info => {
                        if (info.object && info.layer.id.startsWith('cells-')) {
                            const sampleId = info.layer.id.split('-')[1];
                            setHoveredCell({
                                ...info.object,
                                sampleId,
                                x: info.x,
                                y: info.y
                            });
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
                        transform: 'translateX(-50%)'
                    }}>
                        <div>Sample: {hoveredCell.sampleId}</div>
                        <div>Cell Type: {hoveredCell.cell_type}</div>
                    </div>
                )}

                {/* rightCornerControls: Saved regions and drawing controls */}
                <div style={{
                    position: 'absolute',
                    top: 10,
                    right: 10,
                    zIndex: 20,
                    width: 220
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
                        <div style={{ display: "flex", gap: 5 }}>
                            <ColorPicker
                                size='small'
                                value={`rgb(${regionColor.join(',')})`}
                                onChange={color => {
                                    const rgb = color.toRgb();
                                    setRegionColor([rgb.r, rgb.g, rgb.b]);
                                }}
                                style={{ marginBottom: 8 }}
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
                            {regions.length > 0 && (
                                <>
                                    <Select
                                        mode="multiple"
                                        size='small'
                                        placeholder="Select genes"
                                        options={geneList}
                                        value={selectedGenes}
                                        onChange={setSelectedGenes}
                                        filterOption={(input, option) =>
                                            option.label.toLowerCase().indexOf(input.toLowerCase()) >= 0
                                        }
                                        style={{ width: '100%', marginBottom: 10 }}
                                        showSearch
                                    />
                                    <Button onClick={confirmKosaraPlot}>confirm</Button>
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

// CellTypeSettings component with cell count display
const CellTypeSettings = ({
    cellTypes,
    cellData,
    colorMap,
    visibleMap,
    onColorChange,
    onVisibilityChange
}) => {
    const [searchText, setSearchText] = useState('');
    const filteredTypes = useMemo(() =>
        cellTypes?.filter(type =>
            type.toLowerCase().includes(searchText.toLowerCase())
        ) || []
        , [cellTypes, searchText]);
    return (
        <div style={{ maxHeight: 400, overflowY: 'auto' }}>
            <Input.Search
                size='small'
                placeholder="Search cell types"
                value={searchText}
                onChange={e => setSearchText(e.target.value)}
                style={{ marginBottom: 12 }}
            />
            {filteredTypes.map(type => {
                const count = (cellData || []).filter(cell => cell.cell_type === type).length;
                return (
                    <div key={type} style={{
                        display: 'flex',
                        alignItems: 'center',
                        padding: '8px 0',
                        borderBottom: '1px solid #f0f0f0'
                    }}>
                        <Checkbox
                            checked={visibleMap[type] ?? true}
                            onChange={e => onVisibilityChange(type, e.target.checked)}
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
                );
            })}
        </div>
    );
};
