import React, { useRef, useState, useEffect, useMemo } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox, message } from "antd";
import { CloseOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import { booleanPointInPolygon, centroid } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData, cellTypeDir, regions, setRegions }) => {
    const viewerRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [tileSize] = useState(256);
    const [features, setFeatures] = useState({ type: 'FeatureCollection', features: [] });
    const [tempRegion, setTempRegion] = useState({ name: '', data: [] });
    const [selectionMode, setSelectionMode] = useState(false);
    const [visibleCellTypes, setVisibleCellTypes] = useState({});
    const [selectedFeatureIndexes] = useState([]);
    const [searchText, setSearchText] = useState('');
    const [colorMap, setColorMap] = useState({});

    // intialize color map
    useEffect(() => {
        const initialColorMap = {};
        const initialVisibility = {};

        cellTypeDir?.forEach((cellType, index) => {
            const hue = (index * 360) / cellTypeDir.length;
            initialColorMap[cellType] = hslToRgb(hue, 100, 50);
            initialVisibility[cellType] = true;
        });

        setColorMap(initialColorMap);
        setVisibleCellTypes(initialVisibility);
    }, [cellTypeDir]);

    // get image sizes
    useEffect(() => {
        fetch("/get_hires_image_size", {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ sample_id: sampleId })
        })
            .then((response) => response.json())
            .then((data) => setImageSize(data));
    }, [sampleId]);

    // hsl to rgb
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

    // generate random area color
    const generateColor = () => {
        const hue = Math.floor(Math.random() * 360);
        return hslToRgb(hue, 70, 50);
    };

    // visible cell data
    const visibleCellData = useMemo(() => {
        return cellTypeCoordinatesData?.filter(cell =>
            visibleCellTypes[cell.cell_type]
        ) || [];
    }, [cellTypeCoordinatesData, visibleCellTypes]);

    // filter cell types
    const filteredTypes = useMemo(() => {
        return cellTypeDir?.filter(type =>
            type.toLowerCase().includes(searchText.toLowerCase())
        ) || [];
    }, [cellTypeDir, searchText]);

    // cell type counts
    const cellTypeCounts = useMemo(() => {
        const counts = {};
        cellTypeCoordinatesData?.forEach(cell => {
            counts[cell.cell_type] = (counts[cell.cell_type] || 0) + 1;
        });
        return counts;
    }, [cellTypeCoordinatesData]);

    // set region name
    const regionNameChange = (e) => {
        setTempRegion((prev) => ({ ...prev, name: e.target.value }));
    }

    const regionFunction = () => {
        if (!tempRegion.name.trim()) {
            message.error('Please name the region before drawing');
            return;
        }

        if (!selectionMode) {
            const newColor = generateColor();
            setTempRegion(prev => ({ ...prev, color: newColor }));
            setSelectionMode(true);
        } else {
            setSelectionMode(false);
            if (tempRegion.data.length > 0) {
                const lastFeature = features.features[features.features.length - 1];
                const center = centroid(lastFeature.geometry).geometry.coordinates;
                setRegions(prev => [...prev, {
                    ...tempRegion,
                    feature: {
                        ...lastFeature,
                        properties: {
                            color: tempRegion.color,
                            name: tempRegion.name
                        }
                    },
                    center: center
                }]);
            }
            setTempRegion({ name: '', data: [], color: null });
            setFeatures({ ...features, features: [] });
        }
    };

    // filter cells in polygon
    const filterCellsInPolygon = (polygonFeature) => {
        if (!polygonFeature || !cellTypeCoordinatesData) return;
        const polygon = polygonFeature.geometry;
        const filtered = cellTypeCoordinatesData.filter(cell =>
            booleanPointInPolygon([cell.cell_x, cell.cell_y], polygon)
        );
        setTempRegion((prev) => ({ ...prev, data: filtered }));
    };

    const deleteRegion = (index) => {
        setRegions(prev => prev.filter((_, i) => i !== index));
    };

    // tiles layer
    const tileLayer = imageSize.length > 0 && new TileLayer({
        id: 'tif-tiles-layer',
        tileSize,
        extent: [0, 0, imageSize[0], imageSize[1]],
        maxZoom: 0,
        minZoom: 0,
        getTileData: async ({ index }) => {
            try {
                const { x, y } = index;
                const bounds = [
                    [x * tileSize, (y + 1) * tileSize],
                    [x * tileSize, y * tileSize],
                    [(x + 1) * tileSize, y * tileSize],
                    [(x + 1) * tileSize, (y + 1) * tileSize]
                ];
                const response = await fetch(`/get_tile?sample_id=${sampleId}&x=${x}&y=${y}`);
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

                return { image: canvas, bounds };
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

    // editable layer
    const editLayer = useMemo(() => new EditableGeoJsonLayer({
        data: features,
        mode: DrawPolygonMode,
        selectedFeatureIndexes,
        onEdit: ({ updatedData, editType }) => {
            setFeatures(updatedData);
            if (editType === 'addFeature') {
                const lastFeature = updatedData.features[updatedData.features.length - 1];
                filterCellsInPolygon(lastFeature);
            }
        },
        visible: selectionMode,
        getLineColor: tempRegion.color ? [...tempRegion.color, 200] : [255, 0, 0, 200],
        getFillColor: tempRegion.color ? [...tempRegion.color, 50] : [255, 255, 0, 50],
        getEditHandleColor: () => tempRegion.color ? [...tempRegion.color, 200] : [255, 0, 0, 200],
        lineWidthMinPixels: 2
    }), [features, selectionMode, tempRegion.color]);

    const regionsLayer = regions?.length > 0
        ? new EditableGeoJsonLayer({
            id: 'saved-regions',
            data: {
                type: 'FeatureCollection',
                features: regions.map(r => r.feature)
            },
            mode: DrawPolygonMode,
            selectedFeatureIndexes: [],
            pickable: true,
            visible: true,
            getLineColor: d => [...d.properties.color, 200],
            getFillColor: d => [...d.properties.color, 50],
            lineWidthMinPixels: 2
        })
        : null;

    const textLayers = useMemo(() => regions.map(region =>
        new TextLayer({
            id: `text-layer-${region.name}`,
            data: [{
                position: region.center,
                text: region.name
            }],
            getColor: [0, 0, 0, 255],
            getSize: 16,
            getTextAnchor: 'middle',
            getAlignmentBaseline: 'center'
        })
    ), [regions]);

    // layers
    const layers = [
        tileLayer,
        imageSize.length > 0 && new BitmapLayer({
            id: 'mask-layer',
            image: maskUrl,
            bounds: [0, imageSize[1], imageSize[0], 0],
            opacity: 0.05,
        }),
        imageSize.length > 0 && visibleCellData && new ScatterplotLayer({
            id: 'cell-layer',
            data: visibleCellData,
            getPosition: d => [d.cell_x, d.cell_y],
            getRadius: 2,
            getFillColor: d => colorMap[d.cell_type] || [0, 0, 0],
            pickable: true,
        }),
        ...textLayers,
        regionsLayer,
        editLayer,
    ].filter(Boolean);

    return (
        <div ref={viewerRef} style={{ height: '100%', width: '70%', position: 'relative', border: '1px solid #f0f0f0' }}>
            <DeckGL
                layers={layers}
                views={new OrthographicView({ id: 'ortho-view', controller: true })}
                initialViewState={{
                    target: [imageSize[0] / 2, imageSize[1] / 2, 0],
                    zoom: -2,
                    maxZoom: 1,
                    minZoom: -5
                }}
                controller={true}
                getCursor={({ isDragging }) => selectionMode ? 'crosshair' : isDragging ? 'grabbing' : 'grab'}
                style={{ width: '100%', height: '100%' }}
            />

            <Collapse
                size='small'
                defaultActiveKey={[]}
                items={[{
                    key: 'legend',
                    label: <span style={{ fontWeight: 500 }}>Cell Types ({filteredTypes.length})</span>,
                    children: (
                        <div style={{ maxHeight: 400, overflowY: 'auto' }}>
                            <Input.Search
                                placeholder="Search cell types"
                                value={searchText}
                                onChange={e => setSearchText(e.target.value)}
                                style={{ marginBottom: 12 }}
                            />

                            {filteredTypes.map(cellType => {
                                const rgbColor = colorMap[cellType] || [0, 0, 0];
                                return (
                                    <div key={cellType} style={{
                                        display: 'flex',
                                        alignItems: 'center',
                                        padding: '8px 0',
                                        borderBottom: '1px solid #f0f0f0'
                                    }}>
                                        <Checkbox
                                            checked={visibleCellTypes[cellType] ?? true}
                                            onChange={(e) => {
                                                setVisibleCellTypes(prev => ({
                                                    ...prev,
                                                    [cellType]: e.target.checked
                                                }));
                                            }}
                                            style={{ marginRight: 8 }}
                                        />

                                        <ColorPicker
                                            size="small"
                                            value={`rgb(${rgbColor.join(',')})`}
                                            onChange={color => {
                                                const rgb = color.toRgb();
                                                setColorMap(prev => ({
                                                    ...prev,
                                                    [cellType]: [rgb.r, rgb.g, rgb.b]
                                                }));
                                            }}
                                            style={{ marginRight: 12 }}
                                        />

                                        <div style={{ flex: 1 }}>
                                            <div style={{ fontSize: 12 }}>{cellType}</div>
                                            <div style={{ fontSize: 12, color: '#666' }}>
                                                {cellTypeCounts[cellType] || 0} cells
                                            </div>
                                        </div>
                                    </div>
                                );
                            })}
                        </div>
                    )
                }]}
                style={{
                    opacity: 0.7,
                    width: 280,
                    background: '#ffffff',
                    borderRadius: 8,
                    boxShadow: '0 2px 8px rgba(0,0,0,0.1)',
                    position: 'absolute',
                    top: 10,
                    left: 10,
                    zIndex: 1
                }}
                bordered={false}
            />
            <div className="controls" style={{ position: 'absolute', top: 10, right: 10, zIndex: 1 }}>
                <div style={{ display: "flex", justifyContent: "center", alignItems: "center", gap: 5 }}>
                    <span style={{ fontSize: 15, color: '#333' }}>Region: </span>
                    <Input value={tempRegion.name} onChange={regionNameChange} />
                    <Button
                        onClick={regionFunction}
                        style={{
                            opacity: 0.7,
                            background: selectionMode ? '#1890ff' : '#fff',
                        }}
                    >
                        {selectionMode ? 'Confirm' : 'Draw'}
                    </Button>
                </div>

                <Collapse
                    size='small'
                    style={{ background: '#ffffff', marginTop: 5, opacity: 0.7 }}
                    items={[{
                        key: 'regions',
                        label: <span style={{ fontWeight: 500 }}>Regions ({regions.length})</span>,
                        children: (
                            <div style={{ maxHeight: 400, overflowY: 'auto' }}>
                                {regions.map((region, index) => (
                                    <div key={index} style={{
                                        display: 'flex',
                                        justifyContent: 'space-between',
                                        alignItems: 'center',
                                        padding: '8px 0',
                                        borderBottom: '1px solid #f0f0f0'
                                    }}>
                                        <div>
                                            <div style={{ fontSize: 12 }}>{region.name}</div>
                                            <div style={{ fontSize: 12, color: '#666' }}>
                                                {region.data.length} cells
                                            </div>
                                        </div>
                                        <Button
                                            danger
                                            size="small"
                                            onClick={() => deleteRegion(index)}
                                            icon={<CloseOutlined />}
                                        />
                                    </div>
                                ))}
                            </div>
                        )
                    }]}
                />
            </div>
        </div>
    );
};