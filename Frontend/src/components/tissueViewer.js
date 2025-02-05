import React, { useRef, useState, useEffect, useMemo } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox } from "antd";
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer } from '@deck.gl/layers';
import { booleanPointInPolygon } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData, cellTypeDir }) => {
    const viewerRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [tileSize] = useState(256);
    const [features, setFeatures] = useState({ type: 'FeatureCollection', features: [] });
    const [selectionMode, setSelectionMode] = useState(false);
    const [selectedCells, setSelectedCells] = useState([]);
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

    // editable layer
    const editLayer = new EditableGeoJsonLayer({
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
    });

    // filter cells in polygon
    const filterCellsInPolygon = (polygonFeature) => {
        if (!polygonFeature || !cellTypeCoordinatesData) return;
        const polygon = polygonFeature.geometry;
        const filtered = cellTypeCoordinatesData.filter(cell =>
            booleanPointInPolygon([cell.cell_x, cell.cell_y], polygon)
        );
        setSelectedCells(filtered);
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
        editLayer,
    ].filter(Boolean);

    return (
        <div ref={viewerRef} style={{ height: '100%', width: '75%', position: 'relative' }}>
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

            <div className="controls" style={{ position: 'absolute', top: 10, right: 10, zIndex: 1 }}>
                <Button
                    onClick={() => setSelectionMode(!selectionMode)}
                    style={{
                        background: selectionMode ? '#1890ff' : '#fff',
                        marginBottom: 10
                    }}
                >
                    {selectionMode ? 'Exit Selection' : 'Draw Polygon'}
                </Button>

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
                        width: 280,
                        background: '#ffffff',
                        borderRadius: 8,
                        boxShadow: '0 2px 8px rgba(0,0,0,0.1)'
                    }}
                    bordered={false}
                />
            </div>
        </div>
    );
};