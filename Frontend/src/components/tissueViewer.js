import React, { useRef, useState, useEffect } from 'react';
import DeckGL from '@deck.gl/react';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer } from '@deck.gl/layers';
import { booleanPointInPolygon } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import {
    EditableGeoJsonLayer,
    DrawPolygonMode
} from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData }) => {
    const viewerRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [tileSize] = useState(256);
    const [features, setFeatures] = useState({ type: 'FeatureCollection', features: [] });
    const [selectionMode, setSelectionMode] = useState(false);
    const [selectedCells, setSelectedCells] = useState([]);
    const [selectedFeatureIndexes] = useState([]);

    useEffect(() => {
        fetch("/get_hires_image_size", {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ sample_id: sampleId })
        })
            .then((response) => response.json())
            .then((data) => {
                setImageSize(data);
            });
    }, [sampleId]);

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

    const filterCellsInPolygon = (polygonFeature) => {
        if (!polygonFeature || !cellTypeCoordinatesData) return;

        const polygon = polygonFeature.geometry;
        const filtered = cellTypeCoordinatesData.filter(cell => {
            const point = [cell.cell_x, cell.cell_y];
            return booleanPointInPolygon(point, polygon);
        });

        setSelectedCells(filtered);
        console.log('Selected cells:', filtered);
    };

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

                const tileUrl = `/get_tile?sample_id=${sampleId}&x=${x}&y=${y}`;
                const response = await fetch(tileUrl);
                const blob = await response.blob();
                const tiff = await fromBlob(blob);
                const image = await tiff.getImage();
                const rgba = await image.readRasters({ interleave: true });

                const canvas = document.createElement('canvas');
                [canvas.width, canvas.height] = [image.getWidth(), image.getHeight()];
                const ctx = canvas.getContext('2d');

                ctx.translate(canvas.width, 0);
                ctx.scale(-1, 1);

                const imageData = new ImageData(new Uint8ClampedArray(rgba), canvas.width, canvas.height);
                ctx.putImageData(imageData, 0, 0);

                return { image: canvas, bounds };
            } catch (error) {
                console.error('loading tile failed:', error);
                return null;
            }
        },

        renderSubLayers: (props) => {
            const { data, tile, ...rest } = props;
            if (!data) return null;

            return new BitmapLayer(rest, {
                image: data.image,
                opacity: 1,
                bounds: data.bounds
            });
        }
    });

    const layers = [
        tileLayer,

        imageSize.length > 0 && new BitmapLayer({
            id: 'mask-layer',
            image: maskUrl,
            bounds: [0, imageSize[1], imageSize[0], 0],
            opacity: 0.05,
        }),

        imageSize.length > 0 && cellTypeCoordinatesData &&
        new ScatterplotLayer({
            id: 'cell-layer',
            data: cellTypeCoordinatesData,
            getPosition: (d) => [d.cell_x, d.cell_y],
            getRadius: 2,
            getFillColor: (d) => (d.cell_type === 'type1' ? [255, 0, 0] : [0, 0, 255]),
            pickable: true,
        }),
        editLayer,
    ].filter(Boolean);

    const view = new OrthographicView({
        id: 'ortho-view',
        controller: true,
    });

    return (
        <div ref={viewerRef} style={{ height: '100%', width: '50%', position: 'relative' }}>
            <DeckGL
                style={{ width: '100%', height: '100%' }}
                layers={layers}
                views={view}
                initialViewState={{
                    target: [imageSize[0] / 2, imageSize[1] / 2, 0],
                    zoom: -2,
                    maxZoom: 1,
                    minZoom: -5,
                }}
                controller={true}
                getCursor={({ isDragging }) =>
                    selectionMode ? 'crosshair' : isDragging ? 'grabbing' : 'grab'
                }
            />
            <div className='controls' style={{
                position: 'absolute',
                top: 10,
                right: 10,
                zIndex: 1
            }}>
                <button
                    onClick={() => setSelectionMode(!selectionMode)}
                    style={{
                        padding: '8px 12px',
                        background: selectionMode ? '#2196F3' : '#fff',
                        border: '1px solid #ccc',
                        cursor: 'pointer'
                    }}
                >
                    Polygon
                </button>
            </div>
        </div>
    );
};