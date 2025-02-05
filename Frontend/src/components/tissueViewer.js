import React, { useRef, useState, useEffect } from 'react';
import DeckGL from '@deck.gl/react';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer } from '@deck.gl/layers';
import { TileLayer } from '@deck.gl/geo-layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData }) => {
    const viewerRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [tileSize] = useState(256);

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
    ].filter(Boolean);

    const view = new OrthographicView({
        id: 'ortho-view',
        controller: true,
    });

    return (
        <div ref={viewerRef} style={{ height: '100%', width: '50%' }}>
            <DeckGL
                style={{ width: '50%', height: '100%' }}
                layers={layers}
                views={view}
                initialViewState={{
                    target: [imageSize[0] / 2, imageSize[1] / 2, 0],
                    zoom: 1,
                }}
                controller={true}
            />
        </div>
    );
};