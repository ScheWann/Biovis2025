import React, { useRef, useMemo, useState, useEffect } from 'react';
import DeckGL from '@deck.gl/react';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer } from '@deck.gl/layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData }) => {
    const viewerRef = useRef(null);
    const [tifImages, setTifImages] = useState([]);
    const [imageSize, setImageSize] = useState([]);
    const [tileUrls, setTileUrls] = useState([]);

    useEffect(() => {
        fetch("/get_tiles", {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ sample_id: sampleId })
        })
            .then((response) => response.json())
            .then((data) => {
                console.log(data, 'iopiopiop');
                setTileUrls(data);
            });
    }, []);

    // Loading TIF image
    useEffect(() => {
        const loadTifs = async () => {
            const MAX_CONCURRENT_REQUESTS = 500;

            const loadTile = async (tileUrl) => {
                const response = await fetch(tileUrl);
                const blob = await response.blob();
                const tiff = await fromBlob(blob);
                const image = await tiff.getImage();
                const width = image.getWidth();
                const height = image.getHeight();
                const rgba = await image.readRasters({ interleave: true });

                const canvas = document.createElement('canvas');
                canvas.width = width;
                canvas.height = height;
                const ctx = canvas.getContext('2d');
                ctx.translate(width, 0);
                ctx.scale(-1, 1);

                const imageData = new ImageData(new Uint8ClampedArray(rgba), width, height);
                ctx.putImageData(imageData, 0, 0);

                const tileName = tileUrl.split('/').pop();
                const match = tileName.match(/^tile_(\d+)_(\d+)\.tif$/);
                if (match) {
                    const x = parseInt(match[1], 10);
                    const y = parseInt(match[2], 10);

                    const bounds = [
                        [x, y + 256],
                        [x, y],
                        [x + 256, y],
                        [x + 256, y + 256]
                    ];

                    return {
                        image: canvas,
                        bounds: bounds,
                    };
                }
                return null;
            };

            const results = [];
            for (let i = 0; i < tileUrls.length; i += MAX_CONCURRENT_REQUESTS) {
                const chunk = tileUrls.slice(i, i + MAX_CONCURRENT_REQUESTS);
                const chunkResults = await Promise.all(
                    chunk.map(loadTile)
                );
                results.push(...chunkResults.filter(Boolean));
            }

            setTifImages(results);
        };

        if (tileUrls.length > 0) {
            loadTifs();
        }
    }, [tileUrls, imageSize]);

    // Get image size(width, height)
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

    useEffect(() => {
        console.log("tifImages updated:", tifImages);
    }, [tifImages]);

    const layers = [
        // TIF tissue image layer
        ...tifImages.map(({ image, bounds }, index) =>
            new BitmapLayer({
                id: `tif-layer-${index}`,
                image,
                bounds,
                opacity: 1,
            })
        ),

        // Cell boundary mask layer
        imageSize.length > 0 && new BitmapLayer({
            id: 'mask-layer',
            image: maskUrl,
            bounds: [0, imageSize[1], imageSize[0], 0],
            opacity: 0.05,
        }),

        // Nucleus layer
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