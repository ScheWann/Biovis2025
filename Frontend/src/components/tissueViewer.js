import React, { useRef, useMemo, useState, useEffect } from 'react';
import DeckGL from '@deck.gl/react';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer } from '@deck.gl/layers';
import { fromBlob } from 'geotiff';
import maskUrl from '../data/cells_layer.png';
import tifUrl from '../data/wsi.tif';

export const TissueViewer = ({ sampleId, cellTypeCoordinatesData }) => {
    const viewerRef = useRef(null);
    const [tifImage, setTifImage] = useState(null);
    const [imageSize, setImageSize] = useState([]);
    // const [divSize, setDivSize] = useState({ width: 0, height: 0 });

    // useEffect(() => {
    //     if (viewerRef.current) {
    //         setDivSize({
    //             width: viewerRef.current.clientWidth,
    //             height: viewerRef.current.clientHeight
    //         });
    //     }
    // }, []);

    // Loading TIF image
    useEffect(() => {
        const loadTif = async () => {
            const response = await fetch(tifUrl);
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
            const imageData = new ImageData(new Uint8ClampedArray(rgba), width, height);
            ctx.putImageData(imageData, 0, 0);
            setTifImage(canvas);
        };

        loadTif();
    }, []);

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

    // const zoom = useMemo(() => {
    //     if (imageSize.length === 0 || divSize.width === 0 || divSize.height === 0) return 0;
    //     return Math.min(divSize.width / imageSize[0], divSize.height / imageSize[1]);
    // }, [imageSize, divSize]);

    const layers = [
        // TIF tissue image layer
        tifImage && imageSize.length > 0 &&
        new BitmapLayer({
            id: 'tif-layer',
            image: tifImage,
            bounds: [0, imageSize[1], imageSize[0], 0],
            opacity: 1,
        }),

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
            style={{ width: '100%', height: '100%' }}
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