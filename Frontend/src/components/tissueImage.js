import React, { useEffect, useRef, useState } from 'react';
import { RollbackOutlined } from "@ant-design/icons";
import { Button } from "antd";
import { Map, View } from 'ol';

import ImageLayer from 'ol/layer/Image';
import ImageStatic from 'ol/source/ImageStatic';

import 'ol/ol.css';
import hireImage from '../data/tissue_hires_image.png';

export const TissueImage = () => {
    const mapRef = useRef(null);
    const [view, setView] = useState(null);

    useEffect(() => {
        const extent = [0, 0, 4000, 4000];

        const mapView = new View({
            center: [2000, 2000],
            zoom: 1,
            extent,
        });

        const map = new Map({
            target: mapRef.current,
            layers: [
                new ImageLayer({
                    source: new ImageStatic({
                        url: hireImage,
                        imageExtent: extent,
                    }),
                }),
            ],
            view: mapView,
        });

        setView(mapView);

        return () => map.setTarget(null);
    }, []);

    const resetZoom = () => {
        if (view) {
            view.setCenter([2000, 2000]);
            view.setZoom(1);
        }
    };

    return (
        <div style={{ position: 'relative', height: '100%', width: '50%' }}>
            <div ref={mapRef} style={{ height: '100%', width: '100%' }} />
            <Button
                onClick={resetZoom}
                style={{
                    position: 'absolute',
                    top: '10px',
                    right: '10px',
                    zIndex: 10,
                    padding: '10px 15px',
                    color: '#666',
                    borderRadius: '5px',
                    border: '1px solid #d9d9d9',
                    cursor: 'pointer',
                }}
                icon={<RollbackOutlined />}
            >
            </Button>
        </div>
    );
};
