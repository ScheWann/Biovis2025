import React, { useEffect, useRef, useState } from 'react';
import { RollbackOutlined } from "@ant-design/icons";
import { Button, Select } from "antd";
import { Map, View } from 'ol';
import 'ol/ol.css';
import "../styles/tissueImage.css";
import ImageLayer from 'ol/layer/Image';
import ImageStatic from 'ol/source/ImageStatic';
import hireImage from '../data/tissue_hires_image.png';


export const TissueImage = ({kmeansSize, setKmeansSize}) => {
    const mapRef = useRef(null);
    const [view, setView] = useState(null);

    // Kmeans number options
    const kmeansOptions = [
        {
            value: 2,
            label: '2',
        },
        {
            value: 3,
            label: '3',
        },
        {
            value: 4,
            label: '4',
        },
        {
            value: 5,
            label: '5',
        },
        {
            value: 6,
            label: '6',
        },
        {
            value: 7,
            label: '7',
        },
        {
            value: 8,
            label: '8',
        },
        {
            value: 9,
            label: '9',
        },
        {
            value: 10,
            label: '10',
        },
    ];

    const handleChange = (value) => {
        setKmeansSize(value);
    };

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
            <div className="controlButtonGroup">
                <Select
                    value={kmeansSize}
                    style={{
                        width: 120,
                    }}
                    onChange={handleChange}
                    options={kmeansOptions}
                />
                <Button
                    onClick={resetZoom}
                    icon={<RollbackOutlined />}
                >
                </Button>
            </div>
        </div>
    );
};
