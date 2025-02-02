import React, { use, useEffect, useRef, useState } from 'react';
import { RollbackOutlined, SelectOutlined } from "@ant-design/icons";
import { Button, Select } from "antd";
import { Map, View } from 'ol';
import * as d3 from "d3";
import 'ol/ol.css';
import "../styles/tissueImage.css";
import ImageLayer from 'ol/layer/Image';
import ImageStatic from 'ol/source/ImageStatic';
import VectorSource from 'ol/source/Vector';
import Feature from 'ol/Feature';
import Point from 'ol/geom/Point';
import WebGLVectorLayer from 'ol/layer/WebGLVector';
import VectorLayer from 'ol/layer/Vector';
import Draw from 'ol/interaction/Draw';
import Style from 'ol/style/Style';
import Fill from 'ol/style/Fill';
import Stroke from 'ol/style/Stroke';
import hireImage from '../data/cells_layer.png';


export const TissueImage = ({ sampleId, cellTypeCoordinatesData }) => {
    const mapRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [view, setView] = useState(null);
    const [map, setMap] = useState(null);

    // get hire image size(width, height)
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
    }, []);

    // create map
    useEffect(() => {
        if (imageSize.length > 0) {
            const extent = [0, 0, imageSize[0], imageSize[1]];

            const mapView = new View({
                center: [imageSize[0] / 2, imageSize[1] / 2],
                zoom: 1,
                extent,
            });

            const olMap = new Map({
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

            setMap(olMap);
            setView(mapView);

            return () => olMap.setTarget(null);
        }
    }, [imageSize]);

    useEffect(() => {
        if (!map || imageSize.length <= 0) return;
        // if (mode === 'kmeans') {
            const clusterColors = {
                1: 'rgba(255, 0, 0, 0.7)',
                2: 'rgba(0, 255, 0, 0.7)',
                3: 'rgba(0, 0, 255, 0.7)',
                4: 'rgba(255, 255, 0, 0.7)',
                5: 'rgba(0, 255, 255, 0.7)',
                6: 'rgba(255, 0, 255, 0.7)',
                7: 'rgba(128, 0, 128, 0.7)',
                8: 'rgba(255, 165, 0, 0.7)',
                9: 'rgba(128, 128, 128, 0.7)',
                10: 'rgba(0, 0, 0, 0.7)',
            };

            // vector layer
            const vectorSource = new VectorSource();

            cellTypeCoordinatesData.forEach((point) => {
                const { cell_x, cell_y, cell_type } = point;

                const adjustedX = cell_x;
                const adjustedY = imageSize[1] - cell_y;

                const feature = new Feature({
                    geometry: new Point([adjustedX, adjustedY]),
                });

                feature.set('cluster', clusterColors[cell_type]);

                vectorSource.addFeature(feature);
            });

            const webGLVectorLayer = new WebGLVectorLayer({
                source: vectorSource,
                style: {
                    'circle-radius': 1,
                    'circle-fill-color': ['get', 'cluster'],
                },
            });

            // clear old vector layers
            map.getLayers().forEach((layer) => {
                if (layer instanceof WebGLVectorLayer) {
                    map.removeLayer(layer);
                }
            });
            map.addLayer(webGLVectorLayer);
    }, [map, imageSize, cellTypeCoordinatesData]);

    return (
        <div style={{ height: '100%', width: '50%' }}>
            <div style={{ position: 'relative', height: '100%', width: '100%' }}>
                <div ref={mapRef} style={{ height: '100%', width: '100%' }} />
            </div>
        </div>
    );
};