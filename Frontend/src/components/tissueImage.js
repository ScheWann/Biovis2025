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
import hireImage from '../data/tissue_hires_image.png';


const modeOptions = [
    {
        value: 'kmeans',
        label: 'Kmeans',
    },
    {
        value: 'genes',
        label: 'Genes',
    },
]

const kmeansOptions = Array.from({ length: 9 }, (_, i) => ({
    value: i + 2,
    label: `${i + 2}`,
}));


export const TissueImage = ({ mode, setMode, geneName, setGeneName, binSize, kmeansSize, setKmeansSize, tissueData, setTissueData }) => {
    const mapRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [lassoToggleStatus, setLassoToggleStatus] = useState(false);
    const [selectedRegion, setSelectedRegion] = useState(null);
    const [view, setView] = useState(null);
    const [map, setMap] = useState(null);
    const [secondOptions, setSecondOptions] = useState(kmeansOptions);

    const fetchSpecificGeneData = (geneName) => {
        fetch('/get_specific_gene_expression', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ gene_name: geneName, bin_size: binSize })
        })
            .then((response) => response.json())
            .then((data) => {
                setTissueData(data);
            })
    }

    const fetchfullGeneList = () => {
        fetch('/get_full_gene_list', {
            method: 'GET',
        })
            .then((response) => response.json())
            .then((data) => {
                const geneOptions = data.map(item => ({
                    label: item,
                    value: item,
                }));

                setSecondOptions(geneOptions);
            })
    }

    const fetchGeneNameBySearch = async (query = '') => {
        const response = await fetch(`/get_gene_name_search?q=${encodeURIComponent(query || '')}`);
        const data = await response.json();

        const geneOptions = data.map(item => ({
            label: item,
            value: item,
        }));

        setSecondOptions(geneOptions);
    };

    const modeChange = (value) => {
        setMode(value);
        if (value === 'kmeans') {
            setSecondOptions(kmeansOptions);
        } else if (value === 'genes') {
            fetchfullGeneList();
        }
    };

    const handleChange = (value) => {
        if (mode === 'kmeans') {
            setKmeansSize(value);
        } else {
            if (value.length > 0) {
                setGeneName(value);
                fetchSpecificGeneData(value);
            }
        }
    };

    const lassoChange = () => {
        setLassoToggleStatus(!lassoToggleStatus);
        setSelectedRegion(null);
    };

    // get hire image size(width, height)
    useEffect(() => {
        fetch("/get_hires_image_size", {
            method: "GET",
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

    // lasso selection
    useEffect(() => {
        if (map && lassoToggleStatus) {
            const vectorSource = new VectorSource();
            const vectorLayer = new VectorLayer({
                source: vectorSource,
            });

            map.addLayer(vectorLayer);

            const drawInteraction = new Draw({
                source: vectorSource,
                type: 'Polygon',
            });

            map.addInteraction(drawInteraction);

            drawInteraction.on('drawend', (event) => {
                const geometry = event.feature.getGeometry();
                console.log(geometry);
                setSelectedRegion(geometry);

                // The drawn region
                const regionStyle = new Style({
                    fill: new Fill({
                        color: 'rgba(0, 191, 255, 0.3)',
                    }),
                    stroke: new Stroke({
                        color: 'rgba(0, 191, 255, 1)',
                        width: 2,
                    }),
                });

                event.feature.setStyle(regionStyle);
            });

            return () => {
                map.removeLayer(vectorLayer);
                map.removeInteraction(drawInteraction);
            };
        }
    }, [map, lassoToggleStatus]);

    useEffect(() => {
        if (!map || imageSize.length <= 0) return;
        if (mode === 'kmeans') {
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

            tissueData.forEach((point) => {
                const { x, y, cluster } = point;

                const adjustedX = x;
                const adjustedY = imageSize[1] - y;

                const feature = new Feature({
                    geometry: new Point([adjustedX, adjustedY]),
                });

                feature.set('cluster', clusterColors[cluster]);

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
        }
        if (mode === 'genes' && geneName.length > 0) {
            const maxGeneValue = tissueData[0]['max'];
            const minGeneValue = tissueData[0]['min'];
            
            const colorScale = d3.scaleSequential(d3.interpolateOranges).domain([minGeneValue, maxGeneValue]);
            const vectorSource = new VectorSource();

            tissueData.forEach((point) => {
                const { x, y } = point;
                const geneValue = point[geneName];
                const color = colorScale(geneValue)

                const adjustedX = x;
                const adjustedY = imageSize[1] - y;

                const feature = new Feature({
                    geometry: new Point([adjustedX, adjustedY]),
                });

                feature.set('color', color);
                vectorSource.addFeature(feature);
            });

            const webGLVectorLayer = new WebGLVectorLayer({
                source: vectorSource,
                style: {
                    'circle-radius': 1,
                    'circle-fill-color': ['get', 'color'],
                },
            });

            map.getLayers().forEach((layer) => {
                if (layer instanceof WebGLVectorLayer) {
                    map.removeLayer(layer);
                }
            });
            map.addLayer(webGLVectorLayer);
        }
    }, [map, imageSize, mode, tissueData]);

    const resetZoom = () => {
        if (view) {
            view.setCenter([imageSize[0] / 2, imageSize[1] / 2]);
            view.setZoom(1);
        }
    };

    return (
        <div style={{ height: '100%', width: '50%' }}>
            <div style={{ position: 'relative', height: '100%', width: '100%' }}>
                <div ref={mapRef} style={{ height: '100%', width: '100%' }} />
                <div className="controlButtonGroup">
                    <Select
                        value={mode}
                        style={{
                            width: 120,
                        }}
                        onChange={modeChange}
                        options={modeOptions}
                    />
                    <Select
                        mode={mode === 'genes' ? 'tags' : 'default'}
                        showSearch={mode === 'genes'}
                        value={mode === 'genes' ? geneName : kmeansSize}
                        style={{
                            width: 120,
                        }}
                        onChange={handleChange}
                        options={secondOptions}
                        onSearch={mode === 'genes' ? fetchGeneNameBySearch : false}
                    />
                    <Button
                        onClick={lassoChange}
                        icon={<SelectOutlined />}
                    >
                    </Button>
                    <Button
                        onClick={resetZoom}
                        icon={<RollbackOutlined />}
                    >
                    </Button>
                </div>
            </div>
        </div>
    );
};