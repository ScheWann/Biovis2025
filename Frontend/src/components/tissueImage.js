import React, { useEffect, useRef, useState } from 'react';
import { RollbackOutlined, SelectOutlined } from "@ant-design/icons";
import { Button, Select } from "antd";
import { Map, View } from 'ol';
import 'ol/ol.css';
import "../styles/tissueImage.css";
import ImageLayer from 'ol/layer/Image';
import ImageStatic from 'ol/source/ImageStatic';
import VectorImageLayer from 'ol/layer/VectorImage';
import VectorSource from 'ol/source/Vector';
import Feature from 'ol/Feature';
import Point from 'ol/geom/Point';
import { Icon, Style } from 'ol/style';
import CircleStyle from 'ol/style/Circle';
import Fill from 'ol/style/Fill';
import Draw from 'ol/interaction/Draw';
import hireImage from '../data/tissue_hires_image.png';
import scaleFactor from '../data/scalefactors_json.json';

export const TissueImage = ({ positionWithClusterData, kmeansSize, setKmeansSize }) => {
    const mapRef = useRef(null);
    const [imageSize, setImageSize] = useState([]);
    const [lassoToggleStatus, setLassoToggleStatus] = useState(false);
    const [selectedRegion, setSelectedRegion] = useState(null);
    const [view, setView] = useState(null);
    const [map, setMap] = useState(null);

    // Kmeans number options
    const kmeansOptions = Array.from({ length: 9 }, (_, i) => ({
        value: i + 2,
        label: `${i + 2}`,
    }));

    const handleChange = (value) => {
        setKmeansSize(value);
    };

    const lassoChange = () => {
        setLassoToggleStatus(!lassoToggleStatus);
        setSelectedRegion(null);
    };

    useEffect(() => {
        fetch("/get_hires_image_size", {
            method: "GET",
        })
            .then((response) => response.json())
            .then((data) => {
                setImageSize(data);
            });
    }, []);

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
        if (map && lassoToggleStatus) {
            // 创建一个 Draw 交互工具，用于绘制多边形
            const drawInteraction = new Draw({
                source: new VectorSource(), // 用于存储绘制的区域
                type: 'Polygon', // 指定几何类型为多边形
            });

            // 添加到地图中
            map.addInteraction(drawInteraction);

            // 监听绘制完成事件
            drawInteraction.on('drawend', (event) => {
                const geometry = event.feature.getGeometry();
                console.log(geometry);
                setSelectedRegion(geometry); // 保存选中的区域
            });

            // 清理交互工具
            return () => {
                map.removeInteraction(drawInteraction);
            };
        }
    }, [map, lassoToggleStatus]);

    useEffect(() => {
        if (map && imageSize.length > 0) {
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

            positionWithClusterData.forEach((point) => {
                const { x, y, cluster } = point;

                const adjustedX = x * scaleFactor.tissue_hires_scalef;
                const adjustedY = imageSize[1] - y * scaleFactor.tissue_hires_scalef;

                const feature = new Feature({
                    geometry: new Point([adjustedX, adjustedY]),
                });

                const color = clusterColors[cluster] || 'rgba(255, 255, 255, 0.7)';
                feature.setStyle(
                    new Style({
                        image: new CircleStyle({
                            radius: 1,
                            fill: new Fill({ color }),
                        }),
                    })
                );

                vectorSource.addFeature(feature);
            });

            const vectorImageLayer = new VectorImageLayer({
                source: vectorSource,
            });

            // clear old vector layers
            map.getLayers().forEach((layer) => {
                if (layer instanceof VectorImageLayer) {
                    map.removeLayer(layer);
                }
            });
            map.addLayer(vectorImageLayer);
        }
    }, [map, positionWithClusterData, imageSize]);

    const resetZoom = () => {
        if (view) {
            view.setCenter([imageSize[0] / 2, imageSize[1] / 2]);
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
    );
};