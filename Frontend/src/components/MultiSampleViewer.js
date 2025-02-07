import React, { useRef, useState, useEffect, useMemo, useCallback } from 'react';
import DeckGL from '@deck.gl/react';
import { Collapse, Button, Input, ColorPicker, Checkbox, message, Spin } from "antd";
import { CloseOutlined } from '@ant-design/icons';
import { OrthographicView } from '@deck.gl/core';
import { BitmapLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import { booleanPointInPolygon, centroid } from '@turf/turf';
import { TileLayer } from '@deck.gl/geo-layers';
import { EditableGeoJsonLayer, DrawPolygonMode } from '@deck.gl-community/editable-layers';
import { fromBlob } from 'geotiff';

// 工具函数：HSL 转 RGB
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

// 生成随机颜色，每次保存区域后调用以更新下一次绘制区域使用的颜色
const generateRandomColor = () => {
    return hslToRgb(Math.floor(Math.random() * 360), 100, 50);
};

export const MultiSampleViewer = ({
    samples, // 样本数组
    cellTypeCoordinatesData,
    cellTypeDir,
    regions,
    setRegions
}) => {
    const viewerRef = useRef(null);
    const [imageSizes, setImageSizes] = useState({});
    const [tileSize] = useState(256);
    const [features, setFeatures] = useState({});
    const [tempRegions, setTempRegions] = useState({});
    const [selectionMode, setSelectionMode] = useState(null);
    const [visibleCellTypes, setVisibleCellTypes] = useState({});
    const [colorMaps, setColorMaps] = useState({});
    const [hoveredCell, setHoveredCell] = useState(null);
    const [activeSample, setActiveSample] = useState(samples[0]?.id);
    const [visibleSamples, setVisibleSamples] = useState(
        samples.reduce((acc, sample) => ({ ...acc, [sample.id]: true }), {})
    );
    const [regionName, setRegionName] = useState('');
    // 默认区域颜色使用随机颜色；保存区域后会更新为新的随机颜色
    const [regionColor, setRegionColor] = useState(generateRandomColor());

    // 初始化每个样本的颜色映射和细胞类型可见性
    useEffect(() => {
        const initialStates = samples.reduce((acc, sample) => {
            const initialColorMap = {};
            const initialVisibility = {};
            cellTypeDir?.forEach((cellType, index) => {
                const hue = (index * 360) / cellTypeDir.length;
                initialColorMap[cellType] = hslToRgb(hue, 100, 50);
                initialVisibility[cellType] = true;
            });
            return {
                colorMaps: { ...acc.colorMaps, [sample.id]: initialColorMap },
                visibleCellTypes: { ...acc.visibleCellTypes, [sample.id]: initialVisibility }
            };
        }, { colorMaps: {}, visibleCellTypes: {} });
        setColorMaps(initialStates.colorMaps);
        setVisibleCellTypes(initialStates.visibleCellTypes);
    }, [samples, cellTypeDir]);

    // 获取所有样本的图片尺寸
    useEffect(() => {
        const fetchImageSizes = async () => {
            const sizes = {};
            await Promise.all(samples.map(async (sample) => {
                const response = await fetch("/get_hires_image_size", {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ sample_id: sample.id })
                });
                const data = await response.json();
                sizes[sample.id] = data;
            }));
            setImageSizes(sizes);
        };
        if (samples.length) fetchImageSizes();
    }, [samples]);

    // 生成样本瓦片图层
    const generateTileLayers = useCallback(() => {
        return samples.filter(sample => visibleSamples[sample.id]).map(sample => {
            const imageSize = imageSizes[sample.id];
            if (!imageSize || imageSize.length < 2) return null;
            return new TileLayer({
                id: `tif-tiles-${sample.id}`,
                tileSize,
                extent: [0, 0, imageSize[0], imageSize[1]],
                maxZoom: 0,
                minZoom: 0,
                getTileData: async ({ index }) => {
                    try {
                        const { x, y } = index;
                        const response = await fetch(`/get_tile?sample_id=${sample.id}&x=${x}&y=${y}`);
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
                        return {
                            image: canvas,
                            bounds: [
                                [x * tileSize, (y + 1) * tileSize],
                                [x * tileSize, y * tileSize],
                                [(x + 1) * tileSize, y * tileSize],
                                [(x + 1) * tileSize, (y + 1) * tileSize]
                            ]
                        };
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
        }).filter(Boolean);
    }, [samples, imageSizes, tileSize, visibleSamples]);

    // 生成 marker 图层（低透明度，用于背景参考）
    const generateMarkerImageLayers = useCallback(() => {
        return samples.filter(sample => visibleSamples[sample.id]).map(sample => {
            const imageSize = imageSizes[sample.id];
            if (!imageSize || imageSize.length < 2) return null;
            return new BitmapLayer({
                id: `marker-image-${sample.id}`,
                image: `/${sample.id}_cells_layer.png`, // 根据实际路径调整
                bounds: [0, imageSize[1], imageSize[0], 0],
                opacity: 0.05,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [samples, imageSizes, visibleSamples]);

    // 生成细胞图层
    const generateCellLayers = useCallback(() => {
        return samples.filter(sample => visibleSamples[sample.id]).map(sample => {
            const visibleData = cellTypeCoordinatesData?.[sample.id]?.filter(cell =>
                visibleCellTypes[sample.id]?.[cell.cell_type] ?? true
            );
            return visibleData?.length > 0 && new ScatterplotLayer({
                id: `cells-${sample.id}`,
                data: visibleData,
                getPosition: d => [d.cell_x, d.cell_y],
                getRadius: 2,
                getFillColor: d => colorMaps[sample.id]?.[d.cell_type] || [0, 0, 0],
                pickable: true,
                parameters: { depthTest: false }
            });
        }).filter(Boolean);
    }, [samples, cellTypeCoordinatesData, visibleCellTypes, colorMaps, visibleSamples]);

    // 区域编辑逻辑：更新当前绘制区域数据
    const handleRegionUpdate = (sampleId, updatedData) => {
        setFeatures(prev => ({ ...prev, [sampleId]: updatedData }));
        const lastFeature = updatedData.features?.[updatedData.features.length - 1];
        if (lastFeature) {
            setTempRegions(prev => ({
                ...prev,
                [sampleId]: {
                    ...prev[sampleId],
                    data: cellTypeCoordinatesData?.[sampleId]?.filter(cell =>
                        booleanPointInPolygon([cell.cell_x, cell.cell_y], lastFeature.geometry)
                    )
                }
            }));
        }
    };

    // 保存区域，同时保存区域名称和区域颜色，并更新区域颜色为新的随机颜色
    const handleSaveRegion = () => {
        if (!regionName) {
            message.error('请填写区域名称');
            return;
        }
        const sampleId = activeSample;
        const feature = features[sampleId]?.features?.[features[sampleId].features.length - 1];
        if (!feature) {
            message.error('没有可保存的区域');
            return;
        }
        const newRegion = {
            id: `${sampleId}-${regionName}-${Date.now()}`,
            name: regionName,
            color: regionColor,
            feature: {
                type: 'FeatureCollection',
                features: [feature]
            }
        };
        setRegions(prev => [...prev, newRegion]);
        setFeatures(prev => ({ ...prev, [sampleId]: { type: 'FeatureCollection', features: [] } }));
        setTempRegions(prev => ({ ...prev, [sampleId]: null }));
        message.success('区域保存成功');
        setRegionName('');
        setSelectionMode(null);
        // 保存后更新区域颜色为新的随机颜色，保证下次绘制时颜色不同
        setRegionColor(generateRandomColor());
    };

    // 删除区域
    const handleDeleteRegion = (regionId) => {
        setRegions(prev => prev.filter(region => region.id !== regionId));
        message.success('区域删除成功');
    };

    // 生成编辑图层（用于绘制新区域）
    const generateEditLayers = useCallback(() => {
        return samples.map(sample => {
            const tempRegion = tempRegions[sample.id] || {};
            return new EditableGeoJsonLayer({
                id: `edit-layer-${sample.id}`,
                data: features[sample.id] || { type: 'FeatureCollection', features: [] },
                mode: DrawPolygonMode,
                selectedFeatureIndexes: [],
                onEdit: ({ updatedData }) => handleRegionUpdate(sample.id, updatedData),
                visible: selectionMode === sample.id, // 只有当前激活样本的编辑图层可见
                getLineColor: tempRegion.color ? [...tempRegion.color, 200] : [255, 0, 0, 200],
                getFillColor: tempRegion.color ? [...tempRegion.color, 50] : [255, 255, 0, 50],
                lineWidthMinPixels: 2
            });
        });
    }, [samples, features, tempRegions, selectionMode]);

    // 辅助函数：计算区域内的细胞数量
    const getCellCount = (region) => {
        const sampleId = region.id.split('-')[0];
        if (!cellTypeCoordinatesData[sampleId] || region.feature.features.length === 0) {
            return 0;
        }
        const polygon = region.feature.features[0].geometry;
        return cellTypeCoordinatesData[sampleId].filter(cell =>
            booleanPointInPolygon([cell.cell_x, cell.cell_y], polygon)
        ).length;
    };

    const allDataLoaded = useMemo(() =>
        samples.every(sample =>
            imageSizes[sample.id] &&
            cellTypeCoordinatesData[sample.id]
        ),
        [samples, imageSizes, cellTypeCoordinatesData]
    );

    // 合并所有图层，同时增加一个 TextLayer 用于显示已保存区域的名称（在区域中心位置）
    const layers = useMemo(() => {
        const regionLabelLayer = new TextLayer({
            id: 'region-labels',
            data: regions
                .map(region => {
                    if (region.feature.features.length === 0) return null;
                    const center = centroid(region.feature.features[0]);
                    return {
                        position: center.geometry.coordinates,
                        text: region.name
                    };
                })
                .filter(Boolean),
            getPosition: d => d.position,
            getText: d => d.text,
            getColor: () => [0, 0, 0, 255],
            getSize: 16,
            sizeUnits: 'pixels'
        });
        return [
            ...generateTileLayers(),
            ...generateMarkerImageLayers(),
            ...generateCellLayers(),
            ...generateEditLayers(),
            ...regions.map(region =>
                new EditableGeoJsonLayer({
                    id: `region-${region.id}`,
                    data: region.feature,
                    getLineColor: [...region.color, 200],
                    getFillColor: [...region.color, 50]
                })
            ),
            regionLabelLayer
        ].filter(Boolean);
    }, [generateTileLayers, generateCellLayers, generateMarkerImageLayers, generateEditLayers, regions]);

    // 全局绘制模式控制：点击右上角按钮时切换绘制状态
    const toggleDrawingMode = () => {
        if (!selectionMode) {
            setSelectionMode(activeSample);
        } else {
            handleSaveRegion();
        }
    };

    if (!allDataLoaded) {
        return (
            <div style={{
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                height: '100vh'
            }}>
                <Spin tip="加载数据中..." size="large" />
            </div>
        );
    }

    return (
        <div style={{ height: '100vh', display: 'flex' }}>
            {/* 左侧面板：仅显示每个样本的“显示/隐藏样本”和“细胞类型设置” */}
            <div className="controls" style={{ position: 'absolute', top: 10, left: 10, zIndex: 10 }}>
                <Collapse>
                    {samples.map(sample => (
                        <Collapse.Panel key={sample.id} header={sample.name}>
                            <Button
                                block
                                onClick={() => setVisibleSamples(prev => ({
                                    ...prev,
                                    [sample.id]: !prev[sample.id]
                                }))}
                            >
                                {visibleSamples[sample.id] ? '隐藏样本' : '显示样本'}
                            </Button>
                            <Collapse style={{ marginTop: 16 }}>
                                <Collapse.Panel key={`cell-types-${sample.id}`} header="细胞类型设置">
                                    <CellTypeSettings
                                        sampleId={sample.id}
                                        cellTypes={cellTypeDir}
                                        colorMap={colorMaps[sample.id] || {}}
                                        visibleMap={visibleCellTypes[sample.id] || {}}
                                        onColorChange={(type, color) => {
                                            setColorMaps(prev => ({
                                                ...prev,
                                                [sample.id]: { ...prev[sample.id], [type]: color }
                                            }));
                                        }}
                                        onVisibilityChange={(type, visible) => {
                                            setVisibleCellTypes(prev => ({
                                                ...prev,
                                                [sample.id]: { ...prev[sample.id], [type]: visible }
                                            }));
                                        }}
                                    />
                                </Collapse.Panel>
                            </Collapse>
                        </Collapse.Panel>
                    ))}
                </Collapse>
            </div>

            {/* 主视图区域 */}
            <div style={{ flex: 1, position: 'relative' }}>
                <DeckGL
                    layers={layers}
                    views={new OrthographicView({ controller: true })}
                    initialViewState={{
                        target: [5000, 5000, 0],
                        zoom: -2,
                        maxZoom: 1,
                        minZoom: -5
                    }}
                    onHover={info => {
                        if (info.object && info.layer.id.startsWith('cells-')) {
                            const sampleId = info.layer.id.split('-')[1];
                            setHoveredCell({
                                ...info.object,
                                sampleId,
                                x: info.x,
                                y: info.y
                            });
                        } else {
                            setHoveredCell(null);
                        }
                    }}
                    controller={true}
                />

                {/* 悬浮提示 */}
                {hoveredCell && (
                    <div style={{
                        position: 'absolute',
                        left: hoveredCell.x,
                        top: hoveredCell.y - 40,
                        backgroundColor: 'rgba(255, 255, 255, 0.9)',
                        padding: 8,
                        borderRadius: 4,
                        boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
                        transform: 'translateX(-50%)'
                    }}>
                        <div>样本: {hoveredCell.sampleId}</div>
                        <div>细胞类型: {hoveredCell.cell_type}</div>
                    </div>
                )}

                {/* 右上角控制面板：绘制区域控制和全局已保存区域 */}
                <div style={{
                    position: 'absolute',
                    top: 10,
                    right: 10,
                    zIndex: 20,
                    background: 'rgba(255,255,255,0.8)',
                    padding: '10px',
                    borderRadius: '4px',
                    width: 220
                }}>
                    {/* 绘制区域控制部分 */}
                    <div>
                        <Input
                            placeholder="区域名称"
                            value={regionName}
                            onChange={e => setRegionName(e.target.value)}
                            style={{ marginBottom: 8 }}
                        />
                        <ColorPicker
                            value={`rgb(${regionColor.join(',')})`}
                            onChange={color => {
                                const rgb = color.toRgb();
                                setRegionColor([rgb.r, rgb.g, rgb.b]);
                            }}
                            style={{ marginBottom: 8 }}
                        />
                        <Button type="primary" onClick={toggleDrawingMode} block>
                            {selectionMode ? '保存区域' : '开始绘制区域'}
                        </Button>
                    </div>
                    {/* 使用单个 Collapse.Panel 显示“Selected Region (N)” */}
                    <Collapse style={{ marginTop: 10, background: 'rgba(255,255,255,0.9)', borderRadius: 4 }}>
                        <Collapse.Panel header={`Selected Region (${regions.length})`} key="selected-region">
                            {regions.length > 0 ? (
                                regions.map(region => {
                                    const sampleId = region.id.split('-')[0];
                                    const cellCount = getCellCount(region);
                                    return (
                                        <div key={region.id} style={{
                                            display: 'flex',
                                            justifyContent: 'space-between',
                                            alignItems: 'center',
                                            marginBottom: 4,
                                            fontSize: '12px',
                                            color: '#555'
                                        }}>
                                            <div>
                                                ID: {region.id} | Sample: {sampleId} | Cells: {cellCount}
                                            </div>
                                            <CloseOutlined
                                                onClick={(e) => { e.stopPropagation(); handleDeleteRegion(region.id); }}
                                                style={{ fontSize: '12px' }}
                                            />
                                        </div>
                                    );
                                })
                            ) : (
                                <div style={{ fontSize: '12px', color: '#999' }}>No regions selected</div>
                            )}
                        </Collapse.Panel>
                    </Collapse>
                </div>
            </div>
        </div>
    );
};

// 细胞类型设置子组件
const CellTypeSettings = ({
    sampleId,
    cellTypes,
    colorMap,
    visibleMap,
    onColorChange,
    onVisibilityChange
}) => {
    const [searchText, setSearchText] = useState('');
    const filteredTypes = useMemo(() =>
        cellTypes?.filter(type =>
            type.toLowerCase().includes(searchText.toLowerCase())
        ) || []
        , [cellTypes, searchText]);
    return (
        <div style={{ maxHeight: 400, overflowY: 'auto' }}>
            <Input.Search
                placeholder="搜索细胞类型"
                value={searchText}
                onChange={e => setSearchText(e.target.value)}
                style={{ marginBottom: 12 }}
            />
            {filteredTypes.map(type => (
                <div key={type} style={{
                    display: 'flex',
                    alignItems: 'center',
                    padding: '8px 0',
                    borderBottom: '1px solid #f0f0f0'
                }}>
                    <Checkbox
                        checked={visibleMap[type] ?? true}
                        onChange={e => onVisibilityChange(type, e.target.checked)}
                        style={{ marginRight: 8 }}
                    />
                    <ColorPicker
                        size="small"
                        value={`rgb(${(colorMap[type] || [0, 0, 0]).join(',')})`}
                        onChange={color => {
                            const rgb = color.toRgb();
                            onColorChange(type, [rgb.r, rgb.g, rgb.b]);
                        }}
                        style={{ marginRight: 12 }}
                    />
                    <span style={{ fontSize: 14 }}>{type}</span>
                </div>
            ))}
        </div>
    );
};
