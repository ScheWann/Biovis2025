import React, { useEffect, useRef, useState } from 'react';
import { Select, Slider, Button } from 'antd';
import * as d3 from 'd3';

// 根据源端与目标端的坐标及各自厚度，生成一个渐变（tapered）的链接路径。
// 此路径为闭合多边形，从 source 厚度（不变）过渡到 target 厚度（可能较小）。
const taperedLinkPath = (sourceX, sourceY, targetX, targetY, thicknessSource, thicknessTarget) => {
    // 定义四个点：source 顶部、target 顶部、target 底部、source 底部
    const sourceTop = [sourceX, sourceY - thicknessSource / 2];
    const targetTop = [targetX, targetY - thicknessTarget / 2];
    const targetBottom = [targetX, targetY + thicknessTarget / 2];
    const sourceBottom = [sourceX, sourceY + thicknessSource / 2];

    const path = d3.path();
    path.moveTo(...sourceTop);
    // 使用贝塞尔曲线平滑连接
    const midX = (sourceX + targetX) / 2;
    path.bezierCurveTo(midX, sourceY, midX, targetY, ...targetTop);
    path.lineTo(...targetBottom);
    path.bezierCurveTo(midX, targetY, midX, sourceY, ...sourceBottom);
    path.closePath();
    return path.toString();
};

export const NMFGOExpressionViewer = ({ regions, NMFGODataLoading, NMFGOData, setNMFGOData }) => {
    const containerRef = useRef(null);
    const [dimensions, setDimensions] = useState({ width: 600, height: 400 });
    const [regionOptions, setRegionOptions] = useState([]); // 用于 Select 控件的区域选项
    const [selectedRegion, setSelectedRegion] = useState([]); // 选中的 region
    const [selectedComponentCount, setSelectedComponentCount] = useState(null); // 选中的 component 数量
    const [selectedResolution, setSelectedResolution] = useState(null); // 选中的 resolution

    // component count options(2到10)
    const componentCountOptions = [
        { label: "2", value: 2 },
        { label: "3", value: 3 },
        { label: "4", value: 4 },
        { label: "5", value: 5 },
        { label: "6", value: 6 },
        { label: "7", value: 7 },
        { label: "8", value: 8 },
        { label: "9", value: 9 },
        { label: "10", value: 10 },
    ];

    // resolution options(0.25, 0.5, 0.75, 1.0)
    const resolutionOptions = [
        { label: "0.25", value: 0.25 },
        { label: "0.5", value: 0.5 },
        { label: "0.75", value: 0.75 },
        { label: "1.0", value: 1.0 }
    ];

    // 监听容器尺寸变化
    useEffect(() => {
        const container = containerRef.current;
        const resizeObserver = new ResizeObserver(entries => {
            for (let entry of entries) {
                setDimensions({
                    width: entry.contentRect.width,
                    height: entry.contentRect.height,
                });
            }
        });
        if (container) resizeObserver.observe(container);
        return () => resizeObserver.disconnect();
    }, []);

    useEffect(() => {
        regions.forEach(region => {
            setRegionOptions(prevOptions => {
                const existingOption = prevOptions.find(option => option.value === region.id);
                if (!existingOption) {
                    return [...prevOptions, { label: region.name, value: region.id }];
                }
                return prevOptions;
            });
        });
    }, [regions]);

    // 用于计算 component 节点的位置及高度
    const layoutComponentNodes = (components, yOffset, regionHeight) => {
        // 定义区域内可用的布局高度（预留上下各5%边距）
        const maxLayoutHeight = regionHeight * 0.9;
        // 定义最小矩形高度，过多节点时不低于此值
        const minRectHeight = 10;
        const space = maxLayoutHeight / components.length;
        const rectHeight = space < minRectHeight ? minRectHeight : space * 0.6; // 60% 用于节点，其余为间隙
        return components.map((node, index) => {
            const y = yOffset + (index + 0.5) * space;
            return { ...node, x: 350, y, rectHeight };
        });
    };

    // 数据或尺寸变化时重绘 Sankey 图
    useEffect(() => {
        if (!NMFGOData || Object.keys(NMFGOData).length === 0) return;

        d3.select(containerRef.current).selectAll("*").remove();
        const { width, height } = dimensions;
        const svg = d3.select(containerRef.current)
            .append("svg")
            .attr("width", width)
            .attr("height", height);

        // --- 预处理：提取所有 region → cluster 链接中的细胞数量 ---
        const clusterSizeMap = {};
        NMFGOData.links.forEach(link => {
            const sourceParts = link.source.split('_');
            const targetParts = link.target.split('_');
            // 如果 source 节点为 region，target 节点为 cluster，则 link.value 为 cluster 的细胞数量
            const sourceNode = NMFGOData.nodes.find(n => n.id === link.source);
            const targetNode = NMFGOData.nodes.find(n => n.id === link.target);
            if (sourceNode && sourceNode.type === 'region' && targetNode && targetNode.type === 'cluster') {
                if (!clusterSizeMap[link.target]) {
                    clusterSizeMap[link.target] = 0;
                }
                clusterSizeMap[link.target] += link.value;
            }
        });
        const clusterSizes = Object.values(clusterSizeMap);
        const maxSize = d3.max(clusterSizes);
        const minSize = d3.min(clusterSizes);
        const minWidth = 30;
        const maxWidth = 80;
        const clusterWidthScale = d3.scaleLinear()
            .domain([minSize, maxSize])
            .range([minWidth, maxWidth]);

        // --- 按 region 对节点进行分组 ---
        const regionsMap = {};
        NMFGOData.nodes.forEach(node => {
            if (node.type === "region") {
                regionsMap[node.name] = {
                    regionNode: node,
                    clusters: [],
                    components: [],
                };
            }
        });
        // 根据 id 将 cluster 与 component 分配到对应 region（假设 id 格式："cluster_123_1" 或 "comp_123_Component_1"）
        NMFGOData.nodes.forEach(node => {
            if (node.type === "cluster") {
                const parts = node.id.split('_');
                const regionName = parts[1];
                if (regionsMap[regionName]) {
                    regionsMap[regionName].clusters.push(node);
                }
            } else if (node.type === "component") {
                const parts = node.id.split('_');
                const regionName = parts[1];
                if (regionsMap[regionName]) {
                    regionsMap[regionName].components.push(node);
                }
            }
        });

        // --- 对每个 region 内的 cluster 和 component 按标签从小到大排序 ---
        Object.keys(regionsMap).forEach(regionName => {
            const regionGroup = regionsMap[regionName];
            regionGroup.clusters.sort((a, b) => {
                const aLabel = parseInt(a.id.split('_').pop(), 10);
                const bLabel = parseInt(b.id.split('_').pop(), 10);
                return aLabel - bLabel;
            });
            regionGroup.components.sort((a, b) => {
                const aLabel = parseInt(a.id.split('_').pop(), 10);
                const bLabel = parseInt(b.id.split('_').pop(), 10);
                return aLabel - bLabel;
            });
        });

        // --- 根据 region 数量计算区域高度 ---
        const regionNames = Object.keys(regionsMap).sort();
        const regionCount = regionNames.length;
        const regionHeight = height / (regionCount || 1);

        // 固定参数：cluster 节点左对齐 x 坐标；component 节点 x 坐标；默认宽度
        const clusterX = 100;  // cluster 节点左对齐时的 x 坐标（矩形左侧）
        const rectWidth = 60;  // 默认 component 节点宽度
        const clusterRectHeight = 20; // cluster 节点固定高度

        // --- 添加全局列头标签，这里改为在底部显示 ---
        svg.append("text")
            .attr("x", clusterX)
            .attr("y", height)
            .attr("text-anchor", "start")
            .style("font-weight", "bold")
            .style("font-size", "14px")
            .text("Cluster");
        svg.append("text")
            .attr("x", 350)
            .attr("y", height)
            .attr("text-anchor", "middle")
            .style("font-size", "14px")
            .style("font-weight", "bold")
            .text("Component");

        // --- 遍历每个 region 绘制 ---
        regionNames.forEach((regionName, i) => {
            const regionGroup = regionsMap[regionName];
            const yOffset = i * regionHeight;

            // 绘制区域名称
            svg.append("text")
                .attr("x", 10)
                .attr("y", yOffset + regionHeight / 2)
                .attr("dy", ".35em")
                .text(regionName)
                .style("font-weight", "bold")
                .style("font-size", "14px");

            // 绘制 cluster 节点（左对齐绘制），标签取最后一个下划线后的数字
            const clusters = regionGroup.clusters;
            clusters.forEach((node, index) => {
                const clusterSize = clusterSizeMap[node.id] || 0;
                const dynamicWidth = clusterWidthScale(clusterSize);
                node.clusterWidth = dynamicWidth; // 保存宽度
                const y = yOffset + (index + 1) * regionHeight / (clusters.length + 1);
                node.x = clusterX; // 左对齐：x 为矩形左侧
                node.y = y;
                svg.append("rect")
                    .attr("x", node.x)
                    .attr("y", node.y - clusterRectHeight / 2)
                    .attr("width", dynamicWidth)
                    .attr("height", clusterRectHeight)
                    .style("fill", "#69b3a2");
                // 提取最后一个下划线后的数字作为标签
                const label = node.id.split('_').pop();
                svg.append("text")
                    .attr("x", node.x + 4)
                    .attr("y", node.y)
                    .attr("dy", ".35em")
                    .attr("text-anchor", "start")
                    .text(label);
            });

            // 绘制 component 节点 —— 调用 layoutComponentNodes 动态布局；标签取最后一个下划线后的数字
            const components = regionGroup.components;
            const laidOutComponents = layoutComponentNodes(components, yOffset, regionHeight);
            laidOutComponents.forEach(node => {
                svg.append("rect")
                    .attr("x", node.x - rectWidth / 2)
                    .attr("y", node.y - node.rectHeight / 2)
                    .attr("width", rectWidth)
                    .attr("height", node.rectHeight)
                    .style("fill", "#4C9AFF");
                const label = node.id.split('_').pop();
                svg.append("text")
                    .attr("x", node.x)
                    .attr("y", node.y)
                    .attr("dy", ".35em")
                    .attr("text-anchor", "middle")
                    .text(label);
            });

            // --- 绘制链接部分 ---
            // 仅绘制同一区域内 cluster 到 component 的链接
            const linksBySource = {};
            NMFGOData.links.forEach(link => {
                const srcParts = link.source.split('_');
                const tgtParts = link.target.split('_');
                if (srcParts[1] === regionName && tgtParts[1] === regionName) {
                    // 只处理 cluster 到 component 的链接
                    const sourceNode = NMFGOData.nodes.find(n => n.id === link.source);
                    const targetNode = NMFGOData.nodes.find(n => n.id === link.target);
                    if (sourceNode && sourceNode.type === 'cluster' &&
                        targetNode && targetNode.type === 'component') {
                        if (!linksBySource[link.source]) {
                            linksBySource[link.source] = [];
                        }
                        linksBySource[link.source].push(link);
                    }
                }
            });

            clusters.forEach(source => {
                const outgoingLinks = linksBySource[source.id];
                if (!outgoingLinks) return;
                const totalValue = outgoingLinks.reduce((acc, link) => acc + link.value, 0);
                let offset = 0;
                outgoingLinks.forEach(link => {
                    const thicknessOriginal = totalValue !== 0
                        ? (link.value / totalValue) * clusterRectHeight
                        : link.value * clusterRectHeight;
                    const targetCandidate = laidOutComponents.find(n => n.id === link.target);
                    const thicknessTarget = targetCandidate
                        ? (thicknessOriginal > targetCandidate.rectHeight ? targetCandidate.rectHeight : thicknessOriginal)
                        : thicknessOriginal;
                    const startY = source.y - clusterRectHeight / 2 + offset + thicknessOriginal / 2;
                    const targetComponent = laidOutComponents.find(n => n.id === link.target);
                    if (targetComponent) {
                        const targetX = targetComponent.x - rectWidth / 2;
                        const targetY = targetComponent.y;
                        offset += thicknessOriginal;
                        // cluster 左对齐，链接起点 x 为 source.x + source.clusterWidth（即矩形右侧）
                        const pathD = taperedLinkPath(
                            source.x + (source.clusterWidth || rectWidth), startY,
                            targetX, targetY,
                            thicknessOriginal, thicknessTarget
                        );
                        svg.append("path")
                            .attr("d", pathD)
                            .attr("fill", "#999")
                            .attr("opacity", 0.5);
                    }
                });
            });
        });
    }, [dimensions, NMFGOData]);

    const confirm = () => {
        const params = {};
        selectedRegion.forEach(regionValue => {
            const regionData = regions.find(r => r.id === regionValue);
            if (regionData) {
                params[regionData.name] = {
                    sampleid: regionData.sampleId,
                    cell_list: regionData.cellIds,
                };
            }
        });
        fetch('/get_NMF_GO_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ regions: params, n_component: selectedComponentCount, resolution: selectedResolution })
        })
            .then(res => res.json())
            .then(data => {
                setNMFGOData(data.sankey);
            });
    };

    return (
        <div style={{ width: '100%', height: '100%' }}>
            {/* 上方的 Select 控件 */}
            <div style={{ marginBottom: 16, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                <Select
                    size='small'
                    mode="multiple"
                    placeholder="Select Region"
                    options={regionOptions}
                    style={{ width: 200, marginRight: 8 }}
                    onChange={value => setSelectedRegion(value)}
                    allowClear
                />
                <Select
                    size='small'
                    placeholder="Select Component Count"
                    options={componentCountOptions}
                    style={{ width: 120, marginRight: 8 }}
                    onChange={value => setSelectedComponentCount(value)}
                    allowClear
                />
                <Select
                    size='small'
                    placeholder="Select Resolution"
                    options={resolutionOptions}
                    style={{ width: 120 }}
                    onChange={value => setSelectedResolution(value)}
                    allowClear
                />
                <Button size='small' onClick={confirm} style={{ marginLeft: 8 }}>Confirm</Button>
            </div>
            <div ref={containerRef} style={{ width: "100%", height: "calc(100% - 50px)" }} />
        </div>
    );
};

