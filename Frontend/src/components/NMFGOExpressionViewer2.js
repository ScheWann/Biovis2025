import React, { useEffect, useRef, useState } from 'react';
import { Select, Slider, Button } from 'antd';
import * as d3 from 'd3';

const taperedLinkPath = (sourceX, sourceY, targetX, targetY, thicknessSource, thicknessTarget) => {
    const sourceTop = [sourceX, sourceY - thicknessSource / 2];
    const targetTop = [targetX, targetY - thicknessTarget / 2];
    const targetBottom = [targetX, targetY + thicknessTarget / 2];
    const sourceBottom = [sourceX, sourceY + thicknessSource / 2];

    const path = d3.path();
    path.moveTo(...sourceTop);
    const midX = (sourceX + targetX) / 2;
    path.bezierCurveTo(midX, sourceY, midX, targetY, ...targetTop);
    path.lineTo(...targetBottom);
    path.bezierCurveTo(midX, targetY, midX, sourceY, ...sourceBottom);
    path.closePath();
    return path.toString();
};

export const NMFGOExpressionViewer = ({ regions, setNMFclusterCells, NMFGOData, setNMFGOData }) => {
    const containerRef = useRef(null);
    const [dimensions, setDimensions] = useState({ width: 600, height: 400 });
    const [regionOptions, setRegionOptions] = useState([]);
    const [selectedRegion, setSelectedRegion] = useState([]);
    const [selectedComponentCount, setSelectedComponentCount] = useState(null);
    const [selectedResolution, setSelectedResolution] = useState(null);

    const componentCountOptions = Array.from({ length: 9 }, (_, i) => ({ label: `${i + 2}`, value: i + 2 }));
    const resolutionOptions = [0.25, 0.5, 0.75, 1.0].map(v => ({ label: `${v}`, value: v }));

    useEffect(() => {
        const container = containerRef.current;
        const resizeObserver = new ResizeObserver(entries => {
            for (let entry of entries) {
                setDimensions({ width: entry.contentRect.width, height: entry.contentRect.height });
            }
        });
        if (container) resizeObserver.observe(container);
        return () => resizeObserver.disconnect();
    }, []);

    useEffect(() => {
        regions.forEach(region => {
            setRegionOptions(prevOptions => {
                const exists = prevOptions.find(option => option.value === region.id);
                return exists ? prevOptions : [...prevOptions, { label: region.name, value: region.id }];
            });
        });
    }, [regions]);

    const layoutComponentNodes = (components, yOffset, regionHeight, compRightX) => {
        const maxLayoutHeight = regionHeight * 0.9;
        const minRectHeight = 10;
        const space = maxLayoutHeight / components.length;
        const rectHeight = space < minRectHeight ? minRectHeight : space * 0.6;
        return components.map((node, index) => {
            const y = yOffset + (index + 0.5) * space;
            return { ...node, x: compRightX, y, rectHeight };
        });
    };

    useEffect(() => {
        if (!NMFGOData || Object.keys(NMFGOData).length === 0) return;
        const sankeyData = NMFGOData.sankey;
        const GOData = NMFGOData.GO_results;

        d3.select(containerRef.current).selectAll("*").remove();
        const { width, height } = dimensions;
        const svg = d3.select(containerRef.current).append("svg").attr("width", width).attr("height", height);

        // define the x coordinate of the right side of the component area
        const componentRightX = width - 20;

        const clusterSizeMap = {};
        sankeyData.links.forEach(link => {
            const sourceNode = sankeyData.nodes.find(n => n.id === link.source);
            const targetNode = sankeyData.nodes.find(n => n.id === link.target);
            if (sourceNode?.type === 'region' && targetNode?.type === 'cluster') {
                clusterSizeMap[link.target] = (clusterSizeMap[link.target] || 0) + link.value;
            }
        });

        const clusterSizes = Object.values(clusterSizeMap);
        const clusterWidthScale = d3.scaleLinear()
            .domain([d3.min(clusterSizes), d3.max(clusterSizes)])
            .range([30, 80]);

        const regionsMap = {};
        sankeyData.nodes.forEach(node => {
            if (node.type === "region") {
                regionsMap[node.name] = { regionNode: node, clusters: [], components: [] };
            }
        });
        sankeyData.nodes.forEach(node => {
            const parts = node.id.split('_');
            const regionName = parts[1];
            if (node.type === "cluster" && regionsMap[regionName]) {
                regionsMap[regionName].clusters.push(node);
            } else if (node.type === "component" && regionsMap[regionName]) {
                regionsMap[regionName].components.push(node);
            }
        });

        Object.keys(regionsMap).forEach(regionName => {
            const group = regionsMap[regionName];
            group.clusters.sort((a, b) => parseInt(a.id.split('_').pop()) - parseInt(b.id.split('_').pop()));
            group.components.sort((a, b) => parseInt(a.id.split('_').pop()) - parseInt(b.id.split('_').pop()));
        });

        const regionNames = Object.keys(regionsMap).sort();
        const regionHeight = height / (regionNames.length || 1);

        const clusterX = 100;
        const baseCompWidth = 60;
        const clusterRectHeight = 20;

        svg.append("text")
            .attr("x", clusterX)
            .attr("y", height)
            .attr("text-anchor", "start")
            .style("font-weight", "bold")
            .style("font-size", "14px")
            .text("Cluster");
        svg.append("text")
            .attr("x", componentRightX)
            .attr("y", height)
            .attr("text-anchor", "end")
            .style("font-weight", "bold")
            .style("font-size", "14px")
            .text("Component");

        regionNames.forEach((regionName, i) => {
            const regionGroup = regionsMap[regionName];
            const yOffset = i * regionHeight;

            svg.append("text")
                .attr("x", 10)
                .attr("y", yOffset + regionHeight / 2)
                .attr("dy", ".35em")
                .text(regionName)
                .style("font-weight", "bold")
                .style("font-size", "14px");

            regionGroup.clusters.forEach((node, index) => {
                const clusterSize = clusterSizeMap[node.id] || 0;
                const dynamicWidth = clusterWidthScale(clusterSize);
                node.clusterWidth = dynamicWidth;
                const y = yOffset + (index + 1) * regionHeight / (regionGroup.clusters.length + 1);
                node.x = clusterX;
                node.y = y;
                svg.append("rect")
                    .attr("x", node.x)
                    .attr("y", node.y - clusterRectHeight / 2)
                    .attr("width", dynamicWidth)
                    .attr("height", clusterRectHeight)
                    .style("fill", "#69b3a2");
                svg.append("text")
                    .attr("x", node.x + 4)
                    .attr("y", node.y)
                    .attr("dy", ".35em")
                    .attr("text-anchor", "start")
                    .text(node.id.split('_').pop());
            });

            // aligning the component area to the right side of the SVG
            const laidOutComponents = layoutComponentNodes(regionGroup.components, yOffset, regionHeight, componentRightX);
            laidOutComponents.forEach(node => {
                const compKey = node.id.split('_').slice(2).join('_');
                let dynamicCompWidth, fillColor;
                if (GOData && GOData[regionName] && GOData[regionName][compKey]) {
                    const goArray = GOData[regionName][compKey];
                    const effectiveLength = Math.min(goArray.length, 5);
                    dynamicCompWidth = baseCompWidth + effectiveLength * 10;
                    fillColor = "#4C9AFF";
                } else {
                    dynamicCompWidth = 30;
                    fillColor = "#d3d3d3";
                }
                // dynamicCompWidth is the width of the component
                node.dynamicCompWidth = dynamicCompWidth;
                // rect‘s x coordinate is the left side of the component
                const rectX = node.x - dynamicCompWidth;
                svg.append("rect")
                    .attr("x", rectX)
                    .attr("y", node.y - node.rectHeight / 2)
                    .attr("width", dynamicCompWidth)
                    .attr("height", node.rectHeight)
                    .style("fill", fillColor);
                svg.append("text")
                    .attr("x", rectX + dynamicCompWidth / 2)
                    .attr("y", node.y)
                    .attr("dy", ".35em")
                    .attr("text-anchor", "middle")
                    .text(compKey.split('_').pop());
            });

            // drawing the links between clusters and components
            // the target x coordinate of the link is the left side of each component
            const linksBySource = {};
            sankeyData.links.forEach(link => {
                const [srcRegion, tgtRegion] = [link.source.split('_')[1], link.target.split('_')[1]];
                if (srcRegion === regionName && tgtRegion === regionName) {
                    if (!linksBySource[link.source]) linksBySource[link.source] = [];
                    linksBySource[link.source].push(link);
                }
            });

            regionGroup.clusters.forEach(source => {
                const links = linksBySource[source.id];
                if (!links) return;
                const total = links.reduce((sum, l) => sum + l.value, 0);
                let offset = 0;
                links.forEach(link => {
                    // search for the target node in laidOutComponents
                    const target = laidOutComponents.find(n => n.id === link.target);
                    if (!target) return;
                    const thickness = total !== 0 ? (link.value / total) * clusterRectHeight : link.value * clusterRectHeight;
                    const thicknessTarget = Math.min(thickness, target.rectHeight);
                    const startY = source.y - clusterRectHeight / 2 + offset + thickness / 2;
                    offset += thickness;
                    const targetLeftX = target.x - target.dynamicCompWidth;
                    const pathD = taperedLinkPath(
                        source.x + (source.clusterWidth || baseCompWidth), startY,
                        targetLeftX, target.y,
                        thickness, thicknessTarget
                    );
                    svg.append("path")
                        .attr("d", pathD)
                        .attr("fill", "#999")
                        .attr("opacity", 0.5);
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
                setNMFGOData(data);
            });
    };

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div style={{ marginTop: 5, marginBottom: 5, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                <Select
                    size='small'
                    mode="multiple"
                    placeholder="Select Region"
                    options={regionOptions}
                    style={{ width: 200, marginRight: 8 }}
                    onChange={setSelectedRegion}
                    allowClear
                />
                <Select
                    size='small'
                    placeholder="Select Component Count"
                    options={componentCountOptions}
                    style={{ width: 120, marginRight: 8 }}
                    onChange={setSelectedComponentCount}
                    allowClear
                />
                <Select
                    size='small'
                    placeholder="Select Resolution"
                    options={resolutionOptions}
                    style={{ width: 120 }}
                    onChange={setSelectedResolution}
                    allowClear
                />
                <Button size='small' onClick={confirm} style={{ marginLeft: 8 }}>Confirm</Button>
            </div>
            <div ref={containerRef} style={{ width: "100%", height: "calc(100% - 50px)" }} />
        </div>
    );
};