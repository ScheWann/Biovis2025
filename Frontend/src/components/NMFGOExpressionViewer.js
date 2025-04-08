import React, { useRef, useState, useEffect } from 'react';
import { Empty, Spin } from 'antd';
import * as d3 from 'd3';

export const NMFGOExpressionViewer = ({ NMFGOData, NMFGODataLoading }) => {
    const containerRef = useRef(null);
    const svgRef = useRef(null);
    const [dimensions, setDimensions] = useState({ width: 0, height: 0 });

    // observe container size changes
    useEffect(() => {
        const observer = new ResizeObserver(entries => {
            entries.forEach(entry => {
                const { width, height } = entry.contentRect;
                setDimensions({ width, height });
            });
        });
        if (containerRef.current) {
            observer.observe(containerRef.current);
        }
        return () => observer.disconnect();
    }, []);

    useEffect(() => {
        console.log(NMFGODataLoading, 'NMFGODataLoading')
        if (dimensions.width === 0 || dimensions.height === 0 || Object.keys(NMFGOData).length === 0 || !NMFGODataLoading) return;

        const svg = d3.select(svgRef.current);
        // clear previous SVG
        svg.selectAll('*').remove();

        const margin = { top: 40, right: 50, bottom: 60, left: 80 };
        const width = dimensions.width - margin.left - margin.right;
        const height = dimensions.height - margin.top - margin.bottom;

        const clusters = Object.keys(NMFGOData.cluster_means).sort((a, b) => +a - +b);

        // sort the component keys
        const compKeys = Object.keys(NMFGOData.cluster_means[clusters[0]]).sort((a, b) => {
            const numA = parseInt(a.split('_')[1]);
            const numB = parseInt(b.split('_')[1]);
            return numA - numB;
        });

        const matrix = clusters.map(clusterKey => {
            const row = NMFGOData.cluster_means[clusterKey];
            return compKeys.map(compKey => row[compKey]);
        });

        const rows = matrix.length;
        const cols = compKeys.length;

        // XScale
        const xScale = d3.scaleBand()
            .domain(compKeys)
            .range([0, width])
            .padding(0.05);

        // YScale
        const yScale = d3.scaleBand()
            .domain(clusters)
            .range([0, height])
            .padding(0.05);

        // colorScale
        const flatData = matrix.flat();
        const colorScale = d3.scaleSequential()
            .interpolator(d3.interpolateViridis)
            .domain([d3.min(flatData), d3.max(flatData)]);

        const g = svg.append('g')
            .attr('transform', `translate(${margin.left}, ${margin.top})`);

        // draw the heatmap rectangles
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                g.append('rect')
                    .attr('x', xScale(compKeys[j]))
                    .attr('y', yScale(clusters[i]))
                    .attr('width', xScale.bandwidth())
                    .attr('height', yScale.bandwidth())
                    .style('fill', colorScale(matrix[i][j]));
            }
        }

        // x axis
        const xAxis = d3.axisBottom(xScale);
        g.append('g')
            .attr('transform', `translate(0, ${height})`)
            .call(xAxis)
            .selectAll("text")
            .attr("dy", "1em")
            .style("text-anchor", "middle");

        svg.append('text')
            .attr('class', 'x label')
            .attr('text-anchor', 'middle')
            .attr('x', margin.left + width / 2)
            .attr('y', dimensions.height - 20)
            .attr('font-size', '12px')
            .attr('font-weight', 'bold')
            .text('Components');

        // y axis
        const yAxis = d3.axisLeft(yScale)
            .tickFormat(d => `Cluster ${d}`);
        g.append('g')
            .call(yAxis);

        svg.append('text')
            .attr('class', 'y label')
            .attr('text-anchor', 'middle')
            .attr('x', -(margin.top + height / 2))
            .attr('y', 15)
            .attr('transform', 'rotate(-90)')
            .attr('font-size', '12px')
            .attr('font-weight', 'bold')
            .text('Cluster');

        // legend
        const legendWidth = 20;
        const legendHeight = dimensions.height - margin.top - margin.bottom;

        const legend = svg.append('g')
            .attr('class', 'legend')
            .attr('transform', `translate(${dimensions.width - margin.right + 5}, ${margin.top})`);

        // color gradient
        const defs = svg.append('defs');
        const gradient = defs.append('linearGradient')
            .attr('id', 'legend-gradient')
            .attr('x1', '0%')
            .attr('y1', '100%')
            .attr('x2', '0%')
            .attr('y2', '0%');

        const legendStops = 10;
        const legendScale = d3.scaleLinear()
            .domain([0, legendStops - 1])
            .range([d3.min(flatData), d3.max(flatData)]);
        for (let i = 0; i < legendStops; i++) {
            gradient.append('stop')
                .attr('offset', `${(i / (legendStops - 1)) * 100}%`)
                .attr('stop-color', colorScale(legendScale(i)));
        }

        legend.append('rect')
            .attr('width', legendWidth)
            .attr('height', legendHeight)
            .style('fill', 'url(#legend-gradient)');

        // legend axis
        const legendAxisScale = d3.scaleLinear()
            .domain([d3.min(flatData), d3.max(flatData)])
            .range([legendHeight, 0]);
        const legendAxis = d3.axisRight(legendAxisScale)
            .ticks(6);
        legend.append('g')
            .attr('transform', `translate(${legendWidth},0)`)
            .call(legendAxis);
    }, [NMFGOData, dimensions, NMFGODataLoading]);

    return (
        NMFGODataLoading ? (
            <Spin spinning={true} style={{ width: '100%', height: '100%' }} />
        ) : (
            Object.keys(NMFGOData).length > 0 && NMFGOData.cluster_means && Object.keys(NMFGOData.cluster_means).length > 0 ? (
                <div ref={containerRef} style={{ width: '100%', height: '100%' }}>
                    <svg ref={svgRef} width={dimensions.width} height={dimensions.height}></svg>
                </div>
            ) : (
                <Empty
                    image={Empty.PRESENTED_IMAGE_SIMPLE}
                    description="Please select a region first"
                    style={{ width: '100%', height: '100%' }}
                />
            )
        )
    );
};
