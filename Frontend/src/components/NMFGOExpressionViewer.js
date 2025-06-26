import React, { useRef, useState, useEffect } from 'react';
import { Empty, Spin } from 'antd';
import * as d3 from 'd3';

export const NMFGOExpressionViewer = ({ NMFGOData, NMFGODataLoading, setNMFclusterCells }) => {
    const containerRef = useRef(null);
    const svgRef = useRef(null);
    const [activeCluster, setActiveCluster] = useState(null);
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
    }, [NMFGODataLoading]);

    // function to wrap text by chunks
    const wrapByChunks = (text, maxCharsPerLine) => {
        text.each(function () {
            const text = d3.select(this);
            const fullText = text.text();
            const parts = fullText.split(' ');

            const totalLength = fullText.length;
            if (totalLength <= maxCharsPerLine) {
                text.text(null)
                    .append("tspan")
                    .attr("x", 0)
                    .attr("y", 0)
                    .attr("dy", "0em")
                    .attr("alignment-baseline", "middle")
                    .attr("text-anchor", "end")
                    .text(fullText);
            } else {
                let line = [];
                let tspan = text.text(null)
                    .append("tspan")
                    .attr("x", 0)
                    .attr("dy", "0em");
    
                parts.forEach((word, i) => {
                    line.push(word);
                    const currentLine = line.join(' ');
                    if (currentLine.length > maxCharsPerLine || i === parts.length - 1) {
                        if (currentLine.length > maxCharsPerLine) {
                            line.pop(); // remove last word
                            tspan.text(line.join(' '));
                            line = [word];
                        } else {
                            tspan.text(currentLine);
                            line = [];
                        }
                        tspan = text.append("tspan")
                            .attr("x", 0)
                            .attr("dy", "1em");
                    }
                });
            }
        });
    };
    

    useEffect(() => {
        if (dimensions.width === 0 || dimensions.height === 0 || Object.keys(NMFGOData).length === 0 || NMFGODataLoading) return;

        const svg = d3.select(svgRef.current);
        svg.attr("width", dimensions.width).attr("height", dimensions.height);

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
        const xAxis = d3.axisBottom(xScale)
            .tickFormat(d => d.replace('Component_', ''));

        svg.append("text")
            .attr("x", dimensions.width / 2)
            .attr("y", 35)
            .attr("text-anchor", "middle")
            .attr("font-size", "14px")
            .attr("font-weight", "bold")
            .text("NMF Component Expression Heatmap");

        g.append('g')
            .attr('transform', `translate(0, ${height})`)
            .call(xAxis)
            .selectAll("text")
            .attr("dy", "1em")
            .style("text-anchor", "middle")
            .style("cursor", "pointer")
            .on("mouseover", (event, d) => {
                const tooltipData = NMFGOData.GO_results[d];
                if (!tooltipData) return;

                // barchart tooltip
                const tooltip = d3.select('body')
                    .append('div')
                    .attr('class', 'go-tooltip')
                    .style('position', 'absolute')
                    .style('background', 'white')
                    .style('border', '1px solid #ccc')
                    .style('padding', '10px')
                    .style('box-shadow', '0px 0px 5px rgba(0,0,0,0.3)')
                    .style('pointer-events', 'none')
                    .style('z-index', 1000)
                    .style('opacity', 0);

                tooltip.style('left', (event.pageX + 15) + 'px')
                    .style('top', (event.pageY - 40) + 'px');

                // barchart svg
                const barWidth = 500;
                const barHeight = 200;
                const margin = { top: 50, right: 10, bottom: 35, left: 170 };
                const innerWidth = barWidth - margin.left - margin.right;
                const innerHeight = barHeight - margin.top - margin.bottom;

                const svg = tooltip.append('svg')
                    .attr('width', barWidth)
                    .attr('height', barHeight)
                    .append('g')
                    .attr('transform', `translate(${margin.left},${margin.top})`);

                // Go barchart Title
                svg.append("text")
                    .attr("x", innerWidth / 2)
                    .attr("y", -35)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "14px")
                    .attr("font-weight", "bold")
                    .text("GO Term Enrichment");

                const y = d3.scaleBand()
                    .domain(tooltipData.map(d => d.Term))
                    .range([0, innerHeight])
                    .padding(0.1);

                const x = d3.scaleLinear()
                    .domain([0, d3.max(tooltipData, d => d['Combined Score'])])
                    .range([0, innerWidth]);

                const oddsX = d3.scaleLinear()
                    .domain([0, d3.max(tooltipData, d => d['Odds Ratio'])])
                    .range([0, innerWidth]);

                svg.selectAll('.bar')
                    .data(tooltipData)
                    .enter()
                    .append('rect')
                    .attr('class', 'bar')
                    .attr('y', d => y(d.Term))
                    .attr('width', d => x(d['Combined Score']))
                    .attr('height', y.bandwidth())
                    .attr('fill', '#69b3a2');

                // add the y axis
                svg.append('g')
                    .call(d3.axisLeft(y)
                        .tickSize(0)
                        .tickPadding(5)
                        .tickFormat(d => d.replace(/_/g, ' '))
                    )
                    .selectAll("text")
                    .style("font-size", "10px")
                    .call(wrapByChunks, 25);

                // svg.append("text")
                //     .attr("text-anchor", "middle")
                //     .attr("transform", `rotate(-90)`)
                //     .attr("x", -innerHeight / 2)
                //     .attr("y", -margin.left + 15)
                //     .attr("font-size", "12px")
                //     .text("GO Term");

                // add the x axis
                svg.append('g')
                    .attr('transform', `translate(0, ${innerHeight})`)
                    .call(d3.axisBottom(x).ticks(4))
                    .selectAll("text")
                    .style("font-size", "10px");

                svg.append("text")
                    .attr("x", innerWidth - 40)
                    .attr("y", innerHeight + 30)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "10px")
                    .style("font-weight", "bold")
                    .text("Combined Score");

                // add the odds ratio axis
                svg.append('g')
                    .attr('transform', `translate(0, 0)`)
                    .call(d3.axisTop(oddsX).ticks(4))
                    .selectAll("text")
                    .style("font-size", "10px");

                svg.append("text")
                    .attr("x", innerWidth - 30)
                    .attr("y", -margin.top + 30)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "10px")
                    .style("font-weight", "bold")
                    .text("Odds Ratio");

                // add the odds ratio dot
                svg.selectAll(".odds-dot")
                    .data(tooltipData)
                    .enter()
                    .append("circle")
                    .attr("cx", d => oddsX(d["Odds Ratio"]))
                    .attr("cy", d => y(d.Term) + y.bandwidth() / 2)
                    .attr("r", 3)
                    .attr("fill", "#d62728")
                    .attr("stroke", "#000")
                    .attr("stroke-width", 0.5);

                // add the genes label
                svg.selectAll(".genes-label")
                    .data(tooltipData)
                    .enter()
                    .append("text")
                    .attr("x", d => x(d["Combined Score"]) / 2)
                    .attr("y", d => y(d.Term) + y.bandwidth() / 2 + 4) // vertical center
                    .attr("text-anchor", "middle")
                    .attr("font-size", "9px")
                    .attr("fill", "white")
                    .text(d => d.Genes);

                // position the tooltip
                // wait for the tooltip to be added to the DOM
                // before calculating its size and position
                // to avoid flickering
                // and to ensure the tooltip is positioned correctly
                requestAnimationFrame(() => {
                    const tooltipNode = tooltip.node();
                    const rect = tooltipNode.getBoundingClientRect();

                    let left = event.pageX + 15;
                    let top = event.pageY - 40;

                    if (left + rect.width > window.innerWidth) {
                        left = event.pageX - rect.width - 15;
                    }

                    if (top + rect.height > window.innerHeight) {
                        top = window.innerHeight - rect.height - 10;
                    }

                    tooltip
                        .style('left', `${left}px`)
                        .style('top', `${top}px`)
                        .style('opacity', 1);
                });
            })
            .on("mouseout", function () {
                d3.select('.go-tooltip').remove();
            });

        // for each component, check if it has a GO result
        // if not, draw a red cross
        compKeys.forEach(comp => {
            if (!NMFGOData.GO_results[comp]) {
                const xPos = xScale(comp) + xScale.bandwidth() / 2;

                g.append("text")
                    .attr("x", xPos)
                    .attr("y", height + 5)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "14px")
                    .attr("fill", "red")
                    .text("×");
            }
        });
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
            .call(yAxis)
            .selectAll("text")
            .style("font-size", "10px")
            .style("cursor", "pointer")
            .style("fill", d => d === activeCluster ? "#CD853F" : "black")
            .style("font-weight", d => d === activeCluster ? "bold" : "normal")
            .on("click", (event, d) => {
                const clusterIndex = d;
                const isActive = activeCluster === clusterIndex;

                if (isActive) {
                    setActiveCluster(null);
                    setNMFclusterCells([]);
                } else {
                    setActiveCluster(clusterIndex);
                    const cellIds = NMFGOData.cell_ids_by_cluster[clusterIndex] || [];
                    setNMFclusterCells(cellIds);
                }

                event.stopPropagation();
            });

        svg.on("click", () => {
            setActiveCluster(null);
            setNMFclusterCells([]);
        });

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
        const legendWidth = 10;
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
    }, [NMFGOData, dimensions, NMFGODataLoading, activeCluster]);

    return (
        NMFGODataLoading ? (
            <Spin spinning={true} style={{ width: '100%', height: '100%' }} />
        ) : (
            Object.keys(NMFGOData).length > 0 && NMFGOData.cluster_means && Object.keys(NMFGOData.cluster_means).length > 0 ? (
                <div ref={containerRef} style={{ width: '100%', height: '100%' }}>
                    <svg ref={svgRef}></svg>
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
