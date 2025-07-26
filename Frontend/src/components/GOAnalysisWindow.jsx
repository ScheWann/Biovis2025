import React, { useRef, useEffect, useState } from "react";
import { Spin } from "antd";
import { CloseOutlined } from "@ant-design/icons";
import * as d3 from "d3";

export const GOAnalysisWindow = ({
    visible,
    setVisible,
    loading,
    data,
    position, // {x, y} - position where the cluster was clicked
    title = "GO Analysis Results",
}) => {
    const tooltipRef = useRef();
    const svgRef = useRef();
    const [tooltipPosition, setTooltipPosition] = useState({ x: 0, y: 0 });

    // Calculate tooltip position with boundary detection
    useEffect(() => {
        if (!visible || !position || !tooltipRef.current) return;

        const tooltip = tooltipRef.current;
        const tooltipWidth = 500;
        const tooltipHeight = 400;
        const padding = 20;

        // Get viewport dimensions
        const viewportWidth = window.innerWidth;
        const viewportHeight = window.innerHeight;

        let x = position.x + 20; // Default: 20px to the right of click
        let y = position.y - tooltipHeight / 2; // Center vertically on click

        // Boundary detection - horizontal
        if (x + tooltipWidth + padding > viewportWidth) {
            x = position.x - tooltipWidth - 20; // Show to the left instead
        }
        if (x < padding) {
            x = padding; // Ensure it doesn't go off the left edge
        }

        // Boundary detection - vertical
        if (y < padding) {
            y = padding; // Don't go above viewport
        }
        if (y + tooltipHeight + padding > viewportHeight) {
            y = viewportHeight - tooltipHeight - padding; // Don't go below viewport
        }

        setTooltipPosition({ x, y });
    }, [visible, position]);

    // Create D3 chart
    useEffect(() => {
        if (!data || loading || !visible || !svgRef.current) return;

        // Clear previous chart
        d3.select(svgRef.current).selectAll("*").remove();

        const margin = { top: 20, right: 20, bottom: 60, left: 120 };
        const width = 420 - margin.left - margin.right;
        const height = 300 - margin.top - margin.bottom;

        // Sort data by combined_score in descending order and take top 10
        const sortedData = [...data]
            .sort((a, b) => b.combined_score - a.combined_score)
            .slice(0, 10);

        if (sortedData.length === 0) return;

        // Create scales
        const xScale = d3
            .scaleLinear()
            .domain([0, d3.max(sortedData, (d) => d.combined_score)])
            .range([0, width]);

        const yScale = d3
            .scaleBand()
            .domain(sortedData.map((d) => d.name || d.term))
            .range([0, height])
            .padding(0.1);

        // Create SVG
        const svg = d3
            .select(svgRef.current)
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom);

        const g = svg
            .append("g")
            .attr("transform", `translate(${margin.left},${margin.top})`);

        // Create bars
        g.selectAll(".bar")
            .data(sortedData)
            .enter()
            .append("rect")
            .attr("class", "bar")
            .attr("x", 0)
            .attr("y", (d) => yScale(d.name || d.term))
            .attr("width", (d) => xScale(d.combined_score))
            .attr("height", yScale.bandwidth())
            .attr("fill", "#1890ff")
            .attr("opacity", 0.8)
            .on("mouseover", function (event, d) {
                d3.select(this).attr("opacity", 1);

                // Create mini tooltip
                const miniTooltip = d3
                    .select("body")
                    .append("div")
                    .attr("class", "mini-tooltip")
                    .style("position", "absolute")
                    .style("background", "rgba(0, 0, 0, 0.9)")
                    .style("color", "white")
                    .style("padding", "6px 8px")
                    .style("border-radius", "4px")
                    .style("font-size", "11px")
                    .style("pointer-events", "none")
                    .style("z-index", 10001);

                miniTooltip
                    .html(
                        `
                        <strong>${d.name || d.term}</strong><br/>
                        Score: ${d.combined_score.toFixed(3)}
                        `
                    )
                    .style("left", event.pageX + 10 + "px")
                    .style("top", event.pageY - 10 + "px");
            })
            .on("mouseout", function () {
                d3.select(this).attr("opacity", 0.8);
                d3.selectAll(".mini-tooltip").remove();
            });

        // Add value labels on bars
        g.selectAll(".label")
            .data(sortedData)
            .enter()
            .append("text")
            .attr("class", "label")
            .attr("x", (d) => xScale(d.combined_score) + 3)
            .attr("y", (d) => yScale(d.name || d.term) + yScale.bandwidth() / 2)
            .attr("dy", "0.35em")
            .style("font-size", "9px")
            .style("fill", "#333")
            .text((d) => d.combined_score.toFixed(2));

        // Add X axis
        g.append("g")
            .attr("transform", `translate(0,${height})`)
            .call(d3.axisBottom(xScale).ticks(4))
            .selectAll("text")
            .style("font-size", "10px");

        // Add Y axis
        g.append("g")
            .call(d3.axisLeft(yScale))
            .selectAll("text")
            .style("font-size", "9px")
            .call(wrap, margin.left - 10);

        // Add X axis label
        g.append("text")
            .attr(
                "transform",
                `translate(${width / 2}, ${height + margin.bottom - 10})`
            )
            .style("text-anchor", "middle")
            .style("font-size", "11px")
            .style("font-weight", "bold")
            .text("Combined Score");

        // Function to wrap text
        function wrap(text, width) {
            text.each(function () {
                const text = d3.select(this);
                const words = text.text().split(/\s+/).reverse();
                let word;
                let line = [];
                let lineNumber = 0;
                const lineHeight = 1.1;
                const y = text.attr("y");
                const dy = parseFloat(text.attr("dy")) || 0;
                let tspan = text
                    .text(null)
                    .append("tspan")
                    .attr("x", -5)
                    .attr("y", y)
                    .attr("dy", dy + "em");

                while ((word = words.pop())) {
                    line.push(word);
                    tspan.text(line.join(" "));
                    if (tspan.node().getComputedTextLength() > width) {
                        line.pop();
                        tspan.text(line.join(" "));
                        line = [word];
                        tspan = text
                            .append("tspan")
                            .attr("x", -5)
                            .attr("y", y)
                            .attr("dy", ++lineNumber * lineHeight + dy + "em")
                            .text(word);
                    }
                }
            });
        }
    }, [data, loading, visible]);

    if (!visible) return null;

    return (
        <div
            ref={tooltipRef}
            style={{
                position: "fixed",
                left: `${tooltipPosition.x}px`,
                top: `${tooltipPosition.y}px`,
                width: "500px",
                maxHeight: "400px",
                background: "white",
                border: "1px solid #d9d9d9",
                borderRadius: "8px",
                boxShadow: "0 4px 12px rgba(0, 0, 0, 0.15)",
                zIndex: 100000,
                overflow: "hidden",
            }}
        >
            {/* Header */}
            <div
                style={{
                    padding: "5px 5px 5px 10px",
                    borderBottom: "1px solid #f0f0f0",
                    background: "#fafafa",
                    display: "flex",
                    justifyContent: "space-between",
                    alignItems: "center",
                }}
            >
                <div style={{ fontSize: "14px", fontWeight: "600", color: "#262626" }}>
                    {title}
                </div>
                <CloseOutlined
                    style={{ cursor: "pointer", color: "#8c8c8c", fontSize: "12px" }}
                    onClick={() => setVisible(false)}
                />
            </div>

            {/* Content */}
            <div
                style={{
                    padding: "16px",
                    display: "flex",
                    justifyContent: "center",
                    alignItems: "center",
                    minHeight: "300px",
                }}
            >
                {loading ? (
                    <div
                        style={{
                            display: "flex",
                            flexDirection: "column",
                            alignItems: "center",
                            gap: "12px",
                        }}
                    >
                        <Spin size="large" />
                        <div style={{ fontSize: "13px", color: "#666" }}>
                            Loading GO Analysis...
                        </div>
                    </div>
                ) : data && data.length > 0 ? (
                    <svg ref={svgRef}></svg>
                ) : (
                    <div style={{ fontSize: "13px", color: "#999" }}>
                        No GO analysis data available
                    </div>
                )}
            </div>
        </div>
    );
};
