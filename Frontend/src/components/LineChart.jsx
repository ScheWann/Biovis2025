import React, { useRef, useEffect, useState } from "react";
import * as d3 from "d3";


export const LineChart = ({
  data,
  xAccessor,
  yAccessor,
  xLabel = "X Axis",
  yLabel = "Y Axis",
  title = "",
  margin = { top: 40, right: 30, bottom: 50, left: 60 },
  showErrorBands = true,
  yMinAccessor,
  yMaxAccessor,
  errorBandColor = "gray",
  errorBandOpacity = 0.3,
  lineColor = "#1d72b8",
  lineWidth = 2,
}) => {
  const svgRef = useRef();
  const containerRef = useRef();
  const [dimensions, setDimensions] = useState({ width: 400, height: 200 });

  // Detect container size changes
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    // Initial measurement
    const updateDimensions = () => {
      const rect = container.getBoundingClientRect();
      setDimensions({
        width: Math.max(rect.width, 200),
        height: Math.max(rect.height, 150),
      });
    };

    // Measure immediately
    updateDimensions();

    const resizeObserver = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        setDimensions({
          width: Math.max(width, 400),
          height: Math.max(height, 200),
        });
      }
    });

    resizeObserver.observe(container);
    return () => resizeObserver.disconnect();
  }, []);

  // Handle array of numbers or array of objects
  let processedData, xAcc, yAcc;
  if (typeof data[0] === "number") {
    processedData = data.map((y, i) => ({ x: i, y }));
    xAcc = (d) => d.x;
    yAcc = (d) => d.y;
  } else {
    processedData = data;
    xAcc = xAccessor;
    yAcc = yAccessor;
  }

  useEffect(() => {
    if (!processedData || processedData.length === 0) return;
    
    d3.select(svgRef.current).selectAll("*").remove();
    const innerWidth = dimensions.width - margin.left - margin.right;
    const innerHeight = dimensions.height - margin.top - margin.bottom;

    if (innerWidth <= 0 || innerHeight <= 0) return;

    const xValues = processedData.map(xAcc);
    const yValues = processedData.map(yAcc);

    // For error bands, consider ymin and ymax values for the scale
    let yDomain;
    if (showErrorBands && yMinAccessor && yMaxAccessor) {
      const yMinValues = processedData.map(yMinAccessor);
      const yMaxValues = processedData.map(yMaxAccessor);
      yDomain = [
        d3.min([...yValues, ...yMinValues]),
        d3.max([...yValues, ...yMaxValues])
      ];
    } else {
      yDomain = d3.extent(yValues);
    }

    const xScale = d3.scaleLinear().domain(d3.extent(xValues)).nice().range([0, innerWidth]);
    const yScale = d3.scaleLinear().domain(yDomain).nice().range([innerHeight, 0]);

    const line = d3
      .line()
      .x((d) => xScale(xAcc(d)))
      .y((d) => yScale(yAcc(d)))
      .curve(d3.curveMonotoneX);

    // Create area generator for error bands
    const area = showErrorBands && yMinAccessor && yMaxAccessor ? d3.area()
      .x((d) => xScale(xAcc(d)))
      .y0((d) => yScale(yMinAccessor(d)))
      .y1((d) => yScale(yMaxAccessor(d)))
      .curve(d3.curveMonotoneX) : null;

    const svg = d3
      .select(svgRef.current)
      .attr("width", dimensions.width)
      .attr("height", dimensions.height);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Add error band first (so it appears behind the line)
    if (area) {
      g.append("path")
        .datum(processedData)
        .attr("fill", errorBandColor)
        .attr("fill-opacity", errorBandOpacity)
        .attr("d", area);
    }

    // Add axes
    g.append("g")
      .attr("transform", `translate(0,${innerHeight})`)
      .call(d3.axisBottom(xScale));

    g.append("g").call(d3.axisLeft(yScale));

    // Add main line
    g.append("path")
      .datum(processedData)
      .attr("fill", "none")
      .attr("stroke", lineColor)
      .attr("stroke-width", lineWidth)
      .attr("d", line);

    // Add labels
    svg
      .append("text")
      .attr("x", margin.left + innerWidth / 2)
      .attr("y", dimensions.height - 12)
      .attr("text-anchor", "middle")
      .attr("font-size", 12)
      .text(xLabel);

    svg
      .append("text")
      .attr("transform", `rotate(-90)`)
      .attr("x", -dimensions.height / 2)
      .attr("y", 16)
      .attr("text-anchor", "middle")
      .attr("font-size", 12)
      .text(yLabel);

    if (title) {
      svg
        .append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 10)
        .attr("text-anchor", "start")
        .attr("font-size", 14)
        .attr("font-weight", "bold")
        .text(title);
    }
  }, [processedData, xAcc, yAcc, dimensions, xLabel, yLabel, title, margin, showErrorBands, yMinAccessor, yMaxAccessor, errorBandColor, errorBandOpacity, lineColor, lineWidth]);

  return (
    <div ref={containerRef} style={{ width: '100%', height: '100%', minWidth: 0, minHeight: 0 }}>
      <svg ref={svgRef}></svg>
    </div>
  );
};
