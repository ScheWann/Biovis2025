import React, { useRef, useEffect } from "react";
import * as d3 from "d3";

const COLORS = d3.schemeCategory10;

export const ScatterplotUmap = ({
  data, // [{x, y, cluster}]
  width = 400,
  height = 200,
  pointSize = 5,
  clusterAccessor = (d) => d.cluster,
  xAccessor = (d) => d.x,
  yAccessor = (d) => d.y,
  title = "UMAP",
  margin = { top: 30, right: 30, bottom: 30, left: 30 },
}) => {
  const ref = useRef();

  useEffect(() => {
    d3.select(ref.current).selectAll("*").remove();

    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    // Scales
    const xExtent = d3.extent(data, xAccessor);
    const yExtent = d3.extent(data, yAccessor);

    const xScale = d3.scaleLinear().domain(xExtent).range([0, innerWidth]);
    const yScale = d3.scaleLinear().domain(yExtent).range([innerHeight, 0]);

    // Color by cluster
    const clusters = Array.from(new Set(data.map(clusterAccessor)));
    const color = d3.scaleOrdinal().domain(clusters).range(COLORS);

    // SVG
    const svg = d3
      .select(ref.current)
      .attr("width", width)
      .attr("height", height);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Points
    g.selectAll("circle")
      .data(data)
      .enter()
      .append("circle")
      .attr("cx", (d) => xScale(xAccessor(d)))
      .attr("cy", (d) => yScale(yAccessor(d)))
      .attr("r", pointSize)
      .attr("fill", (d) => color(clusterAccessor(d)))
      .attr("opacity", 0.8);

    // Axes
    g.append("g")
      .attr("transform", `translate(0,${innerHeight})`)
      .call(d3.axisBottom(xScale).ticks(5));

    g.append("g").call(d3.axisLeft(yScale).ticks(5));

    // Title
    svg
      .append("text")
      .attr("x", margin.left)
      .attr("y", margin.top - 10)
      .attr("font-size", 16)
      .attr("font-weight", 600)
      .text(title);

    // Legend
    if (clusters.length > 1) {
      const legend = svg
        .append("g")
        .attr(
          "transform",
          `translate(${width - margin.right - 40}, ${margin.top})`
        );
      clusters.forEach((cl, i) => {
        legend
          .append("circle")
          .attr("cx", 0)
          .attr("cy", i * 20)
          .attr("r", 6)
          .attr("fill", color(cl));
        legend
          .append("text")
          .attr("x", 12)
          .attr("y", i * 20 + 4)
          .text(cl)
          .attr("font-size", 12);
      });
    }
  }, [
    data,
    width,
    height,
    pointSize,
    clusterAccessor,
    xAccessor,
    yAccessor,
    title,
    margin,
  ]);

  return <svg ref={ref}></svg>;
};
