import React, { useRef, useEffect } from "react";
import * as d3 from "d3";

// Choose scale based on data type
const getScale = (values, range) => {
  const sample = values[0];
  if (sample instanceof Date) {
    return d3.scaleTime().domain(d3.extent(values)).range(range);
  }
  if (typeof sample === "number") {
    return d3.scaleLinear().domain(d3.extent(values)).nice().range(range);
  }
  return d3.scalePoint().domain(values).range(range).padding(0.5);
};

export const LineChart = ({
  data,
  xAccessor,
  yAccessor,
  width = 450,
  height = 150,
  xLabel = "X Axis",
  yLabel = "Y Axis",
  title = "",
  margin = { top: 40, right: 30, bottom: 50, left: 60 },
}) => {
  const ref = useRef();

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
    d3.select(ref.current).selectAll("*").remove();
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    const xValues = processedData.map(xAcc);
    const yValues = processedData.map(yAcc);

    const xScale = getScale(xValues, [0, innerWidth]);
    const yScale = getScale(yValues, [innerHeight, 0]);

    const line = d3
      .line()
      .x((d) => xScale(xAcc(d)))
      .y((d) => yScale(yAcc(d)));

    const svg = d3
      .select(ref.current)
      .attr("width", width)
      .attr("height", height);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    g.append("g")
      .attr("transform", `translate(0,${innerHeight})`)
      .call(d3.axisBottom(xScale));

    g.append("g").call(d3.axisLeft(yScale));

    g.append("path")
      .datum(processedData)
      .attr("fill", "none")
      .attr("stroke", "#1d72b8")
      .attr("stroke-width", 2)
      .attr("d", line);

    svg
      .append("text")
      .attr("x", margin.left + innerWidth / 2)
      .attr("y", height - 8)
      .attr("text-anchor", "middle")
      .attr("font-size", 14)
      .text(xLabel);

    svg
      .append("text")
      .attr("transform", `rotate(-90)`)
      .attr("x", -height / 2)
      .attr("y", 16)
      .attr("text-anchor", "middle")
      .attr("font-size", 14)
      .text(yLabel);

    if (title) {
      svg
        .append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 15)
        .attr("text-anchor", "start")
        .attr("font-size", 16)
        .attr("font-weight", "bold")
        .text(title);
    }
  }, [processedData, xAcc, yAcc, width, height, xLabel, yLabel, title, margin]);

  return <svg ref={ref}></svg>;
};
