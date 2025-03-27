import React, { useRef, useEffect } from "react";
import useSVGCanvas from "./useSVGCanvas";
import * as d3 from "d3";

const data = {
  senderData: [
    { cell_type: "CD4_Tcell", region: "1", count: 10 },
    { cell_type: "CD4_Tcell", region: "2", count: 20 },
    { cell_type: "CD4_Tcell", region: "3", count: 30 },
    { cell_type: "Macrophage", region: "1", count: 40 },
    { cell_type: "Macrophage", region: "2", count: 50 },
    { cell_type: "Macrophage", region: "3", count: 60 },
    { cell_type: "B_cell", region: "1", count: 70 },
    { cell_type: "B_cell", region: "2", count: 80 },
    { cell_type: "B_cell", region: "3", count: 90 },
    { cell_type: "NK_cell", region: "1", count: 100 },
    { cell_type: "NK_cell", region: "2", count: 110 },
    { cell_type: "NK_cell", region: "3", count: 120 },
  ],
  receiverData: [
    { cell_type: "CD4_Tcell", region: "1", count: 10 },
    { cell_type: "CD4_Tcell", region: "2", count: 20 },
    { cell_type: "CD4_Tcell", region: "3", count: 30 },
    { cell_type: "Macrophage", region: "1", count: 40 },
    { cell_type: "Macrophage", region: "2", count: 50 },
    { cell_type: "Macrophage", region: "3", count: 60 },
    { cell_type: "B_cell", region: "1", count: 70 },
    { cell_type: "B_cell", region: "2", count: 80 },
    { cell_type: "B_cell", region: "3", count: 90 },
    { cell_type: "NK_cell", region: "1", count: 100 },
    { cell_type: "NK_cell", region: "2", count: 110 },
    { cell_type: "NK_cell", region: "3", count: 120 },
  ],
};

export const Cell2CellViewer = ({}) => {
  const d3Container = useRef(null);
  const [svg, height, width] = useSVGCanvas(d3Container);

  const margin = { top: 40, right: 40, bottom: 60, left: 60 };

  useEffect(() => {
    if (data === undefined || svg === undefined) return;

    svg.selectAll("*").remove();

    // Get unique regions and cell types
    const regions = Array.from(
      new Set(data.senderData.map((d) => d.region))
    ).sort();
    const cellTypes = Array.from(
      new Set(data.senderData.map((d) => d.cell_type))
    );

    const processData = (dataArray) => {
      return d3.rollup(
        dataArray,
        (v) => {
          const regionCounts = {};
          regions.forEach((region) => {
            regionCounts[region] =
              v.find((d) => d.region === region)?.count || 0;
          });
          return regionCounts;
        },
        (d) => d.cell_type
      );
    };

    // Process sender data
    const senderData = processData(data.senderData);

    // Process receiver data
    const receiverData = processData(data.receiverData);

    // Calculate maximum stack values
    const calculateMaxStack = (data) => {
      return d3.max(cellTypes, (cellType) => {
        return d3.sum(regions, (region) => data.get(cellType)?.[region] || 0);
      });
    };

    const maxSenderStack = calculateMaxStack(senderData);
    const maxReceiverStack = calculateMaxStack(receiverData);
    const maxStackValue = Math.max(maxSenderStack, maxReceiverStack);

    // Set up scales
    const x = d3
      .scaleBand()
      .domain(cellTypes)
      .range([margin.left, width - margin.right])
      .padding(0.2);

    const y = d3
      .scaleLinear()
      .domain([-maxStackValue * 1.1, maxStackValue * 1.1])
      .range([height - margin.bottom, margin.top]);

    // Color scale
    const color = d3.scaleOrdinal().domain(regions).range(d3.schemeTableau10);

    const stackGenerator = d3
      .stack()
      .keys(regions)
      .value((d, key) => d[key]);

    const senderSeries = stackGenerator(
      Array.from(senderData, ([cell_type, counts]) => ({
        cell_type,
        ...counts,
      }))
    );

    const receiverSeries = stackGenerator(
      Array.from(receiverData, ([cell_type, counts]) => ({
        cell_type,
        ...Object.fromEntries(Object.entries(counts).map(([k, v]) => [k, -v])),
      }))
    );

    // Draw bars
    svg
      .append("g")
      .selectAll("g")
      .data([...senderSeries, ...receiverSeries])
      .join("g")
      .attr("fill", (d) => color(d.key))
      .selectAll("rect")
      .data((d) =>
        d.map((segment) => ({
          ...segment,
          region: d.key,
        }))
      )
      .join("rect")
      .attr("x", (d) => x(d.data.cell_type))
      .attr("y", (d) => y(Math.max(d[0], d[1])))
      .attr("height", (d) => Math.abs(y(d[0]) - y(d[1])))
      .attr("width", x.bandwidth())
      .append("title")
      .text((d) => {
        const value = d[1] - d[0];
        const direction = value > 0 ? "Sender" : "Receiver";
        return `${direction}\nCell Type: ${d.data.cell_type}\nRegion: ${
          d.region
        }\nCount: ${Math.abs(value)}`;
      });

    // Add x-axis
    svg
      .append("g")
      .attr("transform", `translate(0,${y(0)})`)
      .call(d3.axisBottom(x))
      .selectAll("text")
      .attr("transform", "rotate(-45)")
      .style("text-anchor", "end");

    // Add y-axis
    svg
      .append("g")
      .attr("transform", `translate(${margin.left},0)`)
      .call(d3.axisLeft(y).tickFormat((d) => Math.abs(d))); // Show absolute values

    // Add chart title
    svg
      .append("text")
      .attr("x", width / 2)
      .attr("y", margin.top / 2)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .style("font-weight", "bold")
      .text("Cell-to-Cell Interaction");

    // Add legend
    const legend = svg
      .append("g")
      .attr(
        "transform",
        `translate(${width / 2 - (regions.length * 80) / 2},30)`
      ); // Centered horizontally

    regions.forEach((region, i) => {
      const legendItem = legend
        .append("g")
        .attr("transform", `translate(${i * 100},0)`); // Horizontal layout

      legendItem
        .append("rect")
        .attr("width", 15)
        .attr("height", 15)
        .attr("fill", color(region));

      legendItem
        .append("text")
        .attr("x", 20)
        .attr("y", 12)
        .text(`Region ${region}`);
    });

    // Add labels for sender/receiver
    svg
      .append("text")
      .attr("x", 5)
      .attr("y", y(maxStackValue * 0.8))
      .attr("text-anchor", "start")
      .style("font-size", "12px")
      .text("Sender");

    svg
      .append("text")
      .attr("x", 5)
      .attr("y", y(-maxStackValue * 0.8))
      .attr("text-anchor", "start")
      .style("font-size", "12px")
      .text("Receiver");
  }, [svg, data]);

  return (
    <div
      style={{
        height: "100%",
        width: "100%",
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
      }}
      ref={d3Container}
    />
  );
};
