import React, { useRef, useState, useEffect } from "react";
import * as d3 from "d3";
import useSVGCanvas from "./useSVGCanvas";

// Dummy data (cells × genes) reference for python backend
// to run for now, put this at the top of app.js and send this to GeneExpressionViewer
const data = {
  metadata: {
    cell_ids: ["Cell_1", "Cell_2", "Cell_3", "Cell_4", "Cell_5"],
    genes: ["Gene_TP53", "Gene_CDKN1A", "Gene_EGFR", "Gene_KRAS"],
    cell_type_annotations: {
      Cell_1: "CD4_Tcell",
      Cell_2: "CD8_Tcell",
      Cell_3: "B_cell",
      Cell_4: "Macrophage",
      Cell_5: "NK_cell",
    },
  },
  expression_data: [
    { cell_id: "Cell_1", expression: [5.2, 3.1, 8.4, 0.5] },
    { cell_id: "Cell_2", expression: [1.7, 6.3, 2.9, 4.8] },
    { cell_id: "Cell_3", expression: [7.5, 2.4, 0.9, 5.6] },
    { cell_id: "Cell_4", expression: [3.8, 9.1, 4.2, 1.3] },
    { cell_id: "Cell_5", expression: [0.5, 4.8, 5.6, 2.1] },
  ],
};

export const GeneExpressionViewer = ({}) => {
  const d3Container = useRef(null);
  const [svg, height, width] = useSVGCanvas(d3Container);

  const margin = { top: 20, right: 50, bottom: 50, left: 20 };

  useEffect(() => {
    if (data === undefined || svg === undefined) return;
    // Clear previous SVG
    svg.selectAll("*").remove();

    // get the cell types
    const cellTypeMap = data.metadata.cell_type_annotations;

    // sort the data by cell type
    const sortedCells = Object.entries(cellTypeMap)
      .map(([cellId, type]) => ({ cellId, type }))
      .sort((a, b) => a.type.localeCompare(b.type));

    // ordered cell ids and types
    const orderedCellIds = sortedCells.map((d) => d.cellId);
    const cellTypes = sortedCells.map((d) => d.type);

    const genes = data.metadata.genes;
    const expressionMatrix = orderedCellIds.map((cellId) => {
      const cellData = data.expression_data.find((d) => d.cell_id === cellId);
      return cellData ? cellData.expression : Array(genes.length).fill(0);
    });

    const viewWidth = width - margin.left - margin.right;
    const viewHeight = height - margin.top - margin.bottom;

    // create a color scale for heatmap
    // from 0 to positive white to red
    // from 0 to negative white to blue
    const minVal = d3.min(expressionMatrix.flat());
    const maxVal = d3.max(expressionMatrix.flat());
    const extremeVal = Math.max(Math.abs(minVal), Math.abs(maxVal));

    const colorScale = d3
      .scaleLinear()
      .domain([-extremeVal, 0, extremeVal])
      .range(["blue", "white", "red"]);

    // X and Y scales (band scales for genes/cells)
    const xScale = d3
      .scaleBand()
      .domain(genes)
      .range([margin.right, viewWidth]);
    // .padding(0.05);

    const yScale = d3
      .scaleBand()
      .domain(orderedCellIds)
      .range([margin.top, viewHeight]);
    // .padding(0.05);

    // Draw heatmap rectangles
    svg
      .selectAll()
      .data(
        expressionMatrix.flatMap((row, i) =>
          row.map((value, j) => ({
            cell: orderedCellIds[i],
            gene: genes[j],
            value,
          }))
        )
      )
      .enter()
      .append("rect")
      .attr("x", (d) => xScale(d.gene))
      .attr("y", (d) => yScale(d.cell))
      .attr("width", xScale.bandwidth())
      .attr("height", yScale.bandwidth())
      .attr("fill", (d) => colorScale(d.value))
      // .attr("stroke", "#fff")
      // .attr("stroke-width", 0.5)
      .append("title")
      .text((d) => `${d.cell} - ${d.gene}: ${d.value}`);

    // Add X-axis (genes)
    svg
      .append("g")
      .attr("transform", `translate(0,${viewHeight})`)
      .call(d3.axisBottom(xScale))
      .selectAll("text")
      .attr("transform", "rotate(-45)")
      .style("text-anchor", "end");

    // Add Y-axis (cells)
    svg
      .append("g")
      .attr("class", "y-axis")
      .attr("transform", `translate(${margin.right},${0})`)
      .call(d3.axisLeft(yScale));

    // Add VERTICAL color legend
    const legendWidth = 20;
    const legendHeight = 150;
    const legend = svg
      .append("g")
      .attr("transform", `translate(${viewWidth + 10}, 20)`);

    // Gradient definition (3 stops: blue -> white -> red)
    const defs = svg.append("defs");
    const gradient = defs
      .append("linearGradient")
      .attr("id", "legend-gradient")
      .attr("x1", "0%")
      .attr("x2", "0%")
      .attr("y1", "100%") // Bottom (blue)
      .attr("y2", "0%"); // Top (red)

    // Explicit gradient stops
    gradient
      .selectAll("stop")
      .data([
        { offset: "0%", color: "blue" }, // Negative
        { offset: "50%", color: "white" }, // Zero
        { offset: "100%", color: "red" }, // Positive
      ])
      .enter()
      .append("stop")
      .attr("offset", (d) => d.offset)
      .attr("stop-color", (d) => d.color);

    // Draw legend rectangle
    legend
      .append("rect")
      .attr("width", legendWidth)
      .attr("height", legendHeight)
      .attr("fill", "url(#legend-gradient)");

    // Legend axis
    const legendScale = d3
      .scaleLinear()
      .domain([-extremeVal, 0, extremeVal])
      .range([legendHeight, legendHeight / 2, 0]);

    legend
      .append("g")
      .attr("transform", `translate(${legendWidth}, 0)`)
      .call(d3.axisRight(legendScale).ticks(3));

    // Legend title
    legend
      .append("text")
      .attr("x", legendWidth / 2)
      .attr("y", -8)
      .attr("text-anchor", "middle")
      .text("Value");

    // Add titles
    svg
      .append("text")
      .attr("x", viewWidth / 2)
      .attr("y", 13)
      .attr("text-anchor", "middle")
      .text("Gene Expression Heatmap");

    svg
      .append("text")
      .attr("x", viewWidth / 2)
      .attr("y", viewHeight + 50)
      .attr("text-anchor", "middle")
      .text("Genes");

    svg
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("x", -viewHeight / 2)
      .attr("y", 12)
      .attr("text-anchor", "middle")
      .text("Cells");

    // Add divider lines between cell types
    let currentType = cellTypes[0];
    sortedCells.forEach((cell, i) => {
      if (cell.type !== currentType) {
        svg
          .append("line")
          .attr("x1", 0)
          .attr("x2", viewWidth)
          .attr("y1", yScale(cell.cellId))
          .attr("y2", yScale(cell.cellId))
          .attr("stroke", "black")
          .attr("stroke-width", 1);
        currentType = cell.type;
      }
    });

    // Group cells by type and calculate middle positions
    const cellTypeGroups = {};
    sortedCells.forEach((cell, i) => {
      if (!cellTypeGroups[cell.type]) {
        cellTypeGroups[cell.type] = {
          startIndex: i,
          count: 1,
        };
      } else {
        cellTypeGroups[cell.type].count++;
      }
    });

    // Add labels for each cell type group
    Object.entries(cellTypeGroups).forEach(([type, group]) => {
      const middleIndex = group.startIndex + Math.floor(group.count / 2);
      const middleCellId = orderedCellIds[middleIndex];

      svg
        .append("text")
        .attr("x", viewWidth)
        .attr("y", yScale(middleCellId) + yScale.bandwidth() / 2)
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .text(type);
    });
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
    ></div>
  );
};
