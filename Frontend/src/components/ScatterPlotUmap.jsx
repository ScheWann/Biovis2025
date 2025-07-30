import React, { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import { GOAnalysisWindow } from "./GOAnalysisWindow";

const COLORS = d3.schemeCategory10;

export const ScatterplotUmap = ({
  data, // [{x, y, cluster}]
  pointSize = 5,
  clusterAccessor = (d) => d.cluster,
  xAccessor = (d) => d.x,
  yAccessor = (d) => d.y,
  title = "UMAP",
  adata_umap_title,
  margin = { top: 25, right: 10, bottom: 25, left: 25 },
  hoveredCluster,
  setHoveredCluster,
  umapId,
  sampleId,
  setCellName,
}) => {
  const containerRef = useRef();
  const svgRef = useRef();
  const [dimensions, setDimensions] = useState({ width: 400, height: 200 });
  const [clickPosition, setClickPosition] = useState({ x: 0, y: 0 });
  const [currentCellIds, setCurrentCellIds] = useState([]);
  
  // Local GO analysis state for this ScatterPlotUmap instance
  const [GOAnalysisData, setGOAnalysisData] = useState(null);
  const [GOAnalysisLoading, setGOAnalysisLoading] = useState(false);
  const [GOAnalysisVisible, setGOAnalysisVisible] = useState(false);

  const fetchGOAnalysisData = (sampleId, cluster, adata_umap_title) => {
    setGOAnalysisLoading(true);
    setGOAnalysisVisible(true);
    setGOAnalysisData(null); // Clear previous data
    
    const cluster_id = cluster.split(" ")[1]
    fetch("/api/get_go_analysis", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sample_id: sampleId, cluster_id: cluster_id, adata_umap_title: adata_umap_title }),
    })
      .then((res) => res.json())
      .then((data) => {
        setGOAnalysisData(data);
        setGOAnalysisLoading(false);
      })
      .catch((err) => {
        console.error("Failed to fetch GO analysis data", err);
        setGOAnalysisLoading(false);
        setGOAnalysisData(null);
      });
  };

  const fetchPseudotimeData = ( sampleId, cellIds = null ) => {
    fetch("/api/get_pseudotime_data", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ 
        sample_id: sampleId, 
        cell_ids: cellIds 
      }),
    })
      .then((res) => res.json())
      .then((data) => {
        console.log(data, "pseudotime data");
      })
      .catch((err) => {
        console.error("Failed to fetch pseudotime data", err);
      });
  }

  // Detect container size changes
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const resizeObserver = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        setDimensions({
          width: Math.max(width, 200),
          height: Math.max(height, 150),
        });
      }
    });

    resizeObserver.observe(container);
    return () => resizeObserver.disconnect();
  }, []);

  useEffect(() => {
    if (!data || data.length === 0) return;

    d3.select(svgRef.current).selectAll("*").remove();

    const { width, height } = dimensions;
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
      .select(svgRef.current)
      .attr("width", width)
      .attr("height", height);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Group data by cluster and compute hulls
    const clusterGroups = d3.group(data, clusterAccessor);
    clusterGroups.forEach((points, cluster) => {
      if (points.length < 3) return;

      const coords = points.map((d) => [
        xScale(xAccessor(d)),
        yScale(yAccessor(d)),
      ]);

      const hull = d3.polygonHull(coords);

      if (hull && hull.length >= 3) {
        // Create path for the hull with smoother curves
        const line = d3
          .line()
          .curve(d3.curveCardinalClosed.tension(0))
          .x((d) => d[0])
          .y((d) => d[1]);

        g.append("path")
          .datum(hull)
          .attr("d", line)
          .attr("fill", color(cluster))
          .attr("fill-opacity", hoveredCluster?.cluster === cluster ? 0.4 : 0.2)
          .attr("stroke", color(cluster))
          .attr("stroke-width", hoveredCluster?.cluster === cluster ? 3.5 : 2.5)
          .attr(
            "stroke-opacity",
            hoveredCluster?.cluster === cluster ? 0.8 : 0.6
          )
          .style("cursor", "pointer")
          .style("pointer-events", "visibleFill")
          .on("click", function (event) {
            // Capture click position relative to viewport
            setClickPosition({ x: event.clientX, y: event.clientY });
            const cellIds = points
              .map((d) => d.id || d.cell_id)
              .filter(Boolean);
            setCurrentCellIds(cellIds); // Store current cellIds
            fetchGOAnalysisData(sampleId, cluster, adata_umap_title);
          })
          .on("mouseenter", (event) => {
            d3.select(event.currentTarget)
              .attr("fill-opacity", 0.4)
              .attr("stroke-width", 3.5)
              .attr("stroke-opacity", 0.8);

            const cellIds = points
              .map((d) => d.id || d.cell_id)
              .filter(Boolean);

            if (!hoveredCluster || hoveredCluster.cluster !== cluster) {
              setHoveredCluster({
                cluster: cluster,
                cellIds: cellIds,
                points: points,
                umapId: umapId,
                sampleId: sampleId,
              });
            }
          })
          .on("mouseleave", (event) => {
              d3.select(event.currentTarget)
                .attr("fill-opacity", 0.2)
                .attr("stroke-width", 2.5)
                .attr("stroke-opacity", 0.6);

              if (hoveredCluster && hoveredCluster.cluster === cluster) {
                setHoveredCluster(null);
              }
          });
      }
    });

    // Points
    g.selectAll("circle")
      .data(data)
      .enter()
      .append("circle")
      .attr("cx", (d) => xScale(xAccessor(d)))
      .attr("cy", (d) => yScale(yAccessor(d)))
      .attr("r", pointSize)
      .attr("fill", (d) => color(clusterAccessor(d)))
      .attr("opacity", (d) => {
        if (!hoveredCluster) return 0.5;
        return hoveredCluster.cluster === clusterAccessor(d) ? 0.8 : 0.05;
      })
      .attr("stroke", (d) => {
        if (!hoveredCluster) return "none";
        return hoveredCluster.cluster === clusterAccessor(d) ? "#fff" : "none";
      })
      .attr("stroke-width", (d) => {
        if (!hoveredCluster) return 0;
        return hoveredCluster.cluster === clusterAccessor(d) ? 1 : 0;
      })
      .style("cursor", "pointer")
      .on("mouseenter", (event, d) => {
        // Highlight points of the same cluster
        const cluster = clusterAccessor(d);
        const clusterPoints = data.filter(
          (point) => clusterAccessor(point) === cluster
        );
        const cellIds = clusterPoints
          .map((p) => p.id || p.cell_id)
          .filter(Boolean);

        if (!hoveredCluster || hoveredCluster.cluster !== cluster) {
          setHoveredCluster({
            cluster: cluster,
            cellIds: cellIds,
            points: clusterPoints,
            umapId: umapId,
            sampleId: sampleId,
          });
        }
      })
      .on("mouseleave", (event, d) => {
        const cluster = clusterAccessor(d);
        if (hoveredCluster && hoveredCluster.cluster === cluster) {
          setHoveredCluster(null);
        }
      });

    // Axes
    g.append("g")
      .attr("transform", `translate(0,${innerHeight})`)
      .call(d3.axisBottom(xScale).ticks(5));

    g.append("g").call(d3.axisLeft(yScale).ticks(5));

    // Title
    svg
      .append("text")
      .attr("x", margin.left)
      .attr("y", margin.top - 5)
      .attr("font-size", 12)
      .attr("font-weight", 600)
      .text(title)
      .style("cursor", "pointer")
      .on("click", () => {
        // Pseudotime analysis for current cells
        // If no specific cluster selected, use all cells
        const cellIds = currentCellIds.length > 0 
          ? currentCellIds 
          : data.map(d => d.id || d.cell_id).filter(Boolean);
        console.log(cellIds, "cellIds for pseudotime analysis");
        fetchPseudotimeData(sampleId, cellIds);
      });

    // Legend
    const legend = g
      .append("g")
      .attr(
        "transform",
        `translate(${width - margin.right - 70}, ${margin.top - 20})`
      );
    clusters.forEach((cl, i) => {
      legend
        .append("circle")
        .attr("cx", 0)
        .attr("cy", i * 15)
        .attr("r", 4)
        .attr("fill", color(cl))
        .attr(
          "opacity",
          !hoveredCluster || hoveredCluster.cluster === cl ? 1 : 0.3
        );
      legend
        .append("text")
        .attr("x", 8)
        .attr("y", i * 15 + 3)
        .text(cl)
        .attr("font-size", 9)
        .attr(
          "opacity",
          !hoveredCluster || hoveredCluster.cluster === cl ? 1 : 0.3
        );
    });
  }, [
    data,
    dimensions,
    pointSize,
    clusterAccessor,
    xAccessor,
    yAccessor,
    title,
    margin,
    hoveredCluster,
  ]);

  return (
    <div ref={containerRef} style={{ width: "100%", height: "100%" }}>
      <svg ref={svgRef}></svg>
      <GOAnalysisWindow
        visible={GOAnalysisVisible}
        setVisible={setGOAnalysisVisible}
        loading={GOAnalysisLoading}
        data={GOAnalysisData}
        position={clickPosition}
        title="Gene Ontology Analysis"
        setCellName={setCellName}
        cellIds={currentCellIds}
      />
    </div>
  );
};
