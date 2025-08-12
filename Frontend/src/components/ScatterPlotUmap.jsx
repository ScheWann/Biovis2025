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
  setPseudotimeDataSets,
  setPseudotimeLoadingStates,
  setClusterColorMappings,
  hoveredTrajectory,
  coordinatesData,
  cellTypesData,
  setCellTypesData,
  selectedCellTypes,
  setSelectedCellTypes,
  cellTypeColors,
  setCellTypeColors
}) => {
  const containerRef = useRef();
  const svgRef = useRef();
  const pseudotimeDataSetsRef = useRef({});
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

  const fetchPseudotimeData = async (sampleId, cellIds) => {
    // Check if data for this adata_umap_title already exists
    if (pseudotimeDataSetsRef.current[adata_umap_title]) {
      return pseudotimeDataSetsRef.current[adata_umap_title];
    }

    // Parse parameters from adata_umap_title
    // Format: ${formattedName}_${sampleId}_${editNeighbors}_${editNPcas}_${editResolutions}
    let n_neighbors;
    let n_pcas;
    let resolutions;

    // Set loading state for this specific dataset
    setPseudotimeLoadingStates(prevStates => ({
      ...prevStates,
      [adata_umap_title]: true
    }));

    try {
      const parts = adata_umap_title.split('_');
      if (parts.length >= 5) {
        // Extract the last 3 parts as the parameters
        n_neighbors = parseInt(parts[parts.length - 3]) || 10;
        n_pcas = parseInt(parts[parts.length - 2]) || 30;
        resolutions = parseFloat(parts[parts.length - 1]) || 1;
      }
    } catch (error) {
      console.warn("Could not parse parameters from adata_umap_title, using defaults", error);
    }

    try {
      const res = await fetch("/api/get_pseudotime_data", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sample_id: sampleId,
          cell_ids: cellIds,
          adata_umap_title: adata_umap_title,
          n_neighbors: n_neighbors,
          n_pcas: n_pcas,
          resolutions: resolutions,
        }),
      });
      const data = await res.json();

      // Add the new data to the datasets object
      setPseudotimeDataSets(prevDataSets => {
        const newDataSets = {
          ...prevDataSets,
          [adata_umap_title]: data
        };
        pseudotimeDataSetsRef.current = newDataSets;
        return newDataSets;
      });

      // Clear loading state for this specific dataset
      setPseudotimeLoadingStates(prevStates => ({
        ...prevStates,
        [adata_umap_title]: false
      }));

      return data;
    } catch (err) {
      console.error("Failed to fetch pseudotime data", err);
      // Clear loading state for this specific dataset on error
      setPseudotimeLoadingStates(prevStates => ({
        ...prevStates,
        [adata_umap_title]: false
      }));
      throw err;
    }
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

  // Separate useEffect for cluster color mapping to prevent infinite loops
  useEffect(() => {
    if (!data || data.length === 0 || !setClusterColorMappings) return;

    // Color by cluster - sort clusters numerically for consistent ordering
    const clusters = Array.from(new Set(data.map(clusterAccessor))).sort((a, b) => {
      // Extract numeric part from cluster names (e.g., "Cluster 1" -> 1)
      const numA = parseInt(a.toString().replace(/\D/g, '')) || 0;
      const numB = parseInt(b.toString().replace(/\D/g, '')) || 0;
      return numA - numB;
    });
    const color = d3.scaleOrdinal().domain(clusters).range(COLORS);

    // Create cluster color mapping and pass to parent
    const colorMapping = {
      sample_id: sampleId,
      name: adata_umap_title,
      clusters: {}
    };

    clusters.forEach(cluster => {
      // Extract numeric part from cluster name (e.g., "Cluster 4" -> "4")
      const clusterNumber = cluster.toString().replace(/\D/g, '') || cluster;
      colorMapping.clusters[clusterNumber] = color(cluster);
    });

    setClusterColorMappings(prevMappings => {
      // Only update if the mapping has changed
      if (!prevMappings[adata_umap_title] ||
        JSON.stringify(prevMappings[adata_umap_title].clusters) !== JSON.stringify(colorMapping.clusters)) {
        const newMappings = { ...prevMappings };
        newMappings[adata_umap_title] = colorMapping;
        return newMappings;
      }

      return prevMappings;
    });
  }, [data, clusterAccessor, sampleId, adata_umap_title, setClusterColorMappings]);

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

    // Color by cluster - sort clusters numerically for consistent ordering
    const clusters = Array.from(new Set(data.map(clusterAccessor))).sort((a, b) => {
      // Extract numeric part from cluster names (e.g., "Cluster 1" -> 1)
      const numA = parseInt(a.toString().replace(/\D/g, '')) || 0;
      const numB = parseInt(b.toString().replace(/\D/g, '')) || 0;
      return numA - numB;
    });
    const color = d3.scaleOrdinal().domain(clusters).range(COLORS);

    // SVG
    const svg = d3
      .select(svgRef.current)
      .attr("width", width)
      .attr("height", height)
      .on("mouseleave", () => {
        // Clear hover state when mouse leaves the entire SVG area
        if (hoveredCluster && hoveredCluster.umapId === umapId) {
          setHoveredCluster(null);
        }
      });

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Group data by cluster and compute hulls
    const clusterGroups = d3.group(data, clusterAccessor);

    // Check if trajectory is being hovered for this UMAP
    const isTrajectoryHovered = hoveredTrajectory &&
      hoveredTrajectory.adata_umap_title === adata_umap_title &&
      hoveredTrajectory.path &&
      hoveredTrajectory.path.length > 1;

    clusterGroups.forEach((points, cluster) => {
      if (points.length < 3) return;

      const coords = points.map((d) => [
        xScale(xAccessor(d)),
        yScale(yAccessor(d)),
      ]);

      const hull = d3.polygonHull(coords);

      if (hull && hull.length >= 3 && !isTrajectoryHovered) {
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
          .attr("fill-opacity", (hoveredCluster?.cluster === cluster && hoveredCluster?.umapId === umapId) ? 0.4 : 0.2)
          .attr("stroke", color(cluster))
          .attr("stroke-width", (hoveredCluster?.cluster === cluster && hoveredCluster?.umapId === umapId) ? 3.5 : 2.5)
          .attr(
            "stroke-opacity",
            (hoveredCluster?.cluster === cluster && hoveredCluster?.umapId === umapId) ? 0.8 : 0.3
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
            // Stop event propagation to prevent conflicts
            event.stopPropagation();

            d3.select(event.currentTarget)
              .attr("fill-opacity", 0.4)
              .attr("stroke-width", 3.5)
              .attr("stroke-opacity", 0.8);

            const cellIds = points
              .map((d) => d.id || d.cell_id)
              .filter(Boolean);

            if (!hoveredCluster || hoveredCluster.cluster !== cluster || hoveredCluster.umapId !== umapId) {
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
            // Stop event propagation to prevent conflicts
            event.stopPropagation();

            d3.select(event.currentTarget)
              .attr("fill-opacity", 0.2)
              .attr("stroke-width", 2.5)
              .attr("stroke-opacity", 0.3);

            // Use a small delay to prevent flickering when moving between related elements
            setTimeout(() => {
              if (hoveredCluster && hoveredCluster.cluster === cluster && hoveredCluster.umapId === umapId) {
                setHoveredCluster(null);
              }
            }, 10);
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
        if (!hoveredCluster || hoveredCluster.umapId !== umapId) return 0.5;
        return hoveredCluster.cluster === clusterAccessor(d) ? 0.8 : 0.05;
      })
      .attr("stroke", (d) => {
        if (!hoveredCluster || hoveredCluster.umapId !== umapId) return "none";
        return hoveredCluster.cluster === clusterAccessor(d) ? "#fff" : "none";
      })
      .attr("stroke-width", (d) => {
        if (!hoveredCluster || hoveredCluster.umapId !== umapId) return 0;
        return hoveredCluster.cluster === clusterAccessor(d) ? 1 : 0;
      })
      .style("cursor", "pointer")
      .on("mouseenter", (event, d) => {
        // Stop event propagation to prevent conflicts with hull events
        event.stopPropagation();

        // Highlight points of the same cluster
        const cluster = clusterAccessor(d);
        const clusterPoints = data.filter(
          (point) => clusterAccessor(point) === cluster
        );
        const cellIds = clusterPoints
          .map((p) => p.id || p.cell_id)
          .filter(Boolean);

        if (!hoveredCluster || hoveredCluster.cluster !== cluster || hoveredCluster.umapId !== umapId) {
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
        // Stop event propagation to prevent conflicts
        event.stopPropagation();

        const cluster = clusterAccessor(d);
        // Use a small delay to prevent flickering when moving between points in the same cluster
        setTimeout(() => {
          if (hoveredCluster && hoveredCluster.cluster === cluster && hoveredCluster.umapId === umapId) {
            setHoveredCluster(null);
          }
        }, 10);
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
      const legendGroup = legend
        .append("g")
        .style("cursor", "pointer")
        .on("mouseenter", (event) => {
          // Stop event propagation to prevent conflicts
          event.stopPropagation();

          // Find all points in this cluster
          const clusterPoints = data.filter(
            (point) => clusterAccessor(point) === cl
          );
          const cellIds = clusterPoints
            .map((p) => p.id || p.cell_id)
            .filter(Boolean);

          if (!hoveredCluster || hoveredCluster.cluster !== cl || hoveredCluster.umapId !== umapId) {
            setHoveredCluster({
              cluster: cl,
              cellIds: cellIds,
              points: clusterPoints,
              umapId: umapId,
              sampleId: sampleId,
            });
          }
        })
        .on("mouseleave", (event) => {
          // Stop event propagation to prevent conflicts
          event.stopPropagation();

          // Use a small delay to prevent flickering
          setTimeout(() => {
            if (hoveredCluster && hoveredCluster.cluster === cl && hoveredCluster.umapId === umapId) {
              setHoveredCluster(null);
            }
          }, 10);
        });

      legendGroup
        .append("circle")
        .attr("cx", 0)
        .attr("cy", i * 15)
        .attr("r", 4)
        .attr("fill", color(cl))
        .attr(
          "opacity",
          (!hoveredCluster || hoveredCluster.umapId !== umapId) || hoveredCluster.cluster === cl ? 1 : 0.3
        );
      legendGroup
        .append("text")
        .attr("x", 8)
        .attr("y", i * 15 + 3)
        .text(cl)
        .attr("font-size", 9)
        .attr(
          "opacity",
          (!hoveredCluster || hoveredCluster.umapId !== umapId) || hoveredCluster.cluster === cl ? 1 : 0.3
        );
    });

    // Trajectory visualization - draw stars and arrows when trajectory is hovered
    if (hoveredTrajectory &&
      hoveredTrajectory.adata_umap_title === adata_umap_title &&
      hoveredTrajectory.path &&
      hoveredTrajectory.path.length > 1) {

      // Calculate cluster centers
      const clusterCenters = new Map();
      clusters.forEach(cluster => {
        const clusterPoints = data.filter(d => clusterAccessor(d) === cluster);

        if (clusterPoints.length > 0) {
          const centerX = d3.mean(clusterPoints, xAccessor);
          const centerY = d3.mean(clusterPoints, yAccessor);
          clusterCenters.set(cluster, {
            x: xScale(centerX),
            y: yScale(centerY)
          });
        }
      });

      // Create trajectory group
      const trajectoryGroup = g.append("g").attr("class", "trajectory-visualization");

      // Draw stars at cluster centers for clusters in the trajectory path
      hoveredTrajectory.path.forEach((clusterId, index) => {
        // Convert numeric cluster ID to "Cluster X" format to match UMAP data
        const clusterName = `Cluster ${clusterId}`;
        const center = clusterCenters.get(clusterName);

        if (center) {
          // Draw star shape
          const starPoints = [];
          const outerRadius = 12;
          const innerRadius = 6;
          const numPoints = 2;

          for (let i = 0; i < numPoints * 2; i++) {
            const angle = (i * Math.PI) / numPoints;
            const radius = i % 2 === 0 ? outerRadius : innerRadius;
            const x = center.x + Math.cos(angle - Math.PI / 2) * radius;
            const y = center.y + Math.sin(angle - Math.PI / 2) * radius;
            starPoints.push([x, y]);
          }

          trajectoryGroup
            .append("polygon")
            .attr("points", starPoints.map(p => p.join(",")).join(" "))
            .attr("fill", "#4F46E5")
            .attr("stroke", "#312E81")
            .attr("stroke-width", 1)
            .attr("opacity", 0.95);

          // Add cluster label on the star
          trajectoryGroup
            .append("text")
            .attr("x", center.x)
            .attr("y", center.y + 4)
            .attr("text-anchor", "middle")
            .attr("font-size", "10px")
            .attr("font-weight", "bold")
            .attr("fill", "#FFFFFF")
            .attr("stroke", "#312E81")
            .attr("stroke-width", "0.5")
            .style("paint-order", "stroke")
            .text(hoveredTrajectory.path[index]);
        }
      });

      // Draw arrows between consecutive clusters in the trajectory
      for (let i = 0; i < hoveredTrajectory.path.length - 1; i++) {
        const fromClusterId = hoveredTrajectory.path[i];
        const toClusterId = hoveredTrajectory.path[i + 1];
        // Convert numeric cluster IDs to "Cluster X" format
        const fromClusterName = `Cluster ${fromClusterId}`;
        const toClusterName = `Cluster ${toClusterId}`;
        const fromCenter = clusterCenters.get(fromClusterName);
        const toCenter = clusterCenters.get(toClusterName);

        if (fromCenter && toCenter) {
          // Calculate arrow direction
          const dx = toCenter.x - fromCenter.x;
          const dy = toCenter.y - fromCenter.y;
          const length = Math.sqrt(dx * dx + dy * dy);

          if (length > 24) { // Only draw arrow if clusters are far enough apart
            const unitX = dx / length;
            const unitY = dy / length;

            // Adjust start and end points to avoid overlapping with stars
            const startX = fromCenter.x + unitX * 15;
            const startY = fromCenter.y + unitY * 15;
            const endX = toCenter.x - unitX * 15;
            const endY = toCenter.y - unitY * 15;

            // Draw arrow line
            trajectoryGroup
              .append("line")
              .attr("x1", startX)
              .attr("y1", startY)
              .attr("x2", endX)
              .attr("y2", endY)
              .attr("stroke", "#6366F1")
              .attr("stroke-width", 3)
              .attr("opacity", 0.9)
              .attr("marker-end", "url(#arrowhead)");

            // Create arrowhead marker if it doesn't exist
            let defs = svg.select("defs");
            if (defs.empty()) {
              defs = svg.append("defs");
            }

            if (defs.select("#arrowhead").empty()) {
              defs.append("marker")
                .attr("id", "arrowhead")
                .attr("viewBox", "0 -2 4 4")
                .attr("refX", 3)
                .attr("refY", 0)
                .attr("markerWidth", 3)
                .attr("markerHeight", 3)
                .attr("orient", "auto")
                .append("path")
                .attr("d", "M0,-2L4,0L0,2")
                .attr("fill", "#4F46E5")
                .style("filter", "drop-shadow(1px 1px 2px rgba(79, 70, 229, 0.4))");
            }
          }
        }
      }
    }
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
    hoveredTrajectory,
    adata_umap_title,
    sampleId,
  ]);

  return (
    <div
      ref={containerRef}
      style={{ width: "100%", height: "100%" }}
      onMouseLeave={() => {
        // Clear hover state when mouse leaves the entire container
        if (hoveredCluster && hoveredCluster.umapId === umapId) {
          setHoveredCluster(null);
        }
      }}
    >
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
        coordinatesData={coordinatesData}
        cellTypesData={cellTypesData}
        setCellTypesData={setCellTypesData}
        selectedCellTypes={selectedCellTypes}
        setSelectedCellTypes={setSelectedCellTypes}
        cellTypeColors={cellTypeColors}
        setCellTypeColors={setCellTypeColors}
        sampleId={sampleId}
      />
    </div>
  );
};
