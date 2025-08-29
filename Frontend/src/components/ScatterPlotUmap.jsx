import React, { useRef, useEffect, useState, useMemo } from "react";
import * as d3 from "d3";
import { GOAnalysisWindow } from "./GOAnalysisWindow";
import { UmapSettingsPopup } from "./UmapSettingsPopup";
import { COLOR_BREWER2_PALETTE_4 } from "./Utils";
import { Select, Button, Spin } from "antd";

const COLORS = COLOR_BREWER2_PALETTE_4;

export const ScatterplotUmap = ({
  data,
  pointSize = 5,
  clusterAccessor = (d) => d.cluster,
  xAccessor = (d) => d.x,
  yAccessor = (d) => d.y,
  title = "UMAP",
  adata_umap_title,
  margin = { top: 25, right: 10, bottom: 25, left: 25 },
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
  setSelectedCellTypes,
  setCellTypeColors,
  pseudotimeDataSets,
  onUmapDataUpdate,
  onUmapLoadingStart
}) => {
  const containerRef = useRef();
  const svgRef = useRef();

  const [dimensions, setDimensions] = useState({ width: 400, height: 200 });
  const [clickPosition, setClickPosition] = useState({ x: 0, y: 0 });
  const [currentCellIds, setCurrentCellIds] = useState([]);

  // Local hover state management to prevent flickering
  const [localHoveredCluster, setLocalHoveredCluster] = useState(null);
  const hoverTimeoutRef = useRef(null);

  // Local GO analysis state for this ScatterPlotUmap instance
  const [GOAnalysisData, setGOAnalysisData] = useState(null);
  const [GOAnalysisLoading, setGOAnalysisLoading] = useState(false);
  const [GOAnalysisVisible, setGOAnalysisVisible] = useState(false);

  // UMAP settings popup state
  const [umapSettingsVisible, setUmapSettingsVisible] = useState(false);
  const [umapSettingsPosition, setUmapSettingsPosition] = useState({ x: 0, y: 0 });

  // State to store cluster info for GO analysis window
  const [currentClusterInfo, setCurrentClusterInfo] = useState(null);

  // State for direct Slingshot analysis
  const [selectedStartCluster, setSelectedStartCluster] = useState(null);
  const [directSlingshotLoading, setDirectSlingshotLoading] = useState(false);

  // Get available clusters for the selection box
  const availableClusters = useMemo(() => {
    if (!data || data.length === 0) return [];
    
    const clusters = Array.from(new Set(data.map(clusterAccessor))).sort((a, b) => {
      const numA = parseInt(a.toString().replace(/\D/g, '')) || 0;
      const numB = parseInt(b.toString().replace(/\D/g, '')) || 0;
      return numA - numB;
    });
    
    return clusters.map(cluster => ({
      value: cluster,
      label: cluster
    }));
  }, [data, clusterAccessor]);

  // Function to handle direct Slingshot analysis
  const handleDirectSlingshot = async () => {
    if (!selectedStartCluster || !adata_umap_title || !setPseudotimeDataSets || !setPseudotimeLoadingStates) {
      console.warn("Missing required data for direct Slingshot analysis");
      return;
    }

    // Extract cluster number from selected cluster
    const clusterNumber = selectedStartCluster.toString().replace(/\D/g, '');
    if (!clusterNumber) {
      console.warn("No cluster number found in selected cluster");
      return;
    }

    // Parse parameters from adata_umap_title
    let n_neighbors = 15;
    let n_pcas = 30;
    let resolutions = 1;

    try {
      const parts = adata_umap_title.split('_');
      if (parts.length >= 5) {
        // Extract the last 3 parts as the parameters
        n_neighbors = parseInt(parts[parts.length - 3]) || 15;
        n_pcas = parseInt(parts[parts.length - 2]) || 30;
        resolutions = parseFloat(parts[parts.length - 1]) || 1;
      }
    } catch (error) {
      console.warn("Could not parse parameters from adata_umap_title, using defaults", error);
    }

    // Create a unique key for direct slingshot data
    const directSlingshotKey = `${adata_umap_title}_direct_slingshot`;

    // Check if data for this direct slingshot key already exists in cache
    if (pseudotimeDataSets && pseudotimeDataSets[directSlingshotKey]) {
      setPseudotimeLoadingStates(prevStates => ({
        ...prevStates,
        [directSlingshotKey]: true
      }));
      setPseudotimeDataSets(prevDataSets => ({ ...prevDataSets }));
      setTimeout(() => {
        setPseudotimeLoadingStates(prevStates => ({
          ...prevStates,
          [directSlingshotKey]: false
        }));
      }, 0);
      return;
    }

    // Set loading state for this specific dataset
    setDirectSlingshotLoading(true);
    setPseudotimeLoadingStates(prevStates => ({
      ...prevStates,
      [directSlingshotKey]: true
    }));

    try {
      const res = await fetch("/api/get_direct_slingshot_data", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sample_id: sampleId,
          cell_ids: data.map(d => d.id || d.cell_id).filter(Boolean),
          adata_umap_title: adata_umap_title,
          start_cluster: clusterNumber,
          n_neighbors: n_neighbors,
          n_pcas: n_pcas,
          resolutions: resolutions,
        }),
      });
      const responseData = await res.json();

      if (responseData.error) {
        console.error("Direct Slingshot analysis failed:", responseData.error);
        return;
      }

      // Add the new data to the datasets object with the unique key
      setPseudotimeDataSets(prevDataSets => {
        const newDataSets = {
          ...prevDataSets,
          [directSlingshotKey]: responseData
        };
        return newDataSets;
      });
    } catch (err) {
      console.error("Failed to fetch direct Slingshot data", err);
    } finally {
      // Clear loading state for this specific dataset
      setDirectSlingshotLoading(false);
      setPseudotimeLoadingStates(prevStates => ({
        ...prevStates,
        [directSlingshotKey]: false
      }));
    }
  };

  const fetchGOAnalysisData = (sampleId, cluster, adata_umap_title) => {
    setGOAnalysisLoading(true);
    setGOAnalysisVisible(true);
    setGOAnalysisData(null);

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

  // Handler for UMAP settings update
  const handleUmapSettingsUpdate = (newData, newAdataUmapTitle, newSettings, newName) => {
    // Update the data prop by calling a callback from parent
    if (onUmapDataUpdate) {
      onUmapDataUpdate(newData, newAdataUmapTitle, newSettings, newName);
    }
  };

  // Helper function to manage hover state changes
  const updateHoverState = (cluster, points, cellIds) => {
    // Clear any existing timeout
    if (hoverTimeoutRef.current) {
      clearTimeout(hoverTimeoutRef.current);
      hoverTimeoutRef.current = null;
    }

    // If hovering the same cluster, do nothing
    if (localHoveredCluster?.cluster === cluster && localHoveredCluster?.umapId === umapId) {
      return;
    }

    // Update local state immediately
    const newHoverState = cluster ? {
      cluster: cluster,
      cellIds: cellIds,
      points: points,
      umapId: umapId,
      sampleId: sampleId,
    } : null;

    setLocalHoveredCluster(newHoverState);
    setHoveredCluster(newHoverState);
  };

  // Helper function to clear hover state with delay
  const clearHoverState = () => {
    if (hoverTimeoutRef.current) {
      clearTimeout(hoverTimeoutRef.current);
    }
    
    hoverTimeoutRef.current = setTimeout(() => {
      setLocalHoveredCluster(null);
      setHoveredCluster(null);
      hoverTimeoutRef.current = null;
    }, 50); // Small delay to prevent flickering
  };

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
        clearHoverState();
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
          .attr("data-cluster", cluster)
          .attr("fill", color(cluster))
          .attr("fill-opacity", (localHoveredCluster?.cluster === cluster && localHoveredCluster?.umapId === umapId) ? 0.4 : 0.2)
          .attr("stroke", color(cluster))
          .attr("stroke-width", (localHoveredCluster?.cluster === cluster && localHoveredCluster?.umapId === umapId) ? 3.5 : 2.5)
          .attr(
            "stroke-opacity",
            (localHoveredCluster?.cluster === cluster && localHoveredCluster?.umapId === umapId) ? 0.8 : 0.3
          )
          .style("cursor", "pointer")
          .style("pointer-events", "visibleFill")
          .on("click", function (event) {
            setClickPosition({ x: event.clientX, y: event.clientY });
            const cellIds = points
              .map((d) => d.id || d.cell_id)
              .filter(Boolean);
            setCurrentCellIds(cellIds); 

            const clusterNumber = cluster.split(" ")[1];
            setCurrentClusterInfo({
              cluster_name: cluster,
              cluster_number: clusterNumber
            });

            fetchGOAnalysisData(sampleId, cluster, adata_umap_title);
          })
          .on("mouseenter", (event) => {
            event.stopPropagation();

            const cellIds = points
              .map((d) => d.id || d.cell_id)
              .filter(Boolean);

            updateHoverState(cluster, points, cellIds);
          })
          .on("mouseleave", (event) => {
            event.stopPropagation();

            clearHoverState();
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
        if (!localHoveredCluster || localHoveredCluster.umapId !== umapId) return 0.5;
        return localHoveredCluster.cluster === clusterAccessor(d) ? 0.8 : 0.05;
      })
      .attr("stroke", (d) => {
        if (!localHoveredCluster || localHoveredCluster.umapId !== umapId) return "none";
        return localHoveredCluster.cluster === clusterAccessor(d) ? "#fff" : "none";
      })
      .attr("stroke-width", (d) => {
        if (!localHoveredCluster || localHoveredCluster.umapId !== umapId) return 0;
        return localHoveredCluster.cluster === clusterAccessor(d) ? 1 : 0;
      })
      .style("cursor", "pointer")
      .on("click", function (event, d) {
        event.stopPropagation();
        
        setClickPosition({ x: event.clientX, y: event.clientY });
        
        const cluster = clusterAccessor(d);
        const clusterPoints = data.filter(
          (point) => clusterAccessor(point) === cluster
        );
        const cellIds = clusterPoints
          .map((p) => p.id || p.cell_id)
          .filter(Boolean);
        
        setCurrentCellIds(cellIds); 

        const clusterNumber = cluster.split(" ")[1];
        setCurrentClusterInfo({
          cluster_name: cluster,
          cluster_number: clusterNumber
        });

        fetchGOAnalysisData(sampleId, cluster, adata_umap_title);
      })
      .on("mouseenter", (event, d) => {
        event.stopPropagation();

        const cluster = clusterAccessor(d);
        const clusterPoints = data.filter(
          (point) => clusterAccessor(point) === cluster
        );
        const cellIds = clusterPoints
          .map((p) => p.id || p.cell_id)
          .filter(Boolean);

        updateHoverState(cluster, clusterPoints, cellIds);
      })
      .on("mouseleave", (event, d) => {
        event.stopPropagation();

        clearHoverState();
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
      .on("click", (event) => {
        // Open UMAP settings popup
        setUmapSettingsPosition({ x: event.clientX, y: event.clientY });
        setUmapSettingsVisible(true);
      });

    // Legend
    const legend = g
      .append("g")
      .attr(
        "transform",
        `translate(${width - margin.right - 70}, ${margin.top - 20})`
      )
      .style("z-index", 10);

    clusters.forEach((cl, i) => {
      const legendGroup = legend
        .append("g")
        .style("cursor", "pointer")
        .on("mouseenter", (event) => {
          event.stopPropagation();
          const clusterPoints = data.filter(point => clusterAccessor(point) === cl);
          const cellIds = clusterPoints
            .map((p) => p.id || p.cell_id)
            .filter(Boolean);

          updateHoverState(cl, clusterPoints, cellIds);
        })
        .on("mouseleave", (event) => {
          event.stopPropagation();
          clearHoverState();
        });

      legendGroup
        .append("circle")
        .attr("cx", 0)
        .attr("cy", i * 15)
        .attr("r", 4)
        .attr("fill", color(cl))
        .attr(
          "opacity",
          (!localHoveredCluster || localHoveredCluster.umapId !== umapId) || localHoveredCluster.cluster === cl ? 1 : 0.3
        );
      legendGroup
        .append("text")
        .attr("x", 8)
        .attr("y", i * 15 + 3)
        .text(cl)
        .attr("font-size", 9)
        .attr(
          "opacity",
          (!localHoveredCluster || localHoveredCluster.umapId !== umapId) || localHoveredCluster.cluster === cl ? 1 : 0.3
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
    localHoveredCluster,
    hoveredTrajectory,
    adata_umap_title,
    sampleId,
  ]);

  return (
    <div
      ref={containerRef}
      style={{ width: "100%", height: "100%" }}
      onMouseLeave={() => {
        clearHoverState();
      }}
    >
      {/* Direct Slingshot Controls - positioned over the SVG */}
      <div style={{
        position: 'absolute',
        top: '5px',
        left: '25px',
        zIndex: 1000,
        display: 'flex',
        alignItems: 'center',
        gap: '8px',
        backgroundColor: 'rgba(255, 255, 255, 0.9)',
        padding: '4px 8px',
        borderRadius: '4px',
        boxShadow: '0 1px 3px rgba(0,0,0,0.1)'
      }}>
        <Select
          placeholder="Select start cluster"
          value={selectedStartCluster}
          onChange={setSelectedStartCluster}
          style={{ width: '140px' }}
          size="small"
          options={availableClusters}
          allowClear
        />
        <Button
          onClick={handleDirectSlingshot}
          disabled={!selectedStartCluster || directSlingshotLoading}
          loading={directSlingshotLoading}
          type="primary"
          size="small"
        >
          {directSlingshotLoading ? 'Running...' : 'Run Slingshot'}
        </Button>
      </div>

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
        setSelectedCellTypes={setSelectedCellTypes}
        setCellTypeColors={setCellTypeColors}
        sampleId={sampleId}
        clusterInfo={currentClusterInfo}
        adata_umap_title={adata_umap_title}
        setPseudotimeDataSets={setPseudotimeDataSets}
        setPseudotimeLoadingStates={setPseudotimeLoadingStates}
        pseudotimeDataSets={pseudotimeDataSets}
      />
      <UmapSettingsPopup
        visible={umapSettingsVisible}
        setVisible={setUmapSettingsVisible}
        position={umapSettingsPosition}
        onUpdateSettings={handleUmapSettingsUpdate}
        onLoadingStart={onUmapLoadingStart}
        sampleId={sampleId}
        cellIds={data.map(d => d.id || d.cell_id).filter(Boolean)}
        adata_umap_title={adata_umap_title}
        currentTitle={title}
      />
    </div>
  );
};
