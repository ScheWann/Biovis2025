import React, { useRef, useEffect, useState, useCallback } from 'react';
import * as d3 from 'd3';
import { Empty, Spin, Checkbox } from 'antd';

// Custom hook for debounced trajectory hovering
const useDebounceTrajectoryHover = (setHoveredTrajectory, delay = 100) => {
    const timeoutRef = useRef();
    const lastHoverRef = useRef(null);

    const debouncedSetHover = useCallback((trajectoryData) => {
        // Check if the trajectory data is the same to avoid unnecessary updates
        const currentKey = trajectoryData ? `${trajectoryData.adata_umap_title}_${JSON.stringify(trajectoryData.path)}_${trajectoryData.trajectoryIndex}` : null;
        const lastKey = lastHoverRef.current ? `${lastHoverRef.current.adata_umap_title}_${JSON.stringify(lastHoverRef.current.path)}_${lastHoverRef.current.trajectoryIndex}` : null;

        if (currentKey === lastKey) {
            return; // Same trajectory, no need to update
        }

        clearTimeout(timeoutRef.current);
        timeoutRef.current = setTimeout(() => {
            setHoveredTrajectory(trajectoryData);
            lastHoverRef.current = trajectoryData;
        }, delay);
    }, [setHoveredTrajectory, delay]);

    const clearHover = useCallback(() => {
        clearTimeout(timeoutRef.current);
        timeoutRef.current = setTimeout(() => {
            setHoveredTrajectory(null);
            lastHoverRef.current = null;
        }, 50); // Shorter delay for clearing
    }, [setHoveredTrajectory]);

    // Cleanup timeout on unmount
    useEffect(() => {
        return () => {
            if (timeoutRef.current) {
                clearTimeout(timeoutRef.current);
            }
        };
    }, []);

    return { debouncedSetHover, clearHover };
};

const PseudotimeGlyph = ({
    adata_umap_title,
    pseudotimeData,
    pseudotimeLoading,
    isSelected = false,
    onSelectionChange,
    geneExpressionData = null,
    clusterColors = null,
    setHoveredTrajectory,
    trajectoryIndex,
    // sampleId,
    source_title
}) => {
    const containerRef = useRef();
    const svgRef = useRef(null);
    const [dimensions, setDimensions] = useState({ width: 300, height: 300 });
    const [selectedTrajectory, setSelectedTrajectory] = useState(null);

    // Use the debounced hover hook
    const { debouncedSetHover, clearHover } = useDebounceTrajectoryHover(setHoveredTrajectory);

    // Generate a unique ID for this component instance
    const [componentId] = useState(() => `pseudotime-glyph-${Math.random().toString(36).substr(2, 9)}`);

    // Detect container size changes
    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        const updateDimensions = () => {
            const rect = container.getBoundingClientRect();
            setDimensions({
                width: Math.max(rect.width || 400, 300),
                height: Math.max(rect.height || 400, 300),
            });
        };

        // Initial measurement
        updateDimensions();

        const resizeObserver = new ResizeObserver((entries) => {
            for (const entry of entries) {
                const { width, height } = entry.contentRect;
                setDimensions({
                    width: Math.max(width, 300),
                    height: Math.max(height, 300),
                });
            }
        });

        resizeObserver.observe(container);
        return () => resizeObserver.disconnect();
    }, []);

    useEffect(() => {
        if (pseudotimeData && pseudotimeData.length > 0 && dimensions.width > 0 && dimensions.height > 0) {
            createGlyph(pseudotimeData);
        }
    }, [pseudotimeData, dimensions, geneExpressionData, clusterColors]);

    // Cleanup tooltip on unmount
    useEffect(() => {
        return () => {
            d3.select("body").selectAll(`.pseudotime-tooltip-${componentId}`).remove();
        };
    }, [componentId]);

    // Helper function to position tooltip within viewport
    const positionTooltip = (event, tooltip) => {
        const tooltipWidth = 300; // max-width set in CSS
        const tooltipHeight = 100; // estimated height

        let left = event.clientX + 15;
        let top = event.clientY - 10;

        // Check right boundary
        if (left + tooltipWidth > window.innerWidth) {
            left = event.clientX - tooltipWidth - 15;
        }

        // Check bottom boundary
        if (top + tooltipHeight > window.innerHeight) {
            top = event.clientY - tooltipHeight - 15;
        }

        // Check top boundary
        if (top < 0) {
            top = 10;
        }

        // Check left boundary
        if (left < 0) {
            left = 10;
        }

        tooltip.style("left", left + "px")
            .style("top", top + "px");
    };

    const createGlyph = (dataToUse) => {
        const svg = d3.select(svgRef.current);
        svg.selectAll("*").remove();

        // If no data, just clear the SVG and return
        if (!pseudotimeData || pseudotimeData.length === 0) {
            return;
        }

        const { width, height } = dimensions;
        const margin = { top: 10, right: 15, bottom: 15, left: 15 };
        const innerWidth = width - margin.left - margin.right;
        const innerHeight = height - margin.top - margin.bottom;
        const centerX = innerWidth / 2;
        const centerY = innerHeight / 2;

        // Create main group
        const g = svg.append("g")
            .attr("transform", `translate(${margin.left}, ${margin.top})`);

        // Create tooltip (remove existing ones first to avoid duplicates)
        d3.select("body").selectAll(`.pseudotime-tooltip-${componentId}`).remove();
        const tooltip = d3.select("body")
            .append("div")
            .attr("class", `pseudotime-tooltip-${componentId}`)
            .style("position", "fixed")
            .style("visibility", "hidden")
            .style("background", "rgba(0, 0, 0, 0.9)")
            .style("color", "white")
            .style("padding", "10px")
            .style("border-radius", "6px")
            .style("font-size", "12px")
            .style("pointer-events", "none")
            .style("z-index", "9999")
            .style("max-width", "300px")
            .style("box-shadow", "0 4px 8px rgba(0,0,0,0.3)")
            .style("border", "1px solid rgba(255,255,255,0.2)");

        // Draw horizontal dividing line between gene expression (top) and cell trajectories (bottom)
        const axisLength = Math.min(innerWidth, innerHeight) * 0.98;
        g.append("line")
            .attr("x1", centerX - axisLength / 2)
            .attr("y1", centerY)
            .attr("x2", centerX + axisLength / 2)
            .attr("y2", centerY)
            .attr("stroke", "#333")
            .attr("stroke-width", 3)
            .attr("opacity", 0.8);

        // Process trajectory data
        const maxPseudotime = Math.max(...dataToUse.flatMap(traj =>
            traj.pseudotimes.map(pt => parseFloat(pt))
        ));

        // Color scale for different cell states/clusters
        const allClusters = [...new Set(dataToUse.flatMap(traj => traj.path))];

        // Use cluster colors from UMAP if available, otherwise use default colors
        let clusterColorScale;
        if (clusterColors) {
            // Create a color scale using the colors from the UMAP
            // Handle both full cluster names and numeric keys
            const colorRange = allClusters.map(cluster => {
                // Try the full cluster name first
                if (clusterColors[cluster]) {
                    return clusterColors[cluster];
                }
                // If not found, try extracting numeric part and use that as key
                const clusterNumber = cluster.toString().replace(/\D/g, '');
                return clusterColors[clusterNumber] || "#999";
            });
            clusterColorScale = d3.scaleOrdinal()
                .domain(allClusters)
                .range(colorRange);
        } else {
            // Fallback to default color scheme
            clusterColorScale = d3.scaleOrdinal(d3.schemeCategory10)
                .domain(allClusters);
        }

        // Add time point 0 circle at center
        g.append("circle")
            .attr("cx", centerX)
            .attr("cy", centerY)
            .attr("r", 8)
            .attr("fill", "#fff")
            .attr("stroke", "#333")
            .attr("stroke-width", 2)
            .attr("opacity", 0.9);

        // Add time point 0 label
        g.append("text")
            .attr("x", centerX)
            .attr("y", centerY + 4)
            .attr("text-anchor", "middle")
            .attr("font-size", "10px")
            .attr("font-weight", "bold")
            .attr("fill", "#333")
            .text("t0");

        // Create bottom section - macroscopic cell trajectories
        createBottomSection(g, dataToUse, centerX, centerY, axisLength, maxPseudotime, clusterColorScale, tooltip, selectedTrajectory, setSelectedTrajectory);

        // Create top section - gene expression gauge
        createTopSection(g, geneExpressionData, centerX, centerY, axisLength, maxPseudotime, tooltip);

        // Add title
        svg.append("text")
            .attr("x", width / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .attr("font-size", "14px")
            .attr("font-weight", "bold")
            .attr("fill", "#333")
            .text(adata_umap_title);
    };

    const createBottomSection = (g, trajectoryData, centerX, centerY, axisLength, maxPseudotime, clusterColorScale, tooltip, selectedTrajectory, setSelectedTrajectory) => {
        const bottomSection = g.append("g").attr("class", "bottom-section");
        const maxRadius = axisLength / 2 - 15;

        // Draw light brown background for the lower semicircle (soil-like color)
        const arc = d3.arc()
            .innerRadius(0)
            .outerRadius(maxRadius)
            .startAngle(Math.PI / 2)
            .endAngle(Math.PI * 1.5);

        bottomSection.append("path")
            .attr("d", arc)
            .attr("transform", `translate(${centerX}, ${centerY})`)
            .attr("fill", "#D2B48C")
            .attr("opacity", 0.3)
            .attr("stroke", "none")
            .style("cursor", "pointer")
        // .on("click", function (event) {
        //     // Deselect trajectory when clicking on background
        //     event.stopPropagation();
        //     setSelectedTrajectory(null);
        // });

        // Scale for converting pseudotime to radial distance
        const radiusScale = d3.scaleLinear()
            .domain([0, maxPseudotime])
            .range([8, maxRadius]); // Start from time point 0 circle edge (radius 8)

        // Build a tree structure from trajectory data to handle branching
        const buildTrajectoryTree = (trajectoryData) => {
            const nodes = new Map(); // Map of "cluster_pseudotime" -> node info
            const edges = []; // Array of connections between nodes
            const trajectoryPaths = []; // Track which trajectory each edge belongs to

            trajectoryData.forEach((trajectory, trajIndex) => {
                const { path, pseudotimes } = trajectory;

                // Create nodes for each cell state
                path.forEach((cluster, i) => {
                    const pseudotime = parseFloat(pseudotimes[i]);
                    const nodeKey = `${cluster}_${pseudotime}`;

                    if (!nodes.has(nodeKey)) {
                        nodes.set(nodeKey, {
                            cluster: cluster,
                            pseudotime: pseudotime,
                            radius: radiusScale(pseudotime),
                            trajectories: new Set([trajIndex]),
                            isEndpoint: false
                        });
                    } else {
                        // If node already exists, add this trajectory to it
                        nodes.get(nodeKey).trajectories.add(trajIndex);
                    }
                });

                // Mark endpoints
                const lastIndex = path.length - 1;
                const lastNodeKey = `${path[lastIndex]}_${parseFloat(pseudotimes[lastIndex])}`;
                if (nodes.has(lastNodeKey)) {
                    nodes.get(lastNodeKey).isEndpoint = true;
                }

                // Create edges between consecutive nodes
                for (let i = 0; i < path.length - 1; i++) {
                    const fromKey = `${path[i]}_${parseFloat(pseudotimes[i])}`;
                    const toKey = `${path[i + 1]}_${parseFloat(pseudotimes[i + 1])}`;

                    edges.push({
                        from: fromKey,
                        to: toKey,
                        trajectory: trajIndex
                    });
                    trajectoryPaths.push(trajIndex);
                }
            });

            return { nodes, edges, trajectoryPaths };
        };

        // Calculate positions for nodes using a tree layout
        const calculateNodePositions = (nodes, edges) => {
            const nodePositions = new Map();

            // Find root node (pseudotime = 0)
            const rootNode = Array.from(nodes.values()).find(node => node.pseudotime === 0);
            if (rootNode) {
                const rootKey = `${rootNode.cluster}_${rootNode.pseudotime}`;
                nodePositions.set(rootKey, { x: centerX, y: centerY });
            }

            // Build a proper tree structure to handle multiple branching points
            const buildTreeStructure = () => {
                const treeNodes = new Map();
                const childrenMap = new Map(); // parent -> children
                const parentMap = new Map();   // child -> parent

                // Initialize all nodes
                nodes.forEach((node, key) => {
                    treeNodes.set(key, {
                        ...node,
                        children: [],
                        parent: null,
                        branchIndex: -1,  // Which branch from its parent this node belongs to
                        level: 0  // Distance from root
                    });
                    childrenMap.set(key, []);
                });

                // Build parent-child relationships
                edges.forEach(edge => {
                    const parent = edge.from;
                    const child = edge.to;

                    if (!childrenMap.get(parent).includes(child)) {
                        childrenMap.get(parent).push(child);
                        parentMap.set(child, parent);
                        treeNodes.get(child).parent = parent;
                        treeNodes.get(parent).children.push(child);
                    }
                });

                // Calculate levels and assign branch indices
                const assignLevelsAndBranches = () => {
                    const rootKey = Array.from(treeNodes.keys()).find(key =>
                        treeNodes.get(key).pseudotime === 0
                    );

                    if (!rootKey) return;

                    // BFS to assign levels
                    const queue = [{ key: rootKey, level: 0 }];
                    const visited = new Set();

                    while (queue.length > 0) {
                        const { key, level } = queue.shift();
                        if (visited.has(key)) continue;
                        visited.add(key);

                        const node = treeNodes.get(key);
                        node.level = level;

                        // Assign branch indices to children based on their order
                        node.children.forEach((childKey, index) => {
                            const childNode = treeNodes.get(childKey);
                            childNode.branchIndex = index;
                            queue.push({ key: childKey, level: level + 1 });
                        });
                    }
                };

                assignLevelsAndBranches();
                return treeNodes;
            };

            const treeStructure = buildTreeStructure();

            // Define the bottom semicircle range for trajectories
            const minAngle = Math.PI / 6;  // 30 degrees
            const maxAngle = 5 * Math.PI / 6;  // 150 degrees

            // Calculate angle for each node based on its ancestry path
            const calculateNodeAngles = () => {
                const nodeAngles = new Map();

                // Find root and start with straight down angle
                const rootKey = Array.from(treeStructure.keys()).find(key =>
                    treeStructure.get(key).pseudotime === 0
                );

                if (rootKey) {
                    nodeAngles.set(rootKey, Math.PI / 2); // 90 degrees (straight down)
                }

                // Process nodes level by level
                const processedLevels = new Map();
                treeStructure.forEach((node, key) => {
                    if (!processedLevels.has(node.level)) {
                        processedLevels.set(node.level, []);
                    }
                    processedLevels.get(node.level).push({ key, node });
                });

                // Sort levels and process each level
                const sortedLevels = Array.from(processedLevels.entries()).sort((a, b) => a[0] - b[0]);

                sortedLevels.forEach(([level, nodesAtLevel]) => {
                    if (level === 0) return; // Root already processed

                    nodesAtLevel.forEach(({ key, node }) => {
                        const parentKey = node.parent;
                        const parentNode = treeStructure.get(parentKey);
                        const parentAngle = nodeAngles.get(parentKey);

                        if (parentAngle !== undefined && parentNode) {
                            const numSiblings = parentNode.children.length;

                            if (numSiblings === 1) {
                                // Single child: continue in same direction as parent
                                nodeAngles.set(key, parentAngle);
                            } else {
                                // Multiple siblings: spread them out
                                const branchIndex = node.branchIndex;

                                // Calculate angular spread based on number of siblings
                                let spreadAngle;
                                if (numSiblings === 2) {
                                    spreadAngle = Math.PI / 6; // 30 degrees each side
                                } else if (numSiblings === 3) {
                                    spreadAngle = Math.PI / 4; // 45 degrees total spread
                                } else {
                                    spreadAngle = Math.PI / 3; // 60 degrees total spread for more branches
                                }

                                // Calculate angle offset from parent angle
                                const angleStep = (2 * spreadAngle) / (numSiblings - 1);
                                const angleOffset = -spreadAngle + (branchIndex * angleStep);
                                let nodeAngle = parentAngle + angleOffset;

                                // Ensure angle stays within bounds
                                nodeAngle = Math.max(minAngle, Math.min(maxAngle, nodeAngle));
                                nodeAngles.set(key, nodeAngle);
                            }
                        }
                    });
                });

                return nodeAngles;
            };

            const nodeAngles = calculateNodeAngles();

            // Position all nodes based on their calculated angles
            treeStructure.forEach((node, key) => {
                if (node.pseudotime === 0) {
                    // Root is already positioned
                    return;
                }

                const angle = nodeAngles.get(key);
                const radius = radiusScale(node.pseudotime);

                if (angle !== undefined) {
                    const x = centerX + Math.cos(angle) * radius;
                    const y = centerY + Math.sin(angle) * radius;
                    nodePositions.set(key, { x, y });
                }
            });

            return nodePositions;
        };

        const { nodes, edges } = buildTrajectoryTree(trajectoryData);
        const nodePositions = calculateNodePositions(nodes, edges);

        // Use black color for all trajectory paths
        const trajectoryColor = "#333333"; // Black color

        // Draw edges (connections between nodes)
        edges.forEach(edge => {
            const fromPos = nodePositions.get(edge.from);
            const toPos = nodePositions.get(edge.to);

            // Determine if this edge should be grayed out
            const isSelected = selectedTrajectory === null || selectedTrajectory === edge.trajectory;
            const color = isSelected ? trajectoryColor : "#ccc";
            const opacity = isSelected ? 0.8 : 0.4;
            const strokeWidth = isSelected ? 3 : 2;

            if (fromPos && toPos) {
                let x1 = fromPos.x, y1 = fromPos.y;

                // If starting from center (time 0), start from edge of center circle
                const fromNode = nodes.get(edge.from);
                if (fromNode && fromNode.pseudotime === 0) {
                    const angle = Math.atan2(toPos.y - centerY, toPos.x - centerX);
                    x1 = centerX + Math.cos(angle) * 8;
                    y1 = centerY + Math.sin(angle) * 8;
                }

                bottomSection.append("line")
                    .attr("x1", x1)
                    .attr("y1", y1)
                    .attr("x2", toPos.x)
                    .attr("y2", toPos.y)
                    .attr("stroke", color)
                    .attr("stroke-width", strokeWidth)
                    .attr("opacity", opacity)
                    .attr("stroke-linecap", "round")
                    .style("cursor", "pointer")
                    // .on("click", function (event) {
                    //     event.stopPropagation();
                    //     // Toggle selection: if clicking on the same trajectory, deselect it
                    //     if (selectedTrajectory === edge.trajectory) {
                    //         setSelectedTrajectory(null);
                    //     } else {
                    //         setSelectedTrajectory(edge.trajectory);
                    //     }
                    // })
                    .on("mouseover", function (event) {
                        // Only enhance hover effect if this trajectory is not grayed out
                        if (isSelected) {
                            d3.select(this)
                                .attr("stroke-width", strokeWidth + 2)
                                .attr("opacity", 1);
                        }

                        const fromNode = nodes.get(edge.from);
                        const toNode = nodes.get(edge.to);
                        const fromCluster = fromNode ? fromNode.cluster : "unknown";
                        const toCluster = toNode ? toNode.cluster : "unknown";
                        const fromTime = fromNode ? fromNode.pseudotime.toFixed(2) : "0.00";
                        const toTime = toNode ? toNode.pseudotime.toFixed(2) : "0.00";

                        // Get the complete trajectory sequence
                        const trajectory = trajectoryData[edge.trajectory];
                        const trajectorySequence = trajectory.path.map((cluster, idx) =>
                            `${cluster} (t=${parseFloat(trajectory.pseudotimes[idx]).toFixed(2)})`
                        ).join(' → ');

                        // Emit trajectory hover event to parent components with debouncing
                        if (trajectory && source_title) {
                            const currentTrajectory = {
                                path: trajectory.path,
                                adata_umap_title: source_title,
                                // sampleId: sampleId,
                                pseudotimes: trajectory.pseudotimes,
                                trajectoryIndex: trajectoryIndex
                            };

                            debouncedSetHover(currentTrajectory);
                        }

                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Trajectory ${edge.trajectory + 1}</strong></div>
                                <div>Current transition: ${fromCluster} (t=${fromTime}) → ${toCluster} (t=${toTime})</div>
                                <div><strong>Complete sequence:</strong></div>
                                <div style="font-size: 11px; margin-top: 5px;">${trajectorySequence}</div>
                            `);
                        positionTooltip(event, tooltip);
                    })
                    .on("mousemove", function (event) {
                        positionTooltip(event, tooltip);
                    })
                    .on("mouseout", function () {
                        d3.select(this)
                            .attr("stroke-width", strokeWidth)
                            .attr("opacity", opacity);
                        tooltip.style("visibility", "hidden");

                        // Clear trajectory hover event with debouncing
                        clearHover();
                    });
            }
        });

        // Draw nodes
        nodes.forEach((node, key) => {
            const pos = nodePositions.get(key);
            if (!pos || node.radius <= 8) return; // Skip center node

            // Check if this node belongs to the selected trajectory
            const nodeTrajectories = Array.from(node.trajectories);
            const isNodeSelected = selectedTrajectory === null || nodeTrajectories.includes(selectedTrajectory);

            // Use cluster-based color for each cell state
            const baseColor = clusterColorScale(node.cluster);
            const nodeColor = isNodeSelected ? baseColor : "#ccc";
            const nodeOpacity = isNodeSelected ? 0.8 : 0.3;

            // Draw endpoint as star, others as circles
            if (node.isEndpoint) {
                const star = drawStar(bottomSection, pos.x, pos.y, 10, nodeColor, nodeOpacity);
                // Add tooltip functionality to star
                star.style("cursor", "pointer")
                    .on("mouseover", function (event) {
                        // Get all trajectories this node belongs to
                        const trajectoryInfo = Array.from(node.trajectories).map(trajIndex => {
                            const trajectory = trajectoryData[trajIndex];
                            const trajectorySequence = trajectory.path.map((cluster, idx) =>
                                `${cluster} (t=${parseFloat(trajectory.pseudotimes[idx]).toFixed(2)})`
                            ).join(' → ');
                            return `<div><strong>Trajectory ${trajIndex + 1}:</strong> ${trajectorySequence}</div>`;
                        }).join('');

                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Endpoint - Cluster ${node.cluster}</strong></div>
                                <div>Time: ${node.pseudotime.toFixed(2)}</div>
                                <div style="margin-top: 8px;"><strong>Belongs to trajectory(ies):</strong></div>
                                ${trajectoryInfo}
                            `);
                        positionTooltip(event, tooltip);
                    })
                    .on("mousemove", function (event) {
                        positionTooltip(event, tooltip);
                    })
                    .on("mouseout", function () {
                        tooltip.style("visibility", "hidden");
                    });
            } else {
                const circle = bottomSection.append("circle")
                    .attr("cx", pos.x)
                    .attr("cy", pos.y)
                    .attr("r", 6)
                    .attr("fill", nodeColor)
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 2)
                    .attr("opacity", nodeOpacity)
                    .style("cursor", "pointer")
                    .on("mouseover", function (event) {
                        // Get all trajectories this node belongs to
                        const trajectoryInfo = Array.from(node.trajectories).map(trajIndex => {
                            const trajectory = trajectoryData[trajIndex];
                            const trajectorySequence = trajectory.path.map((cluster, idx) =>
                                `${cluster} (t=${parseFloat(trajectory.pseudotimes[idx]).toFixed(2)})`
                            ).join(' → ');
                            return `<div><strong>Trajectory ${trajIndex + 1}:</strong> ${trajectorySequence}</div>`;
                        }).join('');

                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Cluster ${node.cluster}</strong></div>
                                <div>Time: ${node.pseudotime.toFixed(2)}</div>
                                <div style="margin-top: 8px;"><strong>Belongs to trajectory(ies):</strong></div>
                                ${trajectoryInfo}
                            `);
                        positionTooltip(event, tooltip);
                    })
                    .on("mousemove", function (event) {
                        positionTooltip(event, tooltip);
                    })
                    .on("mouseout", function () {
                        tooltip.style("visibility", "hidden");
                    });
            }
        });
    };

    const createTopSection = (g, geneData, centerX, centerY, axisLength, maxPseudotime, tooltip) => {
        const maxRadius = axisLength / 2 - 15;
        const topSection = g.append("g").attr("class", "top-section");

        // Add concentric circles for time progression using trajectory data time range
        // These should always be shown as a time reference
        const trajectoryMaxTime = maxPseudotime;
        const numTimeCircles = 4;
        for (let i = 1; i <= numTimeCircles; i++) {
            const time = (i / numTimeCircles) * trajectoryMaxTime;
            const radius = 8 + (time / trajectoryMaxTime) * (maxRadius - 8);

            // Draw complete circles
            topSection.append("circle")
                .attr("cx", centerX)
                .attr("cy", centerY)
                .attr("r", radius)
                .attr("fill", "none")
                .attr("stroke", "black")
                .attr("stroke-width", 1)
                .attr("opacity", 0.2);

            // Add time labels at the top of each circle
            topSection.append("text")
                .attr("x", centerX)
                .attr("y", centerY - radius - 5)
                .attr("text-anchor", "middle")
                .attr("font-size", "8px")
                .attr("fill", "#666")
                .text(`t${time.toFixed(1)}`);
        }

        // Only proceed with gene expression specific elements if data is provided
        if (!geneData || !Array.isArray(geneData) || geneData.length === 0) {
            return;
        }

        // Draw upper semicircle background for gene expression area
        // Time point scale (radial distance represents time progression)
        // Use the same scaling as concentric circles and trajectory data
        const timeScale = d3.scaleLinear()
            .domain([0, maxPseudotime])
            .range([8, maxRadius]); // Start from time point 0 circle edge (radius 8)

        // Expression scale (angular position - higher expression = more to the right)
        // Upper half only: from left (π) to right (2π) of the upper semicircle
        const expressionScale = d3.scaleLinear()
            .domain([0, 1])
            .range([Math.PI, 2 * Math.PI]); // From left side (180°) to right side (360°) through upper half

        // Gene colors
        const geneColors = ['#ff6b6b', '#4ecdc4', '#45b7d1'];

        geneData.forEach((geneInfo, geneIndex) => {
            const color = geneColors[geneIndex % geneColors.length];

            // Create expression points where angular position is determined by expression level
            const expressionPoints = geneInfo.timePoints.map((timePoint, i) => {
                // Radius is determined by time point (progression outward)
                const radius = timeScale(timePoint);
                // Angle is determined by expression level (higher = more to the right)
                const angle = expressionScale(geneInfo.expressions[i]);

                return {
                    x: centerX + Math.cos(angle) * radius,
                    y: centerY + Math.sin(angle) * radius,
                    timePoint: timePoint,
                    expression: geneInfo.expressions[i],
                    angle: angle,
                    radius: radius
                };
            });

            // Create smooth line generator that passes through all points
            const line = d3.line()
                .x(d => d.x)
                .y(d => d.y)
                .curve(d3.curveCatmullRom.alpha(0.5)); // Smooth curve passing through all points

            // Prepare points for smooth curve (including starting point if needed)
            let curvePoints = [...expressionPoints];

            // Add starting point from time 0 circle edge if first point is not at time 0
            if (expressionPoints.length > 0 && expressionPoints[0].radius > 8) {
                const firstPoint = expressionPoints[0];
                const angle = Math.atan2(firstPoint.y - centerY, firstPoint.x - centerX);
                const startX = centerX + Math.cos(angle) * 8;
                const startY = centerY + Math.sin(angle) * 8;

                curvePoints = [
                    { x: startX, y: startY },
                    ...expressionPoints
                ];
            }

            // Draw smooth curve through all points
            if (curvePoints.length > 1) {
                topSection.append("path")
                    .datum(curvePoints)
                    .attr("d", line)
                    .attr("stroke", color)
                    .attr("stroke-width", 3)
                    .attr("fill", "none")
                    .attr("opacity", 0.6)
                    .attr("stroke-linecap", "round")
                    .attr("stroke-linejoin", "round")
                    .style("cursor", "pointer")
                    .on("mouseover", function (event) {
                        d3.select(this)
                            .attr("stroke-width", 4)
                            .attr("opacity", 1);

                        const minExpression = Math.min(...geneInfo.expressions);
                        const maxExpression = Math.max(...geneInfo.expressions);
                        const avgExpression = (geneInfo.expressions.reduce((a, b) => a + b, 0) / geneInfo.expressions.length).toFixed(2);
                        const timeSpan = `${geneInfo.timePoints[0]} - ${geneInfo.timePoints[geneInfo.timePoints.length - 1]}`;

                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Gene Expression: ${geneInfo.gene}</strong></div>
                                <div>Time span: ${timeSpan}</div>
                                <div>Expression range: ${minExpression.toFixed(2)} - ${maxExpression.toFixed(2)}</div>
                                <div>Average expression: ${avgExpression}</div>
                            `);
                        positionTooltip(event, tooltip);
                    })
                    .on("mousemove", function (event) {
                        positionTooltip(event, tooltip);
                    })
                    .on("mouseout", function () {
                        d3.select(this)
                            .attr("stroke-width", 3)
                            .attr("opacity", 0.6);
                        tooltip.style("visibility", "hidden");
                    });
            }

            // Draw all points
            expressionPoints.forEach((point, i) => {
                topSection.append("circle")
                    .attr("cx", point.x)
                    .attr("cy", point.y)
                    .attr("r", 3)
                    .attr("fill", color)
                    .attr("stroke", "#fff")
                    .attr("opacity", 0.9)
                    .style("cursor", "pointer")
                    .on("mouseover", function (event) {
                        d3.select(this)
                            .attr("r", 5)
                            .attr("opacity", 1);

                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Gene Expression Point</strong></div>
                                <div>Gene: ${geneInfo.gene}</div>
                                <div>Time: ${point.timePoint.toFixed(2)}</div>
                                <div>Expression: ${point.expression.toFixed(3)}</div>
                            `);
                        positionTooltip(event, tooltip);
                    })
                    .on("mousemove", function (event) {
                        positionTooltip(event, tooltip);
                    })
                    .on("mouseout", function () {
                        d3.select(this)
                            .attr("r", 3)
                            .attr("opacity", 0.9);
                        tooltip.style("visibility", "hidden");
                    });
            });
        });

        // Add expression level indicators for upper semicircle
        topSection.append("text")
            .attr("x", centerX - maxRadius * 0.7)
            .attr("y", centerY - 10)
            .attr("text-anchor", "middle")
            .attr("font-size", "10px")
            .attr("fill", "#666")
            .text("Low Expr");

        topSection.append("text")
            .attr("x", centerX + maxRadius * 0.7)
            .attr("y", centerY - 10)
            .attr("text-anchor", "middle")
            .attr("font-size", "10px")
            .attr("fill", "#666")
            .text("High Expr");
    }

    const drawStar = (parent, cx, cy, radius, color, opacity = 0.9) => {
        const starPoints = 5;
        const angle = Math.PI / starPoints;
        let path = "";

        for (let i = 0; i < 2 * starPoints; i++) {
            const r = i % 2 === 0 ? radius : radius * 0.5;
            const x = cx + Math.cos(i * angle) * r;
            const y = cy + Math.sin(i * angle) * r;
            path += i === 0 ? `M ${x} ${y}` : ` L ${x} ${y}`;
        }
        path += " Z";

        return parent.append("path")
            .attr("d", path)
            .attr("fill", color)
            .attr("stroke", "#fff")
            .attr("stroke-width", 1)
            .attr("opacity", opacity);
    };

    return (
        <div ref={containerRef} style={{ width: "100%", height: "100%", position: "relative" }}>
            {/* Checkbox in upper left corner */}
            <div style={{
                position: 'absolute',
                top: '5px',
                left: '5px',
                zIndex: 1000,
                backgroundColor: 'rgba(255, 255, 255, 0.8)',
                borderRadius: '4px',
                padding: '2px'
            }}>
                <Checkbox
                    checked={isSelected}
                    onChange={(e) => onSelectionChange && onSelectionChange(e.target.checked)}
                    size="small"
                />
            </div>
            {pseudotimeLoading && (
                <div style={{
                    position: 'absolute',
                    top: '50%',
                    left: '50%',
                    transform: 'translate(-50%, -50%)'
                }}>
                    <Spin size="large" />
                </div>
            )}
            {!pseudotimeLoading && (!pseudotimeData || pseudotimeData.length === 0) ? (
                <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    width: '100%',
                    height: '100%'
                }}>
                    <Empty
                        description="No pseudotime data available"
                        image={Empty.PRESENTED_IMAGE_SIMPLE}
                    />
                </div>
            ) : (
                <svg
                    ref={svgRef}
                    width={dimensions.width}
                    height={dimensions.height}
                    viewBox={`0 0 ${dimensions.width} ${dimensions.height}`}
                    style={{
                        width: '100%',
                        height: '100%',
                        display: 'block',
                        border: '1px solid #ddd',
                        borderRadius: '5px',
                        backgroundColor: '#f9f9f9'
                    }}
                    preserveAspectRatio="xMidYMid meet"
                />
            )}
        </div>
    );
};

export default PseudotimeGlyph;