import React, { useRef, useEffect, useState } from 'react';
import * as d3 from 'd3';

const PseudotimeGlyph = ({
    sampleId,
    cellIds,
    adata_umap_title,
    early_markers = null,
    width = 400,
    height = 400,
    geneExpressionData = null
}) => {
    console.log(geneExpressionData, "geneExpressionData")
    const svgRef = useRef(null);
    const [trajectoryData, setTrajectoryData] = useState(null);
    const [loading, setLoading] = useState(false);

    // Fetch pseudotime data
    const fetchPseudotimeData = async () => {
        if (!sampleId || !cellIds || cellIds.length === 0) return;

        setLoading(true);
        try {
            const response = await fetch("/api/get_pseudotime_data", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    sample_id: sampleId,
                    cell_ids: cellIds,
                    adata_umap_title: adata_umap_title,
                    early_markers: early_markers,
                    n_neighbors: 15,
                    n_pcas: 30,
                    resolutions: 1,
                }),
            });

            const data = await response.json();
            console.log("Pseudotime data received:", data);
            setTrajectoryData(data);
        } catch (error) {
            console.error("Failed to fetch pseudotime data:", error);
        } finally {
            setLoading(false);
        }
    };

    // useEffect(() => {
    //     fetchPseudotimeData();
    // }, [sampleId, cellIds, adata_umap_title]);

    // Generate mock trajectory data for testing when no real data is available
    const generateMockTrajectoryData = () => {
        return [
            {
                path: [0, 1, 3],
                pseudotimes: ["0.000", "0.200", "0.500", "0.800"]
            },
            {
                path: [0, 1, 4, 5],
                pseudotimes: ["0.000", "0.200", "0.600", "1.000"]
            },
            {
                path: [0, 2, 6],
                pseudotimes: ["0.000", "0.300", "0.700"]
            }
        ];
    };

    useEffect(() => {
        // Use real trajectory data if available, otherwise use mock data for testing
        const dataToUse = trajectoryData && trajectoryData.length > 0 
            ? trajectoryData 
            : generateMockTrajectoryData();
        
        createGlyph(dataToUse);
    }, [trajectoryData, width, height, geneExpressionData]);

    // Cleanup tooltip on unmount
    useEffect(() => {
        return () => {
            d3.select("body").selectAll(".pseudotime-tooltip").remove();
        };
    }, []);

    const createGlyph = (dataToUse) => {
        const svg = d3.select(svgRef.current);
        svg.selectAll("*").remove();

        const margin = { top: 40, right: 40, bottom: 40, left: 40 };
        const innerWidth = width - margin.left - margin.right;
        const innerHeight = height - margin.top - margin.bottom;
        const centerX = innerWidth / 2;
        const centerY = innerHeight / 2;

        // Create main group
        const g = svg.append("g")
            .attr("transform", `translate(${margin.left}, ${margin.top})`);

        // Create tooltip (remove existing ones first to avoid duplicates)
        d3.select("body").selectAll(".pseudotime-tooltip").remove();
        const tooltip = d3.select("body")
            .append("div")
            .attr("class", "pseudotime-tooltip")
            .style("position", "absolute")
            .style("visibility", "hidden")
            .style("background", "rgba(0, 0, 0, 0.8)")
            .style("color", "white")
            .style("padding", "8px")
            .style("border-radius", "4px")
            .style("font-size", "12px")
            .style("pointer-events", "none")
            .style("z-index", "1000");

        // Draw horizontal dividing line between gene expression (top) and cell trajectories (bottom)
        const axisLength = Math.min(innerWidth, innerHeight) * 0.8;
        g.append("line")
            .attr("x1", centerX - axisLength / 2)
            .attr("y1", centerY)
            .attr("x2", centerX + axisLength / 2)
            .attr("y2", centerY)
            .attr("stroke", "#333")
            .attr("stroke-width", 3)
            .attr("opacity", 0.8);

        // Add axis labels
        // g.append("text")
        //     .attr("x", centerX - axisLength / 2 - 10)
        //     .attr("y", centerY - 10)
        //     .attr("text-anchor", "middle")
        //     .attr("font-size", "12px")
        //     .attr("font-weight", "bold")
        //     .attr("fill", "#666")
        //     .text("Gene Expression");

        // g.append("text")
        //     .attr("x", centerX - axisLength / 2 - 10)
        //     .attr("y", centerY + 20)
        //     .attr("text-anchor", "middle")
        //     .attr("font-size", "12px")
        //     .attr("font-weight", "bold")
        //     .attr("fill", "#666")
        //     .text("Cell Trajectories");

        // Process trajectory data
        const maxPseudotime = Math.max(...dataToUse.flatMap(traj =>
            traj.pseudotimes.map(pt => parseFloat(pt))
        ));

        // Color scale for different cell states/clusters
        const allClusters = [...new Set(dataToUse.flatMap(traj => traj.path))];
        const clusterColorScale = d3.scaleOrdinal(d3.schemeCategory10)
            .domain(allClusters);

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
        createBottomSection(g, dataToUse, centerX, centerY, axisLength, maxPseudotime, clusterColorScale, tooltip);

        // Create top section - gene expression gauge
        createTopSection(g, geneExpressionData, centerX, centerY, axisLength, maxPseudotime, tooltip);

        // Add title
        svg.append("text")
            .attr("x", width / 2)
            .attr("y", 20)
            .attr("text-anchor", "middle")
            .attr("font-size", "14px")
            .attr("font-weight", "bold")
            .attr("fill", "#333")
            .text("Pseudotime Analysis Glyph");
    };

    const createBottomSection = (g, trajectoryData, centerX, centerY, axisLength, maxPseudotime, clusterColorScale, tooltip) => {
        const bottomSection = g.append("g").attr("class", "bottom-section");
        const maxRadius = axisLength / 2 - 30;
        
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
            .attr("stroke", "none");
        
        // Extract all unique clusters for legend
        const allClusters = [...new Set(trajectoryData.flatMap(traj => traj.path))];
        
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

                // Draw edges (connections between nodes)
        edges.forEach(edge => {
            const fromPos = nodePositions.get(edge.from);
            const toPos = nodePositions.get(edge.to);
            
            // Use the target node's cluster for edge color
            const toNode = nodes.get(edge.to);
            const color = toNode ? clusterColorScale(toNode.cluster) : "#999";
            
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
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.7)
                    .style("cursor", "pointer")
                    .on("mouseover", function(event) {
                        d3.select(this)
                            .attr("stroke-width", 4)
                            .attr("opacity", 1);
                        
                        const fromNode = nodes.get(edge.from);
                        const toNode = nodes.get(edge.to);
                        const fromCluster = fromNode ? fromNode.cluster : "unknown";
                        const toCluster = toNode ? toNode.cluster : "unknown";
                        const fromTime = fromNode ? fromNode.pseudotime.toFixed(2) : "0.00";
                        const toTime = toNode ? toNode.pseudotime.toFixed(2) : "0.00";
                        
                        tooltip.style("visibility", "visible")
                            .html(`
                                <div><strong>Cell State Transition</strong></div>
                                <div>From: Cluster ${fromCluster} (t=${fromTime})</div>
                                <div>To: Cluster ${toCluster} (t=${toTime})</div>
                                <div>Trajectory: ${edge.trajectory + 1}</div>
                                <div>Cell differentiation pathway</div>
                            `)
                            .style("left", (event.pageX + 10) + "px")
                            .style("top", (event.pageY - 10) + "px");
                    })
                    .on("mousemove", function(event) {
                        tooltip.style("left", (event.pageX + 10) + "px")
                            .style("top", (event.pageY - 10) + "px");
                    })
                    .on("mouseout", function() {
                        d3.select(this)
                            .attr("stroke-width", 2)
                            .attr("opacity", 0.7);
                        tooltip.style("visibility", "hidden");
                    });
            }
        });

        // Draw nodes
        nodes.forEach((node, key) => {
            const pos = nodePositions.get(key);
            if (!pos || node.radius <= 8) return; // Skip center node

            // Use cluster-based color for each cell state
            const nodeColor = clusterColorScale(node.cluster);

            // Draw endpoint as star, others as circles
            if (node.isEndpoint) {
                drawStar(bottomSection, pos.x, pos.y, 10, nodeColor);
            } else {
                bottomSection.append("circle")
                    .attr("cx", pos.x)
                    .attr("cy", pos.y)
                    .attr("r", 6)
                    .attr("fill", nodeColor)
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.8);
            }

            // Add labels
            // const labelOffset = node.radius < 50 ? 15 : 25;
            // bottomSection.append("text")
            //     .attr("x", pos.x)
            //     .attr("y", pos.y + labelOffset)
            //     .attr("text-anchor", "middle")
            //     .attr("font-size", "9px")
            //     .attr("fill", "#666")
            //     .text(`C${node.cluster}`);

            // bottomSection.append("text")
            //     .attr("x", pos.x)
            //     .attr("y", pos.y + labelOffset + 10)
            //     .attr("text-anchor", "middle")
            //     .attr("font-size", "8px")
            //     .attr("fill", "#999")
            //     .text(`${node.pseudotime.toFixed(2)}`);
        });

        // Add cell state color legend
        const legendX = centerX - axisLength / 2 + 20;
        const legendStartY = centerY + axisLength / 2 - 20;
        
        // allClusters.sort((a, b) => a - b).forEach((cluster, index) => {
        //     const legendY = legendStartY - index * 18;
            
        //     bottomSection.append("circle")
        //         .attr("cx", legendX)
        //         .attr("cy", legendY)
        //         .attr("r", 6)
        //         .attr("fill", clusterColorScale(cluster))
        //         .attr("stroke", "#fff")
        //         .attr("stroke-width", 2);

        //     bottomSection.append("text")
        //         .attr("x", legendX + 15)
        //         .attr("y", legendY + 4)
        //         .attr("font-size", "12px")
        //         .attr("font-weight", "bold")
        //         .attr("fill", "#333")
        //         .text(`Cell ${cluster}`);
        // });

        // Add legend title
        // bottomSection.append("text")
        //     .attr("x", legendX)
        //     .attr("y", legendStartY + 20)
        //     .attr("font-size", "12px")
        //     .attr("font-weight", "bold")
        //     .attr("fill", "#333")
        //     .text("Cell States:");

        // Trajectory color scale for labels (separate from cluster colors)
        const trajectoryColorScale = d3.scaleOrdinal(d3.schemeSet2);

        // Add trajectory labels for endpoints
        trajectoryData.forEach((trajectory, trajIndex) => {
            const lastCluster = trajectory.path[trajectory.path.length - 1];
            const lastTime = parseFloat(trajectory.pseudotimes[trajectory.pseudotimes.length - 1]);
            const lastNodeKey = `${lastCluster}_${lastTime}`;
            const lastPos = nodePositions.get(lastNodeKey);
            
            if (lastPos) {
                const angle = Math.atan2(lastPos.y - centerY, lastPos.x - centerX);
                const labelX = centerX + Math.cos(angle) * (maxRadius + 15);
                const labelY = centerY + Math.sin(angle) * (maxRadius + 15);
                
                // bottomSection.append("text")
                //     .attr("x", labelX)
                //     .attr("y", labelY)
                //     .attr("text-anchor", "middle")
                //     .attr("font-size", "10px")
                //     .attr("font-weight", "bold")
                //     .attr("fill", trajectoryColorScale(trajIndex))
                //     .text(`Traj ${trajIndex + 1}`);
            }
        });
    };

    const createTopSection = (g, geneData, centerX, centerY, axisLength, maxPseudotime, tooltip) => {
        const topSection = g.append("g").attr("class", "top-section");
        const maxRadius = axisLength / 2 - 30;

        // Mock gene expression data if not provided
        const mockGeneData = geneData || [
            { gene: "SOX2", timePoints: [0.0, 0.4, 0.7, 1.0], expressions: [0.8, 0.6, 0.3, 0.2] },
            { gene: "NANOG", timePoints: [0.1, 0.4, 0.7, 1.0], expressions: [0.2, 0.5, 0.8, 0.9] },
            { gene: "OCT4", timePoints: [0.0, 0.3, 0.6, 1.0], expressions: [0.9, 0.7, 0.4, 0.1] }
        ];

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

        mockGeneData.forEach((geneInfo, geneIndex) => {
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
                    .attr("stroke-width", 2)
                    .attr("fill", "none")
                    .attr("opacity", 0.6)
                    .style("cursor", "pointer")
                    .on("mouseover", function(event) {
                        d3.select(this)
                            .attr("stroke-width", 4)
                            .attr("opacity", 1);
                        console.log(geneInfo)
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
                                <div>Gene regulation along pseudotime</div>
                            `)
                            .style("left", (event.pageX + 10) + "px")
                            .style("top", (event.pageY - 10) + "px");
                    })
                    .on("mousemove", function(event) {
                        tooltip.style("left", (event.pageX + 10) + "px")
                            .style("top", (event.pageY - 10) + "px");
                    })
                    .on("mouseout", function() {
                        d3.select(this)
                            .attr("stroke-width", 2)
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
                    .attr("opacity", 0.9);

                // Add expression value labels
                // topSection.append("text")
                //     .attr("x", point.x)
                //     .attr("y", point.y - 15)
                //     .attr("text-anchor", "middle")
                //     .attr("font-size", "8px")
                //     .attr("fill", "#333")
                //     .text(point.expression.toFixed(2));

                // Add time point labels
                // topSection.append("text")
                //     .attr("x", point.x)
                //     .attr("y", point.y + 20)
                //     .attr("text-anchor", "middle")
                //     .attr("font-size", "8px")
                //     .attr("fill", "#666")
                //     .text(`t${point.timePoint.toFixed(1)}`);
            });

            // Add gene legend in the top-left area
            const legendX = centerX - axisLength / 2 + 20;
            const legendY = centerY - axisLength / 2 + 20 + geneIndex * 18;
            
            // topSection.append("circle")
            //     .attr("cx", legendX)
            //     .attr("cy", legendY)
            //     .attr("r", 5)
            //     .attr("fill", color);

            // topSection.append("text")
            //     .attr("x", legendX + 15)
            //     .attr("y", legendY + 4)
            //     .attr("font-size", "12px")
            //     .attr("font-weight", "bold")
            //     .attr("fill", "#333")
            //     .text(geneInfo.gene);
        });

        // Add concentric circles for time progression using trajectory data time range
        const trajectoryMaxTime = maxPseudotime; // Use the maxPseudotime calculated from trajectory data
        const numTimeCircles = 4;
        for (let i = 1; i <= numTimeCircles; i++) {
            const time = (i / numTimeCircles) * trajectoryMaxTime;
            const radius = 8 + (time / trajectoryMaxTime) * (maxRadius - 8); // Scale radius from time 0 circle edge to maxRadius
            
            // Draw complete circles instead of semicircles
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
    };

    const drawStar = (parent, cx, cy, radius, color) => {
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

        parent.append("path")
            .attr("d", path)
            .attr("fill", color)
            .attr("stroke", "#fff")
            .attr("stroke-width", 1)
            .attr("opacity", 0.9);
    };

    return (
        <div style={{ width, height }}>
            {loading && (
                <div style={{
                    position: 'absolute',
                    width,
                    height,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    backgroundColor: 'rgba(249, 249, 249, 0.8)',
                    border: '1px solid #ddd',
                    borderRadius: '4px',
                    zIndex: 10
                }}>
                    <div>Loading pseudotime data...</div>
                </div>
            )}
            <svg
                ref={svgRef}
                width={width}
                height={height}
                style={{ border: '1px solid #ddd', borderRadius: '4px' }}
            />
        </div>
    );
};

export default PseudotimeGlyph;