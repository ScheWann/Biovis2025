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
                path: [0, 1, 2, 3],
                pseudotimes: ["0.000", "0.200", "0.500", "0.800"]
            },
            {
                path: [0, 1, 4, 5],
                pseudotimes: ["0.000", "0.200", "0.600", "1.000"]
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
        g.append("text")
            .attr("x", centerX - axisLength / 2 - 10)
            .attr("y", centerY - 10)
            .attr("text-anchor", "middle")
            .attr("font-size", "12px")
            .attr("font-weight", "bold")
            .attr("fill", "#666")
            .text("Gene Expression");

        g.append("text")
            .attr("x", centerX - axisLength / 2 - 10)
            .attr("y", centerY + 20)
            .attr("text-anchor", "middle")
            .attr("font-size", "12px")
            .attr("font-weight", "bold")
            .attr("fill", "#666")
            .text("Cell Trajectories");

        // Process trajectory data
        const maxPseudotime = Math.max(...dataToUse.flatMap(traj =>
            traj.pseudotimes.map(pt => parseFloat(pt))
        ));

        // Color scale for different trajectories
        const colorScale = d3.scaleOrdinal(d3.schemeCategory10);

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
        createBottomSection(g, dataToUse, centerX, centerY, axisLength, maxPseudotime, colorScale);

        // Create top section - gene expression gauge
        createTopSection(g, geneExpressionData, centerX, centerY, axisLength, maxPseudotime);

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

    const createBottomSection = (g, trajectoryData, centerX, centerY, axisLength, maxPseudotime, colorScale) => {
        const bottomSection = g.append("g").attr("class", "bottom-section");
        const maxRadius = axisLength / 2 - 30;
        
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

            // Build a proper tree structure to identify trunk vs branches
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
                        isTrunk: false,
                        branchGroup: -1
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

                // Identify trunk nodes (nodes that are on the path from root to first branching point)
                const markTrunkNodes = () => {
                    let currentNode = Array.from(treeNodes.keys()).find(key => 
                        treeNodes.get(key).pseudotime === 0
                    );
                    
                    while (currentNode) {
                        const node = treeNodes.get(currentNode);
                        node.isTrunk = true;
                        
                        // If this node has multiple children, it's a branching point
                        if (node.children.length > 1) {
                            // Mark children with branch groups
                            node.children.forEach((childKey, index) => {
                                markBranchGroup(childKey, index, treeNodes);
                            });
                            break;
                        } else if (node.children.length === 1) {
                            currentNode = node.children[0];
                        } else {
                            break;
                        }
                    }
                };

                const markBranchGroup = (nodeKey, branchIndex, treeNodes) => {
                    const node = treeNodes.get(nodeKey);
                    if (node) {
                        node.branchGroup = branchIndex;
                        node.children.forEach(childKey => {
                            markBranchGroup(childKey, branchIndex, treeNodes);
                        });
                    }
                };

                markTrunkNodes();
                return treeNodes;
            };

            const treeStructure = buildTreeStructure();

            // Define the bottom semicircle range for trajectories
            const minAngle = Math.PI / 6;  // 30 degrees
            const maxAngle = 5 * Math.PI / 6;  // 150 degrees
            const trunkAngle = Math.PI / 2; // 90 degrees (straight down for trunk)

            // Position nodes based on tree structure
            const positionedNodes = new Set();

            // First, position trunk nodes along a straight line downward
            const trunkNodes = Array.from(treeStructure.entries())
                .filter(([key, node]) => node.isTrunk)
                .sort((a, b) => a[1].pseudotime - b[1].pseudotime);

            trunkNodes.forEach(([key, node]) => {
                if (node.pseudotime === 0) {
                    // Root is already positioned
                    positionedNodes.add(key);
                } else {
                    const radius = radiusScale(node.pseudotime);
                    const x = centerX + Math.cos(trunkAngle) * radius;
                    const y = centerY + Math.sin(trunkAngle) * radius;
                    nodePositions.set(key, { x, y });
                    positionedNodes.add(key);
                }
            });

            // Then, position branch nodes
            const branchGroups = new Map();
            treeStructure.forEach((node, key) => {
                if (!node.isTrunk && node.branchGroup >= 0) {
                    if (!branchGroups.has(node.branchGroup)) {
                        branchGroups.set(node.branchGroup, []);
                    }
                    branchGroups.get(node.branchGroup).push({ key, node });
                }
            });

            // Sort branch groups by branch index
            const sortedBranchGroups = Array.from(branchGroups.entries()).sort((a, b) => a[0] - b[0]);
            
            sortedBranchGroups.forEach(([branchIndex, branchNodes]) => {
                // Find the branching point (parent of first node in this branch)
                const firstBranchNode = branchNodes.find(({ node }) => node.parent && treeStructure.get(node.parent).isTrunk);
                if (!firstBranchNode) return;

                const branchingPoint = nodePositions.get(firstBranchNode.node.parent);
                if (!branchingPoint) return;

                // Calculate branch angle based on branch index
                const numBranches = sortedBranchGroups.length;
                let branchAngle;
                
                if (numBranches === 1) {
                    branchAngle = trunkAngle; // Continue straight if only one branch
                } else if (numBranches === 2) {
                    // For two branches, spread them symmetrically around the trunk
                    const spreadAngle = Math.PI / 4; // 45 degrees spread
                    branchAngle = trunkAngle + (branchIndex === 0 ? -spreadAngle : spreadAngle);
                } else {
                    // For more branches, distribute across available angle range
                    const angleStep = (maxAngle - minAngle) / (numBranches - 1);
                    branchAngle = minAngle + branchIndex * angleStep;
                }

                // Ensure angle is within bounds
                branchAngle = Math.max(minAngle, Math.min(maxAngle, branchAngle));

                // Position all nodes in this branch along the branch angle
                const branchNodesSorted = branchNodes.sort((a, b) => a.node.pseudotime - b.node.pseudotime);
                
                branchNodesSorted.forEach(({ key, node }) => {
                    const radius = radiusScale(node.pseudotime);
                    const x = centerX + Math.cos(branchAngle) * radius;
                    const y = centerY + Math.sin(branchAngle) * radius;
                    nodePositions.set(key, { x, y });
                    positionedNodes.add(key);
                });
            });

            return nodePositions;
        };

        const { nodes, edges } = buildTrajectoryTree(trajectoryData);
        const nodePositions = calculateNodePositions(nodes, edges);

        // Draw edges (connections between nodes)
        edges.forEach(edge => {
            const fromPos = nodePositions.get(edge.from);
            const toPos = nodePositions.get(edge.to);
            const color = colorScale(edge.trajectory);
            
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
                    .attr("opacity", 0.7);
            }
        });

        // Draw nodes
        nodes.forEach((node, key) => {
            const pos = nodePositions.get(key);
            if (!pos || node.radius <= 8) return; // Skip center node

            // Determine color - if multiple trajectories pass through this node, use a neutral color
            const trajArray = Array.from(node.trajectories);
            const nodeColor = trajArray.length === 1 ? colorScale(trajArray[0]) : "#666";

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
            const labelOffset = node.radius < 50 ? 15 : 25;
                bottomSection.append("text")
                .attr("x", pos.x)
                .attr("y", pos.y + labelOffset)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "9px")
                    .attr("fill", "#666")
                .text(`C${node.cluster}`);
                
                bottomSection.append("text")
                .attr("x", pos.x)
                .attr("y", pos.y + labelOffset + 10)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "8px")
                    .attr("fill", "#999")
                .text(`${node.pseudotime.toFixed(2)}`);
        });

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
                
                bottomSection.append("text")
                    .attr("x", labelX)
                    .attr("y", labelY)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "10px")
                    .attr("font-weight", "bold")
                    .attr("fill", colorScale(trajIndex))
                    .text(`Traj ${trajIndex + 1}`);
            }
        });
    };

    const createTopSection = (g, geneData, centerX, centerY, axisLength, maxPseudotime) => {
        const topSection = g.append("g").attr("class", "top-section");
        const maxRadius = axisLength / 2 - 30;

        // Mock gene expression data if not provided
        const mockGeneData = geneData || [
            { gene: "SOX2", timePoints: [0.0, 0.4, 0.7, 1.0], expressions: [0.8, 0.6, 0.3, 0.2] },
            { gene: "NANOG", timePoints: [0.1, 0.4, 0.7, 1.0], expressions: [0.2, 0.5, 0.8, 0.9] },
            { gene: "OCT4", timePoints: [0.0, 0.3, 0.6, 0.9], expressions: [0.9, 0.7, 0.4, 0.1] }
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

            // Draw line from time 0 circle edge to first point (if not at time 0)
            if (expressionPoints.length > 0 && expressionPoints[0].radius > 8) {
                // Calculate starting point on the edge of time 0 circle in direction of first point
                const firstPoint = expressionPoints[0];
                const angle = Math.atan2(firstPoint.y - centerY, firstPoint.x - centerX);
                const startX = centerX + Math.cos(angle) * 8;
                const startY = centerY + Math.sin(angle) * 8;
                
                topSection.append("line")
                    .attr("x1", startX)
                    .attr("y1", startY)
                    .attr("x2", firstPoint.x)
                    .attr("y2", firstPoint.y)
                    .attr("stroke", color)
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.6);
            }

            // Draw connecting lines between consecutive expression points
            for (let i = 0; i < expressionPoints.length - 1; i++) {
                topSection.append("line")
                    .attr("x1", expressionPoints[i].x)
                    .attr("y1", expressionPoints[i].y)
                    .attr("x2", expressionPoints[i + 1].x)
                    .attr("y2", expressionPoints[i + 1].y)
                    .attr("stroke", color)
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.6);
            }

            // Draw all points
            expressionPoints.forEach((point, i) => {
                topSection.append("circle")
                    .attr("cx", point.x)
                    .attr("cy", point.y)
                    .attr("r", 5)
                    .attr("fill", color)
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.9);

                // Add expression value labels
                topSection.append("text")
                    .attr("x", point.x)
                    .attr("y", point.y - 15)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "8px")
                    .attr("fill", "#333")
                    .text(point.expression.toFixed(2));

                // Add time point labels
                topSection.append("text")
                    .attr("x", point.x)
                    .attr("y", point.y + 20)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "8px")
                    .attr("fill", "#666")
                    .text(`t${point.timePoint.toFixed(1)}`);
            });

            // Add gene legend in the top-left area
            const legendX = centerX - axisLength / 2 + 20;
            const legendY = centerY - axisLength / 2 + 20 + geneIndex * 18;
            
            topSection.append("circle")
                .attr("cx", legendX)
                .attr("cy", legendY)
                .attr("r", 5)
                .attr("fill", color);

            topSection.append("text")
                .attr("x", legendX + 15)
                .attr("y", legendY + 4)
                .attr("font-size", "12px")
                .attr("font-weight", "bold")
                .attr("fill", "#333")
                .text(geneInfo.gene);
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