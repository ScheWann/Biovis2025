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
            .range([0, maxRadius]); // Start from 0 (center) instead of 15

        trajectoryData.forEach((trajectory, trajIndex) => {
            const { path, pseudotimes } = trajectory;
            const color = colorScale(trajIndex);

            // Calculate angle for this trajectory to avoid horizontal line
            // Distribute trajectories across bottom half but avoid exact horizontal line (π/6 to 5π/6)
            const numTrajectories = trajectoryData.length;
            const angleSpacing = (2 * Math.PI / 3) / (numTrajectories + 1); // 2π/3 = 120 degrees range
            const baseAngle = Math.PI / 6 + (trajIndex + 1) * angleSpacing; // Start at 30 degrees, avoid 0 and π

            // Draw trajectory path points radiating from center
            const pathPoints = pseudotimes.map((pt, i) => {
                const pseudotime = parseFloat(pt);
                const radius = radiusScale(pseudotime);
                
                return {
                    x: centerX + Math.cos(baseAngle) * radius,
                    y: centerY + Math.sin(baseAngle) * radius,
                    pseudotime: pseudotime,
                    cluster: path[i],
                    isLast: i === pseudotimes.length - 1,
                    timeIndex: i,
                    angle: baseAngle,
                    radius: radius
                };
            });

            // Draw connecting lines between consecutive time points
            for (let i = 0; i < pathPoints.length - 1; i++) {
                bottomSection.append("line")
                    .attr("x1", pathPoints[i].x)
                    .attr("y1", pathPoints[i].y)
                    .attr("x2", pathPoints[i + 1].x)
                    .attr("y2", pathPoints[i + 1].y)
                    .attr("stroke", color)
                    .attr("stroke-width", 2)
                    .attr("opacity", 0.7);
            }

            // Draw all points along the radial trajectory
            pathPoints.forEach((point, i) => {
                // Skip drawing point at exact center (radius = 0) to avoid clutter
                if (point.radius === 0) return;
                
                // Draw a star for the final stage, circle for others
                if (point.isLast) {
                    drawStar(bottomSection, point.x, point.y, 10, color);
                } else {
                    bottomSection.append("circle")
                        .attr("cx", point.x)
                        .attr("cy", point.y)
                        .attr("r", 6)
                        .attr("fill", color)
                        .attr("stroke", "#fff")
                        .attr("stroke-width", 2)
                        .attr("opacity", 0.8);
                }

                // Add cluster and time labels positioned to avoid overlap
                const labelOffset = point.radius < 50 ? 15 : 25;
                bottomSection.append("text")
                    .attr("x", point.x)
                    .attr("y", point.y + labelOffset)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "9px")
                    .attr("fill", "#666")
                    .text(`C${point.cluster}`);
                
                // Add pseudotime value
                bottomSection.append("text")
                    .attr("x", point.x)
                    .attr("y", point.y + labelOffset + 10)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "8px")
                    .attr("fill", "#999")
                    .text(`${point.pseudotime.toFixed(2)}`);
            });

            // Add trajectory label near the outer end
            const lastPoint = pathPoints[pathPoints.length - 1];
            if (lastPoint.radius > 0) {
                const labelX = centerX + Math.cos(baseAngle) * (maxRadius + 15);
                const labelY = centerY + Math.sin(baseAngle) * (maxRadius + 15);
                
                bottomSection.append("text")
                    .attr("x", labelX)
                    .attr("y", labelY)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "10px")
                    .attr("font-weight", "bold")
                    .attr("fill", color)
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
        const timeScale = d3.scaleLinear()
            .domain([0, Math.max(...mockGeneData.flatMap(d => d.timePoints))])
            .range([20, maxRadius - 10]);

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

            // Draw line from center to first point
            if (expressionPoints.length > 0) {
                topSection.append("line")
                    .attr("x1", centerX)
                    .attr("y1", centerY)
                    .attr("x2", expressionPoints[0].x)
                    .attr("y2", expressionPoints[0].y)
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
            const radius = 20 + (time / trajectoryMaxTime) * (maxRadius - 30); // Scale radius based on trajectory time range
            
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