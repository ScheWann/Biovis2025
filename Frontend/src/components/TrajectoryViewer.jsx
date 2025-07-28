import React, { useRef, useEffect, useState } from "react";
import * as d3 from "d3";

// Individual trajectory chart component with error bands
const TrajectoryChart = ({
    data,
    gene,
    width = 300,
    height = 200,
    margin = { top: 30, right: 20, bottom: 50, left: 60 }
}) => {
    const ref = useRef();

    useEffect(() => {
        if (!data || data.length === 0) return;

        // Clear previous chart
        d3.select(ref.current).selectAll("*").remove();

        const innerWidth = width - margin.left - margin.right;
        const innerHeight = height - margin.top - margin.bottom;

        // Create scales
        const xScale = d3.scaleLinear()
            .domain(d3.extent(data, d => d.x))
            .nice()
            .range([0, innerWidth]);

        const yScale = d3.scaleLinear()
            .domain([
                d3.min(data, d => d.ymin),
                d3.max(data, d => d.ymax)
            ])
            .nice()
            .range([innerHeight, 0]);

        // Create line generator
        const line = d3.line()
            .x(d => xScale(d.x))
            .y(d => yScale(d.y))
            .curve(d3.curveMonotoneX);

        // Create area generator for error bands
        const area = d3.area()
            .x(d => xScale(d.x))
            .y0(d => yScale(d.ymin))
            .y1(d => yScale(d.ymax))
            .curve(d3.curveMonotoneX);

        // Create SVG
        const svg = d3.select(ref.current)
            .attr("width", width)
            .attr("height", height);

        const g = svg.append("g")
            .attr("transform", `translate(${margin.left},${margin.top})`);

        // Add error band
        g.append("path")
            .datum(data)
            .attr("fill", "gray")
            .attr("fill-opacity", 0.3)
            .attr("d", area);

        // Add main line
        g.append("path")
            .datum(data)
            .attr("fill", "none")
            .attr("stroke", "#e74c3c")
            .attr("stroke-width", 2)
            .attr("d", line);

        // Add axes
        g.append("g")
            .attr("transform", `translate(0,${innerHeight})`)
            .call(d3.axisBottom(xScale));

        g.append("g")
            .call(d3.axisLeft(yScale));

        // Add labels
        svg.append("text")
            .attr("x", margin.left + innerWidth / 2)
            .attr("y", height - 8)
            .attr("text-anchor", "middle")
            .attr("font-size", 12)
            .text("Distance along Trajectory");

        svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -height / 2)
            .attr("y", 16)
            .attr("text-anchor", "middle")
            .attr("font-size", 12)
            .text("Estimated Expression");

        // Add title
        svg.append("text")
            .attr("x", margin.left)
            .attr("y", margin.top - 10)
            .attr("text-anchor", "start")
            .attr("font-size", 14)
            .attr("font-weight", "bold")
            .text(gene);

    }, [data, gene, width, height, margin]);

    return <svg ref={ref}></svg>;
};

// Main trajectory viewer component
export const TrajectoryViewer = ({ sampleId }) => {
    const [trajectoryData, setTrajectoryData] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);

    const fetchTrajectoryData = async (sampleId) => {
        if (!sampleId) return;

        setLoading(true);
        setError(null);

        try {
            const response = await fetch("http://localhost:5003/api/get_trajectory_data", {
                method: "POST",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify({ sample_id: sampleId }),
            });

            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }

            const data = await response.json();

            if (data.error) {
                throw new Error(data.error);
            }

            setTrajectoryData(data);
        } catch (err) {
            setError(err.message);
            console.error("Error fetching trajectory data:", err);
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        fetchTrajectoryData(sampleId);
    }, [sampleId]);

    if (loading) {
        return (
            <div className="trajectory-viewer" style={{ padding: "20px", textAlign: "center" }}>
                <p>Loading trajectory data...</p>
            </div>
        );
    }

    if (error) {
        return (
            <div className="trajectory-viewer" style={{ padding: "20px", textAlign: "center" }}>
                <p style={{ color: "red" }}>Error: {error}</p>
            </div>
        );
    }

    if (!trajectoryData) {
        return (
            <div className="trajectory-viewer" style={{ padding: "20px", textAlign: "center" }}>
                <p>No trajectory data available for this sample.</p>
            </div>
        );
    }

    const genes = Object.keys(trajectoryData);

    return (
        <div className="trajectory-viewer" style={{ padding: "20px" }}>
            <h3 style={{ marginBottom: "20px", textAlign: "center" }}>
                Gene Expression Trajectory Analysis
            </h3>

            <div
                style={{
                    display: "grid",
                    gridTemplateColumns: "repeat(auto-fit, minmax(320px, 1fr))",
                    gap: "20px",
                    justifyItems: "center"
                }}
            >
                {genes.map((gene) => (
                    <div key={gene} style={{
                        border: "1px solid #ddd",
                        borderRadius: "8px",
                        padding: "10px",
                        backgroundColor: "#f9f9f9"
                    }}>
                        <TrajectoryChart
                            data={trajectoryData[gene].data}
                            gene={gene}
                            width={320}
                            height={220}
                        />
                    </div>
                ))}
            </div>

            <div style={{
                marginTop: "20px",
                padding: "10px",
                backgroundColor: "#f0f0f0",
                borderRadius: "5px",
                fontSize: "12px"
            }}>
                <p><strong>Note:</strong> Gray bands represent confidence intervals. Red lines show estimated gene expression along the trajectory.</p>
            </div>
        </div>
    );
}; 