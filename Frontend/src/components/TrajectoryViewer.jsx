import React, { useRef, useEffect, useState } from "react";
import { LineChart } from "./LineChart";

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
        <div className="trajectory-viewer">
            <div
                style={{
                    display: "grid",
                    gridTemplateColumns: "repeat(auto-fit, minmax(320px, 1fr))",
                    padding: "0px 20px 0px 20px",
                    gap: "20px",
                    justifyItems: "stretch"
                }}
            >
                {genes.map((gene) => (
                    <div key={gene} style={{
                        backgroundColor: "#f9f9f9",
                        height: "240px",
                        display: "flex",
                        flexDirection: "column",
                        width: "100%",
                    }}>
                        <LineChart
                            data={trajectoryData[gene].data}
                            xAccessor={d => d.x}
                            yAccessor={d => d.y}
                            showErrorBands={true}
                            yMinAccessor={d => d.ymin}
                            yMaxAccessor={d => d.ymax}
                            margin={{ top: 30, right: 20, bottom: 50, left: 60 }}
                            xLabel="Distance along Trajectory"
                            yLabel="Estimated Expression"
                            title={gene}
                            lineColor="#e74c3c"
                            lineWidth={2}
                            errorBandColor="gray"
                            errorBandOpacity={0.3}
                        />
                    </div>
                ))}
            </div>
        </div>
    );
}; 