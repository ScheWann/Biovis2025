import React, { useRef, useEffect, useState } from "react";
import { Select, Button, Row, Col, message, Spin, Empty } from "antd";
import { LineChart } from "./LineChart";

const { Option } = Select;

// Main trajectory viewer component
export const TrajectoryViewer = ({ sampleId }) => {
    const [samples, setSamples] = useState([]);
    const [selectedSample, setSelectedSample] = useState(null);
    const [availableGenes, setAvailableGenes] = useState([]);
    const [selectedGenes, setSelectedGenes] = useState([]);
    const [confirmedGenes, setConfirmedGenes] = useState([]);
    const [trajectoryData, setTrajectoryData] = useState({});
    const [loading, setLoading] = useState(false);
    const [geneListLoading, setGeneListLoading] = useState(false);
    const containerRef = useRef();
    const [containerHeight, setContainerHeight] = useState(400);

    // Track container height for dynamic sizing
    useEffect(() => {
        const updateHeight = () => {
            if (containerRef.current) {
                const rect = containerRef.current.getBoundingClientRect();
                setContainerHeight(rect.height);
            }
        };

        updateHeight();
        const resizeObserver = new ResizeObserver(updateHeight);
        if (containerRef.current) {
            resizeObserver.observe(containerRef.current);
        }

        return () => resizeObserver.disconnect();
    }, []);

    // Fetch available samples on component mount
    useEffect(() => {
        fetch("/api/get_samples_option")
            .then((response) => response.json())
            .then((data) => {
                setSamples(data);
                if (sampleId) {
                    setSelectedSample(sampleId);
                }
            })
            .catch((error) => {
                message.error("Get samples failed");
            });
    }, [sampleId]);

    // Fetch gene list when sample changes
    useEffect(() => {
        if (selectedSample) {
            fetchGeneList(selectedSample);
        } else {
            setAvailableGenes([]);
            setSelectedGenes([]);
            setConfirmedGenes([]);
            setTrajectoryData({});
        }
    }, [selectedSample]);

    const fetchGeneList = async (sample_id) => {
        setGeneListLoading(true);

        fetch("/api/get_trajectory_gene_list", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ sample_id }),
        }).then((response) => response.json())
            .then((data) => {
                setAvailableGenes(data);
                setSelectedGenes([]); // Reset selected genes when sample changes
            })
            .catch((error) => {
                setAvailableGenes([]);
            })
            .finally(() => {
                setGeneListLoading(false);
            });
    };

    const fetchTrajectoryData = async (sample_id, genes) => {
        setLoading(true);

        fetch("/api/get_trajectory_data", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ sample_id, selected_genes: genes }),
        }).then((response) => response.json())
            .then((data) => {
                setTrajectoryData(data);
            })
            .catch((error) => {
                setTrajectoryData({});
            })
            .finally(() => {
                setLoading(false);
            });
    };

    const handleConfirm = () => {
        if (!selectedSample) {
            message.warning("Please select a sample first");
            return;
        }
        if (selectedGenes.length === 0) {
            message.warning("Please select at least one gene");
            return;
        }

        setConfirmedGenes([...selectedGenes]);
        fetchTrajectoryData(selectedSample, selectedGenes);
        message.success(`Loading trajectory data for ${selectedGenes.length} gene(s)`);
    };

    // Flatten sample options for Select component
    const sampleOptions = samples.flatMap(group =>
        group.options || []
    );

    // Calculate chart height based on number of confirmed genes
    const getChartHeight = () => {
        if (confirmedGenes.length === 0) return containerHeight;
        // Use the same height calculation for both single and multiple charts
        return containerHeight - 32; // Account for controls
    };

    const chartHeight = getChartHeight();
    const enableScrolling = confirmedGenes.length > 1;


    return (
        <div ref={containerRef} className="trajectory-viewer" style={{ height: "100%", display: "flex", flexDirection: "column" }}>
            {/* Control Panel */}
            <div style={{ display: "flex", justifyContent: "flex-end", alignItems: "center", flexShrink: 0, padding: "0px 10px 8px 10px" }}>
                <Row gutter={16} align="middle">
                    <Col flex="200px">
                        <Select
                            size="small"
                            placeholder="Select Sample"
                            style={{ width: "100%" }}
                            value={selectedSample}
                            onChange={setSelectedSample}
                            showSearch
                            filterOption={(input, option) =>
                                option.children.toLowerCase().indexOf(input.toLowerCase()) >= 0
                            }
                        >
                            {sampleOptions.map(sample => (
                                <Option key={sample.value} value={sample.value}>
                                    {sample.label}
                                </Option>
                            ))}
                        </Select>
                    </Col>
                    <Col flex="200px">
                        <Select
                            size="small"
                            mode="multiple"
                            placeholder="Select Genes"
                            style={{ width: "100%" }}
                            value={selectedGenes}
                            onChange={setSelectedGenes}
                            disabled={!selectedSample || geneListLoading}
                            loading={geneListLoading}
                            showSearch
                            filterOption={(input, option) =>
                                option.children.toLowerCase().indexOf(input.toLowerCase()) >= 0
                            }
                            maxTagCount="responsive"
                        >
                            {availableGenes.map(gene => (
                                <Option key={gene} value={gene}>
                                    {gene}
                                </Option>
                            ))}
                        </Select>
                    </Col>
                    <Col>
                        <Button
                            size="small"
                            type="primary"
                            onClick={handleConfirm}
                            disabled={!selectedSample || selectedGenes.length === 0 || loading}
                            loading={loading}
                        >
                            OK
                        </Button>
                    </Col>
                </Row>
            </div>

            {/* Charts Container */}
            <div
                style={{
                    flex: 1,
                    overflowY: enableScrolling ? "auto" : "hidden",
                    overflowX: "hidden",
                    display: "flex",
                    alignItems: enableScrolling ? "flex-start" : "center",
                    justifyContent: "center",
                }}
            >
                {loading && (
                    <Spin size="large" />
                )}

                {!loading && !selectedSample && (
                    <Empty
                        description="Please select a sample to view trajectory data"
                        image={Empty.PRESENTED_IMAGE_SIMPLE}
                    />
                )}

                {!loading && selectedSample && confirmedGenes.length === 0 && (
                    <Empty
                        description="Select genes and click OK to view trajectory analysis"
                        image={Empty.PRESENTED_IMAGE_SIMPLE}
                    />
                )}

                {!loading && confirmedGenes.length > 0 && Object.keys(trajectoryData).length === 0 && (
                    <Empty
                        description="No trajectory data available for selected genes"
                        image={Empty.PRESENTED_IMAGE_SIMPLE}
                    />
                )}

                {!loading && confirmedGenes.length > 0 && Object.keys(trajectoryData).length > 0 && (
                    <div
                        style={{
                            display: enableScrolling ? "block" : "flex",
                            flexDirection: enableScrolling ? "column" : "row",
                            gap: enableScrolling ? "20px" : "0",
                            height: enableScrolling ? "auto" : "100%",
                            width: "100%",
                        }}
                    >
                        {confirmedGenes.map((gene) => (
                            trajectoryData[gene] && (
                                <div
                                    key={gene}
                                    style={{
                                        backgroundColor: "#f9f9f9",
                                        height: `${chartHeight}px`,
                                        display: "flex",
                                        flexDirection: "column",
                                        width: "100%",
                                        marginBottom: enableScrolling ? "10px" : "0",
                                        flex: enableScrolling ? "none" : "1",
                                        borderRadius: "8px",
                                    }}
                                >
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
                            )
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
}; 