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

    // Calculate chart height - now always the same since we use one chart
    const getChartHeight = () => {
        if (confirmedGenes.length === 0) return containerHeight;
        return containerHeight - 32; // Account for controls
    };

    const chartHeight = getChartHeight();


    return (
        <div ref={containerRef} className="trajectory-viewer" style={{ height: "100%", display: "flex", flexDirection: "column" }}>
            {/* Control Panel */}
            <div style={{ display: "flex", justifyContent: "flex-end", alignItems: "center", flexShrink: 0, padding: "8px 10px 8px 10px" }}>
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
                    overflowY: "hidden",
                    overflowX: "hidden",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                }}
            >
                {loading && (
                    <Spin size="large" />
                )}

                {!loading && confirmedGenes.length === 0 && (
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
                            backgroundColor: "#f9f9f9",
                            height: `${chartHeight}px`,
                            display: "flex",
                            flexDirection: "column",
                            width: "100%",
                            borderRadius: "8px",
                        }}
                    >
                        {confirmedGenes.length === 1 ? (
                            // Single gene: use original approach
                            <LineChart
                                data={trajectoryData[confirmedGenes[0]].data}
                                xAccessor={d => d.x}
                                yAccessor={d => d.y}
                                showErrorBands={true}
                                yMinAccessor={d => d.ymin}
                                yMaxAccessor={d => d.ymax}
                                margin={{ top: 30, right: 20, bottom: 50, left: 60 }}
                                lineColor="#e74c3c"
                                errorBandOpacity={0.3}
                            />
                        ) : (
                            // Multiple genes: combine into single chart
                            <LineChart
                                datasets={confirmedGenes
                                    .filter(gene => trajectoryData[gene])
                                    .map(gene => ({
                                        data: trajectoryData[gene].data,
                                        xAccessor: d => d.x,
                                        yAccessor: d => d.y,
                                        yMinAccessor: d => d.ymin,
                                        yMaxAccessor: d => d.ymax,
                                        label: gene,
                                        lineColor: undefined // Let the chart choose colors automatically
                                    }))
                                }
                                showErrorBands={true}
                                showLegend={true}
                                margin={{ top: 30, right: 20, bottom: 40, left: 60 }}
                                errorBandOpacity={0.3}
                            />
                        )}
                    </div>
                )}
            </div>
        </div>
    );
}; 