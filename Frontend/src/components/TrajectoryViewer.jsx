import React, { useRef, useEffect, useState, useCallback, useMemo } from "react";
import { Select, Button, Row, Col, message, Spin, Empty, Switch } from "antd";
import { LineChart } from "./LineChart";

const { Option } = Select;

// Main trajectory viewer component
export const TrajectoryViewer = ({ sampleId, samples, kosaraDisplayEnabled, onKosaraDisplayToggle, onGeneSelection, onTrajectoryGuidelineChange }) => {
    const [samplesData, setSamplesData] = useState([]);
    const [selectedSample, setSelectedSample] = useState(null);
    const [availableGenes, setAvailableGenes] = useState([]);
    const [selectedGenes, setSelectedGenes] = useState([]);
    const [confirmedGenes, setConfirmedGenes] = useState([]);
    const [trajectoryData, setTrajectoryData] = useState({});
    const [loading, setLoading] = useState(false);
    const [geneListLoading, setGeneListLoading] = useState(false);
    const containerRef = useRef();
    const [containerHeight, setContainerHeight] = useState(400);
    const [isVertical, setIsVertical] = useState(false);

    // Throttle mouse move events to prevent excessive updates
    const lastMouseMoveRef = useRef({ time: 0, position: null, xValue: null });

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

    // Use passed samples and update selected sample when sampleId changes
    useEffect(() => {
        if (samples) {
            setSamplesData(samples);
        }
        if (sampleId) {
            setSelectedSample(sampleId);
        }
    }, [sampleId, samples]);

    // Fetch gene list when sample changes or orientation toggles
    useEffect(() => {
        if (selectedSample) {
            fetchGeneList(selectedSample);
        } else {
            setAvailableGenes([]);
            setSelectedGenes([]);
            setConfirmedGenes([]);
            setTrajectoryData({});
        }
    }, [selectedSample, isVertical]);

    const fetchGeneList = async (sample_id) => {
        setGeneListLoading(true);

        fetch("/api/get_trajectory_gene_list", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ sample_id, is_vertical: isVertical }),
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
            body: JSON.stringify({ sample_id, selected_genes: genes, is_vertical: isVertical }),
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
        
        // Notify parent about gene selection for Kosara display
        if (onGeneSelection && kosaraDisplayEnabled) {
            onGeneSelection([...selectedGenes], selectedSample);
        }
    };

    // Flatten sample options for Select component
    const sampleOptions = samplesData.flatMap(group =>
        group.options || []
    );

    // Calculate chart height - now always the same since we use one chart
    const getChartHeight = () => {
        if (confirmedGenes.length === 0) return containerHeight;
        return containerHeight - 32; // Account for controls
    };

    const chartHeight = getChartHeight();

    // Handle mouse movement over trajectory chart with throttling
    const handleTrajectoryMouseMove = useCallback((normalizedPosition, xValue) => {
        if (!onTrajectoryGuidelineChange || !selectedSample) return;

        const now = Date.now();
        const THROTTLE_MS = 16; // ~60fps
        const lastMove = lastMouseMoveRef.current;

        // Throttle updates to prevent excessive re-renders
        if (now - lastMove.time < THROTTLE_MS) return;

        // Only update if values have changed significantly (prevent floating point drift)
        const positionChanged = Math.abs((lastMove.position || 0) - normalizedPosition) > 0.001;
        const xValueChanged = Math.abs((lastMove.xValue || 0) - xValue) > 0.001;

        if (!positionChanged && !xValueChanged) return;

        // Update our tracking reference
        lastMouseMoveRef.current = { time: now, position: normalizedPosition, xValue: xValue };

        onTrajectoryGuidelineChange({
            sampleId: selectedSample,
            position: normalizedPosition,
            xValue: xValue,
            isVertical: isVertical,
            visible: true
        });
    }, [onTrajectoryGuidelineChange, selectedSample, isVertical]);

    // Handle mouse leave from trajectory chart
    const handleTrajectoryMouseLeave = useCallback(() => {
        if (onTrajectoryGuidelineChange) {
            // Reset our tracking reference
            lastMouseMoveRef.current = { time: 0, position: null, xValue: null };
            onTrajectoryGuidelineChange({
                visible: false
            });
        }
    }, [onTrajectoryGuidelineChange]);

    // Memoize LineChart props to prevent unnecessary re-renders
    const singleGeneChartProps = useMemo(() => ({
        data: trajectoryData[confirmedGenes[0]]?.data,
        xAccessor: d => d.x,
        yAccessor: d => d.y,
        showErrorBands: true,
        yMinAccessor: d => d.ymin,
        yMaxAccessor: d => d.ymax,
        margin: { top: 30, right: 20, bottom: 50, left: 60 },
        lineColor: "#e74c3c",
        errorBandOpacity: 0.3,
        onMouseMove: handleTrajectoryMouseMove,
        onMouseLeave: handleTrajectoryMouseLeave
    }), [trajectoryData, confirmedGenes, handleTrajectoryMouseMove, handleTrajectoryMouseLeave]);

    const multiGeneChartProps = useMemo(() => ({
        datasets: confirmedGenes
            .filter(gene => trajectoryData[gene])
            .map(gene => ({
                data: trajectoryData[gene].data,
                xAccessor: d => d.x,
                yAccessor: d => d.y,
                yMinAccessor: d => d.ymin,
                yMaxAccessor: d => d.ymax,
                label: gene,
                lineColor: undefined
            })),
        showErrorBands: true,
        showLegend: true,
        margin: { top: 30, right: 20, bottom: 40, left: 60 },
        errorBandOpacity: 0.3,
        onMouseMove: handleTrajectoryMouseMove,
        onMouseLeave: handleTrajectoryMouseLeave
    }), [trajectoryData, confirmedGenes, handleTrajectoryMouseMove, handleTrajectoryMouseLeave]);


    return (
        <div ref={containerRef} className="trajectory-viewer" style={{ height: "100%", display: "flex", flexDirection: "column" }}>
            {/* Control Panel */}
            <div style={{
                display: "flex",
                justifyContent: "space-between",
                alignItems: "center",
                flexShrink: 0,
                padding: "8px 10px 8px 10px",
                flexWrap: "wrap",
                gap: "8px"
            }}>
                <div style={{
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "flex-end",
                    gap: "8px",
                    flexWrap: "wrap",
                    flex: 1
                }}>
                    <div style={{ display: "flex", alignItems: "center", gap: "16px" }}>
                        <div style={{ display: "flex", alignItems: "center", gap: "8px" }}>
                            <span style={{ fontSize: "12px", color: "#666" }}>Horizontal</span>
                            <Switch
                                size="small"
                                checked={isVertical}
                                onChange={setIsVertical}
                            />
                            <span style={{ fontSize: "12px", color: "#666" }}>Vertical</span>
                        </div>
                        
                        <div style={{ display: "flex", alignItems: "center", gap: "8px" }}>
                            <span style={{ fontSize: "12px", color: "#666" }}>Kosara Display:</span>
                            <Switch
                                size="small"
                                checked={kosaraDisplayEnabled}
                                onChange={onKosaraDisplayToggle}
                            />
                        </div>
                    </div>
                    <Select
                        size="small"
                        placeholder="Select Sample"
                        style={{ width: "200px", minWidth: "150px" }}
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
                    <Select
                        size="small"
                        mode="multiple"
                        placeholder="Select Genes"
                        style={{ width: "150px", minWidth: "80px" }}
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
                    <Button
                        size="small"
                        type="primary"
                        onClick={handleConfirm}
                        disabled={!selectedSample || selectedGenes.length === 0 || loading}
                        loading={loading}
                        style={{ flexShrink: 0 }}
                    >
                        OK
                    </Button>
                </div>
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
                            <LineChart {...singleGeneChartProps} />
                        ) : (
                            // Multiple genes: combine into single chart
                            <LineChart {...multiGeneChartProps} />
                        )}
                    </div>
                )}
            </div>
        </div>
    );
};