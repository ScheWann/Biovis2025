import React, { useRef, useState, useEffect } from 'react';
import { Box, Button, TextField, Stack, CircularProgress, Alert, Typography, Tooltip, IconButton } from '@mui/material';
import { styled } from '@mui/material/styles';
import axios from 'axios';
import * as d3 from 'd3';
import ZoomInIcon from '@mui/icons-material/ZoomIn';
import ZoomOutIcon from '@mui/icons-material/ZoomOut';
import PanToolIcon from '@mui/icons-material/PanTool';
import FilterListIcon from '@mui/icons-material/FilterList';
import PaletteIcon from '@mui/icons-material/Palette';

const ControlPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    bottom: 0,
    left: 0,
    display: 'flex',
    flexDirection: 'column',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
    '& .MuiTextField-root': {
        width: '60px',
        '& .MuiInputBase-root': {
            height: '25px',
            fontSize: '0.8rem',
        },
        '& .MuiInputLabel-root': {
            fontSize: '0.8rem',
            transform: 'translate(8px, 4px) scale(1)',
        },
    }
}));

const ToolPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    top: theme.spacing(1),
    left: theme.spacing(1),
    display: 'flex',
    flexDirection: 'row',
    gap: theme.spacing(0.1),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(0.1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
}));

const UMAPContainer = styled(Box)(({ theme }) => ({
    width: '100%',
    height: '100%',
    position: 'relative',
    backgroundColor: '#fff',
    borderRadius: theme.shape.borderRadius,
    overflow: 'hidden',
    margin: '0',
}));

const FilterPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(0.5),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
    maxHeight: '200px',
    overflowY: 'auto',
    '& .MuiTypography-root': {
        fontSize: '0.75rem',
        marginBottom: theme.spacing(0.5),
    },
    '& .MuiStack-root': {
        gap: theme.spacing(0.5),
    },
    '& .MuiBox-root': {
        '& .MuiTypography-root': {
            fontSize: '0.7rem',
        },
        '& input[type="checkbox"]': {
            width: '14px',
            height: '14px',
        }
    }
}));

const ColorPickerPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
}));

const defaultColors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
];

export const PseudoTemporalViewer = () => {
    const [samplePercent, setSamplePercent] = useState(1);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [umapData, setUmapData] = useState(null);
    const [selectedCellTypes, setSelectedCellTypes] = useState([]);
    const [showFilters, setShowFilters] = useState(false);
    const svgRef = useRef(null);
    const zoomRef = useRef(null);

    useEffect(() => {
        fetchUMAPData();
    }, [samplePercent]);

    const fetchUMAPData = async () => {
        try {
            setLoading(true);
            setError(null);
            
            console.log('Fetching UMAP data with sample_percent:', samplePercent);
            const response = await axios.get('/get_deaplog_results', {
                params: {
                    sample_percent: samplePercent
                }
            });
            
            console.log('Received UMAP data:', response.data);
            if (!response.data || !response.data.umap_coordinates) {
                throw new Error('Invalid data format received from server');
            }
            
            setUmapData(response.data);
            if (response.data.metadata && response.data.metadata.unique_cell_types) {
                setSelectedCellTypes(response.data.metadata.unique_cell_types);
            } else {
                console.warn('No unique cell types found in metadata');
            }
        } catch (err) {
            console.error('Error fetching UMAP data:', err);
            setError(err.response?.data?.error || 'Failed to fetch UMAP data');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        if (!umapData || !svgRef.current) return;

        d3.select(svgRef.current).selectAll("*").remove();

        const width = svgRef.current.clientWidth;
        const height = svgRef.current.clientHeight;
        const margin = { top: 40, right: 120, bottom: 40, left: 60 };

        const svg = d3.select(svgRef.current)
            .attr("width", width)
            .attr("height", height);

        const xScale = d3.scaleLinear()
            .domain(d3.extent(umapData.umap_coordinates.x))
            .range([margin.left, width - margin.right]);

        const yScale = d3.scaleLinear()
            .domain(d3.extent(umapData.umap_coordinates.y))
            .range([height - margin.bottom, margin.top]);

        const colorScale = d3.scaleOrdinal()
            .domain(umapData.metadata.unique_cell_types)
            .range(d3.schemeCategory10);

        // Create a container for static elements (axes, labels, legend)
        const staticContainer = svg.append("g");

        // Create a container for interactive elements (points)
        const interactiveContainer = svg.append("g");

        // Add static elements
        staticContainer.append("text")
            .attr("x", width / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .text("UMAP1")
            .style("font-size", "14px")
            .style("font-weight", "bold");

        staticContainer.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -height / 2)
            .attr("y", 15)
            .attr("text-anchor", "middle")
            .text("UMAP2")
            .style("font-size", "14px")
            .style("font-weight", "bold");

        staticContainer.append("text")
            .attr("x", width / 2)
            .attr("y", margin.top / 2)
            .attr("text-anchor", "middle")
            .text("Pseudo-temporal Gene Ordering");

        // Add legend
        const legend = staticContainer.append("g")
            .attr("font-family", "sans-serif")
            .attr("font-size", 10)
            .attr("text-anchor", "start")
            .selectAll("g")
            .data(umapData.metadata.unique_cell_types)
            .enter().append("g")
            .attr("transform", (d, i) => {
                const totalHeight = height - margin.top - margin.bottom;
                const itemHeight = totalHeight / umapData.metadata.unique_cell_types.length;
                return `translate(0,${margin.top + (i * itemHeight)})`;
            });

        legend.append("rect")
            .attr("x", width - margin.right + 10)
            .attr("width", 19)
            .attr("height", 19)
            .attr("fill", d => colorScale(d));

        legend.append("text")
            .attr("x", width - margin.right + 35)
            .attr("y", 9.5)
            .attr("dy", "0.32em")
            .text(d => d);

        // Add zoom behavior only to interactive container
        const zoom = d3.zoom()
            .scaleExtent([0.1, 10])
            .on("zoom", (event) => {
                interactiveContainer.attr("transform", event.transform);
            });

        interactiveContainer.call(zoom);

        // Add pan behavior
        let isPanning = false;
        let startX, startY;

        interactiveContainer.on("mousedown", function(event) {
            isPanning = true;
            startX = event.clientX;
            startY = event.clientY;
        });

        interactiveContainer.on("mousemove", function(event) {
            if (!isPanning) return;
            
            const dx = event.clientX - startX;
            const dy = event.clientY - startY;
            
            const currentTransform = d3.zoomTransform(interactiveContainer.node());
            const newTransform = currentTransform.translate(dx / currentTransform.k, dy / currentTransform.k);
            
            interactiveContainer.transition()
                .duration(0)
                .call(zoom.transform, newTransform);
            
            startX = event.clientX;
            startY = event.clientY;
        });

        interactiveContainer.on("mouseup", function() {
            isPanning = false;
        });

        interactiveContainer.on("mouseleave", function() {
            isPanning = false;
        });

        // Add points to interactive container
        const points = interactiveContainer.selectAll("circle")
            .data(umapData.umap_coordinates.x)
            .enter()
            .append("circle")
            .attr("cx", (d, i) => xScale(d))
            .attr("cy", (d, i) => yScale(umapData.umap_coordinates.y[i]))
            .attr("r", 3)
            .attr("fill", (d, i) => {
                const cellType = umapData.cell_types[i];
                return selectedCellTypes.includes(cellType) 
                    ? colorScale(cellType)
                    : "#ddd";
            })
            .attr("opacity", (d, i) => selectedCellTypes.includes(umapData.cell_types[i]) ? 0.6 : 0.1)
            .on("mouseover", function(event, d) {
                const index = event.target.__data__;
                const cellType = umapData.cell_types[index];
                if (selectedCellTypes.includes(cellType)) {
                    d3.select(this)
                        .attr("r", 5)
                        .attr("opacity", 1);

                    const tooltip = d3.select("body")
                        .append("div")
                        .attr("class", "tooltip")
                        .style("opacity", 0)
                        .style("position", "absolute")
                        .style("background-color", "white")
                        .style("border", "1px solid #ddd")
                        .style("border-radius", "4px")
                        .style("padding", "8px")
                        .style("pointer-events", "none")
                        .style("box-shadow", "0 2px 4px rgba(0,0,0,0.1)")
                        .style("font-size", "12px")
                        .style("z-index", "1000");

                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .9);

                    tooltip.html(`Cell Type: ${cellType}<br>Index: ${index}`)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }
            })
            .on("mouseout", function(event, d) {
                const cellType = umapData.cell_types[event.target.__data__];
                if (selectedCellTypes.includes(cellType)) {
                    d3.select(this)
                        .attr("r", 3)
                        .attr("opacity", 0.6);
                    
                    d3.selectAll(".tooltip").remove();
                }
            });

        zoomRef.current = zoom;

    }, [umapData, selectedCellTypes]);

    const handleSampleChange = (event) => {
        const value = Math.min(Math.max(0.1, parseFloat(event.target.value) || 0.1), 100);
        setSamplePercent(value);
    };

    const handleZoomIn = () => {
        if (zoomRef.current) {
            const svg = d3.select(svgRef.current);
            const currentTransform = d3.zoomTransform(svg.node());
            const newScale = currentTransform.k * 1.5;
            svg.transition()
                .duration(300)
                .call(zoomRef.current.transform, d3.zoomIdentity
                    .translate(0, currentTransform.y)
                    .scale(newScale));
        }
    };

    const handleZoomOut = () => {
        if (zoomRef.current) {
            const svg = d3.select(svgRef.current);
            const currentTransform = d3.zoomTransform(svg.node());
            const newScale = currentTransform.k * 0.75;
            svg.transition()
                .duration(300)
                .call(zoomRef.current.transform, d3.zoomIdentity
                    .translate(0, currentTransform.y)
                    .scale(newScale));
        }
    };

    const handleCellTypeToggle = (cellType) => {
        setSelectedCellTypes(prev => 
            prev.includes(cellType)
                ? prev.filter(type => type !== cellType)
                : [...prev, cellType]
        );
    };

    return (
        <Box sx={{ 
            height: '33vh', 
            width: '100%', 
            position: 'relative',
            display: 'flex', 
            flexDirection: 'column',
            borderBottom: '2px solid #e8e8e8',
            gap: 0
        }}>
            <UMAPContainer>
                {loading ? (
                    <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
                        <CircularProgress />
                    </Box>
                ) : error ? (
                    <Alert severity="error" sx={{ m: 1 }}>{error}</Alert>
                ) : (
                    <svg ref={svgRef} style={{ width: '100%', height: '100%' }}></svg>
                )}
            </UMAPContainer>

            <ToolPanel>
                <Tooltip title="Filter Cell Types">
                    <IconButton onClick={() => setShowFilters(!showFilters)} size="small" sx={{ padding: '3px' }}>
                        <FilterListIcon sx={{ fontSize: '1.0rem' }} />
                    </IconButton>
                </Tooltip>
            </ToolPanel>

            {showFilters && umapData && (
                <FilterPanel>
                    <Typography variant="subtitle2">Cell Types</Typography>
                    <Stack spacing={0.5}>
                        {umapData.metadata.unique_cell_types.map(cellType => (
                            <Box key={cellType} sx={{ display: 'flex', alignItems: 'center', gap: 0.1 }}>
                                <input
                                    type="checkbox"
                                    checked={selectedCellTypes.includes(cellType)}
                                    onChange={() => handleCellTypeToggle(cellType)}
                                />
                                <Typography variant="body2">{cellType}</Typography>
                            </Box>
                        ))}
                    </Stack>
                </FilterPanel>
            )}

            <ControlPanel>
                <Typography variant="subtitle2" sx={{ fontSize: '0.6rem', fontWeight: 'small' }}>
                    Sample%
                </Typography>
                <TextField
                    label=" "
                    type="number"
                    value={samplePercent}
                    onChange={handleSampleChange}
                    size="small"
                    inputProps={{
                        min: 0.1, 
                        max: 100, 
                        step: 0.1 
                    }}
                />
            </ControlPanel>
        </Box>
    );
};