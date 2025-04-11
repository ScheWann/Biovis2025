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
    bottom: 20,
    left: 0,
    display: 'flex',
    flexDirection: 'column',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 10,
    '& .MuiTextField-root': {
        width: '75px',
        '& .MuiInputBase-root': {
            height: '20px',
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
    zIndex: 10,
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
    zIndex: 10,
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
    zIndex: 10,
}));

const defaultColors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
];

export const PseudoTemporalViewer = () => {
    const [samplePercent, setSamplePercent] = useState(0.001);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [umapData, setUmapData] = useState(null);
    const [selectedCellTypes, setSelectedCellTypes] = useState([]);
    const [showFilters, setShowFilters] = useState(false);
    const [showArrows, setShowArrows] = useState(true);
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
            
            console.log('Raw response data:', response.data);
            
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

        // Create containers
        const staticContainer = svg.append("g");
        const interactiveContainer = svg.append("g");
        const arrowsContainer = svg.append("g");  // Container for arrows

        // Add static elements
        staticContainer.append("text")
            .attr("x", width / 2)
            .attr("y", height - 25)
            .attr("text-anchor", "middle")
            .text("UMAP1")
            .style("font-size", "12px")
            .style("font-weight", "bold");

        staticContainer.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -height / 2)
            .attr("y", 30)
            .attr("text-anchor", "middle")
            .text("UMAP2")
            .style("font-size", "12px")
            .style("font-weight", "bold");

        staticContainer.append("text")
            .attr("x", width / 2)
            .attr("y", margin.top / 2)
            .attr("text-anchor", "middle")
            .attr("font-size", "14px")
            .attr("font-weight", "bold")
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

        // Add zoom behavior to the entire SVG
        const zoom = d3.zoom()
            .scaleExtent([0.1, 10])
            .on("zoom", (event) => {
                interactiveContainer.attr("transform", event.transform);
                arrowsContainer.attr("transform", event.transform);
            });

        svg.call(zoom);

        // Add pan behavior to the entire SVG
        let isPanning = false;
        let startX, startY;

        svg.on("mousedown", function(event) {
            isPanning = true;
            startX = event.clientX;
            startY = event.clientY;
        });

        svg.on("mousemove", function(event) {
            if (!isPanning) return;
            
            const dx = event.clientX - startX;
            const dy = event.clientY - startY;
            
            const currentTransform = d3.zoomTransform(svg.node());
            const newTransform = currentTransform.translate(dx / currentTransform.k, dy / currentTransform.k);
            
            svg.transition()
                .duration(0)
                .call(zoom.transform, newTransform);
            
            startX = event.clientX;
            startY = event.clientY;
        });

        svg.on("mouseup", function() {
            isPanning = false;
        });

        svg.on("mouseleave", function() {
            isPanning = false;
        });

        // Add points to interactive container
        const points = interactiveContainer.selectAll("circle")
            .data(umapData.umap_coordinates.x)
            .enter()
            .append("circle")
            .attr("cx", (d, i) => xScale(d))
            .attr("cy", (d, i) => yScale(umapData.umap_coordinates.y[i]))
            .attr("r", 2)  // کاهش اندازه نقاط
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
                        .attr("r", 3)
                        .attr("opacity", 1);

                    // Create tooltip div if it doesn't exist
                    let tooltip = d3.select(".tooltip");
                    if (tooltip.empty()) {
                        tooltip = d3.select("body")
                            .append("div")
                            .attr("class", "tooltip")
                            .style("position", "absolute")
                            .style("background-color", "white")
                            .style("border", "1px solid #ddd")
                            .style("border-radius", "4px")
                            .style("padding", "8px")
                            .style("pointer-events", "none")
                            .style("box-shadow", "0 2px 4px rgba(0,0,0,0.1)")
                            .style("font-size", "12px")
                            .style("z-index", "1000");
                    }

                    // Get the first 3 genes for this cell
                    const genes = umapData.gene_names.slice(0, 3);
                    
                    tooltip
                        .style("opacity", 1)
                        .html(`
                            <strong>Cell Type:</strong> ${cellType}<br>
                            <strong>Top 3 Genes:</strong><br>
                            ${genes.map((gene, i) => `${i + 1}. ${gene}`).join('<br>')}
                        `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }
            })
            .on("mousemove", function(event, d) {
                const tooltip = d3.select(".tooltip");
                if (!tooltip.empty()) {
                    tooltip
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }
            })
            .on("mouseout", function(event, d) {
                const cellType = umapData.cell_types[event.target.__data__];
                if (selectedCellTypes.includes(cellType)) {
                    d3.select(this)
                        .attr("r", 2)
                        .attr("opacity", 0.6);
                    
                    d3.select(".tooltip")
                        .style("opacity", 0)
                        .remove();
                }
            });

        // Add arrows for each cell type if showArrows is true
        if (showArrows) {
            // Group data by cell type
            const cellTypeData = {};
            umapData.cell_types.forEach((cellType, i) => {
                if (!selectedCellTypes.includes(cellType)) return;
                
                if (!cellTypeData[cellType]) {
                    cellTypeData[cellType] = {
                        x: [],
                        y: [],
                        pseudotime: []
                    };
                }
                cellTypeData[cellType].x.push(umapData.umap_coordinates.x[i]);
                cellTypeData[cellType].y.push(umapData.umap_coordinates.y[i]);
                cellTypeData[cellType].pseudotime.push(umapData.pseudotime[i]);
            });

            // Calculate arrows for each cell type
            Object.entries(cellTypeData).forEach(([cellType, data]) => {
                // Sort by pseudotime
                const sortedIndices = data.pseudotime
                    .map((time, i) => ({ time, i }))
                    .sort((a, b) => a.time - b.time)
                    .map(item => item.i);

                // Get points for three arrows
                const totalPoints = sortedIndices.length;
                const points = [
                    {
                        start: {
                            x: data.x[sortedIndices[0]],
                            y: data.y[sortedIndices[0]]
                        },
                        end: {
                            x: data.x[sortedIndices[Math.floor(totalPoints * 0.3)]],
                            y: data.y[sortedIndices[Math.floor(totalPoints * 0.3)]]
                        }
                    },
                    {
                        start: {
                            x: data.x[sortedIndices[Math.floor(totalPoints * 0.3)]],
                            y: data.y[sortedIndices[Math.floor(totalPoints * 0.3)]]
                        },
                        end: {
                            x: data.x[sortedIndices[Math.floor(totalPoints * 0.6)]],
                            y: data.y[sortedIndices[Math.floor(totalPoints * 0.6)]]
                        }
                    },
                    {
                        start: {
                            x: data.x[sortedIndices[Math.floor(totalPoints * 0.6)]],
                            y: data.y[sortedIndices[Math.floor(totalPoints * 0.6)]]
                        },
                        end: {
                            x: data.x[sortedIndices[totalPoints - 1]],
                            y: data.y[sortedIndices[totalPoints - 1]]
                        }
                    }
                ];

                // Draw three arrows for each cell type
                points.forEach(({ start, end }) => {
                    // Calculate arrow properties
                    const dx = end.x - start.x;
                    const dy = end.y - start.y;
                    const angle = Math.atan2(dy, dx);

                    // Draw arrow line
                    arrowsContainer.append("line")
                        .attr("x1", xScale(start.x))
                        .attr("y1", yScale(start.y))
                        .attr("x2", xScale(end.x))
                        .attr("y2", yScale(end.y))
                        .attr("stroke", "black")
                        .attr("stroke-width", 1)
                        .attr("opacity", 0.8);

                    // Draw arrow head
                    arrowsContainer.append("path")
                        .attr("d", `M ${xScale(end.x)} ${yScale(end.y)} 
                                  L ${xScale(end.x - 0.3 * Math.cos(angle - Math.PI / 6))} 
                                    ${yScale(end.y - 0.3 * Math.sin(angle - Math.PI / 6))}
                                  L ${xScale(end.x - 0.3 * Math.cos(angle + Math.PI / 6))} 
                                    ${yScale(end.y - 0.3 * Math.sin(angle + Math.PI / 6))}
                                  Z`)
                        .attr("fill", "black")
                        .attr("opacity", 0.8);
                });
            });
        }

        zoomRef.current = zoom;

    }, [umapData, selectedCellTypes, showArrows]);

    const handleSampleChange = (event) => {
        const value = parseFloat(event.target.value) || 0.001;
        const validValues = [0.001, 0.01, 0.1, 1];
        const currentIndex = validValues.indexOf(samplePercent);
        let newValue;
        
        if (value > samplePercent) {
            // افزایش مقدار
            newValue = currentIndex < validValues.length - 1 ? validValues[currentIndex + 1] : validValues[validValues.length - 1];
        } else if (value < samplePercent) {
            // کاهش مقدار
            newValue = currentIndex > 0 ? validValues[currentIndex - 1] : validValues[0];
        } else {
            newValue = value;
        }
        
        setSamplePercent(newValue);
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
                <Tooltip title="Toggle Arrows">
                    <IconButton onClick={() => setShowArrows(!showArrows)} size="small" sx={{ padding: '3px' }}>
                        <PanToolIcon sx={{ fontSize: '1.0rem' }} />
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