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
    bottom: theme.spacing(1),
    right: theme.spacing(1),
    display: 'flex',
    flexDirection: 'column',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
}));

const ToolPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    top: theme.spacing(1),
    right: theme.spacing(1),
    display: 'flex',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(0.5),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
}));

const UMAPContainer = styled(Box)(({ theme }) => ({
    width: '80%',
    height: '100%',
    position: 'relative',
    backgroundColor: '#fff',
    borderRadius: theme.shape.borderRadius,
    overflow: 'hidden',
    margin: '0 auto',
}));

const FilterPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    top: theme.spacing(1),
    left: theme.spacing(1),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
    maxHeight: '200px',
    overflowY: 'auto',
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

// تعریف رنگ‌های پیش‌فرض
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
    const [showColorPicker, setShowColorPicker] = useState(false);
    const [customColors, setCustomColors] = useState({});
    const svgRef = useRef(null);
    const zoomRef = useRef(null);

    useEffect(() => {
        fetchUMAPData();
    }, [samplePercent]);

    const fetchUMAPData = async () => {
        try {
            setLoading(true);
            setError(null);
            
            const response = await axios.get('/get_deaplog_results', {
                params: {
                    sample_percent: samplePercent
                }
            });
            
            setUmapData(response.data);
            // تنظیم تمام انواع سلول‌ها به عنوان انتخاب شده
            setSelectedCellTypes(response.data.metadata.unique_cell_types);
        } catch (err) {
            console.error('Error fetching UMAP data:', err);
            setError(err.response?.data?.error || 'Failed to fetch UMAP data');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        if (!umapData || !svgRef.current) return;

        // پاک کردن SVG قبلی
        d3.select(svgRef.current).selectAll("*").remove();

        const width = svgRef.current.clientWidth;
        const height = svgRef.current.clientHeight;
        const margin = { top: 20, right: 20, bottom: 30, left: 40 };

        const svg = d3.select(svgRef.current)
            .attr("width", width)
            .attr("height", height);

        // ایجاد مقیاس‌ها
        const xScale = d3.scaleLinear()
            .domain(d3.extent(umapData.umap_coordinates.x))
            .range([margin.left, width - margin.right]);

        const yScale = d3.scaleLinear()
            .domain(d3.extent(umapData.umap_coordinates.y))
            .range([height - margin.bottom, margin.top]);

        // ایجاد رنگ‌ها برای انواع سلول‌ها
        const colorScale = d3.scaleOrdinal()
            .domain(umapData.metadata.unique_cell_types)
            .range(d3.schemeCategory10);

        const zoom = d3.zoom()
            .scaleExtent([0.1, 10])
            .on("zoom", (event) => {
                const transform = event.transform;
                transform.x = 0; 
                g.attr("transform", transform);
            });

        svg.call(zoom);

        const g = svg.append("g");

        const points = g.selectAll("circle")
            .data(umapData.umap_coordinates.x)
            .enter()
            .append("circle")
            .attr("cx", (d, i) => xScale(d))
            .attr("cy", (d, i) => yScale(umapData.umap_coordinates.y[i]))
            .attr("r", 3)
            .attr("fill", (d, i) => {
                const cellType = umapData.cell_types[i];
                return selectedCellTypes.includes(cellType) 
                    ? (customColors[cellType] || colorScale(cellType))
                    : "#ddd";
            })
            .attr("opacity", (d, i) => selectedCellTypes.includes(umapData.cell_types[i]) ? 0.6 : 0.1)
            .on("mouseover", function(event, d) {
                const cellType = umapData.cell_types[event.target.__data__];
                if (selectedCellTypes.includes(cellType)) {
                    d3.select(this)
                        .attr("r", 5)
                        .attr("opacity", 1);
                }
            })
            .on("mouseout", function(event, d) {
                const cellType = umapData.cell_types[event.target.__data__];
                if (selectedCellTypes.includes(cellType)) {
                    d3.select(this)
                        .attr("r", 3)
                        .attr("opacity", 0.6);
                }
            });

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

        points.on("mouseover", function(event, d) {
            const cellType = umapData.cell_types[event.target.__data__];
            if (selectedCellTypes.includes(cellType)) {
                const cellInfo = {
                    type: cellType,
                    index: event.target.__data__,
                    x: umapData.umap_coordinates.x[event.target.__data__].toFixed(2),
                    y: umapData.umap_coordinates.y[event.target.__data__].toFixed(2)
                };

                if (umapData.gene_expression) {
                    const topGenes = Object.entries(umapData.gene_expression)
                        .map(([gene, values]) => ({
                            gene,
                            value: values[event.target.__data__]
                        }))
                        .sort((a, b) => b.value - a.value)
                        .slice(0, 3);

                    cellInfo.topGenes = topGenes;
                }

                tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                let tooltipHtml = `
                    <div style="font-weight: bold; margin-bottom: 4px;">Cell Information</div>
                    <div>Type: ${cellInfo.type}</div>
                    <div>Index: ${cellInfo.index}</div>
                    <div>UMAP1: ${cellInfo.x}</div>
                    <div>UMAP2: ${cellInfo.y}</div>
                `;

                if (cellInfo.topGenes) {
                    tooltipHtml += `
                        <div style="margin-top: 8px; font-weight: bold;">Top Expressed Genes:</div>
                        ${cellInfo.topGenes.map(g => `
                            <div>${g.gene}: ${g.value.toFixed(2)}</div>
                        `).join('')}
                    `;
                }

                tooltip.html(tooltipHtml)
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 28) + "px");
            }
        })
        .on("mouseout", function(event, d) {
            tooltip.transition()
                .duration(500)
                .style("opacity", 0);
        });

        g.append("g")
            .attr("transform", `translate(0,${height - margin.bottom})`)
            .call(d3.axisBottom(xScale));

        g.append("g")
            .attr("transform", `translate(${margin.left},0)`)
            .call(d3.axisLeft(yScale));

        g.append("text")
            .attr("x", width / 2)
            .attr("y", margin.top / 2)
            .attr("text-anchor", "middle")
            .text("Pseudo-Temporal UMAP");

        const legend = g.append("g")
            .attr("font-family", "sans-serif")
            .attr("font-size", 10)
            .attr("text-anchor", "start")
            .selectAll("g")
            .data(umapData.metadata.unique_cell_types)
            .enter().append("g")
            .attr("transform", (d, i) => `translate(0,${i * 20})`);

        legend.append("rect")
            .attr("x", width - margin.right + 10)
            .attr("width", 19)
            .attr("height", 19)
            .attr("fill", d => customColors[d] || colorScale(d));

        legend.append("text")
            .attr("x", width - margin.right + 35)
            .attr("y", 9.5)
            .attr("dy", "0.32em")
            .text(d => d);

        zoomRef.current = zoom;

    }, [umapData, selectedCellTypes, customColors]);

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

    const handleColorChange = (cellType, color) => {
        setCustomColors(prev => ({
            ...prev,
            [cellType]: color
        }));
    };

    const getDefaultColor = (index) => {
        return defaultColors[index % defaultColors.length];
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
                <Tooltip title="Zoom In">
                    <IconButton onClick={handleZoomIn} size="small">
                        <ZoomInIcon />
                    </IconButton>
                </Tooltip>
                <Tooltip title="Zoom Out">
                    <IconButton onClick={handleZoomOut} size="small">
                        <ZoomOutIcon />
                    </IconButton>
                </Tooltip>
                <Tooltip title="Pan Mode">
                    <IconButton onClick={() => setShowFilters(!showFilters)} size="small">
                        <PanToolIcon />
                    </IconButton>
                </Tooltip>
                <Tooltip title="Filter Cell Types">
                    <IconButton onClick={() => setShowFilters(!showFilters)} size="small">
                        <FilterListIcon />
                    </IconButton>
                </Tooltip>
                <Tooltip title="Customize Colors">
                    <IconButton onClick={() => setShowColorPicker(!showColorPicker)} size="small">
                        <PaletteIcon />
                    </IconButton>
                </Tooltip>
            </ToolPanel>

            {showFilters && umapData && (
                <FilterPanel>
                    <Typography variant="subtitle2" sx={{ mb: 1 }}>Cell Types</Typography>
                    <Stack spacing={1}>
                        {umapData.metadata.unique_cell_types.map(cellType => (
                            <Box key={cellType} sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
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

            {showColorPicker && umapData && (
                <ColorPickerPanel>
                    <Typography variant="subtitle2" sx={{ mb: 1 }}>Customize Colors</Typography>
                    <Stack spacing={1}>
                        {umapData.metadata.unique_cell_types.map((cellType, index) => (
                            <Box key={cellType} sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                                <Typography variant="body2">{cellType}</Typography>
                                <input
                                    type="color"
                                    value={customColors[cellType] || getDefaultColor(index)}
                                    onChange={(e) => handleColorChange(cellType, e.target.value)}
                                />
                            </Box>
                        ))}
                    </Stack>
                </ColorPickerPanel>
            )}

            <ControlPanel>
                <TextField
                    label="Sample %"
                    type="number"
                    value={samplePercent}
                    onChange={handleSampleChange}
                    size="small"
                    inputProps={{ 
                        min: 0.1, 
                        max: 100, 
                        step: 0.1 
                    }}
                    sx={{ width: 80 }}
                />
            </ControlPanel>
        </Box>
    );
};