async function drawFlowerPlot() {
    // --- Configuration ---
    const config = {
        region1CsvPath: 'region1_NMF.csv',
        region2CsvPath: 'region2_NMF.csv',
        goResult1Path: 'sorted_gene1_counts.json', // Use the sorted gene file
        goResult2Path: 'sorted_gene2_counts.json', // Use the sorted gene file
        geneColorPath: 'gene_colors.json',
        maxClustersPerRegion: 10,
        numComponentsPerCluster: 10,
        topNGenes: 5,
        regionColors: { region1: '#ddd25d', region2: '#2ca02c' }, // Set region1 color same as region2
        unknownGeneColor: 'lightgray',
        componentStrokeColor: '#aaa',
        componentStrokeWidth: 0.5,
        geneStrokeColor: 'white',
        geneStrokeWidth: 0.3,
        petalPaddingAngle: 0.02, // radians, padding between petals
        svgWidth: 1200,
        svgHeight: 700,
        margin: { top: 40, right: 20, bottom: 20, left: 20 },

        // --- Scaling Factor ---
        globalScaleFactor: 1.5,

        // Apply scaling factor to size/distance parameters
        get flowerCenterXOffset() { return 300 * this.globalScaleFactor; },
        get regionCircleRadius() { return 30 * this.globalScaleFactor; },
        get clusterEllipseRx() { return 20 * this.globalScaleFactor; },
        get clusterEllipseRy() { return 35 * this.globalScaleFactor; },
        get componentEllipseBaseRx() { return 3 * this.globalScaleFactor; },
        get componentEllipseMaxRy() { return 15 * this.globalScaleFactor; },
        get componentEllipseMinRy() { return 1 * this.globalScaleFactor; },
        get geneMaxRadiusRelativeToComponentFactor() { return 0.4; }, // Factor itself doesn't scale
        get geneMinAbsoluteRadius() { return 1 * this.globalScaleFactor; },

        // --- Original values before scaling (kept for reference if needed) ---
        // flowerCenterXOffset: 300,
        // regionCircleRadius: 30,
        // clusterEllipseRx: 20,
        // clusterEllipseRy: 35,
        // componentEllipseBaseRx: 3,
        // componentEllipseMaxRy: 15,
        // componentEllipseMinRy: 1,
        // geneMaxRadiusRelativeToComponentFactor: 0.4,
        // geneMinAbsoluteRadius: 1,
    };

    const chartWidth = config.svgWidth - config.margin.left - config.margin.right;
    const chartHeight = config.svgHeight - config.margin.top - config.margin.bottom;
    const flowerCenterY = chartHeight / 2;
    const flowerCenters = {
        region1: { x: chartWidth / 2 - config.flowerCenterXOffset / 2, y: flowerCenterY },
        region2: { x: chartWidth / 2 + config.flowerCenterXOffset / 2, y: flowerCenterY }
    };

    // --- Tooltip Setup ---
    const tooltip = d3.select("body").append("div")
        .attr("class", "tooltip");

    // --- Helper: Process Data ---
    async function loadAndProcessData() {
        const [csvData1, csvData2, sortedGenes1, sortedGenes2, geneColors] = await Promise.all([
            d3.text(config.region1CsvPath),
            d3.text(config.region2CsvPath),
            d3.json(config.goResult1Path),
            d3.json(config.goResult2Path),
            d3.json(config.geneColorPath)
        ]).catch(handleError);

        console.log("Loaded geneColors object:", geneColors);

        const parseWideCsv = (csvText, regionId, sortedGenes) => {
            return d3.csvParse(csvText, (d, i, columns) => {
                const clusterId = d[columns[0]];
                const components = [];
                for (let k = 1; k <= config.numComponentsPerCluster; k++) {
                    const colName = columns[k];
                    const componentId = `Comp${k}`;
                    const value = +d[colName] || 0;
                    
                    const jsonComponentKey = `component_${k}`;
                    let genesData = [];
                    const componentGenesObject = sortedGenes && sortedGenes[jsonComponentKey] ? sortedGenes[jsonComponentKey] : null;

                    let topNRawGenes = [];
                    if (componentGenesObject) {
                        const genesArray = Object.entries(componentGenesObject).map(([gene, count]) => ({ gene, count }));
                        genesArray.sort((a, b) => b.count - a.count);
                        topNRawGenes = genesArray.slice(0, config.topNGenes);
                    }

                    if (topNRawGenes.length > 0) {
                        const maxCount = topNRawGenes[0].count;
                        if (maxCount > 0) {
                            genesData = topNRawGenes.map(g => ({
                                gene: g.gene,
                                count: g.count,
                                normalized_count: g.count / maxCount
                            }));
                        }
                    }

                    if (k < columns.length) {
                        components.push({
                            component_id: componentId,
                            value: value,
                            genes: genesData
                        });
                    } else {
                        console.warn(`Expected component column ${k} not found for cluster ${clusterId}`);
                        components.push({
                            component_id: componentId,
                            value: 0,
                            genes: []
                        });
                    }
                }
                return { 
                    cluster_id: clusterId, 
                    region: regionId,
                    components: components
                };
            });
        };

        const parsedCsv1 = parseWideCsv(csvData1, 'region1', sortedGenes1);
        const parsedCsv2 = parseWideCsv(csvData2, 'region2', sortedGenes2);

        const limitedData1 = parsedCsv1.slice(0, config.maxClustersPerRegion);
        const limitedData2 = parsedCsv2.slice(0, config.maxClustersPerRegion);

        if (limitedData1.length === 0 && limitedData2.length === 0) throw new Error("No cluster data after parsing/limiting CSVs.");

        const hierarchy = { 
            region1: { 
                clusters: limitedData1.sort((a, b) => a.cluster_id.localeCompare(b.cluster_id, undefined, {numeric: true}))
            },
            region2: {
                clusters: limitedData2.sort((a, b) => a.cluster_id.localeCompare(b.cluster_id, undefined, {numeric: true}))
            }
        };

        console.log("Processed Hierarchy (Both Regions, Wide CSV):", hierarchy);
        return { hierarchy, geneColors };
    }

    // --- Helper: Error Handling ---
    function handleError(error) {
        console.error("Error:", error);
        d3.select("#petal-chart").html(`<p style="color:red;">Failed to load or process data: ${error.message}. Check console.</p>`);
        throw error;
    }

    // --- Helper: Petal Arc Generator ---
    // Creates a basic arc path that can be scaled for length
    // NO LONGER USED for main drawing - keeping for potential future use or remove later
    const petalArc = d3.arc()
        .innerRadius(0) // Start from center
        .startAngle(-Math.PI / 16) // Small angle width for petal look
        .endAngle(Math.PI / 16);
        // outerRadius will be set dynamically based on scale
        
    // --- Scales ---
    // Scales need domains set after data processing
    // Re-enable component scale for component ellipse ry
    let componentValueScale = d3.scaleLinear().range([config.componentEllipseMinRy, config.componentEllipseMaxRy]); // Range is now scaled via config getters
    // geneProportionScale is no longer used directly in drawing
    // let geneProportionScale = d3.scaleLinear().range([5, 25]); 

    // --- Gene Size Scale ---
    // Gene size scale removed, now calculated relative to component Ry

    // --- Color Scales ---
    // clusterColorScale removed, using fixed color
    let componentColorScale;

    // --- Main Drawing Logic ---
    try {
        const { hierarchy, geneColors } = await loadAndProcessData();
        console.log("Processed Hierarchy:", hierarchy);

        // Set scale domains based on actual data
        const allComponentValues = Object.values(hierarchy)
                                    .flatMap(r => r.clusters)
                                    .flatMap(c => c.components.map(p => p.value));
        componentValueScale.domain([0, d3.max(allComponentValues) || 1]); // Use 0 as min
        
        // Define color scales based on data
        const componentIds = Array.from({length: config.numComponentsPerCluster}, (_, i) => `Comp${i+1}`);
        componentColorScale = d3.scaleOrdinal(d3.schemeCategory10).domain(componentIds); // Standard categorical for components

        // Gene size scale removed

        // --- SVG Setup ---
        d3.select("#petal-chart svg").remove(); // Clear previous
        const svg = d3.select("#petal-chart")
            .append("svg")
            .attr("width", config.svgWidth)
            .attr("height", config.svgHeight)
            .append("g")
            .attr("transform", `translate(${config.margin.left},${config.margin.top})`);

        // --- Draw Flowers (Both Regions) ---
        for (const regionId in hierarchy) {
            const regionData = hierarchy[regionId];
            if (!regionData || !regionData.clusters || regionData.clusters.length === 0) continue; // Skip if no data for region

            const center = flowerCenters[regionId];
            const regionColor = config.regionColors[regionId];
            const numClusters = regionData.clusters.length;
            const clusterAngleStep = (2 * Math.PI) / numClusters;

            // Main group for the whole visualization, centered
            const mainGroup = svg.append("g")
                 .attr("transform", `translate(${center.x},${center.y})`);
            
            // Re-draw Central Region Circle
            mainGroup.append("circle")
                .attr("r", config.regionCircleRadius)
                .style("fill", regionColor)
                .style("opacity", 0.7);

            // Add text label in the center of the region circle
            mainGroup.append("text")
                .attr("x", 0)
                .attr("y", 0)
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "middle")
                .style("font-size", "10px")
                .style("font-weight", "bold")
                .style("fill", "#333") // Darker color for visibility
                .text(regionId === 'region1' ? 'Reg1' : 'Reg2');

            // Add Region Label (Positioned at the top)
            mainGroup.append("text")
                .attr("x", 0)
                .attr("y", -center.y + config.margin.top / 2)
                .attr("text-anchor", "middle")
                .style("font-size", "16px")
                .style("font-weight", "bold")
                .text(regionId === 'region1' ? 'Region 1' : 'Region 2');

            // Draw Cluster Ellipses and their Components
            regionData.clusters.forEach((cluster, clusterIndex) => {
                const clusterAngle = clusterIndex * clusterAngleStep;

                // Calculate position for the cluster ellipse center to touch the central circle
                const clusterDistance = config.regionCircleRadius + config.clusterEllipseRy; // Distance from center to cluster ellipse center
                const clusterX = Math.cos(clusterAngle - Math.PI / 2) * clusterDistance;
                const clusterY = Math.sin(clusterAngle - Math.PI / 2) * clusterDistance;

                // Create a group for the cluster ellipse and its components
                const clusterGroup = mainGroup.append("g")
                    // Translate to position AND rotate to point outwards
                    .attr("transform", `translate(${clusterX}, ${clusterY}) rotate(${clusterAngle * 180 / Math.PI})`);

                    // Draw the main Cluster Ellipse (fixed size)
                    const clusterFillColor = 'lightskyblue'; // Fixed light blue color
                    clusterGroup.append("ellipse")
                        // Ellipse is now drawn at (0,0) within the rotated group
                        .attr("rx", config.clusterEllipseRx)
                        .attr("ry", config.clusterEllipseRy)
                        .style("fill", clusterFillColor)
                        .style("opacity", 0.5) // Make it semi-transparent
                        .style("stroke", config.componentStrokeColor)
                        .style("stroke-width", config.componentStrokeWidth)
                        .on("mouseover", (event) => {
                            tooltip.transition().duration(200).style("opacity", .9);
                            tooltip.html(`Region: ${regionId}<br>Cluster: ${cluster.cluster_id}`)
                                .style("left", (event.pageX + 5) + "px")
                                .style("top", (event.pageY - 28) + "px");
                        })
                        .on("mouseout", () => {
                            tooltip.transition().duration(500).style("opacity", 0);
                        });

                // Draw Component Ellipses orbiting the cluster ellipse center
                const numComponents = cluster.components.length;
                const componentAngleStep = (2 * Math.PI) / numComponents; // Distribute components evenly

                cluster.components.forEach((component, compIndex) => {
                    if (component.value <= 0) return; // Skip components with no value

                    const componentAngle = compIndex * componentAngleStep;
                    // 1. Calculate position on the boundary of the cluster ellipse
                    const rx = config.clusterEllipseRx;
                    const ry = config.clusterEllipseRy;
                    const cos_a = Math.cos(componentAngle - Math.PI / 2);
                    const sin_a = Math.sin(componentAngle - Math.PI / 2);
                    const r_boundary = 1 / Math.sqrt(Math.pow(cos_a / rx, 2) + Math.pow(sin_a / ry, 2));
                    const boundaryX = r_boundary * cos_a;
                    const boundaryY = r_boundary * sin_a;

                    // 2. Calculate offset inwards (use component rx for simplicity)
                    const offset = config.componentEllipseBaseRx;
                    const length = Math.sqrt(boundaryX * boundaryX + boundaryY * boundaryY); // == r_boundary
                    const unitX = length > 0 ? boundaryX / length : 0;
                    const unitY = length > 0 ? boundaryY / length : 0;

                    // 3. Final position slightly inside the boundary
                    const compX = boundaryX - unitX * offset;
                    const compY = boundaryY - unitY * offset;

                    const componentRy = componentValueScale(component.value);

                    const componentColor = componentColorScale(component.component_id);
                    // Append component ellipse directly to the cluster group
                    clusterGroup.append("ellipse") 
                        .attr("cx", compX)
                        .attr("cy", compY)
                        .attr("rx", config.componentEllipseBaseRx)
                        .attr("ry", componentRy)
                        .style("fill", componentColor)
                        .style("opacity", 0.9)
                        .style("stroke", config.componentStrokeColor)
                        .style("stroke-width", config.componentStrokeWidth * 0.5)
                        // Rotate component ellipse to point outwards from the cluster center
                        .attr("transform", `rotate(${componentAngle * 180 / Math.PI}, ${compX}, ${compY})`)
                        .on("mouseover", (event) => {
                            tooltip.transition().duration(200).style("opacity", .9);
                            // Updated tooltip content
                            tooltip.html(`Cluster: ${cluster.cluster_id}<br>Comp: ${component.component_id}<br>Value: ${component.value.toFixed(3)}`) // Show more precision
                                .style("left", (event.pageX + 5) + "px")
                                .style("top", (event.pageY - 28) + "px");
                        })
                        .on("mouseout", () => {
                            tooltip.transition().duration(500).style("opacity", 0);
                        });

                    // Draw Gene Circles orbiting the component center (compX, compY)
                    const genes = component.genes;
                    if (genes && genes.length > 0) {
                        const numGenes = genes.length;
                        const geneAngleStep = (2 * Math.PI) / numGenes;
                        const componentGeneGroup = clusterGroup.append("g")
                            .attr("class", "component-gene-group") // Add class for easier selection
                            .attr("transform", `translate(${compX}, ${compY})`); // Center gene group on component

                        console.log(`Drawing genes for Comp: ${component.component_id}, Cluster: ${cluster.cluster_id}`); // Log entry into loop
                        genes.forEach((geneData, geneIndex) => {
                            const geneAngle = geneIndex * geneAngleStep;
                            // Calculate position on the boundary of the component ellipse
                            const compRx = config.componentEllipseBaseRx;
                            const compRy = componentRy; // Use the calculated Ry for this component
                            const angleRadGene = geneAngle - Math.PI / 2; // Angle for positioning (0 is top)
                            const cos_a_gene = Math.cos(angleRadGene);
                            const sin_a_gene = Math.sin(angleRadGene);
                            let geneX = 0;
                            let geneY = 0;
                            if (compRx > 0 && compRy > 0) { // Avoid division by zero if rx/ry are 0
                                // Formula for radius to ellipse boundary at a given angle
                                const r_boundary_gene = 1 / Math.sqrt(Math.pow(cos_a_gene / compRx, 2) + Math.pow(sin_a_gene / compRy, 2));
                                geneX = r_boundary_gene * cos_a_gene;
                                geneY = r_boundary_gene * sin_a_gene;
                            }

                            // Calculate gene circle radius relative to component Ry
                            const targetRadius = componentRy * geneData.normalized_count * config.geneMaxRadiusRelativeToComponentFactor;
                            const geneCircleR = Math.max(config.geneMinAbsoluteRadius, targetRadius);

                            const geneColor = geneColors[geneData.gene] || config.unknownGeneColor;

                            if (!geneColors[geneData.gene]) {
                                console.warn(`Color not found for gene: '${geneData.gene}'. Using default.`);
                            }

                            console.log(`  Gene: ${geneData.gene}, NormCount: ${geneData.normalized_count.toFixed(3)}, Radius: ${geneCircleR.toFixed(3)}, X: ${geneX.toFixed(1)}, Y: ${geneY.toFixed(1)}, Color: ${geneColor}`); // Log calculated values

                            componentGeneGroup.append("circle")
                                .attr("cx", geneX)
                                .attr("cy", geneY)
                                .attr("r", geneCircleR)
                                // Optional: Rotate gene circle to point outwards (relative to its angle)
                                // .attr("transform", `rotate(${geneAngle * 180 / Math.PI}, ${geneX}, ${geneY})`) 
                                .style("fill", geneColor)
                                .style("stroke", "white")
                                .style("stroke-width", 0.3)
                                .on("mouseover", (event) => {
                                    tooltip.transition().duration(200).style("opacity", .9);
                                    tooltip.html(`Gene: ${geneData.gene}<br>Count: ${geneData.count}<br>Norm. Count: ${geneData.normalized_count.toFixed(2)}<br>Comp: ${component.component_id}<br>Cluster: ${cluster.cluster_id}`)
                                        .style("left", (event.pageX + 5) + "px")
                                        .style("top", (event.pageY - 28) + "px");
                                })
                                .on("mouseout", () => {
                                    tooltip.transition().duration(500).style("opacity", 0);
                                });
                        });
                    }
                });
            }); // 
        } // 
    
        // --- Draw Legends ---
        drawFlowerLegend(config, config.regionColors, geneColors); 
        console.log("Visualization complete (Circle/Ellipse version - Both Regions).");

    } catch (error) {
        handleError(error); // Log and display error
    }
}

// --- Legend Drawing Function ---
// Updated to accept config and only show region/component legend
function drawFlowerLegend(config, regionColors, geneColorData) {
    const container = d3.select("#legend-container");
    container.html(""); // Clear previous legend
    const itemHeight = 18;
    const padding = 4;
    const symbolSize = itemHeight - padding * 2;
    const textOffset = symbolSize + padding * 2;
    const legendWidth = 200;

    // Region Legend (Only Region 1)
    container.append("h4").text("Regions");
    const regionLegendSvg = container.append("svg")
        // Adjust height for multiple regions
        .attr("height", Object.keys(regionColors).length * itemHeight)
        .attr("width", legendWidth);
        
    // Draw legend items for all defined regions
    Object.entries(regionColors).forEach(([name, color], i) => {
        const g = regionLegendSvg.append("g").attr("transform", `translate(0, ${i * itemHeight})`);
        g.append("rect") // Using rect for region
            .attr("x", padding)
            .attr("y", padding)
            .attr("width", symbolSize)
            .attr("height", symbolSize)
            .style("fill", color)
            // Remove border handling as region1 is no longer white
            .style("stroke", "none"); 
        g.append("text")
            .attr("x", textOffset)
            .attr("y", itemHeight / 2)
            .attr("dy", "0.35em")
            .style("font-size", "11px")
            .text(name === 'region1' ? 'Region 1' : 'Region 2');
    });

    container.append("h4").text("Cluster");
    const clusterLegendSvg = container.append("svg")
        .attr("height", itemHeight)
        .attr("width", legendWidth);
    const gCluster = clusterLegendSvg.append("g");
    gCluster.append("ellipse")
        .attr("cx", padding + config.clusterEllipseRx / 1.5) // Adjust centering for smaller legend ellipse
        .attr("cy", itemHeight / 2)
        .attr("rx", config.clusterEllipseRx / 1.5) // Smaller representation
        .attr("ry", config.clusterEllipseRy / 1.5)
        .style("fill", "lightskyblue") // Show the actual cluster color
        .style("opacity", 0.5)
        .style("stroke", config.componentStrokeColor)
        .style("stroke-width", config.componentStrokeWidth);
    gCluster.append("text")
        // Adjust text position based on smaller ellipse
        .attr("x", config.clusterEllipseRx / 1.5 * 2 + padding * 2) 
        .attr("y", itemHeight / 2)
        .attr("dy", "0.35em")
        .style("font-size", "11px")
        .text("Cluster Ellipse (fixed color)");

    container.append("h4").text("Component");
    const compLegendSvg = container.append("svg")
        .attr("height", itemHeight)
        .attr("width", legendWidth);
    const gComp = compLegendSvg.append("g");
     gComp.append("ellipse")
            .attr("cx", padding + config.componentEllipseBaseRx)
            .attr("cy", itemHeight / 2)
            .attr("rx", config.componentEllipseBaseRx)
            .attr("ry", (config.componentEllipseMinRy + config.componentEllipseMaxRy) / 2.5) // Example height
            .style("fill", d3.schemeCategory10[0]) // Show first color from the scale
            .style("opacity", 0.9)
            .style("stroke", config.componentStrokeColor)
            .style("stroke-width", config.componentStrokeWidth * 0.5);
     gComp.append("text")
            .attr("x", textOffset)
            .attr("y", itemHeight / 2)
            .attr("dy", "0.35em")
            .style("font-size", "11px")
            .text("Ellipse (color by Comp ID, height ~ value)");

    // --- Gene Legend ---
    container.append("h4").text("Gene");

    // Gene Circle Example
    const geneExampleSvg = container.append("svg")
        .attr("height", itemHeight)
        .attr("width", legendWidth);
    const gGeneCircle = geneExampleSvg.append("g");
    gGeneCircle.append("circle")
        .attr("cx", padding + config.geneMinAbsoluteRadius + 1.5)
        .attr("cy", itemHeight / 2)
        .attr("r", config.geneMinAbsoluteRadius + 1.5)
        .style("fill", "#ddd") // Generic color
        .style("stroke", "white")
        .style("stroke-width", 0.3);
     gGeneCircle.append("text")
        .attr("x", textOffset + (config.geneMinAbsoluteRadius + 1.5) - symbolSize)
        .attr("y", itemHeight / 2)
        .attr("dy", "0.35em")
        .style("font-size", "11px")
        .text("Circle (size ~ norm. count relative to component)");

    // Gene Colors
    container.append("h5").style("margin-top", "8px").style("margin-bottom", "2px").text("Colors:");
    const genes = typeof geneColorData === 'object' && geneColorData !== null ? Object.keys(geneColorData).sort() : [];
    const maxGenesInLegend = 15; // Limit number of genes shown in legend if too many
    const genesToShow = genes.slice(0, maxGenesInLegend);

    const geneLegendSvg = container.append("svg")
        .attr("height", genesToShow.length * itemHeight + (genes.length > maxGenesInLegend ? itemHeight : 0))
        .attr("width", legendWidth);

    genesToShow.forEach((gene, i) => {
        const color = geneColorData[gene];
        const g = geneLegendSvg.append("g").attr("transform", `translate(0, ${i * itemHeight})`);
        g.append("rect")
            .attr("x", padding)
            .attr("y", padding)
            .attr("width", symbolSize)
            .attr("height", symbolSize)
            .style("fill", color || config.unknownGeneColor);
        g.append("text")
            .attr("x", textOffset)
            .attr("y", itemHeight / 2)
            .attr("dy", "0.35em")
            .style("font-size", "11px")
            .text(gene);
    });
    if (genes.length > maxGenesInLegend) {
        geneLegendSvg.append("text")
           .attr("x", padding)
           .attr("y", genesToShow.length * itemHeight + itemHeight / 2)
            .attr("dy", "0.35em")
           .style("font-size", "10px")
           .style("font-style", "italic")
           .text(`...and ${genes.length - maxGenesInLegend} more`);
    }
}

// Run the visualization
drawFlowerPlot();