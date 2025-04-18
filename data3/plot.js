// --- Configuration ---
// Region 1
const NMF_FILENAME_R1 = 'region1_NMF.csv';
const GENE_COUNTS_FILENAME_R1 = 'sorted_gene1_counts.json';
const CLUSTER_COLOR_R1 = 'blue';
// Region 2
const NMF_FILENAME_R2 = 'region2_NMF.csv'; // New
const GENE_COUNTS_FILENAME_R2 = 'sorted_gene2_counts.json'; // New
const CLUSTER_COLOR_R2 = 'green'; // New

const NUM_CLUSTERS = 8; // Per region
const NUM_COMPONENTS = 10;
const NUM_TOP_GENES_OVERALL = 10; // For distinct legend colors (across both regions)
const NUM_TOP_GENES_PER_COMPONENT = 5; // Genes considered within each component
const DEFAULT_COMPONENT_COLOR = 'lightgrey';

// --- Sizing Configuration ---
const VISUAL_SCALE = 50; // Factor to make conceptual units (0-1) visually significant
const CLUSTER_TARGET_RADIUS = 1 * VISUAL_SCALE;   // Target radius 1, scaled
const COMPONENT_MAX_RADIUS_PIXELS = 1 * VISUAL_SCALE; // Max pixel radius when NMF = 1

const CLUSTER_GRID_COLS = 4; 
const CLUSTER_SPACING = 250; 

// --- D3 Setup ---
const margin = { top: 50, right: 200, bottom: 50, left: 50 }; 
// Calculate height needed for two regions
const rowsPerRegion = Math.ceil(NUM_CLUSTERS / CLUSTER_GRID_COLS);
const chartHeightNeeded = (rowsPerRegion * 2) * CLUSTER_SPACING; // Height for two sets of rows

const svgWidth = CLUSTER_GRID_COLS * CLUSTER_SPACING + margin.left + margin.right;
const svgHeight = chartHeightNeeded + margin.top + margin.bottom; 

const svg = d3.select("#chart").append("svg")
    .attr("width", svgWidth)
    .attr("height", svgHeight)
  .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

const chartWidth = svgWidth - margin.left - margin.right;
// const chartHeight = svgHeight - margin.top - margin.bottom; // Not strictly needed now

const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip");

// --- Data Loading (Load all 4 files) --- 
Promise.all([
    d3.csv(NMF_FILENAME_R1),
    d3.json(GENE_COUNTS_FILENAME_R1),
    d3.csv(NMF_FILENAME_R2),       // Load R2 NMF
    d3.json(GENE_COUNTS_FILENAME_R2) // Load R2 Genes
]).then(([nmfDataCsvR1, geneDataJsonR1, nmfDataCsvR2, geneDataJsonR2]) => {

    console.log("R1 NMF Data Loaded:", nmfDataCsvR1);
    console.log("R1 Gene Data Loaded:", geneDataJsonR1);
    console.log("R2 NMF Data Loaded:", nmfDataCsvR2);
    console.log("R2 Gene Data Loaded:", geneDataJsonR2);

    // --- Data Processing ---
    let maxNmfValue = -Infinity;
    let minNmfValue = Infinity;

    // Helper function to process NMF data
    function processNmf(nmfDataCsv) {
        const lookup = {};
        nmfDataCsv.forEach(row => {
            const firstColumnKey = Object.keys(row)[0]; 
            if (firstColumnKey === undefined || row[firstColumnKey] === undefined || row[firstColumnKey] === null) return; 
            const clusterKey = String(row[firstColumnKey]).trim(); 
            if (!clusterKey) return; 
            lookup[clusterKey] = {};
            for (let i = 1; i <= NUM_COMPONENTS; i++) {
                const compKey = `Component_${i}`;
                const value = +row[compKey]; 
                if (!isNaN(value)) {
                    const clampedValue = Math.max(0, Math.min(1, value));
                    lookup[clusterKey][compKey] = clampedValue; 
                    maxNmfValue = Math.max(maxNmfValue, clampedValue);
                    minNmfValue = Math.min(minNmfValue, clampedValue);
                } else {
                     lookup[clusterKey][compKey] = 0; 
                     minNmfValue = Math.min(minNmfValue, 0);
                }
            }
        });
        return lookup;
    }

    // Process both NMF datasets and find combined range
    const nmfLookupR1 = processNmf(nmfDataCsvR1);
    const nmfLookupR2 = processNmf(nmfDataCsvR2);

    console.log("Processed R1 NMF Lookup:", nmfLookupR1);
    console.log("Processed R2 NMF Lookup:", nmfLookupR2);
    if (minNmfValue === Infinity) minNmfValue = 0;
    if (maxNmfValue === -Infinity) maxNmfValue = 1;
    if (minNmfValue >= maxNmfValue) maxNmfValue = minNmfValue + 1e-6; 
    console.log(`Combined NMF Value Range: [${minNmfValue}, ${maxNmfValue}]`);
        
    // --- Create Bubble RADIUS Scale (Based on COMBINED Min/Max NMF) ---
    const radiusScale = d3.scaleLinear()
        .domain([0, 1]) // NMF values are clamped 0-1
        .range([0, COMPONENT_MAX_RADIUS_PIXELS]); // Output SVG radius range

    // --- Identify Top Overall Genes & Create Color Scale (COMBINED from R1 & R2) ---
    const overallGeneCounts = {};
    function countTopGenes(geneData) {
         for (let i = 1; i <= NUM_COMPONENTS; i++) {
            const compKeyJson = `component_${i}`; 
            if (geneData && geneData[compKeyJson]) {
                const componentGenes = Object.entries(geneData[compKeyJson])
                                          .slice(0, NUM_TOP_GENES_PER_COMPONENT);
                componentGenes.forEach(([gene, count]) => {
                    overallGeneCounts[gene] = (overallGeneCounts[gene] || 0) + count;
                });
            }
        }
    }
    countTopGenes(geneDataJsonR1);
    countTopGenes(geneDataJsonR2);

    const sortedOverallGenes = Object.entries(overallGeneCounts)
                                  .sort((a, b) => b[1] - a[1]) 
                                  .slice(0, NUM_TOP_GENES_OVERALL) 
                                  .map(d => d[0]); 
    console.log("Top Overall Genes (Combined R1 & R2):", sortedOverallGenes);

    const colorScale = d3.scaleOrdinal(d3.schemeCategory10)
        .domain(sortedOverallGenes);
        
    // --- Pre-process Gene Data for Pie Charts (for R1 & R2) ---
    function processGeneData(geneData) {
        const processedData = {};
        for (let i = 1; i <= NUM_COMPONENTS; i++) {
            const compKeyJson = `component_${i}`; 
            processedData[compKeyJson] = []; 
            if (geneData && geneData[compKeyJson] && Object.keys(geneData[compKeyJson]).length > 0) {
                const topGenes = Object.entries(geneData[compKeyJson])
                                    .slice(0, NUM_TOP_GENES_PER_COMPONENT); 
                processedData[compKeyJson] = topGenes.map(([gene, count]) => ({ gene: gene, count: count }));
            }
        }
        return processedData;
    }
    const processedGeneDataR1 = processGeneData(geneDataJsonR1);
    const processedGeneDataR2 = processGeneData(geneDataJsonR2);
    console.log("Processed R1 Gene Data for Pies:", processedGeneDataR1);
    console.log("Processed R2 Gene Data for Pies:", processedGeneDataR2);

    // --- Calculate Layout Positions for Clusters --- 
    const clusterPositionsR1 = [];
    const clusterPositionsR2 = [];
    const verticalOffsetR2 = rowsPerRegion * CLUSTER_SPACING; // Start R2 below R1

    for (let i = 0; i < NUM_CLUSTERS; i++) { 
         // Region 1 positions
         clusterPositionsR1.push({
            id: String(i), 
            region: 1,
            x: (i % CLUSTER_GRID_COLS) * CLUSTER_SPACING + CLUSTER_SPACING / 2, 
            y: Math.floor(i / CLUSTER_GRID_COLS) * CLUSTER_SPACING + CLUSTER_SPACING / 2
         });
         // Region 2 positions (offset vertically)
         clusterPositionsR2.push({
            id: String(i),
            region: 2,
            x: (i % CLUSTER_GRID_COLS) * CLUSTER_SPACING + CLUSTER_SPACING / 2, 
            y: Math.floor(i / CLUSTER_GRID_COLS) * CLUSTER_SPACING + CLUSTER_SPACING / 2 + verticalOffsetR2
         });
    }

    // --- Drawing Function (to avoid duplicating code) ---
    function drawClustersAndComponents(clusterData, nmfLookup, processedGeneData, clusterColor, regionId) {
        const clusters = svg.append("g")
            .attr("class", `clusters region-${regionId}`)
            .selectAll("g")
            .data(clusterData)
            .join("g") 
            .attr("transform", d => `translate(${d.x},${d.y})`);

        clusters.append("circle")
            .attr("r", CLUSTER_TARGET_RADIUS)
            .attr("fill", clusterColor) // Use region-specific color
            .attr("opacity", 0.7);
            
        clusters.append("text")
            .attr("y", CLUSTER_TARGET_RADIUS + 12) 
            .attr("text-anchor", "middle")
            .attr("font-weight", "bold")
            .style("font-size", "12px")
            .text(d => `R${regionId}-C${d.id}`); // Add region prefix

        clusters.each(function(currentClusterData) {
            const clusterGroup = d3.select(this); 
            const clusterKeyNmf = currentClusterData.id; 

            for (let i = 0; i < NUM_COMPONENTS; i++) {
                const angle = (2 * Math.PI / NUM_COMPONENTS) * i - (Math.PI / 2); 
                const compKeyNmf = `Component_${i + 1}`; 
                const compKeyJson = `component_${i + 1}`; 
                
                const nmfValue = nmfLookup[clusterKeyNmf]?.[compKeyNmf] ?? 0; 
                const componentRadius = radiusScale(nmfValue); 

                if (componentRadius < 1) continue; 

                const dynamicOffset = CLUSTER_TARGET_RADIUS + componentRadius; 
                const cx = dynamicOffset * Math.cos(angle);
                const cy = dynamicOffset * Math.sin(angle);
                const pieData = processedGeneData[compKeyJson];

                if (pieData && pieData.length > 0) {
                    const pieGroup = clusterGroup.append("g")
                        .attr("transform", `translate(${cx},${cy})`);
                    const pie = d3.pie().value(d => d.count).sort(null);
                    const arc = d3.arc().innerRadius(0).outerRadius(componentRadius);

                    pieGroup.selectAll('path')
                        .data(pie(pieData))
                        .join('path')
                            .attr('d', arc)
                            .attr('fill', d => colorScale.domain().includes(d.data.gene) ? colorScale(d.data.gene) : DEFAULT_COMPONENT_COLOR)
                            .attr("stroke", "white")
                            .style("stroke-width", "0.5px")
                            .style("cursor", "pointer")
                            .on("mouseover", function(event, d) {
                                d3.select(this).attr("stroke", "black").attr("stroke-width", 1.5);
                                tooltip.transition().duration(200).style("opacity", .9);
                                tooltip.html(`Region ${regionId} | Cluster: ${clusterKeyNmf} | Comp: ${compKeyNmf}<br>Gene: <strong>${d.data.gene}</strong><br>Count: ${d.data.count}<br><span style='font-size:10px;'>(Comp NMF: ${nmfValue.toFixed(3)})</span>`) 
                                    .style("left", (event.pageX + 10) + "px") 
                                    .style("top", (event.pageY - 28) + "px");
                            })
                            .on("mouseout", function(d) { 
                                d3.select(this).attr("stroke", "white").style("stroke-width", "0.5px");
                                tooltip.transition().duration(500).style("opacity", 0);
                            });
                } else {
                     clusterGroup.append("circle") // Fallback circle
                        .attr("cx", cx)
                        .attr("cy", cy)
                        .attr("r", componentRadius) 
                        .attr("fill", DEFAULT_COMPONENT_COLOR)
                        .attr("opacity", 0.7)
                        .on("mouseover", function(event, d) { 
                            tooltip.transition().duration(200).style("opacity", .9);
                            tooltip.html(`Region ${regionId} | Cluster: ${clusterKeyNmf}<br>Component: ${compKeyNmf}<br>NMF Value: ${nmfValue.toFixed(3)}<br>---No Gene Data---`) 
                                .style("left", (event.pageX + 10) + "px") 
                                .style("top", (event.pageY - 28) + "px");
                        })
                        .on("mouseout", function(d) { 
                            tooltip.transition().duration(500).style("opacity", 0);
                        });
                }
            }
        });
    }

    // --- Draw Both Regions ---
    drawClustersAndComponents(clusterPositionsR1, nmfLookupR1, processedGeneDataR1, CLUSTER_COLOR_R1, 1);
    drawClustersAndComponents(clusterPositionsR2, nmfLookupR2, processedGeneDataR2, CLUSTER_COLOR_R2, 2);
    
    // --- Drawing Legends --- 
    const legendX = chartWidth + 40; 
    let legendY = 0; 
    const legendGroup = svg.append("g").attr("class", "legends")
                           .attr("transform", `translate(${legendX}, ${legendY})`);

    // Cluster Legends
    const clusterLegendR1 = legendGroup.append("g").attr("class", "legend-item");
    clusterLegendR1.append("circle").attr("cx", CLUSTER_TARGET_RADIUS).attr("cy", CLUSTER_TARGET_RADIUS).attr("r", CLUSTER_TARGET_RADIUS).style("fill", CLUSTER_COLOR_R1).attr("opacity", 0.7);
    clusterLegendR1.append("text").attr("x", CLUSTER_TARGET_RADIUS * 2 + 10).attr("y", CLUSTER_TARGET_RADIUS).text("Region 1 Cluster (R=1)").attr("alignment-baseline","middle");
    legendY += CLUSTER_TARGET_RADIUS * 2 + 10;
    
    const clusterLegendR2 = legendGroup.append("g").attr("class", "legend-item").attr("transform", `translate(0, ${legendY})`); // Offset R2 legend
    clusterLegendR2.append("circle").attr("cx", CLUSTER_TARGET_RADIUS).attr("cy", CLUSTER_TARGET_RADIUS).attr("r", CLUSTER_TARGET_RADIUS).style("fill", CLUSTER_COLOR_R2).attr("opacity", 0.7);
    clusterLegendR2.append("text").attr("x", CLUSTER_TARGET_RADIUS * 2 + 10).attr("y", CLUSTER_TARGET_RADIUS).text("Region 2 Cluster (R=1)").attr("alignment-baseline","middle");
    legendY += CLUSTER_TARGET_RADIUS * 2 + 20;
    
    // Gene Legend (uses combined top genes)
    const geneLegendGroup = legendGroup.append("g").attr("transform", `translate(0, ${legendY})`);
    geneLegendGroup.append("text").attr("x", 0).attr("y", 0).text("Top Genes (Both Regions)").attr("class", "legend-title");
    legendY += 25;
    sortedOverallGenes.forEach((gene, i) => {
        const itemY = legendY + i * 20;
        const geneItemGroup = legendGroup.append("g").attr("class", "legend-item").attr("transform", `translate(0, ${itemY})`);
        geneItemGroup.append("circle").attr("cx", 0).attr("cy", 0).attr("r", 6).style("fill", colorScale(gene));
        geneItemGroup.append("text").attr("x", 15).attr("y", 0).text(gene);
    });
    legendY += sortedOverallGenes.length * 20 + 15; 

    // Component Size Legend (based on combined NMF range)
    const sizeLegendGroup = legendGroup.append("g").attr("transform", `translate(0, ${legendY})`);
    sizeLegendGroup.append("text").attr("x", 0).attr("y", 0).text("Component Size (Radius = NMF)").attr("class", "legend-title");
    legendY += 25;
    const legendNmfValues = [0, 0.5, 1.0]; 
    let sizeItemY = legendY;
    legendNmfValues.forEach((value, i) => {
        const radius = radiusScale(value); 
        if (value === 0 && radius < 0.1) return; 
        const itemGroup = legendGroup.append("g").attr("class", "legend-item").attr("transform", `translate(0, ${sizeItemY})`);
        itemGroup.append("circle").attr("cx", radius + 5).attr("cy", radius).attr("r", radius).style("fill", DEFAULT_COMPONENT_COLOR).attr("opacity", 0.6);
        itemGroup.append("text").attr("x", radius * 2 + 15).attr("y", radius).text(`NMF = ${value.toFixed(1)}`).style("font-size", "10px").attr("alignment-baseline","middle");
        sizeItemY += radius * 2 + 15; 
    });

}).catch(error => {
    console.error('Error processing data:', error);
    svg.append("text")
       .attr("x", chartWidth / 2)
       .attr("y", svgHeight / 2) // Center error in the whole SVG
       .attr("text-anchor", "middle")
       .style("fill", "red")
       .text(`Error loading data. Check console and file paths (region1_NMF..., sorted_gene1..., region2_NMF..., sorted_gene2...).`);
}); 