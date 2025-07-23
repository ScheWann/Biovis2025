import { ScatterplotUmap } from "./ScatterPlotUmap";

function generateDummyData(numPoints = 100, numClusters = 5) {
  // Cluster centers
  const centers = Array.from({ length: numClusters }, () => ({
    x: (Math.random() - 0.5) * 10,
    y: (Math.random() - 0.5) * 10,
  }));

  return Array.from({ length: numPoints }, () => {
    const clIdx = Math.floor(Math.random() * numClusters);
    const center = centers[clIdx];
    return {
      x: center.x + (Math.random() - 0.5) * 2,
      y: center.y + (Math.random() - 0.5) * 2,
      cluster: `Cluster ${clIdx + 1}`,
    };
  });
}

export const UmapComponent = () => {
  const data = generateDummyData(100, 3);
  return (
    <div>
      <ScatterplotUmap
        data={data}
        xAccessor={(d) => d.x}
        yAccessor={(d) => d.y}
        clusterAccessor={(d) => d.cluster}
        title="UMAP Scatter Plot"
        pointSize={4}
      />
    </div>
  );
};
