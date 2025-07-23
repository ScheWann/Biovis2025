import { ScatterplotUmap } from "./ScatterPlotUmap";
import { Spin } from "antd";

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

export const UmapComponent = ({ umapData, umapLoading, title = "UMAP Scatter Plot" }) => {
  return (
    <div style={{ height: '100%', width: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', position: 'relative' }}>
      {umapLoading ? (
        <div style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          justifyContent: 'center',
          gap: '12px',
          width: '100%',
          height: '100%',
          padding: 5
        }}>
          <Spin size="large"/>
          <div style={{ fontSize: '12px', color: '#999' }}>
            Generating {title}...
          </div>
        </div>
      ) : umapData ? (
        <ScatterplotUmap
          data={umapData}
          xAccessor={(d) => d.x}
          yAccessor={(d) => d.y}
          clusterAccessor={(d) => d.cluster}
          title={title}
          pointSize={3}
        />
      ) : null}
    </div>
  );
};