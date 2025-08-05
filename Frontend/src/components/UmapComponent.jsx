import { ScatterplotUmap } from "./ScatterPlotUmap";
import { Spin } from "antd";

export const UmapComponent = ({ 
  umapData, 
  umapLoading, 
  title = "UMAP Scatter Plot", 
  adata_umap_title, 
  hoveredCluster, 
  setHoveredCluster, 
  umapId, 
  sampleId, 
  setCellName,
  setPseudotimeData,
  setPseudotimeLoading
}) => {
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
          adata_umap_title={adata_umap_title}
          pointSize={3}
          hoveredCluster={hoveredCluster}
          setHoveredCluster={setHoveredCluster}
          umapId={umapId}
          sampleId={sampleId}
          setCellName={setCellName}
          setPseudotimeData={setPseudotimeData}
          setPseudotimeLoading={setPseudotimeLoading}
        />
      ) : null}
    </div>
  );
};