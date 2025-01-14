import { useEffect, useState } from 'react';
import './App.css';
import { TissueImage } from './components/tissueImage';
import { Umap } from './components/umap';


function App() {
  const [kmeansSize, setKmeansSize] = useState(8);
  const [binSize, setBinSize] = useState('008');
  const [tissueData, setTissueData] = useState([]);
  const [umapPositionWithClusterData, setUmapPositionWithClusterData] = useState([]);
  const [geneName, setGeneName] = useState([]);
  const [mode, setMode] = useState('kmeans');

  const fetchPositions_with_clusters_data = () => {
    fetch('/get_um_positions_with_clusters', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ kmeans: kmeansSize, bin_size: binSize })
    })
      .then((response) => response.json())
      .then((data) => {
        setTissueData(data);
      })
  };

  const fetchUmapPositions_with_clusters_data = () => {
    fetch('/get_umap_positions', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ kmeans: kmeansSize, bin_size: binSize })
    })
      .then((response) => response.json())
      .then((data) => {
        setUmapPositionWithClusterData(data);
      })
  };

  useEffect(() => {
    fetchPositions_with_clusters_data();
    fetchUmapPositions_with_clusters_data();
  }, [kmeansSize, binSize]);

  return (
    <div className="App">
      <div className='content'>
        <TissueImage
          mode={mode}
          setMode={setMode}
          geneName={geneName}
          setGeneName={setGeneName}
          binSize={binSize}
          kmeansSize={kmeansSize}
          setKmeansSize={setKmeansSize}
          tissueData={tissueData}
          setTissueData={setTissueData}
        />
        <Umap
          kmeansSize={kmeansSize}
          umapPositionWithClusterData={umapPositionWithClusterData}
        />
      </div>
    </div>
  );
}

export default App;
