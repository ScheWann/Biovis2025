import { useEffect, useState } from 'react';
import './App.css';
import { TissueImage } from './components/tissueImage';

function App() {
  const [kmeansSize, setKmeansSize] = useState(8);
  const [positionWithClusterData, setPositionWithClusterData] = useState([]);

  const fetchPositions_with_clusters_data = () => {
    fetch('/get_um_008_positions_with_clusters', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ kmeans: kmeansSize })
    })
      .then((response) => response.json())
      .then((data) => {
        setPositionWithClusterData(data);
      })
  };

  useEffect(() => {
    fetchPositions_with_clusters_data();
  }, [kmeansSize]);

  return (
    <div className="App">
      <TissueImage
        kmeansSize={kmeansSize}
        setKmeansSize={setKmeansSize}
        positionWithClusterData={positionWithClusterData}
      />
    </div>
  );
}

export default App;
