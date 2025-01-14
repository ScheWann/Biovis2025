import { useEffect, useState } from 'react';
import './App.css';
import { Select } from 'antd';
import { TissueImage } from './components/tissueImage';
import { Umap } from './components/umap';

function App() {
  const [kmeansSize, setKmeansSize] = useState(8);
  const [binSize, setBinSize] = useState('008');
  const [positionWithClusterData, setPositionWithClusterData] = useState([]);
  const [umapPositionWithClusterData, setUmapPositionWithClusterData] = useState([]);
  const [geneName, setGeneName] = useState('');
  const [geneNameList, setGeneNameList] = useState([]);

  const geneNameSearch = value => {
    fetchGeneNameBySearch(value);
  }

  const geneNameChange = value => {
      setGeneName(value);
  };

  const fetchGeneNameBySearch = async (query = '') => {
    const response = await fetch(`/get_gene_name_search?q=${encodeURIComponent(query || '')}`);
    const data = await response.json();
  
    const geneOptions = data.map(item => ({
      label: item,
      value: item,
    }));
  
    setGeneNameList(geneOptions);
  };

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
        setPositionWithClusterData(data);
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
    fetchGeneNameBySearch();
    fetchPositions_with_clusters_data();
    fetchUmapPositions_with_clusters_data();
  }, [kmeansSize, binSize]);

  return (
    <div className="App">
      <div className='content'>
        <TissueImage
          kmeansSize={kmeansSize}
          setKmeansSize={setKmeansSize}
          positionWithClusterData={positionWithClusterData}
        />
        <Select
          showSearch
          value={geneName}
          size="small"
          style={{
            width: "10%",
            marginRight: 20,
          }}
          onChange={geneNameChange}
          onSearch={geneNameSearch}
          options={geneNameList}
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
