import { useEffect, useState } from 'react';
import './App.css';
import { TissueViewer } from './components/tissueViewer';
// import { Umap } from './components/umap';


function App() {
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState([]);
  const [sampleId, setSampleId] = useState("skin_TXK6Z4X_A1");
  const [cellTypeDir, setCellTypeDir] = useState([]);

  const fetch_Cell_Type_With_Coordinates_Data = () => {
    fetch('/get_cell_type_coordinates', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ sample_id: sampleId })
    })
      .then((response) => response.json())
      .then((data) => {
        setCellTypeCoordinatesData(data);
      })
  };

  const fetct_Cell_Type_Directory = () => {
    fetch('/get_unique_cell_types', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ sample_id: sampleId })
    })
      .then((response) => response.json())
      .then((data) => {
        setCellTypeDir(data);
      })
  }

  useEffect(() => {
    fetch_Cell_Type_With_Coordinates_Data();
    fetct_Cell_Type_Directory()
  }, [sampleId]);

  return (
    <div className="App">
      <div className='content'>
        <TissueViewer
          sampleId={sampleId}
          cellTypeDir={cellTypeDir}
          cellTypeCoordinatesData={cellTypeCoordinatesData}
        />
        <div style={{ height: '100%', width: '30%' }}>123</div>
      </div>
    </div>
  );
}

export default App;
