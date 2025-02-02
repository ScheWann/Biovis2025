import { useEffect, useState } from 'react';
import './App.css';
import { TissueImage } from './components/tissueImage';
import { Umap } from './components/umap';


function App() {
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState([]);
  const [sampleId, setSampleId] = useState("skin_TXK6Z4X_A1");

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

  useEffect(() => {
    fetch_Cell_Type_With_Coordinates_Data();
  }, [sampleId]);

  return (
    <div className="App">
      <div className='content'>
        <TissueImage
          sampleId={sampleId}
          cellTypeCoordinatesData={cellTypeCoordinatesData}
        />
      </div>
    </div>
  );
}

export default App;
