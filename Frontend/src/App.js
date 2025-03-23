import { useEffect, useState } from 'react';
import { Select, Spin, message } from 'antd';
import './App.css';
import { MultiSampleViewer } from './components/MultiSampleViewer';
import { Cell2CellViewer } from './components/Cell2CellViewer';
import { GeneExpressionViewer } from './components/GeneExpressionViewer';
import { PseudoTemporalViewer } from './components/PseudoTemporalViewer';


function App() {
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState({});
  const [samples, setSamples] = useState([]);
  const [selectedSamples, setSelectedSamples] = useState([]);
  const [cellTypeDir, setCellTypeDir] = useState({});
  const [regions, setRegions] = useState([]);
  const [loading, setLoading] = useState(false);

  // get all aviailable samples
  const fetchAvailableSamples = async () => {
    try {
      const response = await fetch('/get_available_samples');
      const data = await response.json();
      setSamples(data.map(sample => ({
        id: sample.id,
        name: sample.name || `${sample.id}`
      })));
    } catch (error) {
      message.error('Get samples failed');
      console.error(error);
    }
  };

  // get cell information(cell_type, cell_x, cell_y, id) for selected samples
  const fetchCellTypeData = async (sampleIds) => {
    setLoading(true);
    try {
      const responses = await Promise.all(
        sampleIds.map(sampleId =>
          fetch('/get_cell_type_coordinates', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: sampleId })
          })
            .then(res => {
              if (!res.ok) throw new Error(`HTTP error! status: ${res.status}`);
              return res.json();
            })
        )
      );

      // save all samples data when merge data
      const newData = sampleIds.reduce((acc, sampleId, index) => {
        acc[sampleId] = responses[index];
        return acc;
      }, {});

      setCellTypeCoordinatesData(prev => ({ ...prev, ...newData }));
    } catch (error) {
      message.error(`Error fetching cell data: ${error.message}`);
    } finally {
      setLoading(false);
    }
  };

  // get cell type directory for each selected sample
  const fetchCellTypeDirectory = async (sampleIds) => {
    try {
      const cellTypeMap = {};
      await Promise.all(
        sampleIds.map(sampleId =>
          fetch('/get_unique_cell_types', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: sampleId })
          })
            .then(res => res.json())
            .then(data => {
              cellTypeMap[sampleId] = data;
            })
        )
      );
      setCellTypeDir(cellTypeMap);
    } catch (error) {
      message.error('Get cell type directory failed');
      console.error(error);
    }
  };

  // initial loading (get all available sample list)
  useEffect(() => {
    fetchAvailableSamples();
  }, []);

  // loading data for selected samples
  useEffect(() => {
    if (selectedSamples.length > 0) {
      fetchCellTypeData(selectedSamples);
      fetchCellTypeDirectory(selectedSamples);
    }
  }, [selectedSamples]);

  return (
    <div className="App">
      <div className='main'>
        {/* select samples */}
        <Select
          size='small'
          mode="multiple"
          placeholder="Select samples"
          value={selectedSamples}
          onChange={setSelectedSamples}
          options={samples.map(sample => ({
            label: sample.name,
            value: sample.id
          }))}
          style={{ width: '100%', margin: 8 }}
          maxTagCount="responsive"
          loading={loading}
        />

        {/* all views */}
        <div className='content'>
          {loading ? (
            <Spin spinning={true} size="large" style={{ height: '100%', width: '100%', display: 'flex', justifyContent: 'center', alignItems: 'center'}} />
          ) : (
            selectedSamples.length > 0 ? (
              <>
                <MultiSampleViewer
                  samples={samples.filter(s => selectedSamples.includes(s.id))}
                  cellTypeCoordinatesData={cellTypeCoordinatesData}
                  cellTypeDir={cellTypeDir}
                  regions={regions}
                  setRegions={setRegions}
                />

                <div className='auxiliaryViews'>
                  <PseudoTemporalViewer />
                  <Cell2CellViewer />
                  <GeneExpressionViewer />
                </div>
              </>
            ) : (
              <div style={{
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                height: '100%',
                width: '100%',
                color: '#999'
              }}>
                Please select at least one sample to view
              </div>
            )
          )}
        </div>
      </div>
    </div>
  );
}

export default App;