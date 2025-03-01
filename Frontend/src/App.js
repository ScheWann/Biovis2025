import { useEffect, useState } from 'react';
import './App.css';
import { MultiSampleViewer } from './components/MultiSampleViewer';
import { Select, Spin, message } from 'antd';

function App() {
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState({});
  const [samples, setSamples] = useState([]);
  const [selectedSamples, setSelectedSamples] = useState([]);
  const [cellTypeDir, setCellTypeDir] = useState([]);
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

  // get cell type data for selected samples
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

  // get cell type directory for selected samples
  const fetchCellTypeDirectory = async (sampleIds) => {
    try {
      const uniqueTypes = new Set();
      await Promise.all(
        sampleIds.map(sampleId =>
          fetch('/get_unique_cell_types', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: sampleId })
          })
            .then(res => res.json())
            .then(data => {
              data.forEach(type => uniqueTypes.add(type));
            })
        )
      );
      setCellTypeDir(Array.from(uniqueTypes));
    } catch (error) {
      message.error('Get cell type directory failed');
      console.error(error);
    }
  };

  // initial loading
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
      <div className="content">
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

        {/* main view */}
        <Spin spinning={loading}>
          {selectedSamples.length > 0 ? (
            <MultiSampleViewer
              samples={samples.filter(s => selectedSamples.includes(s.id))}
              cellTypeCoordinatesData={cellTypeCoordinatesData}
              cellTypeDir={cellTypeDir}
              regions={regions}
              setRegions={setRegions}
            />
          ) : (
            <div style={{
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              height: '80vh',
              color: '#999'
            }}>
              Please select at least one sample to view
            </div>
          )}
        </Spin>
      </div>
    </div>
  );
}

export default App;