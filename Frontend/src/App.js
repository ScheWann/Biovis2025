import { useEffect, useState } from 'react';
import { Select, Spin, message, Button, Splitter } from 'antd';
import './App.css';
import { MultiSampleViewer } from './components/MultiSampleViewer';
import { NMFGOExpressionViewer } from './components/NMFGOExpressionViewer';
import { Cell2CellViewer } from './components/Cell2CellViewer';
import { Cell2CellViewer2 } from './components/Cell2CellViewer2';
import { GeneExpressionViewer } from './components/GeneExpressionViewer';
import { PseudoTemporalViewer } from './components/PseudoTemporalViewer';


function App() {
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState({}); // each sample's cell type directory(e.g. {"skin_TXK6Z4X_A1": [{"cell_type": "cd19+cd20+ b","cell_x": 3526, "cell_y": 3780, "id": "ID_1}, ...}])'
  const [selectOptions, setSelectOptions] = useState([]); // Available sample Option(e.g. [{value: 'skin_TXK6Z4X_A1', label: 'skin_TXK6Z4X_A1'}, ...])
  const [samples, setSamples] = useState([]); // Select the sample to be displayed (the beginning input of the entire project)[{id: 'sample_id', name: 'sample_id'}, ...]
  const [selectedSamples, setSelectedSamples] = useState([]); // Select the sample to analyze(e.g. '[skin_TXK6Z4X_D1]')
  const [cellTypeDir, setCellTypeDir] = useState({}); // Each sample's cell type directory based on users' sample selection(e.g. {"skin_TXK6Z4X_A1": ["cd19+cd20+ b", "cd19+cd20+ t"], ...})
  const [regions, setRegions] = useState([]); // Selected regions in the multiSampleViewer to perform analysis (e.g. [{name: 'region1', sampleId: 'skin_TXK6Z4X_A1', cellIds: ['ID_1', 'ID_2']}, ...])
  const [loading, setLoading] = useState(false); // Full page loading
  const [analyzedRegion, setAnalyzedRegion] = useState(null); // Selected in the multiSampleViewer to perform analysis
  // NMF data and its loading settings
  const [NMFGOData, setNMFGOData] = useState({});
  const [NMFGODataLoading, setNMFGODataLoading] = useState(false);
  const [NMFclusterCells, setNMFclusterCells] = useState([]);
  // cell2cell data and its loading settings
  const [cell2cellData, setCell2cellData] = useState({});
  const [cell2cellDataLoading, setCell2cellDataLoading] = useState(false);
  const [selectedRegionGeneExpressionData, setSelectedRegionGeneExpressionData] = useState({});

  // get all aviailable samples
  const fetchAvailableSamples = () => {
    fetch('/get_available_samples')
      .then(response => response.json())
      .then(data => {
        setSelectOptions(data);
      })
      .catch(error => {
        message.error('Get samples failed');
        console.error(error);
      });
  };

  // get cell information(cell_type, cell_x, cell_y, id) for selected samples
  const fetchCellTypeData = (sampleIds) => {
    setLoading(true);
    fetch('/get_cell_type_coordinates', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sample_ids: sampleIds })
    })
      .then(res => res.json())
      .then(data => {
        setCellTypeCoordinatesData(data);
        setLoading(false);
      })
  };

  // get cell type directory for each selected sample
  const fetchCellTypeDirectory = (sampleIds) => {
    fetch('/get_unique_cell_types', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sample_ids: sampleIds })
    })
      .then(res => res.json())
      .then(data => {
        setCellTypeDir(data);
      }).catch(error => {
        message.error('Get cell type directory failed');
        console.error(error);
      });
  };

  // confirm selected samples
  const confirmSamples = () => {
    if (selectedSamples.length === 0) {
      message.warning('Please select at least one sample');
    } else {
      fetchCellTypeData(selectedSamples);
      fetchCellTypeDirectory(selectedSamples);
      setSamples(selectedSamples.map(sample => ({ id: sample, name: sample })));
    }
  };

  // initial loading (get all available sample list)
  useEffect(() => {
    fetchAvailableSamples();
  }, []);

  return (
    <div className="App">
      <div className='main'>
        {/* select samples */}
        <div className="selectSamples">
          <Select
            size='small'
            mode="multiple"
            placeholder="Select samples"
            value={selectedSamples}
            onChange={setSelectedSamples}
            options={selectOptions}
            style={{ width: '100%', margin: 8 }}
            maxTagCount="responsive"
            loading={loading}
          />
          <Button size='small' onClick={confirmSamples}>Confirm</Button>
        </div>

        {/* all views */}
        <div className="content" style={{ position: "relative" }}>
          {loading && (
            <div style={{
              position: "absolute",
              top: 0,
              left: 0,
              width: "100%",
              height: "100%",
              display: "flex",
              justifyContent: "center",
              alignItems: "center",
              background: "rgba(0, 0, 0, 0.5)",
              zIndex: 20
            }}>
              <Spin spinning={true} size="large" />
            </div>
          )}

          {samples.length > 0 ? (
            <Splitter lazy style={{ width: "100%", height: "100%" }}>
              <Splitter.Panel defaultSize="70%" min="50%" max="80%">
                <MultiSampleViewer
                  setLoading={setLoading}
                  samples={samples}
                  cellTypeCoordinatesData={cellTypeCoordinatesData}
                  cellTypeDir={cellTypeDir}
                  regions={regions}
                  setRegions={setRegions}
                  analyzedRegion={analyzedRegion}
                  setAnalyzedRegion={setAnalyzedRegion}
                  setNMFGOData={setNMFGOData}
                  setNMFGODataLoading={setNMFGODataLoading}
                  NMFclusterCells={NMFclusterCells}
                  setSelectedRegionGeneExpressionData={setSelectedRegionGeneExpressionData}
                />
              </Splitter.Panel>
              <Splitter.Panel defaultSize="30%" min="20%" max="50%">
                <Splitter lazy layout='vertical'>
                  <Splitter.Panel defaultSize="33%" min="20%" max="45%" style={{ borderBottom: "1px solid #e8e8e8" }}>
                    {/* <GeneExpressionViewer
                      data={selectedRegionGeneExpressionData}
                    /> */}
                    <NMFGOExpressionViewer 
                      NMFGOData={NMFGOData}
                      NMFGODataLoading={NMFGODataLoading}
                      setNMFclusterCells={setNMFclusterCells}
                    />
                  </Splitter.Panel>
                  <Splitter.Panel defaultSize="33%" min="20%" max="45%" style={{ borderBottom: "1px solid #e8e8e8" }}>
                    <PseudoTemporalViewer />
                  </Splitter.Panel>
                  <Splitter.Panel defaultSize="33%" min="20%" max="45%">
                    {/* <Cell2CellViewer /> */}
                    <Cell2CellViewer2
                      regions={regions}
                      analyzedRegion={analyzedRegion}
                      cell2cellData={cell2cellData}
                      setCell2cellData={setCell2cellData}
                      cell2cellDataLoading={cell2cellDataLoading}
                    />
                  </Splitter.Panel>
                </Splitter>
              </Splitter.Panel>
            </Splitter>
          ) : (
            <div style={{
              display: "flex",
              justifyContent: "center",
              alignItems: "center",
              height: "100%",
              width: "100%",
              color: "#999"
            }}>
              Please select at least one sample to view
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

export default App;