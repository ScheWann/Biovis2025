import { useEffect, useState } from 'react';
import { Select, Spin, message, Button, Splitter, Modal, Form, Input, Upload, ConfigProvider } from 'antd';
import './App.css';
// import { MultiSampleViewer } from './components/MultiSampleViewer';
import { SampleViewer } from './components/SampleViewer';
import { PlusOutlined, InboxOutlined, PaperClipOutlined } from '@ant-design/icons';
import '@ant-design/v5-patch-for-react-19';

// Custom theme configuration
const customTheme = {
  token: {
    colorPrimary: '#1890ff',
    colorPrimaryHover: '#40a9ff',
    colorPrimaryActive: '#096dd9',
  },
};

function App() {
  const [coordinatesData, setCoordinatesData] = useState({}); // each sample's cell type directory(e.g. {"skin_TXK6Z4X_A1": [{"cell_type": "cd19+cd20+ b","cell_x": 3526, "cell_y": 3780, "id": "ID_1}, ...}])'
  const [selectOptions, setSelectOptions] = useState([]); // Available sample Option(e.g. [{value: 'skin_TXK6Z4X_A1', label: 'skin_TXK6Z4X_A1'}, ...])
  const [selectedSamples, setSelectedSamples] = useState([]); // Confirmed sample to be displayed(e.g. [{id: 'sample_id', name: 'sample_id'}, ...])
  const [tempSamples, setTempSamples] = useState([]); // The sample identified in the selector
  const [sampleDataLoading, setSampleDataLoading] = useState(false); // Sample Data Loading
  const [uploadFormVisible, setUploadFormVisible] = useState(false); // Upload form visibility

  useEffect(() => {
    fetchSamplesOption();
  }, []);

  // get all aviailable sample options
  const fetchSamplesOption = () => {
    fetch('/api/get_samples_option')
      .then(response => response.json())
      .then(data => {
        setSelectOptions(data);
      })
      .catch(error => {
        message.error('Get samples failed');
      });
  };

  // get cell coordinates for selected samples(cell or dot)
  const fetchCoordinates = (sampleIds) => {
    setSampleDataLoading(true);
    fetch('/api/get_coordinates', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sample_ids: sampleIds })
    })
      .then(res => res.json())
      .then(data => {
        setCoordinatesData(data);
        setSampleDataLoading(false);
      })
  };

  // confirm selected samples
  const confirmSamples = () => {
    if (tempSamples.length === 0) {
      message.warning('Please select at least one sample');
    } else {
      fetchCoordinates(tempSamples);
      setSelectedSamples(tempSamples.map(sample => ({ id: sample, name: sample })));
    }
  };

  const handleUploadSTData = async (values) => {
    const { name, description, folder } = values;
    const formData = new FormData();
    formData.append('name', name);
    formData.append('description', description || '');

    // Only include relevant files
    folder.forEach(fileObj => {
      const file = fileObj.originFileObj
      const path = fileObj.originFileObj.webkitRelativePath;
      const segments = path.split('/');
      const relativePath = segments.slice(1).join('/');

      if (
        path.endsWith('binned_outputs/square_002um/filtered_feature_bc_matrix.h5') ||
        path.endsWith('binned_outputs/square_008um/filtered_feature_bc_matrix.h5') ||
        path.endsWith('binned_outputs/square_016um/filtered_feature_bc_matrix.h5') ||
        (segments.length > 2 && segments[1] === 'spatial')
      ) {
        formData.append('files', file, relativePath);
      }
    });

    try {
      const response = await fetch('/upload_spaceranger', {
        method: 'POST',
        body: formData,
      });
      if (response.ok) {
        message.success('Upload successful!');
        setUploadFormVisible(false);
      } else {
        message.error('Upload failed.');
      }
    } catch (err) {
      message.error('Upload error: ' + err.message);
    }
  };

  return (
    <ConfigProvider theme={customTheme}>
      <div className="App">
        <div className='main'>
          {/* select samples */}
          <div className="selectSamples">
            <Select
              size='small'
              mode="multiple"
              placeholder="Select samples"
              value={tempSamples}
              onChange={setTempSamples}
              options={selectOptions}
              style={{ width: '100%', marginTop: 8, marginBottom: 8 }}
              maxTagCount="responsive"
              loading={sampleDataLoading}
            />
            <Button size='small' onClick={() => setUploadFormVisible(true)} icon={<PlusOutlined />} />
            <Button size='small' onClick={confirmSamples}>Confirm</Button>
          </div>

          {/* Upload Data Form Modal */}
          <Modal
            title="Upload Data"
            open={uploadFormVisible}
            onCancel={() => setUploadFormVisible(false)}
            footer={null}
            destroyOnHidden
          >
            <Form
              layout="vertical"
              onFinish={handleUploadSTData}
            >
              <Form.Item
                label="Name"
                name="name"
                rules={[{ required: true, message: 'Please input a name!' }]}
              >
                <Input placeholder="Custom name" />
              </Form.Item>
              <Form.Item
                label="Description"
                name="description"
              >
                <Input.TextArea placeholder="Description (optional)" rows={2} />
              </Form.Item>
              <Form.Item
                label="Upload Folder"
                name="folder"
                valuePropName="fileList"
                getValueFromEvent={e => Array.isArray(e) ? e : e && e.fileList}
                rules={[{ required: true, message: 'Please upload a spaceranger output folder!' }]}
              >
                <Upload.Dragger
                  directory
                  multiple
                  showUploadList={true}
                  beforeUpload={(file) => {
                    const path = file.webkitRelativePath || file.name;
                    const matrixH5Pattern = /binned_outputs\/square_(002|008|016)um\/filtered_feature_bc_matrix\.h5$/;
                    const spatialPattern = /\/spatial\//;
                    if (matrixH5Pattern.test(path)) {
                      return false;
                    }
                    if (spatialPattern.test(path) && !/\/\./.test(path)) {
                      return false;
                    }
                    return Upload.LIST_IGNORE;
                  }}
                  itemRender={(originNode, file) => (
                    <div className='ant-upload-list-item-name'>
                      <PaperClipOutlined style={{ marginRight: 6 }} />
                      {file.name}
                    </div>
                  )}
                >
                  <p className="ant-upload-drag-icon">
                    <InboxOutlined />
                  </p>
                  <p className="ant-upload-hint">Click or drag folder to this area to upload</p>
                </Upload.Dragger>
              </Form.Item>
              <Form.Item>
                <Button type="primary" htmlType="submit" block>Upload</Button>
              </Form.Item>
            </Form>
          </Modal>

          {/* all views */}
          <div className="content" style={{ position: "relative" }}>
            {sampleDataLoading && (
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
                zIndex: 1000
              }}>
                <Spin spinning={true} size="large" />
              </div>
            )}

            {selectedSamples.length > 0 ? (
              <Splitter lazy style={{ width: "100%", height: "100%" }}>
                <Splitter.Panel defaultSize="70%" min="50%" max="80%">
                  <SampleViewer
                    selectedSamples={selectedSamples}
                    coordinatesData={coordinatesData}
                  />
                </Splitter.Panel>
                <Splitter.Panel defaultSize="30%" min="20%" max="50%">
                  <Splitter lazy layout='vertical'>
                    <Splitter.Panel defaultSize="33%" min="20%" max="45%" style={{ borderBottom: "1px solid #e8e8e8" }}>
                      Gene Expression
                    </Splitter.Panel>
                    <Splitter.Panel defaultSize="33%" min="20%" max="45%" style={{ borderBottom: "1px solid #e8e8e8" }}>
                      UMAP
                    </Splitter.Panel>
                    <Splitter.Panel defaultSize="33%" min="20%" max="45%">
                      Glyphs
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
    </ConfigProvider>
  );
}

export default App;