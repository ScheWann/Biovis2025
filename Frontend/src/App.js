import { useEffect, useState } from 'react';
import './App.css';
import { MultiSampleViewer } from './components/MultiSampleViewer'; // 修改为使用新的多样本组件
import { Select, Spin, message } from 'antd';

function App() {
  // 状态管理
  const [cellTypeCoordinatesData, setCellTypeCoordinatesData] = useState({});
  const [samples, setSamples] = useState([]); // 改为存储样本列表
  const [selectedSamples, setSelectedSamples] = useState([]); // 当前选中的样本
  const [cellTypeDir, setCellTypeDir] = useState([]);
  const [regions, setRegions] = useState([]);
  const [loading, setLoading] = useState(false);

  // 获取所有可用样本
  const fetchAvailableSamples = async () => {
    try {
      const response = await fetch('/get_available_samples');
      const data = await response.json();
      setSamples(data.map(sample => ({
        id: sample.id,
        name: sample.name || `样本 ${sample.id}` // 默认名称
      })));
    } catch (error) {
      message.error('获取样本列表失败');
      console.error(error);
    }
  };

  // 获取细胞类型数据
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
            .then(data => {
              console.log(`Loaded data for ${sampleId}:`, data.slice(0, 5)); // 打印前5条数据
              return data;
            })
        )
      );

      // 合并数据时确保保留所有样本数据
      const newData = sampleIds.reduce((acc, sampleId, index) => {
        acc[sampleId] = responses[index];
        return acc;
      }, {});

      setCellTypeCoordinatesData(prev => ({ ...prev, ...newData }));
    } catch (error) {
      message.error(`获取细胞数据失败: ${error.message}`);
      console.error('Error fetching cell data:', error);
    } finally {
      setLoading(false);
    }
  };

  // 获取细胞类型目录
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
      message.error('获取细胞类型目录失败');
      console.error(error);
    }
  };

  // 初始化加载
  useEffect(() => {
    fetchAvailableSamples();
  }, []);

  // 当选中样本变化时加载数据
  useEffect(() => {
    if (selectedSamples.length > 0) {
      fetchCellTypeData(selectedSamples);
      fetchCellTypeDirectory(selectedSamples);
    }
  }, [selectedSamples]);

  return (
    <div className="App">
      <div className="content">
        {/* 样本选择器 */}
        <div>
          <Select
            size='small'
            mode="multiple"
            placeholder="选择要查看的样本"
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
        </div>

        {/* 主视图区域 */}
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
              请从上方选择至少一个样本进行查看
            </div>
          )}
        </Spin>
      </div>
    </div>
  );
}

export default App;