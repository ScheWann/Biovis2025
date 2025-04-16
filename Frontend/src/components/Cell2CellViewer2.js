import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Button, Select } from 'antd';
import * as d3 from 'd3';

// StackBarChart 组件：使用 d3 绘制堆积条形图，同时支持 svg 自适应父容器
const StackBarChart = ({ cell2cellData }) => {
    const svgRef = useRef();
    const containerRef = useRef();
    const [dimensions, setDimensions] = useState({ width: 800, height: 400 });

    // 使用 ResizeObserver 监听父容器尺寸变化
    const resizeObserver = useRef(null);

    const observeResize = useCallback(() => {
        if (containerRef.current) {
            const { clientWidth, clientHeight } = containerRef.current;
            setDimensions({ width: clientWidth, height: clientHeight });
        }
    }, []);

    useEffect(() => {
        observeResize();
        resizeObserver.current = new ResizeObserver(observeResize);
        if (containerRef.current) {
            resizeObserver.current.observe(containerRef.current);
        }
        return () => {
            if (resizeObserver.current && containerRef.current) {
                resizeObserver.current.unobserve(containerRef.current);
            }
        };
    }, [observeResize]);

    useEffect(() => {
        if (!cell2cellData || Object.keys(cell2cellData).length === 0) return;

        // 获取所有区域key并排序
        const regionKeys = Object.keys(cell2cellData).sort();
        
        // 1. 统计每个ID在不同区域的收发次数
        const counts = {};
        regionKeys.forEach(regionKey => {
            cell2cellData[regionKey].forEach(interaction => {
                const receiverId = interaction.Receiver;
                const senderId = interaction.Sender;

                // 初始化数据结构
                if (!counts[receiverId]) {
                    counts[receiverId] = { 
                        receiver: Object.fromEntries(regionKeys.map(k => [k, 0])),
                        sender: Object.fromEntries(regionKeys.map(k => [k, 0]))
                    };
                }
                if (!counts[senderId]) {
                    counts[senderId] = { 
                        receiver: Object.fromEntries(regionKeys.map(k => [k, 0])),
                        sender: Object.fromEntries(regionKeys.map(k => [k, 0]))
                    };
                }

                // 累加计数
                counts[receiverId].receiver[regionKey]++;
                counts[senderId].sender[regionKey]++;
            });
        });

        // 2. 构建数据集
        const data = Object.keys(counts).map(id => ({
            id,
            receiver: counts[id].receiver,
            sender: counts[id].sender
        }));

        // 3. 创建颜色比例尺
        const colorScale = d3.scaleOrdinal()
            .domain(regionKeys)
            .range(d3.schemeCategory10);

        // 清空SVG
        d3.select(svgRef.current).selectAll('*').remove();

        // 4. 计算尺寸
        const margin = { top: 20, right: 20, bottom: 50, left: 50 };
        const innerWidth = dimensions.width - margin.left - margin.right;
        const innerHeight = dimensions.height - margin.top - margin.bottom;

        // 创建SVG容器
        const svg = d3.select(svgRef.current)
            .attr('width', dimensions.width)
            .attr('height', dimensions.height)
            .attr('viewBox', `0 0 ${dimensions.width} ${dimensions.height}`);

        const chartGroup = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // 5. 创建比例尺
        const xScale = d3.scaleBand()
            .domain(data.map(d => d.id))
            .range([0, innerWidth])
            .padding(0.2);

        // 创建堆叠生成器
        const stackGenerator = d3.stack()
            .keys(regionKeys);

        // 生成收发数据的堆叠结构
        const receiverStack = stackGenerator
            .value((d, key) => d.receiver[key] || 0)
            (data);

        const senderStack = stackGenerator
            .value((d, key) => -(d.sender[key] || 0)) // 发送方转为负值
            (data);

        // 计算Y轴范围
        const maxReceiver = d3.max(receiverStack.flat(2));
        const maxSender = d3.max(senderStack.flat(2).map(Math.abs));
        const yScale = d3.scaleLinear()
            .domain([-maxSender, maxReceiver])
            .nice()
            .range([innerHeight, 0]);

        // 6. 绘制坐标轴
        chartGroup.append('g')
            .attr('transform', `translate(0,${yScale(0)})`)
            .call(d3.axisBottom(xScale).tickFormat(() => ''));

        chartGroup.append('g')
            .call(d3.axisLeft(yScale));

        // 7. 绘制接收方堆叠条
        chartGroup.selectAll('.receiver-layer')
            .data(receiverStack)
            .enter().append('g')
            .attr('class', 'receiver-layer')
            .attr('fill', d => colorScale(d.key))
            .selectAll('rect')
            .data(d => d)
            .enter().append('rect')
            .attr('x', d => xScale(d.data.id))
            .attr('y', d => yScale(d[1]))
            .attr('height', d => yScale(d[0]) - yScale(d[1]))
            .attr('width', xScale.bandwidth());

        // 8. 绘制发送方堆叠条
        chartGroup.selectAll('.sender-layer')
            .data(senderStack)
            .enter().append('g')
            .attr('class', 'sender-layer')
            .attr('fill', d => colorScale(d.key))
            .selectAll('rect')
            .data(d => d)
            .enter().append('rect')
            .attr('x', d => xScale(d.data.id))
            .attr('y', d => yScale(d[0]))
            .attr('height', d => yScale(d[1]) - yScale(d[0]))
            .attr('width', xScale.bandwidth());

        // 9. 添加图例
        const legend = svg.append('g')
            .attr('font-family', 'sans-serif')
            .attr('font-size', 10)
            .attr('text-anchor', 'start')
            .selectAll('g')
            .data(regionKeys)
            .enter().append('g')
            .attr('transform', (d, i) => `translate(30,${i * 20})`);

        legend.append('rect')
            .attr('x', dimensions.width - 120)
            .attr('width', 19)
            .attr('height', 19)
            .attr('fill', colorScale);

        legend.append('text')
            .attr('x', dimensions.width - 96)
            .attr('y', 9.5)
            .attr('dy', '0.32em')
            .text(d => `${d}`);

    }, [cell2cellData, dimensions]);

    return (
        <div ref={containerRef} style={{ width: '100%', height: '100%' }}>
            <svg ref={svgRef}></svg>
        </div>
    );
};

// Cell2CellViewer2 组件：参数选择与数据请求
export const Cell2CellViewer2 = ({ regions, cell2cellData, setCell2cellData }) => {
    const [regionOptions, setRegionOptions] = useState([]);
    const [selectedRegion, setSelectedRegion] = useState([]);
    const [cellTypeList, setCellTypeList] = useState([]);
    const [genelist, setGenelist] = useState([]);
    const [receiver, setReceiver] = useState(null);
    const [sender, setSender] = useState(null);
    const [receiverGene, setReceiverGenes] = useState(null);
    const [senderGene, setSenderGenes] = useState(null);

    // 根据样本名称请求细胞类型列表
    const fetchCellTypeList = (sample_name) => {
        fetch('/get_cell_types', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name })
        })
            .then(res => res.json())
            .then(data => {
                setCellTypeList(data);
            });
    };

    // 根据样本名称请求基因列表
    const fetchGeneList = (sample_name) => {
        fetch('/get_cell2cell_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name })
        })
            .then(res => res.json())
            .then(data => {
                setGenelist(data);
            });
    };

    // 根据选中的区域构建参数
    const getSelectedRegionsInfo = () => {
        const params = {};
        selectedRegion.forEach(regionId => {
            const regionInfo = regions.find(r => r.id === regionId);
            if (regionInfo) {
                params[regionInfo.name] = {
                    sampleid: regionInfo.sampleId,
                    cell_list: regionInfo.cellIds,
                };
            }
        });
        return params;
    };

    // 确认参数并请求 cell2cell 数据
    const confirmCell2cellparameters = () => {
        if (!receiver || !sender || !receiverGene || !senderGene) {
            alert('Please select all parameters');
            return;
        }
        const selectedRegions = getSelectedRegionsInfo();
        const params = {
            regions: selectedRegions,
            receiver,
            sender,
            receiverGene,
            senderGene,
        };
        console.log('Selected parameters:', params);
        fetch('/get_cell_cell_interaction_data', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params)
        })
            .then(res => res.json())
            .then(data => {
                setCell2cellData(data);
            });
    };

    // 初始化区域下拉选项
    useEffect(() => {
        regions.forEach(region => {
            setRegionOptions(prevOptions => {
                const exists = prevOptions.find(option => option.value === region.id);
                return exists ? prevOptions : [...prevOptions, { label: region.name, value: region.id }];
            });
        });
    }, [regions]);

    // 当选中区域发生变化时，请求细胞类型和基因数据
    useEffect(() => {
        if (selectedRegion.length > 0) {
            const selectedRegions = regions.filter(region => selectedRegion.includes(region.id));
            const uniqueSampleIds = [...new Set(selectedRegions.map(r => r.sampleId))];
            if (uniqueSampleIds.length > 0) {
                fetchCellTypeList(uniqueSampleIds[0]);
                fetchGeneList(uniqueSampleIds[0]);
            }
        }
    }, [selectedRegion, regions]);

    // 简单的筛选函数（不区分大小写）
    const filterSelect = (input, option) =>
        option?.label?.toLowerCase().includes(input.toLowerCase());

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div style={{ margin: 5, fontSize: 14, fontWeight: 'bold' }}>
                Cell to Cell Interaction Viewer
            </div>
            {/* 参数选择区 */}
            <div
                style={{
                    width: '100%',
                    height: '10%',
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    gap: 5,
                    marginTop: 10,
                }}
            >
                <Select
                    size="small"
                    mode="multiple"
                    placeholder="Select Region"
                    options={regionOptions}
                    style={{ width: 200 }}
                    onChange={setSelectedRegion}
                    allowClear
                />
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size="small"
                    value={receiver}
                    options={cellTypeList}
                    onChange={setReceiver}
                    placeholder="Select Receiver"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size="small"
                    value={sender}
                    options={cellTypeList}
                    onChange={setSender}
                    placeholder="Select Sender"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size="small"
                    value={receiverGene}
                    options={genelist}
                    onChange={setReceiverGenes}
                    placeholder="Select Receiver Gene"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size="small"
                    value={senderGene}
                    options={genelist}
                    onChange={setSenderGenes}
                    placeholder="Select Sender Gene"
                    filterOption={filterSelect}
                />
                <Button style={{ width: 100 }} size="small" onClick={confirmCell2cellparameters}>
                    Confirm
                </Button>
            </div>

            {/* 绘制图表 */}
            <div style={{ marginTop: 20, width: '100%', height: '80%' }}>
                {cell2cellData && Object.keys(cell2cellData).length > 0 ? (
                    <StackBarChart cell2cellData={cell2cellData} />
                ) : (
                    <div style={{ textAlign: 'center', padding: '20px' }}>No Data to Display</div>
                )}
            </div>
        </div>
    );
};
