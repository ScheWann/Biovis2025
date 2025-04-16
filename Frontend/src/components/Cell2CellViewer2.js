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
        // 如果 cell2cellData 为空则不绘制
        if (!cell2cellData || Object.keys(cell2cellData).length === 0) {
            return;
        }

        // 统计每个 id 作为 receiver 和 sender 出现的次数
        const counts = {};
        Object.keys(cell2cellData).forEach(regionKey => {
            cell2cellData[regionKey].forEach(interaction => {
                const receiverId = interaction.Receiver;
                const senderId = interaction.Sender;

                if (!counts[receiverId]) {
                    counts[receiverId] = { receiver: 0, sender: 0 };
                }
                if (!counts[senderId]) {
                    counts[senderId] = { receiver: 0, sender: 0 };
                }
                counts[receiverId].receiver += 1;
                counts[senderId].sender += 1;
            });
        });

        // 构建数据数组，每个元素包含 id, receiverCount 和 senderCount (senderCount 转换为负数)
        const data = Object.keys(counts).map(id => ({
            id,
            receiver: counts[id].receiver,
            // sender 转为负数用于绘制在下方
            sender: -counts[id].sender,
        }));

        // 清空上一次的 svg 内容
        d3.select(svgRef.current).selectAll('*').remove();

        // 计算图表内边距与尺寸
        const margin = { top: 20, right: 20, bottom: 50, left: 50 };
        const innerWidth = dimensions.width - margin.left - margin.right;
        const innerHeight = dimensions.height - margin.top - margin.bottom;

        // 创建 svg，并设定 viewBox 实现响应式
        const svg = d3.select(svgRef.current)
            .attr('width', dimensions.width)
            .attr('height', dimensions.height)
            .attr('viewBox', `0 0 ${dimensions.width} ${dimensions.height}`);

        // 添加一个分组用于绘制图表内容
        const chartGroup = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // X 轴
        const xScale = d3.scaleBand()
            .domain(data.map(d => d.id))
            .range([0, innerWidth])
            .padding(0.2);

        // Y 轴
        const maxReceiver = d3.max(data, d => d.receiver);
        const maxSender = d3.max(data, d => -d.sender);
        const yScale = d3.scaleLinear()
            .domain([-maxSender, maxReceiver])
            .range([innerHeight, 0]);

        // 添加 X 轴：位于 y= yScale(0) 处
        chartGroup.append('g')
            .attr('transform', `translate(0, ${yScale(0)})`)
            .call(d3.axisBottom(xScale))
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-0.5em")
            .attr("dy", "0.2em")
            .attr("transform", "rotate(-40)");

        // 添加 Y 轴
        chartGroup.append('g').call(d3.axisLeft(yScale));

        // 绘制 receiver bar（正方向）
        chartGroup.selectAll('.bar-receiver')
            .data(data)
            .enter()
            .append('rect')
            .attr('class', 'bar-receiver')
            .attr('x', d => xScale(d.id))
            .attr('y', d => yScale(d.receiver))
            .attr('width', xScale.bandwidth())
            .attr('height', d => yScale(0) - yScale(d.receiver))
            .attr('fill', '#4CAF50');

        // 绘制 sender bar（负方向）
        chartGroup.selectAll('.bar-sender')
            .data(data)
            .enter()
            .append('rect')
            .attr('class', 'bar-sender')
            .attr('x', d => xScale(d.id))
            .attr('y', yScale(0))
            .attr('width', xScale.bandwidth())
            .attr('height', d => yScale(d.sender) - yScale(0))
            .attr('fill', '#F44336');

        // 可添加数值标签
        /* chartGroup.selectAll('.label-receiver')
              .data(data)
              .enter()
              .append('text')
              .attr('class', 'label-receiver')
              .attr('x', d => xScale(d.id) + xScale.bandwidth() / 2)
              .attr('y', d => yScale(d.receiver) - 5)
              .attr('text-anchor', 'middle')
              .text(d => d.receiver);
    
           chartGroup.selectAll('.label-sender')
              .data(data)
              .enter()
              .append('text')
              .attr('class', 'label-sender')
              .attr('x', d => xScale(d.id) + xScale.bandwidth() / 2)
              .attr('y', d => yScale(d.sender) + 15)
              .attr('text-anchor', 'middle')
              .text(d => -d.sender);
        */
    }, [cell2cellData, dimensions]);

    return (
        // 外层容器用 ref 获得尺寸信息
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
