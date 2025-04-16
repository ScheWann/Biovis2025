import React, { useState, useEffect, useRef } from 'react';
import { Button, Select } from 'antd';
import * as d3 from 'd3';

// StackBarChart 组件：使用 d3 绘制堆积条形图
const StackBarChart = ({ cell2cellData, width = 800, height = 400 }) => {
    const svgRef = useRef();

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
            // 取负数绘制 sender 部分
            sender: -counts[id].sender,
        }));

        // 清空上一次的 svg 内容
        d3.select(svgRef.current).selectAll('*').remove();

        // 设置图表边距
        const margin = { top: 20, right: 20, bottom: 50, left: 50 };
        const innerWidth = width - margin.left - margin.right;
        const innerHeight = height - margin.top - margin.bottom;

        // 创建 svg
        const svg = d3
            .select(svgRef.current)
            .attr('width', width)
            .attr('height', height);

        // 添加一个分组用于绘制图表内容
        const chartGroup = svg
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // X 轴: 每个 id 显示为分类
        const xScale = d3
            .scaleBand()
            .domain(data.map(d => d.id))
            .range([0, innerWidth])
            .padding(0.2);

        // Y 轴: domain 根据 receiver 和 sender 的最大绝对值设置
        const maxReceiver = d3.max(data, d => d.receiver);
        const maxSender = d3.max(data, d => -d.sender); // sender 为负数
        const yScale = d3
            .scaleLinear()
            .domain([-maxSender, maxReceiver])
            .range([innerHeight, 0]);

        // 添加 X 轴：在 y= yScale(0) 处绘制
        const xAxis = d3.axisBottom(xScale);
        chartGroup
            .append('g')
            .attr('transform', `translate(0, ${yScale(0)})`)
            .call(xAxis)
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-0.5em")
            .attr("dy", "0.2em")
            .attr("transform", "rotate(-40)");

        // 添加 Y 轴
        const yAxis = d3.axisLeft(yScale);
        chartGroup.append('g').call(yAxis);

        // 为 receiver 绘制 bar（正方向）
        chartGroup
            .selectAll('.bar-receiver')
            .data(data)
            .enter()
            .append('rect')
            .attr('class', 'bar-receiver')
            .attr('x', d => xScale(d.id))
            .attr('y', d => yScale(d.receiver))
            .attr('width', xScale.bandwidth())
            .attr('height', d => yScale(0) - yScale(d.receiver))
            .attr('fill', '#4CAF50'); // receiver 使用绿色

        // 为 sender 绘制 bar（负方向）
        chartGroup
            .selectAll('.bar-sender')
            .data(data)
            .enter()
            .append('rect')
            .attr('class', 'bar-sender')
            .attr('x', d => xScale(d.id))
            .attr('y', yScale(0))
            .attr('width', xScale.bandwidth())
            .attr('height', d => yScale(d.sender) - yScale(0))
            .attr('fill', '#F44336'); // sender 使用红色

        // 如果需要添加数值标签，可取消下面代码的注释
        /*
        chartGroup.selectAll('.label-receiver')
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

    }, [cell2cellData, width, height]);

    return <svg ref={svgRef}></svg>;
};

// Cell2CellViewer2 组件：数据请求和参数选择，包含 StackBarChart 部分
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
                // 假设返回数据的格式是 [{ label: '细胞类型1', value: 'ID_XXX' }, ...]
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
                // 假设返回数据的格式是 [{ label: '基因1', value: 'Gene_1' }, ...]
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

    // 当选中区域发生变化时，根据第一个样本 id 请求细胞类型和基因数据
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

            {/* 当 cell2cellData 存在时绘制图表，否则显示提示 */}
            <div style={{ marginTop: 20 }}>
                {cell2cellData && Object.keys(cell2cellData).length > 0 ? (
                    <StackBarChart cell2cellData={cell2cellData} width={800} height={400} />
                ) : (
                    <div style={{ textAlign: 'center', padding: '20px' }}>No Data to Display</div>
                )}
            </div>
        </div>
    );
};
