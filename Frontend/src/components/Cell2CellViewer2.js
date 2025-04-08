import React, { useRef, useState, useEffect } from 'react';
import { Button, Select, Spin, Empty } from 'antd';

export const Cell2CellViewer2 = ({ regions, analyzedRegion, cell2cellData, setCell2cellData, cell2cellDataLoading }) => {
    const [receiver, setReceiver] = useState(null);
    const [sender, setSender] = useState(null);
    const [receiverGenes, setReceiverGenes] = useState([]);
    const [senderGenes, setSenderGenes] = useState([]);
    const [selectedCellIds, setSelectedCellIds] = useState([]);

    const fetchCellTypeList = (sample_name) => {
        fetch('/get_cell_types', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name: sample_name })
        })
            .then(res => res.json())
            .then(data => {
                console.log(data, "cell type list");
                setReceiver(data.receiver);
                setSender(data.sender);
            });
    }

    const fetchGeneList = (sample_name) => {
        fetch('/get_cell2cell_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name: sample_name })
        })
            .then(res => res.json())
            .then(data => {
                console.log(data, "gene list");
                setReceiverGenes(data);
                setSenderGenes(data);
            });

    }
    // Method to find and store cellIds based on analyzedRegion
    const findCellIdsByRegion = () => {
        const region = regions.find(region => region.name === analyzedRegion);
        if (region) {
            // Found the matching region
            const { sampleId, cellIds } = region;
            fetchCellTypeList(sampleId);
            fetchGeneList(sampleId);
            setSelectedCellIds(cellIds);
        }
    };

    // Call findCellIdsByRegion once when the component is mounted or analyzedRegion changes
    useEffect(() => {
        if (analyzedRegion) {
            findCellIdsByRegion();
        }
    }, [analyzedRegion, regions]);

    return (
        cell2cellDataLoading ? (
            <Spin spinning={true} style={{ width: '100%', height: '100%' }} />
        ) : (
            Object.keys(cell2cellData).length > 0 ? (
                <div style={{ width: '100%', height: '100%' }}>
                    <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', gap: 5 }}>
                        <Select>
                        </Select>
                    </div>
                </div>
            ) : (
                <Empty
                    image={Empty.PRESENTED_IMAGE_SIMPLE}
                    description="Please set receiver and sender first"
                    style={{ width: '100%', height: '100%' }}
                />
            )
        )
    );
};