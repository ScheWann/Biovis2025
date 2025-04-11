import React, { useState, useEffect } from 'react';
import { Button, Select, Spin, Empty } from 'antd';

export const Cell2CellViewer2 = ({ regions, analyzedRegion, cell2cellData, setCell2cellData, cell2cellDataLoading }) => {
    const [cellTypeList, setCellTypeList] = useState([]);
    const [genelist, setGenelist] = useState([]);
    const [receiver, setReceiver] = useState(null);
    const [sender, setSender] = useState(null);
    const [receiverGene, setReceiverGenes] = useState(null);
    const [senderGene, setSenderGenes] = useState(null);
    const [selectedCellIds, setSelectedCellIds] = useState([]);

    const fetchCellTypeList = (sample_name) => {
        fetch('/get_cell_types', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name: sample_name })
        })
            .then(res => res.json())
            .then(data => {
                setCellTypeList(data);
            });
    };

    const fetchGeneList = (sample_name) => {
        fetch('/get_cell2cell_gene_list', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_name: sample_name })
        })
            .then(res => res.json())
            .then(data => {
                setGenelist(data);
            });
    };

    const getSelectedRegionsInfo = () => {
        return regions
            .filter(region => analyzedRegion[region.name])
            .map(region => ({
                sampleId: region.sampleId,
                cellIds: region.cellIds
            }));
    };

    // Method to find and store cellIds based on analyzedRegion
    const findCellIdsByRegion = () => {
        const selectedRegions = getSelectedRegionsInfo();
        const uniqueSampleIds = [...new Set(selectedRegions.map(r => r.sampleId))];

        const allCellIds = selectedRegions.flatMap(r => r.cellIds);
        setSelectedCellIds(allCellIds);

        // now it only focus on one sample
        if (uniqueSampleIds.length > 0) {
            fetchCellTypeList(uniqueSampleIds[0]);
            fetchGeneList(uniqueSampleIds[0]);
        }
    };

    const confirmCell2cellparameters = () => {
        if (!receiver || !sender || !receiverGene || !senderGene) {
            alert('Please select all parameters');
            return;
        }

        const selectedRegions = getSelectedRegionsInfo();
        const allCellIds = selectedRegions.flatMap(r => r.cellIds);
        const sampleIds = [...new Set(selectedRegions.map(r => r.sampleId))];

        const params = {
            sample_ids: sampleIds,
            receiver,
            sender,
            receiverGene,
            senderGene,
            cellIds: allCellIds
        };

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

    useEffect(() => {
        if (Object.values(analyzedRegion).some(v => v)) {
            findCellIdsByRegion();
        }
    }, [analyzedRegion, regions]);

    // A simple filter function (case-insensitive) for the Selects.
    const filterSelect = (input, option) =>
        option?.label?.toLowerCase().includes(input.toLowerCase());

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div style={{ margin: 5, fontSize: 14, fontWeight: 'bold' }}>Cell to Cell interaction Viewer</div>
            {/* Selectors and Button */}
            <div style={{ width: '100%', height: '10%', display: 'flex', justifyContent: 'center', alignItems: 'center', gap: 5, marginTop: 10 }}>
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size='small'
                    value={receiver}
                    options={cellTypeList}
                    onChange={setReceiver}
                    disabled={!analyzedRegion || cell2cellDataLoading}
                    placeholder="Select Receiver"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size='small'
                    value={sender}
                    options={cellTypeList}
                    onChange={setSender}
                    disabled={!analyzedRegion || cell2cellDataLoading}
                    placeholder="Select Sender"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size='small'
                    value={receiverGene}
                    options={genelist}
                    onChange={setReceiverGenes}
                    disabled={!analyzedRegion || cell2cellDataLoading}
                    placeholder="Select Receiver Gene"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size='small'
                    value={senderGene}
                    options={genelist}
                    onChange={setSenderGenes}
                    disabled={!analyzedRegion || cell2cellDataLoading}
                    placeholder="Select Sender Gene"
                    filterOption={filterSelect}
                />
                <Button
                    style={{ width: 120 }}
                    size='small'
                    disabled={!analyzedRegion || cell2cellDataLoading}
                    onClick={confirmCell2cellparameters}
                >
                    Confirm
                </Button>
            </div>

            {cell2cellDataLoading ? (
                <Spin spinning={true} style={{ width: '100%', height: '100%' }} />
            ) : (
                <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                    {!analyzedRegion && (
                        <Empty
                            image={Empty.PRESENTED_IMAGE_SIMPLE}
                            description="Please set receiver and sender first"
                            style={{ width: '100%', height: '100%' }}
                        />
                    )}

                    {/* Show content only if cell2cellData is available */}
                    {Object.keys(cell2cellData).length > 0 && analyzedRegion && (
                        <div>
                            <p>Content goes here...</p>
                        </div>
                    )}
                </div>
            )}
        </div>
    );
};
