import React, { useState, useEffect } from 'react';
import { Button, Select, Spin, Empty } from 'antd';

export const Cell2CellViewer2 = ({ regions, cell2cellData, setCell2cellData }) => {
    const [regionOptions, setRegionOptions] = useState([]);
    const [selectedRegion, setSelectedRegion] = useState([]);
    const [cellTypeList, setCellTypeList] = useState([]);
    const [genelist, setGenelist] = useState([]);
    const [receiver, setReceiver] = useState(null);
    const [sender, setSender] = useState(null);
    const [receiverGene, setReceiverGenes] = useState(null);
    const [senderGene, setSenderGenes] = useState(null);

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
        const params = {};
        selectedRegion.forEach(region => {
            const regionInfo = regions.find(r => r.id === region);
            if (regionInfo) {
                params[regionInfo.name] = {
                    sampleid: regionInfo.sampleId,
                    cell_list: regionInfo.cellIds,
                };
            }
        })
        return params;
    };

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

    useEffect(() => {
        regions.forEach(region => {
            setRegionOptions(prevOptions => {
                const exists = prevOptions.find(option => option.value === region.id);
                return exists ? prevOptions : [...prevOptions, { label: region.name, value: region.id }];
            });
        });
    }, [regions]);

    useEffect(() => {
        if(selectedRegion.length > 0) {
            const selectedRegions = regions.filter(region => selectedRegion.includes(region.id));
            const uniqueSampleIds = [...new Set(selectedRegions.map(r => r.sampleId))];

            if (uniqueSampleIds.length > 0) {
                fetchCellTypeList(uniqueSampleIds[0]);
                fetchGeneList(uniqueSampleIds[0]);
            }
        }
    }, [selectedRegion]);

    // A simple filter function (case-insensitive) for the Selects.
    const filterSelect = (input, option) =>
        option?.label?.toLowerCase().includes(input.toLowerCase());

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div style={{ margin: 5, fontSize: 14, fontWeight: 'bold' }}>Cell to Cell interaction Viewer</div>
            {/* Selectors and Button */}
            <div style={{ width: '100%', height: '10%', display: 'flex', justifyContent: 'center', alignItems: 'center', gap: 5, marginTop: 10 }}>
                <Select
                    size='small'
                    mode="multiple"
                    placeholder="Select Region"
                    options={regionOptions}
                    style={{ width: 200}}
                    onChange={setSelectedRegion}
                    allowClear
                />
                <Select
                    showSearch
                    style={{ width: 120 }}
                    size='small'
                    value={receiver}
                    options={cellTypeList}
                    onChange={setReceiver}
                    placeholder="Select Receiver"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size='small'
                    value={sender}
                    options={cellTypeList}
                    onChange={setSender}
                    placeholder="Select Sender"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size='small'
                    value={receiverGene}
                    options={genelist}
                    onChange={setReceiverGenes}
                    placeholder="Select Receiver Gene"
                    filterOption={filterSelect}
                />
                <Select
                    showSearch
                    style={{ width: 100 }}
                    size='small'
                    value={senderGene}
                    options={genelist}
                    onChange={setSenderGenes}
                    placeholder="Select Sender Gene"
                    filterOption={filterSelect}
                />
                <Button
                    style={{ width: 100 }}
                    size='small'
                    onClick={confirmCell2cellparameters}
                >
                    Confirm
                </Button>
            </div>

        </div>
    );
};
