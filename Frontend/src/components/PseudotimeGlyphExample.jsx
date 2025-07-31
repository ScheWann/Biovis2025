import React, { useState } from 'react';
import PseudotimeGlyph from './PseudotimeGlyph';
import { Button, Card, Space, InputNumber, Select } from 'antd';

const { Option } = Select;

export const PseudotimeGlyphExample = ({ 
    initialSampleId = 'sample1',
    initialCellIds = [1, 2, 3, 4, 5],
    initialUmapTitle = 'test_umap',
    initialDimensions = { width: 400, height: 400 },
    customGeneData = null,
    customEarlyMarkers = null,
    title = null
}) => {
    const [sampleId, setSampleId] = useState(initialSampleId);
    const [cellIds] = useState(initialCellIds);
    const [adata_umap_title] = useState(initialUmapTitle);
    const [dimensions, setDimensions] = useState(initialDimensions);

    // Use custom gene data if provided, otherwise use default example data
    const exampleGeneData = customGeneData || [
        { gene: "SOX2", timePoints: [0.0, 0.4, 0.7, 1.0], expressions: [0.8, 0.6, 0.3, 0.2] },
        { gene: "NANOG", timePoints: [0.1, 0.4, 0.7, 1.0], expressions: [0.2, 0.5, 0.8, 0.9] },
        { gene: "OCT4", timePoints: [0.0, 0.3, 0.6, 1.0], expressions: [0.9, 0.7, 0.4, 0.1] }
    ];

    const earlyMarkers = customEarlyMarkers || ['NANOG', 'SOX2', 'OCT4'];

    return (
        <div style={{ padding: '10px' }}>
            {/* <Card
                title={title}
                size="small"
            > */}
                {/* <div style={{ marginBottom: '15px' }}>
                    <Space>
                        <span>Width:</span>
                        <InputNumber
                            value={dimensions.width}
                            onChange={(value) => setDimensions(prev => ({ ...prev, width: value }))}
                            min={300}
                            max={800}
                            step={50}
                        />
                        <span>Height:</span>
                        <InputNumber
                            value={dimensions.height}
                            onChange={(value) => setDimensions(prev => ({ ...prev, height: value }))}
                            min={300}
                            max={800}
                            step={50}
                        />
                    </Space>
                </div> */}

                <div style={{
                    marginTop: '10px',
                    display: 'flex',
                    gap: '20px',
                    flexWrap: 'wrap',
                    justifyContent: 'center'
                }}>
                    <div style={{ textAlign: 'center' }}>
                        <h4>Component A</h4>
                        <p style={{ fontSize: '12px', color: '#666', marginBottom: '10px' }}>
                            Click trajectories to see independent selection
                        </p>
                        <PseudotimeGlyph
                            sampleId={sampleId}
                            cellIds={cellIds}
                            adata_umap_title={adata_umap_title}
                            width={dimensions.width}
                            height={dimensions.height}
                            geneExpressionData={exampleGeneData}
                            early_markers={earlyMarkers}
                        />
                    </div>
                    <div style={{ textAlign: 'center' }}>
                        <h4>Component B</h4>
                        <p style={{ fontSize: '12px', color: '#666', marginBottom: '10px' }}>
                            Selection state is independent from Component A
                        </p>
                        <PseudotimeGlyph
                            sampleId={sampleId + "_copy"}
                            cellIds={cellIds}
                            adata_umap_title={adata_umap_title}
                            width={dimensions.width}
                            height={dimensions.height}
                            geneExpressionData={exampleGeneData}
                            early_markers={earlyMarkers}
                        />
                    </div>
                </div>
            {/* </Card> */}
        </div>
    );
};