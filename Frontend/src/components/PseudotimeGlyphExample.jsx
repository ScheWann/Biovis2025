import React, { useState } from 'react';
import PseudotimeGlyph from './PseudotimeGlyph';
import { Button, Card, Space, InputNumber, Select } from 'antd';

const { Option } = Select;

export const PseudotimeGlyphExample = () => {
    const [sampleId, setSampleId] = useState('sample1');
    const [cellIds] = useState([1, 2, 3, 4, 5]); // Example cell IDs
    const [adata_umap_title] = useState('test_umap');
    const [dimensions, setDimensions] = useState({ width: 500, height: 500 });

    // Example gene expression data (this would come from your backend in real usage)
    const exampleGeneData = [
        {
            gene: "NANOG",
            timePoints: [0.0, 0.3, 0.6, 1.0],
            expressions: [0.9, 0.7, 0.4, 0.1]
        },
        // {
        //     gene: "SOX2",
        //     timePoints: [0.0, 0.4, 0.7, 1.0],
        //     expressions: [0.8, 0.6, 0.3, 0.2]
        // },
        {
            gene: "OCT4",
            timePoints: [0.0, 0.2, 0.5, 0.8],
            expressions: [0.95, 0.8, 0.5, 0.15]
        }
    ];

    return (
        <div style={{ padding: '20px' }}>
            <Card
                title="Pseudotime Glyph Visualization"
                style={{ marginBottom: '20px' }}
                extra={
                    <Space>
                        <span>Sample:</span>
                        <Select value={sampleId} onChange={setSampleId} style={{ width: 120 }}>
                            <Option value="sample1">Sample 1</Option>
                            <Option value="sample2">Sample 2</Option>
                            <Option value="sample3">Sample 3</Option>
                        </Select>
                    </Space>
                }
            >
                <div style={{ marginBottom: '15px' }}>
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
                </div>

                <div style={{
                    display: 'flex',
                    justifyContent: 'center',
                    marginTop: '20px',
                    padding: '20px',
                    backgroundColor: '#fafafa',
                    borderRadius: '8px'
                }}>
                    <PseudotimeGlyph
                        sampleId={sampleId}
                        cellIds={cellIds}
                        adata_umap_title={adata_umap_title}
                        width={dimensions.width}
                        height={dimensions.height}
                        geneExpressionData={exampleGeneData}
                        early_markers={['NANOG', 'SOX2', 'OCT4']}
                    />
                </div>

                <div style={{ marginTop: '20px', padding: '15px', backgroundColor: '#f0f8ff', borderRadius: '6px' }}>
                    <h4>Glyph Components:</h4>
                    <ul>
                        <li><strong>Central Time Axis:</strong> Vertical line dividing top (gene expression) and bottom (cell trajectories)</li>
                        <li><strong>Bottom Section:</strong> Macroscopic cell trajectory paths with colored circles for intermediate stages and stars for final stages</li>
                        <li><strong>Top Section:</strong> Gauge-like visualization showing gene expression changes over pseudotime</li>
                        <li><strong>Radial Direction:</strong> Distance from center represents pseudotime progression</li>
                        <li><strong>Color Coding:</strong> Different colors for different trajectories and genes</li>
                    </ul>
                </div>
            </Card>
        </div>
    );
};