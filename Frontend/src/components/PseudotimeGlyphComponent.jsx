import PseudotimeGlyph from './PseudotimeGlyph';


export const PseudotimeGlyphComponent = ({ 
    adata_umap_title,
    pseudotimeData,
    pseudotimeLoading
}) => {
    // If pseudotimeData is not an array or is empty, show loading or empty state
    if (pseudotimeLoading) {
        return (
            <div style={{ padding: '10px', textAlign: 'center' }}>
                Loading pseudotime data...
            </div>
        );
    }

    if (!pseudotimeData || !Array.isArray(pseudotimeData) || pseudotimeData.length === 0) {
        return (
            <div style={{ padding: '10px', textAlign: 'center' }}>
                No pseudotime data available
            </div>
        );
    }

    // Calculate responsive dimensions
    const numGlyphs = pseudotimeData.length;
    const maxPerRow = 3;
    const glyphsPerRow = Math.min(numGlyphs, maxPerRow);
    const numRows = Math.ceil(numGlyphs / maxPerRow);
    
    // Calculate width percentage for each glyph (with gap consideration)
    const gapSize = 20; // gap in px
    const padding = 20; // total horizontal padding
    const widthPercentage = `calc(${100 / glyphsPerRow}% - ${(gapSize * (glyphsPerRow - 1)) / glyphsPerRow}px)`;
    
    // Calculate height for each row
    const heightPercentage = `calc(${100 / numRows}% - ${(gapSize * (numRows - 1)) / numRows}px)`;

    return (
        <div style={{ 
            padding: '10px', 
            width: '100%', 
            height: '100%',
            boxSizing: 'border-box'
        }}>
            <div style={{
                display: 'grid',
                gridTemplateColumns: `repeat(${glyphsPerRow}, 1fr)`,
                gridTemplateRows: `repeat(${numRows}, 1fr)`,
                gap: `${gapSize}px`,
                width: '100%',
                height: `calc(100% - 20px)`, // Account for padding
                justifyItems: 'center',
                alignItems: 'center'
            }}>
                {pseudotimeData.map((trajectoryData, index) => (
                    <div 
                        key={index} 
                        style={{ 
                            width: '100%',
                            height: '100%',
                            minWidth: '200px',
                            minHeight: '200px',
                            textAlign: 'center'
                        }}
                    >
                        <PseudotimeGlyph
                            adata_umap_title={`${adata_umap_title} - Trajectory ${index + 1}`}
                            pseudotimeData={[trajectoryData]}
                            pseudotimeLoading={pseudotimeLoading}
                        />
                    </div>
                ))}
            </div>
        </div>
    );
};