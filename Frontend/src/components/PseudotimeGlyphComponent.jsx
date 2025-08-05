import PseudotimeGlyph from './PseudotimeGlyph';


export const PseudotimeGlyphComponent = ({ 
    adata_umap_title,
    pseudotimeDataSets,
    pseudotimeLoadingStates
}) => {
    // Convert pseudotimeDataSets object to array of all trajectory data
    const allPseudotimeData = [];
    Object.entries(pseudotimeDataSets).forEach(([title, dataArray]) => {
        if (Array.isArray(dataArray)) {
            dataArray.forEach((trajectoryData, index) => {
                allPseudotimeData.push({
                    ...trajectoryData,
                    source_title: title,
                    display_title: `${title} - Trajectory ${index + 1}`,
                    isLoading: pseudotimeLoadingStates[title] || false
                });
            });
        }
    });

    // Add loading placeholders for datasets that are currently being loaded
    Object.entries(pseudotimeLoadingStates).forEach(([title, isLoading]) => {
        if (isLoading && !pseudotimeDataSets[title]) {
            allPseudotimeData.push({
                source_title: title,
                display_title: `${title} - Loading...`,
                isLoading: true,
                isPlaceholder: true
            });
        }
    });

    // Check if there's any global loading happening and no data exists yet
    const anyLoading = Object.values(pseudotimeLoadingStates).some(loading => loading);
    const hasNoData = Object.keys(pseudotimeDataSets).length === 0;

    // If allPseudotimeData is not an array or is empty, show loading or empty state
    if (anyLoading && hasNoData) {
        return (
            <div style={{ padding: '10px', textAlign: 'center' }}>
                Loading pseudotime data...
            </div>
        );
    }

    if (!allPseudotimeData || allPseudotimeData.length === 0) {
        return (
            <div style={{ padding: '10px', textAlign: 'center' }}>
                No pseudotime data available
            </div>
        );
    }

    // Calculate responsive dimensions
    const numGlyphs = allPseudotimeData.length;
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
            {/* Summary header */}
            <div style={{ 
                marginBottom: '10px', 
                textAlign: 'center', 
                fontSize: '12px', 
                color: '#666',
                fontWeight: 'bold'
            }}>
                {Object.keys(pseudotimeDataSets).length} dataset(s), {allPseudotimeData.length} trajectory(s)
            </div>
            
            <div style={{
                display: 'grid',
                gridTemplateColumns: `repeat(${glyphsPerRow}, 1fr)`,
                gridTemplateRows: `repeat(${numRows}, 1fr)`,
                gap: `${gapSize}px`,
                width: '100%',
                height: `calc(100% - 40px)`, // Account for padding and header
                justifyItems: 'center',
                alignItems: 'center'
            }}>
                {allPseudotimeData.map((trajectoryData, index) => (
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
                        {trajectoryData.isPlaceholder ? (
                            <div style={{
                                width: '100%',
                                height: '100%',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center',
                                border: '2px dashed #ccc',
                                borderRadius: '8px',
                                color: '#666',
                                fontSize: '14px'
                            }}>
                                Loading {trajectoryData.source_title}...
                            </div>
                        ) : (
                            <PseudotimeGlyph
                                adata_umap_title={trajectoryData.display_title || `${adata_umap_title} - Trajectory ${index + 1}`}
                                pseudotimeData={[trajectoryData]}
                                pseudotimeLoading={trajectoryData.isLoading}
                            />
                        )}
                    </div>
                ))}
                        </div>
        </div>
    );
};