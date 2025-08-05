import PseudotimeGlyph from './PseudotimeGlyph';


export const PseudotimeGlyphComponent = ({ 
    adata_umap_title,
    pseudotimeData,
    pseudotimeLoading
}) => {
    return (
        <div style={{ padding: '10px' }}>
                <div style={{
                    marginTop: '10px',
                    display: 'flex',
                    gap: '20px',
                    flexWrap: 'wrap',
                    justifyContent: 'center'
                }}>
                    <div style={{ textAlign: 'center' }}>
                        <PseudotimeGlyph
                            adata_umap_title={adata_umap_title}
                            pseudotimeData={pseudotimeData}
                            pseudotimeLoading={pseudotimeLoading}
                        />
                    </div>
                </div>
        </div>
    );
};