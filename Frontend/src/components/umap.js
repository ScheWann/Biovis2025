import React, { useState } from 'react';
import { DeckGL, ScatterplotLayer } from 'deck.gl';


const INITIAL_VIEW_STATE = {
    longitude: 0,
    latitude: 0,
    zoom: 4,
    pitch: 0,
    bearing: 0,
};

export const Umap = ({ umapPositionWithClusterData }) => {
    const [hoverInfo, setHoverInfo] = useState(null);

    const scatterplotLayer = new ScatterplotLayer({
        id: 'scatterplot-layer',
        data: umapPositionWithClusterData,
        getPosition: d => [d['UMAP-1'], d['UMAP-2']],
        getFillColor: d => {
            const colors = [
                [255, 99, 71],
                [135, 206, 250],
                [144, 238, 144],
                [238, 130, 238],
            ];
            return colors[d.cluster % colors.length];
        },
        getRadius: 1000,
        pickable: true,
        onHover: info => setHoverInfo(info),
    });


    return (
        <div style={{ height: '100%',  width: '50%', position: 'relative' }}>
            <DeckGL
                initialViewState={INITIAL_VIEW_STATE}
                controller={true}
                layers={[scatterplotLayer]}
            >
                {hoverInfo && hoverInfo.object && (
                    <div
                        style={{
                            position: 'absolute',
                            zIndex: 1,
                            pointerEvents: 'none',
                            left: hoverInfo.x,
                            top: hoverInfo.y,
                            backgroundColor: 'white',
                            padding: '5px',
                            borderRadius: '3px',
                            fontSize: '12px',
                            border: '1px solid gray',
                        }}
                    >
                        <div>{`Barcode: ${hoverInfo.object.barcode}`}</div>
                        <div>{`Cluster: ${hoverInfo.object.clusterluster}`}</div>
                        <div>{`UMAP-1: ${hoverInfo.object['UMAP-1'].toFixed(2)}`}</div>
                        <div>{`UMAP-2: ${hoverInfo.object['UMAP-2'].toFixed(2)}`}</div>
                    </div>
                )}
            </DeckGL>
        </div>
    );
}