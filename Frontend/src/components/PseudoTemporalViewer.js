import React, { useRef, useState, useEffect } from 'react';
import { Box, Button, TextField, Stack, CircularProgress, Alert, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import axios from 'axios';

const ControlPanel = styled(Box)(({ theme }) => ({
    position: 'absolute',
    bottom: theme.spacing(1),
    right: theme.spacing(1),
    display: 'flex',
    flexDirection: 'column',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(1),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
}));

const NavigationButtons = styled(Box)(({ theme }) => ({
    position: 'absolute',
    bottom: theme.spacing(1),
    left: theme.spacing(1),
    display: 'flex',
    gap: theme.spacing(0.5),
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    padding: theme.spacing(0.5),
    borderRadius: theme.shape.borderRadius,
    boxShadow: theme.shadows[2],
    zIndex: 1000,
    '& .MuiButton-root': {
        padding: '2px 4px',
        minWidth: '60px',
        fontSize: '0.75rem'
    }
}));

const ResultsPanel = styled(Box)(({ theme }) => ({
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    padding: theme.spacing(1),
    gap: theme.spacing(1),
    overflow: 'auto',
}));

const ImageGrid = styled(Box)(({ theme }) => ({
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
    gap: theme.spacing(2),
    marginTop: theme.spacing(2),
}));

const ImageContainer = styled(Box)(({ theme }) => ({
    position: 'relative',
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    '& img': {
        width: '100%',
        height: '100%',
        objectFit: 'contain',
    },
}));

export const PseudoTemporalViewer = () => {
    const [samplePercent, setSamplePercent] = useState(1);
    const [currentStep, setCurrentStep] = useState(0);
    const [results, setResults] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [availableImages, setAvailableImages] = useState([]);
    const [imageTitles, setImageTitles] = useState([]);
    const [currentImageIndex, setCurrentImageIndex] = useState(0);

    useEffect(() => {
        // Fetch available images when component mounts
        const fetchImages = async () => {
            try {
                const response = await axios.get('/test-image');
                const images = response.data.images || [];
                setAvailableImages(images);
                
                // Create titles for each image
                const titles = images.map(img => {
                    if (img.startsWith('gene_expression_') && !img.includes('heatmap')) {
                        return `Gene Expression: ${img.replace('gene_expression_', '').replace('.png', '')}`;
                    } else if (img === 'umap_cell_types.png') {
                        return 'UMAP: Cell Types';
                    } else if (img === 'umap_gene_expression.png') {
                        return 'UMAP: Gene Expression';
                    } else if (img === 'gene_expression_heatmap.png') {
                        return 'Gene Expression Heatmap';
                    }
                    return img;
                });
                setImageTitles(titles);
            } catch (error) {
                console.error('Error fetching images:', error);
                setAvailableImages([]);
                setImageTitles([]);
            }
        };
        fetchImages();
    }, []);

    const fetchResults = async () => {
        try {
            setLoading(true);
            setError(null);
            
            const response = await axios.get('/get_deaplog_results', {
                params: {
                    sample_percent: samplePercent,
                    step: currentStep
                }
            });
            
            console.log('Received results:', response.data);
            setResults(response.data);
        } catch (err) {
            console.error('Error fetching results:', err);
            setError(err.response?.data?.error || 'Failed to fetch results');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        fetchResults();
    }, [samplePercent, currentStep]);

    const handleSampleChange = (event) => {
        const value = Math.min(Math.max(0.1, parseFloat(event.target.value) || 0.1), 100);
        setSamplePercent(value);
    };

    const handleNext = () => {
        if (availableImages.length > 0) {
            setCurrentImageIndex((prev) => (prev + 1) % availableImages.length);
        }
    };

    const handlePrevious = () => {
        if (availableImages.length > 0) {
            setCurrentImageIndex((prev) => (prev - 1 + availableImages.length) % availableImages.length);
        }
    };

    const currentImage = availableImages.length > 0 ? availableImages[currentImageIndex] : null;
    const currentTitle = imageTitles.length > 0 ? imageTitles[currentImageIndex] : '';
    const imageExists = currentImage && availableImages.includes(currentImage);

    const renderResults = () => {
        if (!results) return null;

        return (
            <Box sx={{ p: 2 }}>
                <Typography variant="h6">Analysis Results</Typography>
                {/* <Box sx={{ mt: 2 }}>
                    <Typography variant="body1">
                        {currentGene} ({currentImageIndex + 1} of {testGenes.length})
                    </Typography>
                </Box> */}
                
                <Box sx={{ 
                    display: 'flex', 
                    justifyContent: 'center', 
                    alignItems: 'center',
                    height: 'calc(33vh - 100px)',
                    mt: 2 
                }}>
                    <ImageContainer>
                        <Typography variant="subtitle2" sx={{ mb: 1, textAlign: 'center' }}>
                            {currentTitle} {!imageExists && '(Image not found)'}
                        </Typography>
                        {imageExists && (
                            <img 
                                src={`/figures/${currentImage}?t=${Date.now()}`}
                                alt={currentTitle}
                                loading="lazy"
                                style={{ 
                                    width: '100%',
                                    height: '100%',
                                    objectFit: 'contain'
                                }}
                                onError={(e) => {
                                    console.error('Failed to load image:', e.target.src);
                                    e.target.style.display = 'none';
                                }}
                            />
                        )}
                    </ImageContainer>
                </Box>
            </Box>
        );
    };

    return (
        <Box sx={{ 
            height: '33vh', 
            width: '100%', 
            position: 'relative',
            display: 'flex', 
            flexDirection: 'column',
            borderBottom: '2px solid #e8e8e8',
            gap: 0
        }}>
            <ResultsPanel>
                {loading ? (
                    <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
                        <CircularProgress />
                    </Box>
                ) : error ? (
                    <Alert severity="error" sx={{ m: 1 }}>{error}</Alert>
                ) : (
                    renderResults()
                )}
            </ResultsPanel>

            <ControlPanel>
                <TextField
                    label="Sample %"
                    type="number"
                    value={samplePercent}
                    onChange={handleSampleChange}
                    size="small"
                    inputProps={{ 
                        min: 0.1, 
                        max: 100, 
                        step: 0.1 
                    }}
                    sx={{ width: 80 }}
                />
            </ControlPanel>

            <NavigationButtons>
                <Button 
                    variant="contained" 
                    onClick={handlePrevious}
                    disabled={loading}
                    sx={{ m: 0 }}
                >
                    Previous
                </Button>
                <Button 
                    variant="contained" 
                    onClick={handleNext}
                    disabled={loading}
                    sx={{ m: 0 }}
                >
                    Next
                </Button>
            </NavigationButtons>
        </Box>
    );
};