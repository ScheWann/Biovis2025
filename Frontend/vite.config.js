import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
    plugins: [react()],
    server: {
        port: 3003, // or your preferred port
        proxy: {
            '/api': 'http://localhost:5003', // adjust if you use a backend
        },
    },
});