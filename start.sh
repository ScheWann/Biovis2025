#!/bin/bash

set -e

cd Backend
echo "Activating virtual environment..."
source venv/bin/activate

echo "Running the Backend..."
python server.py &


cd ..

cd Frontend
echo "Running the Frontend..."
npm start &

echo "Project started successfully!"

wait