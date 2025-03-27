#!/bin/bash

set -e


REPO_URL="https://github.com/ScheWann/Spacia.git"
BACKEND="Backend"
FRONTEND="Frontend"
PYTHON_ENV="venv"

echo "Cloning repository from $REPO_URL..."
git clone "$REPO_URL"

cd "$BACKEND" || { echo "Error: Cannot enter directory $BACKEND"; exit 1; }

echo "Creating virtual environment..."
python3 -m venv "$PYTHON_ENV"

echo "Activating virtual environment..."
source "$PYTHON_ENV/bin/activate"

if [ -f "requirements.txt" ]; then
    echo "Installing dependencies from requirements.txt..."
    pip install --upgrade pip
    pip install -r requirements.txt
else
    echo "Warning: requirements.txt not found, skipping dependency installation."
fi

cd ..

if [ -d "$FRONTEND" ]; then
    cd "$FRONTEND"
    echo "Installing frontend dependencies with npm..."
    npm install

else
    echo "Error: Frontend directory $FRONTEND not found."
    exit 1
fi

echo "Project setup completed successfully!"
