#!/bin/bash
# TTPicker startup script
echo ""
echo "  TTPicker — Gene Therapy Analysis Tool"
echo "  ======================================"
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "  ERROR: Python 3 not found. Please install Python 3.8+"
    exit 1
fi

# Install dependencies if needed
echo "  Checking dependencies..."
pip install flask requests reportlab python-dotenv --break-system-packages -q

# Load .env if present
if [ -f .env ]; then
    echo "  Found .env — loading API key..."
    export $(grep -v '^#' .env | xargs)
fi

echo "  Starting server..."
echo ""
echo "  Open your browser at:  http://localhost:5000"
echo ""
cd "$(dirname "$0")"
python3 app.py
