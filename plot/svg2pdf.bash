#!/bin/bash

# Exit on errors and undefined variables
set -euo pipefail

# Check if Inkscape is installed
if ! command -v inkscape &> /dev/null; then
    echo "Error: inkscape is not installed. Please install it first."
    exit 1
fi

# Convert all SVG files to PDF
for file in *.svg; do
    output="${file%.svg}.pdf"
    inkscape "$file" --export-pdf="$output"
    echo "Converted: $file → $output"
done

echo "Conversion complete!"
