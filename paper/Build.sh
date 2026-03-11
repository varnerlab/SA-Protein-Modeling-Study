#!/bin/bash
# Build the paper (requires pdflatex and bibtex)
set -e
cd "$(dirname "$0")"

pdflatex -interaction=nonstopmode Paper_v1.tex
bibtex Paper_v1
pdflatex -interaction=nonstopmode Paper_v1.tex
pdflatex -interaction=nonstopmode Paper_v1.tex

echo "Done: Paper_v1.pdf"
