#!/bin/bash
# Build PNAS main paper and SI appendix
set -e

echo "=== Building main paper ==="
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex

echo ""
echo "=== Building SI Appendix ==="
pdflatex -interaction=nonstopmode si-appendix.tex
bibtex si-appendix
pdflatex -interaction=nonstopmode si-appendix.tex
pdflatex -interaction=nonstopmode si-appendix.tex

echo ""
echo "=== Done ==="
echo "Main paper: main.pdf"
echo "SI Appendix: si-appendix.pdf"
