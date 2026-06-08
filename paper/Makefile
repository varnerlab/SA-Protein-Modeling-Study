# Makefile for the arXiv (NeurIPS-style) version of the paper.
#
# Builds Paper_v1.pdf from Paper_v1.tex, the section files under sections/,
# and the References_v1.bib bibliography. latexmk drives pdflatex and bibtex,
# rerunning as many times as needed to settle cross-references and citations.
#
# Usage:
#   make            build Paper_v1.pdf (default)
#   make clean      remove auxiliary files, keep the PDF
#   make cleanall   remove auxiliary files and the PDF
#
# Requires: latexmk, pdflatex, bibtex.

MAIN     := Paper_v1
BIB      := References_v1.bib
STYLE    := neurips_2026.sty
SECTIONS := $(wildcard sections/*.tex)
LATEXMK  := latexmk -pdf -interaction=nonstopmode

.PHONY: all clean cleanall

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex $(BIB) $(STYLE) $(SECTIONS)
	$(LATEXMK) $(MAIN).tex

clean:
	latexmk -c $(MAIN).tex

cleanall:
	latexmk -C $(MAIN).tex
