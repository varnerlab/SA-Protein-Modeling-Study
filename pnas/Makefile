# Makefile for the PNAS version of the paper.
#
# Builds two documents that share references.bib:
#   main.pdf         the research article (pnas-new class)
#   si-appendix.pdf  the supporting information appendix
# latexmk drives pdflatex and bibtex, rerunning as needed to settle
# cross-references and citations.
#
# Usage:
#   make            build both main.pdf and si-appendix.pdf (default)
#   make main       build only main.pdf
#   make si         build only si-appendix.pdf
#   make clean      remove auxiliary files, keep the PDFs
#   make cleanall   remove auxiliary files and the PDFs
#
# Requires: latexmk, pdflatex, bibtex.

BIB     := references.bib
CLASS   := pnas-new.cls pnasresearcharticle.sty
LATEXMK := latexmk -pdf -interaction=nonstopmode

.PHONY: all main si clean cleanall

all: main.pdf si-appendix.pdf

main: main.pdf
si:   si-appendix.pdf

main.pdf: main.tex $(BIB) $(CLASS)
	$(LATEXMK) main.tex

si-appendix.pdf: si-appendix.tex $(BIB)
	$(LATEXMK) si-appendix.tex

clean:
	latexmk -c main.tex si-appendix.tex

cleanall:
	latexmk -C main.tex si-appendix.tex
