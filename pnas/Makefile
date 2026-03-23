.PHONY: all main si clean

all: main si

main: main.pdf

si: si-appendix.pdf

main.pdf: main.tex references.bib
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

si-appendix.pdf: si-appendix.tex references.bib
	pdflatex -interaction=nonstopmode si-appendix.tex
	bibtex si-appendix
	pdflatex -interaction=nonstopmode si-appendix.tex
	pdflatex -interaction=nonstopmode si-appendix.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.synctex.gz *.fls *.fdb_latexmk
