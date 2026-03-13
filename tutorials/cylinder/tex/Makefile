TEX  = pdflatex -shell-escape -interaction=nonstopmode
BIB  = bibtex
SRC  = main

all: $(SRC).pdf

$(SRC).pdf: $(SRC).tex bib.bib
	$(TEX) $(SRC)
	$(BIB) $(SRC) || true
	$(TEX) $(SRC)
	$(TEX) $(SRC)

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pyg *.fls *.fdb_latexmk
	rm -f *.markdown.in *.markdown.lua
	rm -rf _minted-main _markdown_main

.PHONY: all clean
