
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTSPDE.ps 
	$(RM) *.aux *.log

allpdf: OTSPDE.pdf 
	$(RM) *.aux *.log

OTSPDE.ps:	OTSPDE.tex \
                                OTSPDE.bib \
                                macros.tex
	$(TEX) OTSPDE.tex
	$(BIBTEX) OTSPDE.aux
	$(TEX) OTSPDE.tex
	$(TEX) OTSPDE.tex
	$(DVIPS) OTSPDE.dvi -o OTSPDE.ps

OTSPDE.pdf:  OTSPDE.tex \
                                  OTSPDE.bib \
                                  macros.tex
	$(PDFLATEX) OTSPDE.tex
	$(BIBTEX) OTSPDE.aux
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
