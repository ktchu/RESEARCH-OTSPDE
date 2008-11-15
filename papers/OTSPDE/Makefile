
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

OTSPDE.ps:	OTSPDE.tex 
	$(TEX) OTSPDE.tex
	$(TEX) OTSPDE.tex
	$(TEX) OTSPDE.tex
	$(DVIPS) OTSPDE.dvi -o OTSPDE.ps

OTSPDE.pdf:  OTSPDE.tex 
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
