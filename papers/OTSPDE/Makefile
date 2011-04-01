
RM = rm -f
TEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
DVIPS = dvips -t letter
PS2PDF = ps2pdf

all: 	OTSPDE.pdf
	$(RM) *.aux *.log

OTSPDE.pdf: OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex
#	$(BIBTEX) OTSPDE.aux
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex

OTSPDE.ps:	OTSPDE.tex 
	$(TEX) OTSPDE.tex
#	$(BIBTEX) OTSPDE.aux
	$(TEX) OTSPDE.tex
	$(TEX) OTSPDE.tex
	$(DVIPS) OTSPDE.dvi -o OTSPDE.ps

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
