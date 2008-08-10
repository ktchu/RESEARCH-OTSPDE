
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTSPDE-Systems.ps 
	$(RM) *.aux *.log

allpdf: OTSPDE-Systems.pdf 
	$(RM) *.aux *.log

OTSPDE-Systems.ps:	OTSPDE-Systems.tex \
                    OTSPDE-Systems.bib \
                    macros.tex
	$(TEX) OTSPDE-Systems.tex
	$(BIBTEX) OTSPDE-Systems.aux
	$(TEX) OTSPDE-Systems.tex
	$(TEX) OTSPDE-Systems.tex
	$(DVIPS) OTSPDE-Systems.dvi -o OTSPDE-Systems.ps

OTSPDE-Systems.pdf:  OTSPDE-Systems.tex \
                     OTSPDE-Systems.bib \
                     macros.tex
	$(PDFLATEX) OTSPDE-Systems.tex
	$(BIBTEX) OTSPDE-Systems.aux
	$(PDFLATEX) OTSPDE-Systems.tex
	$(PDFLATEX) OTSPDE-Systems.tex
	$(PDFLATEX) OTSPDE-Systems.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
