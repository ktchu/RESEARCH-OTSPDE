
RM = rm -f
TEX = latex
BIBTEX = bibtex
#DVIPS = dvips -t letter
DVIPS = dvips 
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTSPDE-WCE2009.ps 
	$(RM) *.aux *.log

allpdf: OTSPDE-WCE2009.pdf 
	$(RM) *.aux *.log

OTSPDE-WCE2009.ps:  OTSPDE-WCE2009.tex \
                    OTSPDE-WCE2009.bib
	$(TEX) OTSPDE-WCE2009.tex
	$(BIBTEX) OTSPDE-WCE2009.aux
	$(TEX) OTSPDE-WCE2009.tex
	$(TEX) OTSPDE-WCE2009.tex
	$(DVIPS) OTSPDE-WCE2009.dvi -o OTSPDE-WCE2009.ps

OTSPDE-WCE2009.pdf:  OTSPDE-WCE2009.tex \
                    OTSPDE-WCE2009.bib
	$(PDFLATEX) OTSPDE-WCE2009.tex
	$(BIBTEX) OTSPDE-WCE2009.aux
	$(PDFLATEX) OTSPDE-WCE2009.tex
	$(PDFLATEX) OTSPDE-WCE2009.tex
	$(PDFLATEX) OTSPDE-WCE2009.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
