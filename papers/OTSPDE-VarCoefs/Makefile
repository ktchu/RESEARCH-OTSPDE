
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTSPDE-VarCoefs.ps 
	$(RM) *.aux *.log

allpdf: OTSPDE-VarCoefs.pdf 
	$(RM) *.aux *.log

OTSPDE-VarCoefs.ps:  OTSPDE-VarCoefs.tex \
                     OTSPDE-VarCoefs.bib \
                     macros.tex
	$(TEX) OTSPDE-VarCoefs.tex
	$(BIBTEX) OTSPDE-VarCoefs.aux
	$(TEX) OTSPDE-VarCoefs.tex
	$(TEX) OTSPDE-VarCoefs.tex
	$(DVIPS) OTSPDE-VarCoefs.dvi -o OTSPDE-VarCoefs.ps

OTSPDE-VarCoefs.pdf:  OTSPDE-VarCoefs.tex \
                      OTSPDE-VarCoefs.bib \
                      macros.tex
	$(PDFLATEX) OTSPDE-VarCoefs.tex
	$(BIBTEX) OTSPDE-VarCoefs.aux
	$(PDFLATEX) OTSPDE-VarCoefs.tex
	$(PDFLATEX) OTSPDE-VarCoefs.tex
	$(PDFLATEX) OTSPDE-VarCoefs.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
