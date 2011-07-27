
RM = rm -f
TEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex

all: 	OTSPDE.pdf
	$(RM) *.aux *.log

OTSPDE.pdf: OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex
	#$(BIBTEX) OTSPDE.aux
	$(PDFLATEX) OTSPDE.tex
	$(PDFLATEX) OTSPDE.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
