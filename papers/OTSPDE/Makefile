
RM = rm -f
TEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex

all: 	OTSPDE.pdf OTSPDE-supplemental.pdf
	$(RM) *.aux *.log

OTSPDE.pdf: OTSPDE.tex
	$(PDFLATEX) $^
	#$(BIBTEX) OTSPDE.aux
	$(PDFLATEX) $^
	$(PDFLATEX) $^

OTSPDE-supplemental.pdf: OTSPDE-supplemental.tex
	$(PDFLATEX) $^
	$(PDFLATEX) $^
	$(PDFLATEX) $^

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
