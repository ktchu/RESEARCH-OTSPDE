
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	cover_letter.ps 
	$(RM) *.aux *.log

allpdf: cover_letter.pdf 
	$(RM) *.aux *.log

cover_letter.ps:	cover_letter.tex 
	$(TEX) cover_letter.tex
	$(DVIPS) cover_letter.dvi -o cover_letter.ps

cover_letter.pdf:  cover_letter.tex
	$(PDFLATEX) cover_letter.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
