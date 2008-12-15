
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTS-MaxwellsEqns.ps 
	$(RM) *.aux *.log

allpdf: OTS-MaxwellsEqns.pdf 
	$(RM) *.aux *.log

OTS-MaxwellsEqns.ps:	OTS-MaxwellsEqns.tex \
                      OTS-MaxwellsEqns.bib \
                      macros.tex
	$(TEX) OTS-MaxwellsEqns.tex
	$(BIBTEX) OTS-MaxwellsEqns.aux
	$(TEX) OTS-MaxwellsEqns.tex
	$(TEX) OTS-MaxwellsEqns.tex
	$(DVIPS) OTS-MaxwellsEqns.dvi -o OTS-MaxwellsEqns.ps

OTS-MaxwellsEqns.pdf:  OTS-MaxwellsEqns.tex \
                       OTS-MaxwellsEqns.bib \
                       macros.tex
	$(PDFLATEX) OTS-MaxwellsEqns.tex
	$(BIBTEX) OTS-MaxwellsEqns.aux
	$(PDFLATEX) OTS-MaxwellsEqns.tex
	$(PDFLATEX) OTS-MaxwellsEqns.tex
	$(PDFLATEX) OTS-MaxwellsEqns.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
