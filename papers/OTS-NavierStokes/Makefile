
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	OTS-NavierStokes.ps 
	$(RM) *.aux *.log

allpdf: OTS-NavierStokes.pdf 
	$(RM) *.aux *.log

OTS-NavierStokes.ps:	OTS-NavierStokes.tex \
                      OTS-NavierStokes.bib \
                      macros.tex
	$(TEX) OTS-NavierStokes.tex
	$(BIBTEX) OTS-NavierStokes.aux
	$(TEX) OTS-NavierStokes.tex
	$(TEX) OTS-NavierStokes.tex
	$(DVIPS) OTS-NavierStokes.dvi -o OTS-NavierStokes.ps

OTS-NavierStokes.pdf:  OTS-NavierStokes.tex \
                       OTS-NavierStokes.bib \
                       macros.tex
	$(PDFLATEX) OTS-NavierStokes.tex
	$(BIBTEX) OTS-NavierStokes.aux
	$(PDFLATEX) OTS-NavierStokes.tex
	$(PDFLATEX) OTS-NavierStokes.tex
	$(PDFLATEX) OTS-NavierStokes.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
