
RM = rm -f
TEX = latex
BIBTEX = bibtex
DVIPS = dvips -t letter
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	optimal_time_steps_for_pdes.ps 
	$(RM) *.aux *.log

allpdf: optimal_time_steps_for_pdes.pdf 
	$(RM) *.aux *.log

optimal_time_steps_for_pdes.ps:	optimal_time_steps_for_pdes.tex \
                                optimal_time_steps_for_pdes.bib \
                                macros.tex
	$(TEX) optimal_time_steps_for_pdes.tex
	$(BIBTEX) optimal_time_steps_for_pdes.aux
	$(TEX) optimal_time_steps_for_pdes.tex
	$(TEX) optimal_time_steps_for_pdes.tex
	$(DVIPS) optimal_time_steps_for_pdes.dvi -o optimal_time_steps_for_pdes.ps

optimal_time_steps_for_pdes.pdf:  optimal_time_steps_for_pdes.tex \
                                  optimal_time_steps_for_pdes.bib \
                                  macros.tex
	$(PDFLATEX) optimal_time_steps_for_pdes.tex
	$(BIBTEX) optimal_time_steps_for_pdes.aux
	$(PDFLATEX) optimal_time_steps_for_pdes.tex
	$(PDFLATEX) optimal_time_steps_for_pdes.tex
	$(PDFLATEX) optimal_time_steps_for_pdes.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
