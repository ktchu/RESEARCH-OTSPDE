
RM = rm -f
TEX = latex
PDFLATEX = pdflatex

.SUFFIXES:	.ps .pdf

all: 	cover_letter.pdf responses_to_reviewers.pdf
	$(RM) *.aux *.log

cover_letter.pdf:	cover_letter.tex 
	$(PDFLATEX) cover_letter.tex

responses_to_reviewers.pdf:	responses_to_reviewers.tex
	$(PDFLATEX) responses_to_reviewers.tex

clean:
	$(RM) *.dvi *.aux *.log *.ps *.pdf *.bbl *.blg
