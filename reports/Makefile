# Makefile to use knitr for package vignettes
# Original source from Yihui Xie:
# http://yihui.github.com/knitr/demo/vignette/

# For more about Makefiles and R, see:
# http://www.stat.auckland.ac.nz/~stat782/downloads/make-tutorial.pdf

.PHONY: clean

# Put all PDF targets here, separated by spaces.
PDFS= outline-opencyto.pdf paper-opencyto.pdf NOTES.pdf

all: $(PDFS) 

%.pdf: %.Rnw
	Rscript -e "library(knitr); knit2pdf('$*.Rnw')"

# We create a header file named 'header.tex' that contains two directives.
# First, we set the margins of the LaTeX document generated with the geometry package.
# Second, we turn off the page numbering.
# See the following link for more: http://en.wikibooks.org/wiki/LaTeX/Page_Layout 

%.pdf: %.md
	echo '\usepackage[top=0.75in, bottom=0.75in, left=1in, right=1in]{geometry}' > header.tex
	echo '\pagestyle{empty}' >> header.tex
	echo '\usepackage{setspace, natbib}' >> header.tex
	echo '\doublespacing' >> header.tex
	pandoc -s $*.md -t latex -H header.tex -o $*.pdf
	rm header.tex

clean:
	@echo "Cleaning Up"
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary demo* *.synctex.gz
