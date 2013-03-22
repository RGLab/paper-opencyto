# Makefile to use knitr for package vignettes
# Original source from Yihui Xie:
# http://yihui.github.com/knitr/demo/vignette/

# For more about Makefiles and R, see:
# http://www.stat.auckland.ac.nz/~stat782/downloads/make-tutorial.pdf

# Put all PDF targets here, separated by spaces.
PDFS= quantiles-summary-HVTN065.pdf
#PDFS= quantiles-summary-HVTN065.pdf gates-HVTN065.pdf summary-HVTN065.pdf gates-HVTN080.pdf summary-HVTN080.pdf

all: $(PDFS) 

clean:
	rm -f *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary demo*

%.pdf: %.Rnw
	R -e "library(knitr); knit2pdf('$*.Rnw')"
