# Makefile to use knitr for package vignettes
# Original source from Yihui Xie:
# http://yihui.github.com/knitr/demo/vignette/

# For more about Makefiles and R, see:
# http://www.stat.auckland.ac.nz/~stat782/downloads/make-tutorial.pdf

# Put all PDF targets here, separated by spaces.
REPORTS= classification-summary-HVTN065.html

all: $(REPORTS) 

clean:
	rm -f classification-summary-HVTN065.md classification-summary-HVTN065.Rmd
	rm -f cytokines-HVTN065.md cytokines-HVTN065.Rmd

%.html: %.r
	R -e "library(knitr); spin('$*.r')"

