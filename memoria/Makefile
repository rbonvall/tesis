SHELL=/bin/bash
LATEX=pdflatex

# Which document file should I construct?
DOCUMENT=$(shell grep -H ^[^%].*\\documentclass *.tex | cut -d: -f1 | cut -d. -f1)

DEPENDENCIES=$(DOCUMENT).tex chapters/*.tex

ifeq ($(strip $(wildcard *.bib) ),)
  HAVE_BIB=no
else
  HAVE_BIB=yes
endif

ifeq '$(HAVE_BIB)' 'yes'
  DEPENDENCIES+= $(strip $(wildcard *.bib) )
endif

IMGDIR = img

.PHONY: all bib clean distclean help

all: $(DOCUMENT).pdf

$(DOCUMENT).pdf: $(DEPENDENCIES)
	@echo "   Making document for first time..."
	@$(LATEX) -halt-on-error -interaction=nonstopmode $(DOCUMENT).tex &> .my_log || (cat .my_log && rm .my_log && exit 1)
	@rm .my_log
	@if [[ "${HAVE_BIB}" == "yes" ]]; then \
	    echo "   Compiling references..." && \
	    (bibtex $(DOCUMENT) &> .my_log || (cat .my_log && rm .my_log && exit 1)) && \
	    rm .my_log && \
	    echo "   Re-making document including bibliography..." && \
	    $(LATEX) -halt-on-error -interaction=nonstopmode $(DOCUMENT).tex &> .my_log || (cat .my_log && rm .my_log && exit 1) && \
	    rm .my_log; \
	fi
	@echo "   Re-making document for satisfying references..."
	@$(LATEX) $(DOCUMENT) &> /dev/null

clean:
	@echo "   Deleting auxiliary files..."
	-@rm -f $(DOCUMENT).{aux,toc,log,tmp,dvi,idx,ilg,ind,bbl,blg,out,lof,lot} .my_log
	-@rm -f chapters/*.aux

distclean: clean
	@echo "   Deleting all generated files..."
	@-rm -f $(DOCUMENT).pdf

