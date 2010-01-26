SHELL=/bin/bash
LATEX=pdflatex

# Which document file should I construct?
NOMBRE=$(shell grep -H ^[^%].*\\documentclass *.tex | cut -d: -f1 | cut -d. -f1)

DEPENDENCIES=$(NOMBRE).tex chapters/*.tex

ifeq ($(strip $(wildcard $(NOMBRE).bib) ),)
HAVE_BIB=no
else
HAVE_BIB=yes
endif

ifeq '$(HAVE_BIB)' 'yes'
DEPENDENCIES+= $(strip $(wildcard $(NOMBRE).bib) )
endif

IMGDIR = img

.PHONY: all bib clean distclean help

all: $(NOMBRE).pdf

bib: clean
	@echo "   Making document to get needed references..."
	@$(LATEX) -halt-on-error -interaction=nonstopmode $(NOMBRE).tex &> .my_log || (cat .my_log && rm .my_log && exit 1)
	@rm .my_log
	@echo "   Generating bibliography..."
	@bibtool -q -x $(NOMBRE).aux $(BIBDIR)/*.bib > $(NOMBRE).bib 

$(NOMBRE).pdf: $(DEPENDENCIES)
	@echo "   Making document for first time..."
	@$(LATEX) -halt-on-error -interaction=nonstopmode $(NOMBRE).tex &> .my_log || (cat .my_log && rm .my_log && exit 1)
	@rm .my_log
	@if [[ "${HAVE_BIB}" == "yes" ]]; then \
	echo "   Compiling references..." && \
	(bibtex $(NOMBRE) &> .my_log || (cat .my_log && rm .my_log && exit 1)) && \
	rm .my_log && \
	echo "   Re-making document including bibliography..." && \
	$(LATEX) -halt-on-error -interaction=nonstopmode $(NOMBRE).tex &> .my_log || (cat .my_log && rm .my_log && exit 1) && \
	rm .my_log; \
	fi
	@echo "   Re-making document for satisfying references..."
	@$(LATEX) $(NOMBRE) &> /dev/null

clean:
	@echo "   Deleting auxiliary files..."
	-@rm -f $(NOMBRE).{aux,toc,log,tmp,dvi,idx,ilg,ind,bbl,blg,out,lof,lot} .my_log

distclean: clean
	@echo "   Deleting all generated files..."
	@-rm -f $(NOMBRE).pdf

