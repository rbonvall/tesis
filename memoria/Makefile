SHELL = /bin/bash
LATEX = pdflatex
TEXFLAGS = --halt-on-error --interaction nonstopmode &> $$T
ON_SUCCESS = (rm $$T; echo done)
ON_ERROR = (cat $$T; rm $$T; exit 1)

DOCUMENT := $(basename $(wildcard *.tex))
DEPENDENCIES := $(DOCUMENT).tex chapters/*.tex

ifeq ($(strip $(wildcard *.bib) ),)
  HAVE_BIB = no
else
  HAVE_BIB = yes
endif
ifeq '$(HAVE_BIB)' 'yes'
  DEPENDENCIES += $(strip $(wildcard *.bib) )
endif

.PHONY: all bib clean distclean help

all: $(DOCUMENT).pdf

$(DOCUMENT).pdf: $(DEPENDENCIES)
	@echo -n "Making $@ for first time... "
	@T=$$(mktemp); \
	$(LATEX) $(TEXFLAGS) $< && $(ON_SUCCESS) || $(ON_ERROR)
	@if [[ "${HAVE_BIB}" == "yes" ]]; then \
	    echo "Compiling references..." && \
	    (bibtex $(DOCUMENT) &> .my_log || (cat .my_log && rm .my_log && exit 1)) && \
	    rm .my_log && \
	    echo -n "Re-making $@ including bibliography..." && \
	    T=$$(mktemp); \
	    $(LATEX) $(TEXFLAGS) $< && $(ON_SUCCESS) || $(ON_ERROR) \
	fi
	@echo -n "Re-making $@ for satisfying references... "
	@T=$$(mktemp); \
	$(LATEX) $(TEXFLAGS) $< && $(ON_SUCCESS) || $(ON_ERROR)

clean:
	@echo "Deleting auxiliary files..."
	-@rm -f $(DOCUMENT).{aux,toc,log,tmp,dvi,idx,ilg,ind,bbl,blg,out,lof,lot} .my_log
	-@rm -f chapters/*.aux

distclean: clean
	@echo "Deleting all generated files..."
	@-rm -f $(DOCUMENT).pdf

