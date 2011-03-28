TEXFILES = $(wildcard *.tex)
PDFFILES = $(TEXFILES:.tex=.pdf)

pdf: $(PDFFILES)

%.pdf: %.tex
	@rubber --pdf $<
clean:
	@rubber --clean --pdf $(TEXFILES:.tex=)
	@rm -f *.vrb

.PHONY: pdf clean
