
all: docx html
docx: README.docx
#pdf: README.pdf
html: README.html

%.docx: %.md
	pandoc -f markdown -t docx -o $@ $<

#%.pdf: %.tex
#	latex README
#	bibtex README
#	latex README
#	latex README
#	dvips README.dvi
#	ps2pdf README.ps

#%.tex: %.md
#	pandoc -f markdown -t latex -o $@ $<

%.html: %.md
	pandoc -f markdown -t html -o $@ $<

clean:
	rm -f README.html
	rm -f README.docx
#	rm -f README.pdf
#	rm -f README.ps
#	rm -f README.dvi

