# type "make" command in Unix to create asme2ej.pdf file
all:
	latex writeup
	bibtex writeup
	latex writeup
	latex writeup
	dvips -o writeup.ps writeup
	ps2pdf writeup.ps writeup.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
