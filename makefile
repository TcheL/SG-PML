file := SG-PML
sufs := aux \
        bbl bcf blg \
		fdb_latexmk fls \
		idx ind ilg \
		listing loc lof log lol los lot ltx \
		nav nlo nls \
		out \
		toc \
		run.xml \
		snm synctex\(busy\) synctex.gz \
		vrb \
		xdv

all : tex2pdf view

cc : clean clear

tex2pdf :
	xelatex $(file)

full :
	latexmk -pdfxe -interaction=nonstopmode $(file)

view :
	evince $(file).pdf &

clean :
	-rm -f $(foreach suf,$(sufs),$(file).$(suf))
	-rm -f $(foreach suf,$(sufs),Figure/*.$(suf))

clear :
	-rm -f $(file).pdf

# vim:ft=make:noet:ts=4
