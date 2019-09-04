xlx := latexmk -pdfxe -interaction=nonstopmode
file := SG-PML

all : tex2pdf backup view

cc : clean clear

tex2pdf :
	$(xlx) $(file).tex

view :
	evince $(file).pdf &

backup : $(file).tex $(file).pdf
	tar -zpcv -f Backup.tar.gz $(file).tex $(file).pdf

clean :
	-rm -f $(addprefix $(file), .aux .blg .bbl .log .synctex.gz .tex.bak \
	  .fls .xdv .fdb_latexmk)

clear :
	-rm -f $(file).pdf

# vim:ft=make:noet
