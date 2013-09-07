### Makefile to prepare YourCast package

# prepare the package for release
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

# Rules
all: $(PKGSRC).pdf clean

$(PKGSRC).tex: vignettes/$(PKGSRC).Rnw
	cd vignettes;\
	echo "Sweave(\"YourCast.Rnw\", debug=TRUE, eval=TRUE)" | R --slave

### Rd2tex here
Rd2tex: man/yourcast.Rd man/yourprep.Rd man/lifetable.Rd man/plot.yourcast.Rd man/lifetable.Rd man/array.yourcast.Rd man/histograph.Rd man/summary.yourcast.Rd
	echo 'library(tools); Rd2latex(parse_Rd("man/yourcast.Rd"), out="vignettes/man-yourcast.tex"); Rd2latex(parse_Rd("man/yourprep.Rd"), out="vignettes/man-yourprep.tex"); Rd2latex(parse_Rd("man/plot.yourcast.Rd"), out="vignettes/man-plot.yourcast.tex"); Rd2latex(parse_Rd("man/lifetable.Rd"), out="vignettes/man-lifetable.tex"); Rd2latex(parse_Rd("man/array.yourcast.Rd"), out="vignettes/man-array.yourcast.tex"); Rd2latex(parse_Rd("man/histograph.Rd"), out="vignettes/man-histograph.tex"); Rd2latex(parse_Rd("man/summary.yourcast.Rd"), out="vignettes/man-summary.yourcast.tex")' | R --slave

$(PKGSRC).pdf: $(PKGSRC).tex Rd2tex
	cd vignettes;\
	pdflatex $(PKGSRC);\
	bibtex $(PKGSRC);\
	pdflatex $(PKGSRC);\
	pdflatex $(PKGSRC);\
	pdflatex $(PKGSRC)

build:	
	cd ..;\
	R CMD build $(PKGSRC) --resave-data

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran --timings                 

clean:
	cd vignettes;\
	rm -rf *.aux *.tex *.bbl *.blg *.bcf *.log *.out *.rel *.toc *.idx