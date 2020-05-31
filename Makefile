check:
	Rscript -e 'devtools::load_all()'
	Rscript -e 'devtools::check()'

build:
	Rscript -e 'devtools::build()'

test:
	Rscript -e 'devtools::test()'

install:
	Rscript -e 'devtools::install()'

manual: docs/sveval-manual.pdf

docs/sveval-manual.pdf: man/*.Rd
	rm -f docs/sveval-manual.pdf
	R CMD Rd2pdf . --no-preview -o docs/sveval-manual.pdf
