ukbbgwas_%.tar.gz:
	Rscript -e "devtools::document(path = \".\")"
	Rscript -e "devtools::build()"
