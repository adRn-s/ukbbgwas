# UK BioBank GWAS Results

Explore SNVs associations with traits of UK BioBank. GWAS study results where
made publicly available courtesy of Neale\'s Lab on August 1st, 2018. For more
information please visit: <http://www.nealelab.is/uk-biobank>

## INSTALL

In R: `devtools::install_github("adRn-s/ukbbgwas")`.

## DEV(s)

### Dependencies

Mandatory and optional dependencies can all be handled with
[bioconda](https://bioconda.github.io/). Here's a complete `YAML`
specification if you choose to do so:

```
dependencies:
 - r-essentials
 - r-devtools
 - r-curl
 - r-data.table
 - r-tidyr
 - bioconductor-biomart
 - bioconductor-genomicranges
 - bioconductor-biovizbase
 - bioconductor-ggbio
 - bioconductor-gviz
 - bioconductor-sushi
```

## USAGE

See vignette: `browseVignette(package = "ukbbgwas")`.



<!--
## TODO

- add `plot*()` functions to vignette.
- built in tutorial following vignette  : https://education.rstudio.com/blog/2020/09/delivering-learnr-tutorials-in-a-package/
- write tests (`testthat`) and/ or runnable examples (`biocCheck`)
- add proper `CITATION` (incl. Neale's URL)
- add 2019 new results, i.e.
  [biomarkers](https://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results).
  See [changelog](https://github.com/Nealelab/UK_Biobank_GWAS) and new file
  manifest.
- conectar con prsice, ldpred2, bayes jerarquico, etc
- ver analisis de sumHer y agregarles
- retomar MR
- tests de burden
- agregar otros GWAS y graficos de meta-analisis como forestplot : https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
-->
