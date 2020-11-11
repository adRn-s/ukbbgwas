#' Generate a hg19 table with biomaRt attributes
#'
#' `generateHG19()` retrieves metadata columns using biomaRt.
#'
#' User won't need this function because we already provide the returned object
#' built-in the package. It only needs to run `data(hg19)`.
#'
#' @param columntributes A character vector indicating biomaRt attributes to
#'     retrieve. See:
#'     <https://www.ensembl.org/info/data/biomart/biomart_r_package.html> for
#'     help.
#'
#' @return A dataframe with each of the 'columntributes' and corresponding
#'     genomic coordinates.

generateM37 <- function () {
    #biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    #                 host = "grch37.ensembl.org",
    #                 path = "/biomart/martservice",
    biomaRt::useEnsembl(biomart = "ensembl", GRCh = 37, mirror = "eastusa",
                      dataset = "hsapiens_gene_ensembl")
}

generateHG19 <- function(columntributes = c("hgnc_symbol", "gene_biotype", "description")) {
    hg19 <- NA
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        warning("Please run BiocManager::install(\"biomaRt\") and retry.")
    }
    else {
        atributosDeseados <- c("chromosome_name", "start_position", "end_position",
                               columntributes)
        #if (!all(sapply(c("mart37", "hg19"), `%in%`, ls()))) {
        checkOnline(); mart37 <- generateM37()
        hg19 <- biomaRt::getBM(attributes = atributosDeseados, mart = mart37)
        #}
    }
    hg19
}
