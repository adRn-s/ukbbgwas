#' Retrieving Phenotype Description
#'
#' `getPhenoDescrip()` retrieves phenotype description for a paired phenotype
#' and sex.
#'
#' @param pheno Phenotype code factor.
#'
#' @param sex Sex factor.
#'
#' @return A character vector of length 3 with sex, phenotype.code and phenotype
#'     description; ready to use as tittle in a manhattan plot.

getPhenoDescrip <- function(pheno, sex) {
    sex <- as.character(sex)
    pheno <- as.character(pheno)
    cmd.zcat <- sprintf("zcat %s/phenotypes.%s.tsv.gz",
                        getOption("ukbbgwas.path"), sex)
    cmd.grep <- paste("| grep -w", pheno)
    dt <- data.table::fread(cmd = paste(cmd.zcat, cmd.grep), sep = '\t')
    if (nrow(dt) == 0) {
        descrip <- NA
        warning(sprintf("Phenotype \"%s\" is not available for \"%s\".", pheno, sex))
    } else { descrip <- as.character(dt[, 2]) }
    c(sex, pheno, descrip)
}
