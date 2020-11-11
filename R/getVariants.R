#' Searching through variants files.
#'
#' `getVariants()` looks for matching variants.
#'
#' @param vectorstr A character vector corresponding to SNP rs identification
#'     number(s) to retrieve.
#'
#' @return Data frame with variant information. All columns are inherited from
#'     Neale's Lab published file (compressed TSV format).
#'
#' @examples
#' \dontrun{getVariants(c("rs236918", "rs2369"))}
#'
#' @export

getVariants <- function(vectorstr) {
    cmd.zcat <- paste("zcat", file.path(getOption("ukbbgwas.path"),
                                        "variants.tsv.bgz"))
    grep.args <- paste(vectorstr, collapse = " -e ")
    cmd.grep <- paste(cmd.zcat, "| grep -w -e ", grep.args)
    variante <- data.table::fread(cmd = cmd.grep)
    cmd.head <- paste(cmd.zcat, "| head")
    names(variante) <- names(data.table::fread(cmd = cmd.head, header = TRUE))
    if (nrow(variante) < length(vectorstr)) {
        warning("Not all rs IDs were found.")
    }
    variante
}
