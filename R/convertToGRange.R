#' Converting data to Genomic Range
#'
#' `convertToGRange()` converts a dataframe with search results into a
#' genomic range with metadata columns.
#'
#' @param df A dataframe as returned by findAtPhenotypes().
#'
#' @param cnames A character vector with either two or three strings, each
#'     indicating the corresponding column name for chromosome, start and end
#'     positions. If no end position column is present, indicate only first two
#'     strings and we'll use "start + 1" for end position.
#'
#' @return A GRange object.
#'
#' @examples
#' \dontrun{convertToGRange(findAtPhenotypes(dfgroup, "pval", 1e-5))}
#'
#' @export

convertToGRange <- function(df, cnames = c("chromosome", "position")) {
    if (!requireNamespace("biovizBase", quietly = TRUE)) {
        warning("Please run BiocManager::install(\"biovizBase\") and retry.")
        return <- df
    } else if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        warning("Please run BiocManager::install(\"GenomicRanges\") and retry.")
        return <- df
    } else {
        # tidyr::separate(df, 3, c("chr", "pos", "ref", "alt"), ":"))
        if ("variant" %in% names(df)) { df <- splitVariantColumn(df) }
        if (length(cnames) == 2) {
            df <- cbind(df, final = df[ , cnames[2]] + 1)
            cnames <- c(cnames, "final")
        }
        return <- biovizBase::transformDfToGr(df, seqnames = cnames[1],
                                              start = cnames[2], end = cnames[3])
        GenomicRanges::mcols(return)[ , cnames] <- NULL
    }
    return
}
