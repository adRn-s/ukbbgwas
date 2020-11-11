#' Splitting variant column in data
#'
#' `splitVariantColumn()` splits variant column into "chromosome", "position",
#' "ref" and "alt" columns.
#'
#' This functions is auxilliary for convertToGRange().
#'
#' @param df A dataframe as returned by findAtPhenotypes()
#'
#' @return A dataframe with variant column removed and four new columns added.

splitVariantColumn <- function(df) {
    newCols <- c("chromosome", "position", "ref", "alt")
    if (all(newCols %in% names(df))) {
        warning(sprintf("All columns (%s) already present.",
                        paste0(newCols, collapse = ", ")))
    } else {
        stopifnot(names(df)[3] == "variant")
        N <- nrow(df)
        matriz <- matrix(unlist(strsplit(as.character(df[ , 3]), ":")),
                         byrow = TRUE, nrow = N)
        df[ , 3] <- NULL
        nombres <- c(newCols, names(df))
        df <- cbind(matriz, df)
        rm(matriz)
        names(df) <- nombres
        df$position <- as.numeric(as.character(df$position))
    }
    df
}
