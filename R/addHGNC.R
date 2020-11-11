#' Adding hg19's attributes as new columns to dataframe
#'
#' `addHGNC()` takes a dataframe and adds matching HGNC symbols, gene biotype
#' and gene description columns to it.
#'
#' @param df A dataframe with genomic coordinates as columns.
#'
#' @param coords A character vector indicating the column names of corresponding
#'     and required genomic coordinates.
#'
#' @return A dataframe with new columns.

addHGNC <- function(df, coords = c("chromosome", "position")) {
    data(hg19)  # from generateHG19 (biomaRt)
    newColumns <- c("hgnc_symbol", "gene_biotype", "description")
    if (!all(newColumns %in% names(df))) {
        nombres <- names(df)
        N <- nrow(df)
        for (j in seq_len(length(newColumns))) {
            df <- cbind(Vn = rep(NA, N), df)
        }
        names(df) <- c(newColumns, nombres)
    }
    for (i in seq_len(N)) {
        target_chr <- df[i, coords[1]]
        target_pos <- df[i, coords[2]]
        rtado <- hg19[with(hg19, chromosome_name == target_chr &
                                 start_position <= target_pos &
                                 end_position >= target_pos), ]
        if (nrow(rtado) > 0) {
            if ("hgnc_symbol" %in% names(rtado)) {
                # skip pseudogenes and other(s):
                rtado <- rtado[rtado$hgnc_symbol != "", ]
            }
            for (k in newColumns) {
                df[i, k] <- paste(rtado[ , k], collapse = ", ")
            }
        }
    }
    df
}
