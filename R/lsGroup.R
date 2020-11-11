#' Listing downloaded files.
#'
#' `lsGroup()` lists downloaded files in a corresponding group.
#'
#' @param groupName one of the previously defined groups with `defGroup()`.
#'
#' @return A dataframe with each of the phenotype-sex pairs available in the
#'     group.
#'
#' @examples
#' \dontrun{lsGroup("MyResearch")}
#' lsGroup("example")
#'
#' @export

lsGroup <- function(groupName) {
    if (groupName == "example") {
        directorio <- system.file("extdata", "phenotypes", package = "ukbbgwas")
    } else {
        directorio <- file.path(getOption("ukbbgwas.path"), groupName)
    }
    my_files <- list.files(directorio)
    patron <- ".*\\.(female|both_sexes|male)\\.tsv\\.bgz"
    df.group.raw <- strsplit(my_files[grepl(patron, my_files)], "\\.")
    N <- length(df.group.raw)
    df <- data.frame(group = rep(groupName, N),
                     phenotype = rep(NA, N),
                     sex = rep(NA, N))
    for (i in seq_len(N)) {
        df$phenotype[i] <- df.group.raw[[i]][1]
        sex.index <- ifelse(length(df.group.raw[[i]]) == 4, 2, 4)
        df$sex[i] <- df.group.raw[[i]][sex.index]
    }
    df
}
