#' Retrieving phenotype file path.
#'
#' `getPhenotypePath()` retrieves the phenotype file path with its full name and
#' the previously defined group sub-folder.
#'
#' @param  phenotypeName A string with the phenotype name.
#'
#' @param pairedSex A string with the corresponding sex to the specified
#'     phenotype.
#'
#' @param groupName A string indicating the previously defined group where
#'     phenotype was chosen and downloaded.
#'
#' @return A string with the absolute path to desired phenotype file.
#'
#' @examples
#' \dontrun{getPhenotypePath("Z85", "both_sexes", "MyResearch")}

getPhenotypePath <- function(phenotypeName, pairedSex, groupName) {
    directorio <- file.path(getOption("ukbbgwas.path"), groupName)
    my_file.fenotipo <- list.files(directorio)
    # my_file.fenotipo <- my_file.fenotipo[grepl("\\.bgz", my_file.fenotipo)]
    my_file.fenotipo <- my_file.fenotipo[grepl(phenotypeName, my_file.fenotipo)]
    my_file.fenotipo <- my_file.fenotipo[grepl(paste0("\\.", pairedSex),
                                               my_file.fenotipo)]
    file.path(directorio, my_file.fenotipo)
}
