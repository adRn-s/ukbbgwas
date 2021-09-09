#' Downloading phenotype files
#'
#' defGroup() downloads specified phenotype files into a chosen directory named
#' as a group. This folder will reside inside ~/.ukbbgwas/ by default. Select
#' the download directory by issuing: options(ukbbgwas.path = "/path/new") or
#' soft-linking the new path to the default location using your shell:
#' "$ ln -s /path/new ~/.ukbbgwas".
#'
#' If the directory specified at "ukbbgwas.path" doesn't exists or is empty,
#' this function will create it for downloading first 9 core files that are
#' essential. They are listed as the first 9 rows of the manifest, see
#' data(neale.manifest)[1:9, ]. Most of these files are described as answers in
#' the FAQ: <https://www.nealelab.is/uk-biobank/faq>. Please take a look at
#' 'README' file in the `ukbbgwas.path` specified location.
#'
#' After running this function a "neale.manifest" dataframe will appear in the
#' Global Environment if not already present.
#'
#' Defining groups with too many files risks using all your available random
#' acces memory (RAM) when running searches with findAtPhenotypes() or further
#' on downstream analysis. In such cases your system will seem unresponsive or
#' either it will halt entirely. Avoid exceeding your limits by splitting groups
#' intelligently. As a rule of thumb, consider that every 2-3 files uses 1 Gb of
#' available RAM later on.
#'
#' @param group A character string naming all phenotypes. A sub-folder will be
#'     created (as needed) at the ukbbgwas.path download directory.
#'
#' @param phenotypes A character vector with phenotype codes as they appear in
#'     neale.manifest$Phenotype.Code.
#'
#' @param sexes A character vector of the same length as `phenotypes`, with a
#'     1:1 correspondence specifying one of each sex using the corresponding
#'     character string: "male", "female" or "both_sexes".
#'
#' @param index A boolean value indicating if an "index.tsv" file with
#'     phenotype's code and descriptions should be placed inside download
#'     folder.
#'
#' @return None
#'
#' @examples
#' \dontrun{defGroup("Gestational", c("O13", "O14"), rep("female", 2))}
#' \dontrun{defGroup("Posdoc",
#'                   c(c("O13", "O14", "4041", "HYPTENSPREG", "3849"),
#'                     c("IV_ENDOCRIN_NUTRIT", "23105_irnt", "E66", "E4_OBESITY",
#'                       "I9_HYPTENS", "I10")),
#'                   c(rep("female", 5),
#'                     rep("both_sexes", 6)))}
#'
#' @export

defGroup <- function(group, phenotypes, sexes, index = TRUE, clean = FALSE) {
    data(neale.manifest)
    directorio <- getOption("ukbbgwas.path")
    estaElDir <- !(dir.exists(directorio))
    estanLos9 <- length(list.files(directorio)) <= 22-1  #2nd url is 404...
    if (estaElDir || !estanLos9) {
      dir.create(directorio, recursive = TRUE)
      manifest_select <- neale.manifest[c(1,3:22), "AWS.File"]
      N <- nrow(manifest_select)
      for (i in seq_len(N)) {
          my_file <- file.path(directorio, strsplit(as.character(manifest_select[i]), "/")[[1]][6])
          checkOnline()
          utils::download.file(url = as.character(manifest_select[i]), destfile = my_file,
                               method = "wget", extra = "-c")
      }
    }
    N <- length(phenotypes)
    stopifnot(N == length(sexes))
    if (group == "example") {
        warning("The \"example\" group name is reserved for the vignette. You won't be able to list files (and therefore search through them) if such group is created.")
        group <- readline(prompt = "Please choose a new name: ")
    }
    fpath <- file.path(directorio, group)
    if (!(dir.exists(fpath))) { dir.create(fpath) }
    #else if (clean == TRUE) {
    #    #' @param clean A boolean value, if TRUE and group already exists, it
    #    #'     will clean ' all files not listed in the paired vectors
    #    #'     phenotypes/ sexes.
    #    arhivosActuales <- list.files(fpath)
    #}
    for (k in seq_len(N)) {
        # si el par fenotipo/sex no existe en phenotypes.gz, skip.
        cmd.zcat <- sprintf("zcat %s/phenotypes.%s.tsv.gz",
                            directorio, sexes[k])
        cmd.grep <- paste("| grep -w", phenotypes[k])
        dt <- data.table::fread(cmd = paste(cmd.zcat, cmd.grep))
        if (nrow(dt) == 1) {
            filtro <- c(neale.manifest$Phenotype.Code == phenotypes[k] &
                        neale.manifest$Sex == sexes[k])
            filtrado <- subset(neale.manifest, filtro)
            direccion <- sub("wget ", "", filtrado$wget.command, fixed = TRUE)
            archivo <- file.path(fpath, as.character(filtrado$File))
            checkOnline()
            utils::download.file(url = direccion, destfile = archivo,
                                 method = "wget", extra = "-c",
                                 progress = TRUE)
        } else if (nrow(dt) == 0) {
            warning(sprintf("Phenotype \"%s\" is not available for \"%s\".",
                            phenotypes[k], sexes[k]))
        }
    }
    if (index == TRUE) {
        elegidos <- neale.manifest[neale.manifest$Phenotype.Code %in% phenotypes, c(1,2)]
        write.table(unique(elegidos), file = file.path(fpath, "index.tsv"),
                    sep = "\t", row.names = FALSE)
    }
}
