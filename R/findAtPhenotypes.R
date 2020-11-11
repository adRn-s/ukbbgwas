#' Searching in all phenotypes from a given group.
#'
#' `findAtPhenotypes()` is this package main function. It searches through all
#' phenotypes in a previously defined group. There are three search modes, either
#' by a threshold p-value, by a given vector of coordinates in the format
#' `chr:position` or by a region from any chromosome.
#'
#' @param df.group A dataframe with all the phenotypes in a given group, just
#'     like the one returned from `lsGroup()`.
#'
#' @param search.mode A string keyword enabling one of the available modes,
#'     either "pval", "coords" or "region".
#'
#' @param search.args A number to use as maximun in case of "pval" mode. Or a
#'     vector with coordinates in the format `chr:position:ref:alt` (i.e.
#'     `c("X:349813:T:C", "3:843028:A:G")`) in case of "coords" mode. Last, but
#'     not least, the "region" mode would take coordinates too but only for a
#'     single genomic interval using a vector of three integers, namely:
#'     chromosome, start and end position.
#'
#' @param p.threshold Optionally, filter out matching hits with p-values above a
#'     given number (i.e. > 5e-10). This argument only makes sense when mode is
#'     not "pval".
#'
#' @param write Optionally, set to TRUE and write to a file in "ukbbgwas.path"
#'     (package option) instead of returning a dataframe (that must fit in RAM).
#'     The filename will start with a timestamp (year, month, day, hour, minute)
#'     followed by "_search_", the corresponding phenotype group, search.mode in
#'     use and ending with the ".tsv" extension.
#'
#' @param lowconfilter Optionally, set to TRUE to filter out variants that
#'     are low confidence under the heuristic(s) applied by Neale's lab.
#'     (Case/control phenotypes: expected_case_minor_AC < 25 or minor_AF <
#'     0.001. Categorical phenotypes with less than 5 categories:
#'     expected_min_category_minor_AC < 25 or minor_AF < 0.001 Quantitative
#'     phenotypes: minor_AF < 0.001.)
#'
#' @return A dataframe containing variants that match the search criteria.
#'
#' @examples
#' \dontrun{findAtPhenotypes(lsGroup("MyResearch"), "pval", 5e-20)}
#' findAtPhenotypes(lsGroup("example"), "pval", 5e-3)
#' findAtPhenotypes(lsGroup("example"), "coords", "1:69487:G:A")
#' \dontrun{findAtPhenotypes(lsGroup("MyResearch"), "coords", with(getVariants("rs236918"), variant))}
#' findAtPhenotypes(lsGroup("example"), "region", c(1, 69200, 69700))
#' findAtPhenotypes(lsGroup("example"), "region", c(1, 69200, 69700), 5e-10)
#'
#' @export

findAtPhenotypes <- function(df.group, search.mode, search.args,
                             p.threshold = NA, write = FALSE, lowconfilter = FALSE)
  {
  phenotypeName <- as.character(df.group$phenotype[1])
  pairedSex <- as.character(df.group$sex[1])
  groupName <- as.character(df.group$group[1])
  cmd.zcat <- paste("zcat", getPhenotypePath(
                                phenotypeName, pairedSex, groupName))
  nombres <- c("variant", "minor_allele", "minor_AF",
               "low_confidence_variant", "n_complete_samples", "AC", "ytx",
               "beta", "se", "tstat", "pval")
  if (write == FALSE) {
    orig <- as.data.frame(matrix(ncol = 11))
    names(orig) <- nombres
    orig <- cbind(phenotype = NA, sex = NA, orig)
  } else {
    file.name <- paste0(format(Sys.time(), "%Y%m%d%H%M"),
                        "_search_", groupName, "_", search.mode, ".tsv")
    file.handler <- file.path(getOption("ukbbgwas.path"), file.name)
    write(paste0(c("phenotype", "sex", nombres), collapse = "\t"), file.handler)
  }
  # PVAL
  if (search.mode == "pval") {
    cmd.awk <- paste0("awk -F '\t' -v pv=", search.args,
                      " '{ if($NF <= pv) print $0 }'")
  # COORDS
  } else if (search.mode == "coords") {
    search.args <- paste0("(", paste(search.args, collapse = "|"), ")")
    if (is.na(p.threshold)) {
      cmd.awk <- paste0("awk -F '\t' -v id=\"", search.args,
                        "\" '{ if($1 ~ id) print $0 }'")
    } else {
      cmd.awk <- paste0("awk -F '\t' -v pt=", p.threshold,
                        " -v id=\"", search.args,
                        "\" '{ if($1 ~ id && $NF <= pt) print $0 }'")
    }
  # REGION
  } else if (search.mode == "region") {
      # TODO integrar con bioconductor (GenomicRanges)
      if (is.na(p.threshold)) {
          cmd.awk <- paste0("awk -F '[:\t]' -v chr=", search.args[1],
              " -v start=", search.args[2], " -v end=", search.args[3],
              " '{ if($1 == chr && $2 >= start && $2 <= end) print $0 }'")
      } else {
          cmd.awk <- paste0("awk -F '[:\t]' -v chr=", search.args[1], " -v start=",
      search.args[2], " -v end=", search.args[3], " -v pt=", p.threshold,
      " '{ if($1 == chr && $2 >= start && $2 <= end && $NF <= pt) print $0 }'")
      }
  }
  # SEARCHING
  M <- nrow(df.group)
  for (i in 1:M) {
    phenotypeName <- as.character(df.group$phenotype[i])
    pairedSex <- as.character(df.group$sex[i])
    groupName <- as.character(df.group$group[i])
    cmd.zcat <- paste("zcat", getPhenotypePath(phenotypeName, pairedSex,
                                               groupName))
    solo.nombres <- data.table::fread(cmd = paste(cmd.zcat, "| head -1"))
    message(paste(paste0("(", i, "/", M, ")"), "Looking through:",
                  groupName, phenotypeName, pairedSex, "..."))
    stopifnot(length(phenotypeName) == 1, length(pairedSex) == 1)
    rtado.parcial <- data.table::fread(cmd = paste(cmd.zcat, "|", cmd.awk))
    if (lowconfilter == TRUE) {
        rtado.parcial <- rtado.parcial[rtado.parcial$low_confidence_variant == "FALSE", ]
    }
    Q <- nrow(rtado.parcial)
    if (Q > 0) {
      message(paste("\t Found", Q, "variants."))
      names(rtado.parcial) <- names(solo.nombres)
      columnas <- c("expected_min_category_minor_AC", "expected_case_minor_AC")
      for (columna in columnas) {
          if (columna %in% names(rtado.parcial)) {
              rtado.parcial[[columna]] <- NULL
          }
      }
      rtado.parcial <- cbind(phenotype = phenotypeName, sex = pairedSex, rtado.parcial)
      if (write == FALSE) {
        orig <- rbind(orig, rtado.parcial)
      } else {
          message("\t Writing to file, please wait.")
          for (j in 1:Q) {
            write(paste0(as.character(rtado.parcial[j, ]), collapse = "\t"),
                  file.handler, append = TRUE)
          }
      }
    }
  }
  # RETURN
  if (write == FALSE) {
      M <- nrow(orig)
      if (M < 2) {
          message("None found.")
          return(NA)
      } else {
          return(orig[2:M, ])
      }
  } else {
      message(paste("Results written to file:", as.character(file.handler)))
  }
}
