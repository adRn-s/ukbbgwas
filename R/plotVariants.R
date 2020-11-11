#' Plotting variants in skyscrapper style.
#'
#' `plotVariants()` plots `-log10(p-values)` along genomic coordinates. This
#' kind of plot is named manhattan, grandlinear or skyscrapper.
#'
#' PNG & PDF files will be saved in path specified by "ukbbgwas.path" option.
#' Each file corresponds to a pair of Phenotype and Sex (as defined when
#' creating the group).
#'
#' Heuristic for cutoff is six times the mean `-log10(p-values)`. This cutoff is
#' finally rounded so that it has no decimal positions.
#'
#' @param data A dataframe as returned by `findAtPhenotypes(mode = "pval")`.
#'
#' @param group Character string corresponding to group name where
#'     findAtPhenotypes was run to generate the data. This will be used as a
#'     subfolder inside 'ukbbgwas.path' for saving images.
#'
#' @param writeFiles One of "PNG", "PDF", "RDATA", "ALL" or "NONE" indicating if
#'     you want the files to be written to disk or not; and in which format.
#'
#' @param color A character vector with colors to pass into
#'     `ggbio::plotGrandLinear()` function.
#'
#' @param cutoff.color A color to use in cutoff line and value.
#'
#' @param cutoff A numeric vector with cutoff values for each pair of phenotypes
#'     - sexes. Pairs are ordered alphabetically. If this vector is shorter of
#'     larger than the number of pairs, an heuristic cutoff value will be used
#'     in each plot. In such case, writeLoci will be disabled.
#'
#' @param writeLoci A boolean value indicating if you want a TSV file
#'     accompanying each plot with the corresponding HGNC symbol and chromosome
#'     position of variants above cutoff value. Beware choosing low cutoff
#'     values because too many positions may slow the procedure considerably. If
#'     no group name was passed as argument and/ or writeFiles is "none"; these
#'     TSV files will be written in "ukbbgwas.path" root directory.
#'
#' @return A vector of figures you may retrieve with print(). Only when
#'     "writeFiles" option is set to "none".
#'
#' @examples
#' \dontrun{plotVariants(data, group = "example", writeFiles = "all", cutoff = TRUE)  # will trigger heuristic}
#' \dontrun{plots <- plotVariants(data, writeFiles = "none"); print(plots[[1]])}
#'
#' @export

plotVariants <- function(data, group = NA, writeFiles = "ALL", color = c("#FF7F0EFF", "#2CA02CFF"),
                         cutoff.color = "purple", cutoff = NULL, writeLoci = FALSE) {
    if (!requireNamespace("ggbio", quietly = TRUE)) {
        warning("Please run BiocManager::install(\"ggbio\") and retry.")
    } else {
        writeFiles <- toupper(writeFiles)
        valoresPosibles <- c("PNG", "PDF", "RDATA", "ALL", "NONE")
        if (all(writeFiles != valoresPosibles)) {
            stop(paste("Unrecognized writeFiles parameter value:", writeFiles))
        }
        if (is.na(group) & writeFiles != "NONE") {
            warning("Please define a group name in order to use writeFiles parameter.")
        }
        #prepareHGNC <- function(symbols, hg19) {
        #    N <- length(symbols)
        #    symbols <- strsplit(symbols, ", ")
        #    rtados <- data.frame()
        #    for (i in seq_len(N)) {
        #        rtados <- rbind(rtados, hg19[hg19$hgnc_symbol == symbols[[i]][1], ])
        #    }
        #    rtados
        #}
        data <- splitVariantColumn(data)
        recorrido <- unique(data.frame(phen = data$phenotype, sex = data$sex))
        recorrido <- with(recorrido, recorrido[order(phen, sex), ])
        tengoValoresCutoff <- !is.null(cutoff)
        smartCutoff <- FALSE
        if (tengoValoresCutoff) {
            if (length(cutoff) == nrow(recorrido)) {
                recorrido <- cbind(recorrido, cutoff)
            } else { smartCutoff <- TRUE }
        }
        N <- nrow(recorrido)
        if (writeFiles == "NONE") {
            plots <- list()
        }
        for (k in seq_len(N)) {
            titulo <- getPhenoDescrip(recorrido$phen[k], recorrido$sex[k])
            nombre.archivo <- paste0(titulo[2], "_", titulo[1])
            subsetiado <- data$phenotype == recorrido$phen[k] & data$sex == recorrido$sex[k]
            data2 <- data[subsetiado, ]
            if (tengoValoresCutoff) {
                if (smartCutoff) {
                    cutoff <- round(mean(unlist(-log10(data[subsetiado, "pval"]))) * 6, 0)
                    if (writeLoci) {
                        writeLoci <- FALSE
                        warning("HGNC symbols not written to tsv file because we're using \
                                 heuristic cutoff values.")
                    }
                } else { cutoff <- recorrido$cutoff[k] }
                if (writeLoci) {
                    if (is.na(group)) {
                        archivo.tsv <- file.path(getOption("ukbbgwas.path"),
                                                 paste0(nombre.archivo, ".tsv"))
                    } else {
                        archivo.tsv <- file.path(getOption("ukbbgwas.path"), group,
                                                 paste0(nombre.archivo, ".tsv"))
                    }
                    data3 <- addHGNC(data2[unlist(-log10(data[subsetiado, "pval"])) > cutoff, ])
                    data3 <- unique(data3[ , c("phenotype", "sex", "chromosome",
                                               "hgnc_symbol", "description")])
                    write.table(data3, sep = '\t', row.names = FALSE, file = archivo.tsv)
                    rm(data3)
                }
            }
           #  data <- data[!is.na(data$hgnc_symbol) & data$hgnc_symbol != "", ]
            gr.snp <- convertToGRange(data2)
#hgnc.df <- prepareHGNC(symbols, hg19)
#gr.hgnc <- convertToGRange(hgnc.df, c("chromosome_name", "start_position", "end_position"))
#names(gr.hgnc) <- hgnc.df$hgnc_symbol  # symbols asi vuelven todos los nombres
            plot <- ggbio::plotGrandLinear(gr.snp, ggplot2::aes(y = -log10(pval)), color = color,
                cutoff = cutoff, cutoff.size = 0.5, cutoff.color = cutoff.color) +
                ggplot2::ggtitle(sprintf("%s, %s: %s", titulo[1], titulo[2], titulo[3])) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
                ggplot2::annotate("text", label = cutoff, x = 0, y = cutoff + 1.5,
                size = 3, color = cutoff.color)
            if (writeFiles == "NONE") {
                plots[[k]] <- plot
            } else {
                if ((writeFiles == "PDF" | writeFiles == "ALL") & !is.na(group)) {
                    pdf(file.path(getOption("ukbbgwas.path"), group, paste0(nombre.archivo, ".pdf")))
                    print(plot); dev.off()
                }
                if ((writeFiles == "PNG" | writeFiles == "ALL") & !is.na(group)) {
                    png(file.path(getOption("ukbbgwas.path"), group,
                                  paste0(nombre.archivo, ".png")),
                        width = 5, height = 4, units = "in", res = 96)
                    print(plot); dev.off()
                }
                if ((writeFiles == "RDATA" | writeFiles == "ALL") & !is.na(group)) {
                    save(plot, file = file.path(getOption("ukbbgwas.path"),
                                                group, paste0(nombre.archivo, ".RData")))
                }
            }
            rm(data2, plot)
        }
    }
    if (writeFiles == "NONE") { return(plots) }
    else { return(NULL) }
}
