#' Plotting variants in tracklayer style.
#'
#' `plotRegion()` plots `-log10(p-values)` along genomic coordinates, as a track
#' of a genome browser. Ideogram with chromosome position and track with
#' transcript RNA are added.
#'
#' PNG & PDF files will be saved in path specified by "ukbbgwas.path" option.
#'
#' @param data A dataframe as returned by `findAtPhenotypes(mode = "pval")`.
#'
#' @param start An integer value indicating where the loci and plot should
#'     start.
#'
#' @param end An integer value indicating where the loci and plot should end.
#'
#' @param writeFiles One of "PNG", "PDF", "RDATA", "ALL" or "NONE" indicating if
#'     you want the files to be written to disk or not; and in which format.
#'
#' @param fname A character string to name file(s) to write.
#'
#' @param group Character string corresponding to group name where
#'     findAtPhenotypes was run to generate the data. This will be used as a
#'     subfolder inside 'ukbbgwas.path' for saving images.
#'
#' @return A Gviz plot object, if saved to a variable; retrieving with print()
#'     is unsupported by upstream. Only when "writeFiles" option is set to
#'     "none".
#'
#' @examples
#' \dontrun{plotRegion(myFoundByPval, gene_start - 1000, gene_end + 1000, fname = "myLoci")}
#' \dontrun{p <- plotRegion(myFoundByPval, gene_start, gene_end, "none")}
#'
#' @export

# el cromosoma está dado por el retorno de findAtPhenotypes() ...¿ok?
plotRegion <- function(data, start, end,
                       writeFiles = "ALL", fname = NA, group = NA) {
    dependencies <- c("Gviz", "tidyr", "biovizBase",
                      "GenomicRanges", "GenomeInfoDb", "rtracklayer")
    checkDeps <- sapply(dependencies, requireNamespace, quietly = TRUE)
    if (writeFiles != "NONE" & is.na(fname)) {
        stop("Please provide 'fname' argument when using 'writeFiles'.")
    }
    if (!all(checkDeps)) {
        warning("Please run BiocManager::install() through \
                dependencies: ", paste0(dependencies[checkDeps],
                                        collapse = ", "))
    } else {
        data <- data[ , c("phenotype", "sex", "variant", "pval")]
        data$pval <- -log10(data$pval)
        data <- tidyr::pivot_wider(data, names_from = c("phenotype", "sex"),
                                   values_from = "pval")
        columnas <- c("chr", "pos", "ref", "alt")
        data <- tidyr::separate(data, 1, columnas, ":")
        data$pos <- as.numeric(data$pos)
        data <- cbind(data, final = data$pos + 1)
        gr <- biovizBase::transformDfToGr(data, seqnames = "chr",
                                          start = "pos", end = "final")
        GenomicRanges::mcols(gr)[ , c(columnas, "final")] <- NULL
        GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
        GenomeInfoDb::genome(gr) <- "hg19"
        #options(ucscChromosomeNames=FALSE)
        chr <- as.character(unique(GenomicRanges::seqnames(gr)))
        session <- rtracklayer::browserSession()
        rtracklayer::genome(session) <- "hg19"
        ucscTrack <- Gviz::UcscTrack(genome = "hg19", track = "ensGene",
                               chromosome = chr, strand = "strand",
                               from = start, to = end,
                               trackType = "GeneRegionTrack",
                               rstarts = "exonStarts", rends = "exonEnds",
                               gene = "name", symbol = "name2",
                               transcript = "name",
                               name = "Ensembl Gene Model")
        iTrack <- Gviz::IdeogramTrack(genome = "hg19", chromosome = chr)
        gTrack <- Gviz::GenomeAxisTrack(littleTicks = TRUE)
        aTrack <- Gviz::AnnotationTrack(gr, name = "Variants")
        dTrack <- Gviz::DataTrack(gr, name = "-log10(P-Values)",
                                  groups = names(GenomicRanges::mcols(gr)),
                                  legend = TRUE, type = c("p", "g"))
        if (writeFiles != "NONE") {
            rutArchivo <- function(fname, group, ext) {
                if (is.na(group)) {
                    ruta <- file.path(getOption("ukbbgwas.path"), paste0(fname, ext))
                } else {
                    ruta <- file.path(getOption("ukbbgwas.path"), group, paste0(fname, ext))
                }
                ruta
            }
            if ((writeFiles == "PDF" | writeFiles == "ALL") & !is.na(fname)) {
                pdf(rutArchivo(fname, group, ".pdf"))
                Gviz::plotTracks(list(iTrack, gTrack, ucscTrack, aTrack, dTrack))
                dev.off()
            }
            if ((writeFiles == "PNG" | writeFiles == "ALL") & !is.na(fname)) {
                png(rutArchivo(fname, group, ".png"),
                    width = 5, height = 4, units = "in", res = 96)
                Gviz::plotTracks(list(iTrack, gTrack, ucscTrack, aTrack, dTrack))
                dev.off()
            }
            if ((writeFiles == "RDATA" | writeFiles == "ALL") & !is.na(fname)) {
                save(Gviz::plotTracks(list(
                               iTrack, gTrack, ucscTrack, aTrack, dTrack)),
                     file = rutArchivo(fname, group, ".RData"))
            }
        }
    }
    if (writeFiles == "NONE") {
        Gviz::plotTracks(list(iTrack, gTrack, ucscTrack, aTrack, dTrack))
    } else { return(NULL) }
}
