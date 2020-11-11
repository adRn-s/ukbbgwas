.onLoad <- function(libname, pkgname) {
    op <- options()
    op.ukbbgwas <- list(
        ukbbgwas.path = "~/.ukbbgwas"
    )
    toset <- !(names(op.ukbbgwas) %in% names(op))
    if(any(toset)) options(op.ukbbgwas[toset])
    invisible()
}
