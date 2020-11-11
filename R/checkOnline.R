#' Checking Internet Connectivity.
#'
#' `checkOnline()` checks if you're online before downloading any data files.
#'
#' @return NULL

checkOnline <- function() {
    if (requireNamespace("curl", quietly = TRUE)) {
        offline <- is.null(curl::nslookup("r-project.org", error = FALSE))
        if (offline) {
            warning("Couldn't resolve NS lookup.
  Please ensure you have a realiable Internet connection.
  If access to the Internet is only available to you from
  behind a proxy, please configure it properly in your OS
  through environmental variables. See \"Setting Proxies\"
  at help(\"download.file\") page.")
        }
    }
    else {
        warning("Please install 'curl' package to check Internet
  availability and smooth download experience.")
    }
}
