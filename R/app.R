#' Launch the interactive app
#'
#' Opens a Shiny app for analytical power, sample size, power curves,
#' simulation checks and concentration profiles. Needs the 'shiny' and 'bslib'
#' packages.
#'
#' @param ... Passed to [shiny::runApp()], e.g. `port` or `launch.browser`.
#' @return Called for its side effect.
#' @examples
#' if (interactive()) cqtc_app()
#' @export
cqtc_app <- function(...) {
  need <- c("shiny", "bslib")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(sprintf("The app needs: %s. Install with install.packages(c(%s)).",
                 paste(miss, collapse = ", "),
                 paste0('"', miss, '"', collapse = ", ")), call. = FALSE)
  }
  shiny::runApp(system.file("shiny", package = "cqtcpower"), ...)
}
