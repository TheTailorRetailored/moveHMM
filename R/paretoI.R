#' Pareto Type I distribution helpers
#'
#' Internal density, distribution and random-generation functions used by the
#' moveHMM step-length machinery. The distribution has minimum value `xmin`
#' and shape parameter `mu`.

dparetoI <- function(x, xmin, mu, log = FALSE) {
    if (length(xmin) != 1L || !is.finite(xmin) || xmin <= 0)
        stop("xmin must be a finite positive number")
    if (length(mu) != 1L || !is.finite(mu) || mu <= 0)
        stop("mu must be a finite positive number")

    log_density <- ifelse(
        x >= xmin,
        log(mu) - log(xmin) - (mu + 1) * log(x / xmin),
        -Inf
    )
    if (log) log_density else exp(log_density)
}

pparetoI <- function(q, xmin, mu, lower.tail = TRUE, log.p = FALSE) {
    if (length(xmin) != 1L || !is.finite(xmin) || xmin <= 0)
        stop("xmin must be a finite positive number")
    if (length(mu) != 1L || !is.finite(mu) || mu <= 0)
        stop("mu must be a finite positive number")

    probability <- ifelse(q < xmin, 0, 1 - (xmin / q)^mu)
    if (!lower.tail) probability <- 1 - probability
    if (log.p) log(probability) else probability
}

rparetoI <- function(n, xmin, mu) {
    if (length(n) > 1L) n <- length(n)
    if (length(xmin) != 1L || !is.finite(xmin) || xmin <= 0)
        stop("xmin must be a finite positive number")
    if (length(mu) != 1L || !is.finite(mu) || mu <= 0)
        stop("mu must be a finite positive number")

    xmin * (1 - stats::runif(n))^(-1 / mu)
}
