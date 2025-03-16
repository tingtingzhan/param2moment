

#' @title Moments of Normal Distribution
#' 
#' @description
#' Moments of \href{https://en.wikipedia.org/wiki/Normal_distribution}{normal distribution}, parameter nomenclature follows
#' \link[stats]{dnorm} function.
#' 
#' @param mean \link[base]{numeric} scalar or \link[base]{vector}, 
#' mean parameter \eqn{\mu}
#' 
#' @param sd \link[base]{numeric} scalar or \link[base]{vector}, 
#' standard deviation \eqn{\sigma}
#' 
#' @returns
#' Function [moment_norm()] returns a \linkS4class{moment} object.
#' 
#' @examples
#' moment_norm(mean = 1.2, sd = .7)
#' 
#' @export
moment_norm <- function(mean = 0, sd = 1) {
  do.call(what = new, args = c(list(Class = 'moment', location = mean, scale = sd), moment_int(distname = 'norm', mu = 0, raw2 = 1, raw3 = 0, raw4 = 3)))
}

