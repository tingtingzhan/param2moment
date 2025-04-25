

#' @title Moments of Normal Distribution
#' 
#' @description
#' Moments of \href{https://en.wikipedia.org/wiki/Normal_distribution}{normal distribution}.
#' 
#' @param mean,sd \link[base]{numeric} scalars or \link[base]{vector}s, 
#' see function \link[stats]{dnorm}
#' 
#' @returns
#' Function [moment_norm()] returns a \linkS4class{moment} object.
#' 
#' @keywords internal
#' @export
moment_norm <- function(mean = 0, sd = 1) {
  new(Class = 'moment', 
      distname = 'norm', 
      location = mean, scale = sd, 
      raw1 = 0, raw2 = 1, raw3 = 0, raw4 = 3)
}

