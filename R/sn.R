

#' @title Moments of Skew-Normal Distribution
#' 
#' @description
#' Moments of \href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{skew-normal distribution}, parameter nomenclature follows
#' \link[sn]{dsn} function.
#' 
#' @param xi \link[base]{numeric} scalar or \link[base]{vector}, 
#' location parameter \eqn{\xi}
#' 
#' @param omega \link[base]{numeric} scalar or \link[base]{vector}, 
#' scale parameter \eqn{\omega}
#' 
#' @param alpha \link[base]{numeric} scalar or \link[base]{vector}, 
#' slant parameter \eqn{\alpha}
#' 
#' @returns
#' Function [moment_sn()] returns a \linkS4class{moment} object.
#' 
#' @examples
#' xi = 2; omega = 1.3; alpha = 3
#' moment_sn(xi, omega, alpha)
#' curve(sn::dsn(x, xi = 2, omega = 1.3, alpha = 3), from = 0, to = 6)
#' 
#' @importFrom methods new
#' @export
moment_sn <- function(xi = 0, omega = 1, alpha = 0) {
  do.call(what = new, args = c(list(Class = 'moment', location = xi, scale = omega), moment_sn_(alpha = alpha)))
}

moment_sn_ <- function(alpha = 0) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(2/pi)
  mu <- b * delta
  moment_int(
    distname = 'sn', 
    mu = mu,
    raw2 = 1,
    raw3 = - pi/2 *mu^3 + 3*mu,
    raw4 = 3)
}



#' @title Solve Skew-Normal Parameters from Moments
#' 
#' @description
#' Solve skew-normal parameters from mean, standard deviation and skewness.
#' 
#' @param mean \link[base]{numeric} scalar, mean \eqn{\mu}, default value 0
#' 
#' @param sd \link[base]{numeric} scalar, standard deviation \eqn{\sigma}, default value 1
#' 
#' @param skewness \link[base]{numeric} scalar
#' 
#' @details
#' Function [moment2sn()] solves the 
#' location \eqn{\xi}, scale \eqn{\omega} and slant \eqn{\alpha} parameters 
#' of skew-normal distribution,
#' from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1) and 
#' skewness.  
#' 
#' @returns
#' Function [moment2sn()] returns a \link[base]{length}-3 
#' \link[base]{numeric} \link[base]{vector} \eqn{(\xi, \omega, \alpha)}.
#' 
#' @examples
#' moment2sn(skewness = .3)
#' 
#' @importFrom stats optim
#' @export
moment2sn <- function(mean = 0, sd = 1, skewness) {
  optim(par = c(xi = 0, omega = 1, alpha = 0), fn = function(x) {
    mm <- moment_sn_(alpha = x[3L])
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm)) - c(mean, sd, skewness))
  })$par
}


