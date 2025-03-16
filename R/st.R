

#' @title Moments of Skew-\eqn{t} Distribution
#' 
#' @description
#' Moments of skew-\eqn{t} distribution, parameter nomenclature follows
#' \link[sn]{dst} function.
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
#' @param nu \link[base]{numeric} scalar or \link[base]{vector}, 
#' degree of freedom \eqn{\nu}
#' 
#' @returns
#' Function [moment_st()] returns a \linkS4class{moment} object.
#' 
#' @references
#' Raw moments of skew-\eqn{t}: \url{https://arxiv.org/abs/0911.2342}
#' 
#' @examples
#' xi = 2; omega = 1.3; alpha = 3; nu = 6
#' curve(sn::dst(x, xi = xi, omega = omega, alpha = alpha, nu = nu), from = 0, to = 6)
#' moment_st(xi, omega, alpha, nu)
#' 
#' @export
moment_st <- function(xi = 0, omega = 1, alpha = 0, nu = Inf) {
  do.call(what = new, args = c(list(Class = 'moment', location = xi, scale = omega), moment_st_(alpha = alpha, nu = nu)))
}
moment_st_ <- function(alpha = 0, nu = Inf) { # xi = 0, omega = 1, 
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(nu/pi) * gamma(nu/2 - 1/2) / gamma(nu/2) # equation (29); https://arxiv.org/pdf/0911.2342.pdf
  mu <- b * delta
  moment_int(
    distname = 'st', 
    mu = mu,
    raw2 = nu/(nu-2),
    raw3 = mu * (3-delta^2) * nu/(nu-3),
    raw4 = 3*nu^2/(nu-2)/(nu-4)
  )
}



# @title Solve \eqn{t}-Distribution from Moments
# 
# @description
# A short description...
# 



#' @title Solve Skew-\eqn{t} Parameters from Moments
#' 
#' @description
#' Solve skew-\eqn{t} parameters from mean, standard deviation, skewness and kurtosis.
#' 
#' @param mean \link[base]{numeric} scalar, mean \eqn{\mu}, default value 0
#' 
#' @param sd \link[base]{numeric} scalar, standard deviation \eqn{\sigma}, default value 1
#' 
#' @param skewness \link[base]{numeric} scalar
#' 
#' @param kurtosis \link[base]{numeric} scalar
#' 
#' @details
#' Function [moment2st()] solves the 
#' location \eqn{\xi}, scale \eqn{\omega}, slant \eqn{\alpha} 
#' and degree of freedom \eqn{\nu} parameters of skew-\eqn{t} distribution,
#' from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1), 
#' skewness and kurtosis.  
#' 
#' @returns
#' Function [moment2st()] returns a \link[base]{length}-4 \link[base]{numeric} \link[base]{vector} 
#' \eqn{(\xi, \omega, \alpha, \nu)}.
#' 
#' @examples
#' moment2st(skewness = .2, kurtosis = .3)
#' 
#' @name moment2st
#' @importFrom stats optim
#' @export
moment2st <- function(mean = 0, sd = 1, skewness, kurtosis) {
  optim(par = c(
    xi = 0, omega = 1, alpha = 0, 
    nu = 10 # do not allow Inf starting value
  ), fn = function(x) {
    mm <- moment_st_(alpha = x[3L], nu = x[4L])
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm), kurtosis_moment_(mm)) - c(mean, sd, skewness, kurtosis))
  })$par
}


#' @rdname moment2st
#' 
#' @details
#' An educational and demonstration function [moment2t_demo] solves 
#' \eqn{(\omega, \nu)} parameters of \eqn{t}-distribution,
#' from user-specified \eqn{\sigma} and kurtosis.
#' This is a non-skewed distribution, thus 
#' the location parameter \eqn{\xi=\mu=0}, and the slant parameter \eqn{\alpha=0}.
#' 
#' @returns
#' Function [moment2t_demo()] returns a \link[base]{length}-2 
#' \link[base]{numeric} \link[base]{vector} \eqn{(\omega, \nu)}.
#' 
#' @examples
#' moment2t_demo(kurtosis = .3)
#' 
#' @importFrom stats optim
#' @export
moment2t_demo <- function(sd = 1, kurtosis) {
  optim(par = c(
    omega = 1, nu = 10 # do not allow Inf starting value
  ), fn = function(x) {
    mm <- moment_st_(alpha = 0, nu = x[2L])
    crossprod(c(sd_moment_(mm, scale = x[1L]), kurtosis_moment_(mm)) - c(sd, kurtosis))
  })$par
}
