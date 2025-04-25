

#' @title Moments of Skew-Normal Distribution
#' 
#' @param xi,omega,alpha \link[base]{numeric} scalars or \link[base]{vector}s, 
#' location \eqn{\xi}, scale \eqn{\omega} and slant \eqn{\alpha},
#' see function \link[sn]{dsn}.
#' 
#' @returns
#' Function [moment_sn()] returns a \linkS4class{moment} object.
#' 
#' @keywords internal
#' @export
moment_sn <- function(xi = 0, omega = 1, alpha = 0) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(2/pi)
  r1 <- b * delta
  new(Class = 'moment', 
      distname = 'sn', 
      location = xi, scale = omega,
      raw1 = r1,
      raw2 = 1,
      raw3 = - pi/2 *r1^3 + 3*r1,
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
#' @keywords internal
#' @importFrom stats optim
#' @export
moment2sn <- function(mean = 0, sd = 1, skewness) {
  opt <- optim(par = c(xi = 0, omega = 1, alpha = 0), fn = \(x) {
    alpha <- x[3L]
    delta <- alpha / sqrt(1 + alpha^2)
    b <- sqrt(2/pi)
    r1 <- b * delta
    m <- moment_init(raw1 = r1, raw2 = 1, raw3 = - pi/2 *r1^3 + 3*r1, raw4 = 3)
    
    z1 <- mean_moment_(m, location = x[1L], scale = x[2L])
    z2 <- sd_moment_(m, scale = x[2L])
    z3 <- skewness_moment_(m)
    crossprod(c(z1, z2, z3) - c(mean, sd, skewness))
  })
  return(opt$par)
}






