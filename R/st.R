

#' @title Moments of Skew-\eqn{t} Distribution
#' 
#' @param xi,omega,alpha,nu \link[base]{numeric} scalars or \link[base]{vector}s, 
#' location \eqn{\xi}, scale \eqn{\omega}, slant \eqn{\alpha} and degree of freedom \eqn{\nu}
#' of function \link[sn]{dst}.
#' 
#' @returns
#' Function [moment_st()] returns a \linkS4class{moment} object.
#' 
#' @references
#' Raw moments of skew-\eqn{t}: \url{https://arxiv.org/abs/0911.2342}
#' 
#' @keywords internal
#' @export
moment_st <- function(xi = 0, omega = 1, alpha = 0, nu = Inf) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(nu/pi) * gamma(nu/2 - 1/2) / gamma(nu/2) # equation (29); https://arxiv.org/pdf/0911.2342.pdf
  r1 <- b * delta
  new(Class = 'moment', 
      distname = 'st', 
      location = xi, scale = omega,
      raw1 = r1,
      raw2 = nu/(nu-2),
      raw3 = r1 * (3-delta^2) * nu/(nu-3),
      raw4 = 3*nu^2/(nu-2)/(nu-4))
}








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
#' @param alpha \link[base]{numeric} scalar, 
#' constrained slant parameter \eqn{\alpha}, currently only `alpha=0` allowed.
#' 
#' @details
#' Function [moment2st()] solves the 
#' location \eqn{\xi}, scale \eqn{\omega}, slant \eqn{\alpha} 
#' and degree of freedom \eqn{\nu} parameters of skew-\eqn{t} distribution,
#' from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1), 
#' skewness and kurtosis.  
#' 
#' `alpha=0` solves 
#' \eqn{(\omega, \nu)} parameters of \eqn{t}-distribution,
#' from user-specified \eqn{\sigma} and kurtosis.
#' This is a non-skewed distribution, thus 
#' the location parameter \eqn{\xi=\mu=0}, and the slant parameter \eqn{\alpha=0}.
#' 
#' 
#' 
#' @returns
#' Function [moment2st()] returns a \link[base]{length}-4 \link[base]{numeric} \link[base]{vector} 
#' \eqn{(\xi, \omega, \alpha, \nu)}.
#' 
#' @keywords internal
#' @importFrom stats optim
#' @export
moment2st <- function(
    mean = 0, sd = 1, skewness, kurtosis,
    alpha,
    ...
) {
  
  # starting value `nu = 10`
  # ?stats::optim does not allow Inf starting value
  
  opt <- if (!missing(alpha)) {
    
    if (!identical(alpha, 0)) stop('use `alpha = 0` to constrain slant parameter alpha at 0')
    if (!missing(skewness) && (skewness != 0)) stop('skewness must be zero, for non-slant t-distribution')
    
    if (mean != 0) {
      target <- c(mean, sd, kurtosis)
      optim(par = c(xi = 0, omega = 1, nu = 10), fn = \(x) {
        nu <- x[3L]
        m <- moment_init(raw1 = 0, raw2 = nu/(nu-2), raw3 = 0, raw4 = 3*nu^2/(nu-2)/(nu-4))
        z <- c(
          mean = mean_moment_(m, location = x[1L], scale = x[2L]),
          sd = sd_moment_(m, scale = x[2L]),
          kurtosis = kurtosis_moment_(m)
        )
        crossprod(z - target)
      })
    } else {
      target <- c(sd, kurtosis)
      optim(par = c(omega = 1, nu = 10), fn = \(x) {
        nu <- x[2L]
        m <- moment_init(raw1 = 0, raw2 = nu/(nu-2), raw3 = 0, raw4 = 3*nu^2/(nu-2)/(nu-4))
        z <- c(
          sd = sd_moment_(m, scale = x[1L]),
          kurtosis = kurtosis_moment_(m)
        )
        crossprod(z - target)
      })
    }
    
  } else {
  
    target <- c(mean, sd, skewness, kurtosis)
  
    optim(par = c(xi = 0, omega = 1, alpha = 0, nu = 10), fn = \(x) {
      alpha <- x[3L]
      nu <- x[4L]
      delta <- alpha / sqrt(1 + alpha^2)
      b <- sqrt(nu/pi) * gamma(nu/2 - 1/2) / gamma(nu/2) # equation (29); https://arxiv.org/pdf/0911.2342.pdf
      r1 <- b * delta
      m <- moment_init(
        raw1 = r1,
        raw2 = nu/(nu-2),
        raw3 = r1 * (3-delta^2) * nu/(nu-3),
        raw4 = 3*nu^2/(nu-2)/(nu-4)
      )
      z <- c(
        mean = mean_moment_(m, location = x[1L], scale = x[2L]),
        sd = sd_moment_(m, scale = x[2L]),
        skewness = skewness_moment_(m),
        kurtosis = kurtosis_moment_(m)
      )
      crossprod(z - target)
    })
    
  }
  
  return(opt$par)
  
}



