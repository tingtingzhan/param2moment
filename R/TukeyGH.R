



#' @title Moments of Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description
#' Moments of Tukey \eqn{g}-&-\eqn{h} distribution.
#' 
#' @param A \link[base]{numeric} scalar or \link[base]{vector},
#' location parameter \eqn{A}
#' 
#' @param B \link[base]{numeric} scalar or \link[base]{vector},
#' scale parameter \eqn{B}
#' 
#' @param g \link[base]{numeric} scalar or \link[base]{vector},
#' skewness parameter \eqn{g}
#' 
#' @param h \link[base]{numeric} scalar or \link[base]{vector},
#' elongation parameter \eqn{h}
#' 
#' @returns
#' Function [moment_GH()] returns a \linkS4class{moment} object.
#' 
#' @references 
#' Raw moments of Tukey \eqn{g}-&-\eqn{h} distribution: \doi{10.1002/9781118150702.ch11}
#' 
#' @examples
#' A = 3; B = 1.5; g = .7; h = .01
#' moment_GH(A = A, B = B, g = 0, h = h)
#' moment_GH(A = A, B = B, g = g, h = 0)
#' moment_GH(A = A, B = B, g = g, h = h)
#' 
#' @importFrom methods new
#' @export
moment_GH <- function(A = 0, B = 1, g = 0, h = 0) {
  do.call(what = new, args = c(list(Class = 'moment', location = A, scale = B), moment_GH_(g = g, h = h)))
}

moment_GH_ <- function(g = 0, h = 0) {
  tmp <- data.frame(g, h) # recycling
  g <- tmp[[1L]]
  h <- tmp[[2L]]
  
  g0 <- (g == 0)
  
  # 1st-4th raw moment E(Y^n), when `g = 0`
  mu <- r3 <- rep(0, times = length(g))
  r2 <- 1 / (1-2*h) ^ (3/2) # (45a), page 502
  r4 <- 3 / (1-4*h) ^ (5/2) # (45b), page 502
  
  if (any(!g0)) {
    # {r}aw moment E(Y^n), when `g != 0`
    r_g <- function(n) { # equation (47), page 503
      tmp <- lapply(0:n, FUN = function(i) {
        (-1)^i * choose(n,i) * exp((n-i)^2 * g^2 / 2 / (1-n*h))
      })
      suppressWarnings(Reduce(f = `+`, tmp) / g^n / sqrt(1-n*h)) # warnings for `h > 1/n`
    }
    
    mu[!g0] <- r_g(1L)[!g0]
    r2[!g0] <- r_g(2L)[!g0]
    r3[!g0] <- r_g(3L)[!g0]
    r4[!g0] <- r_g(4L)[!g0]
  }
  
  moment_int(distname = 'GH', mu = mu, raw2 = r2, raw3 = r3, raw4 = r4)
  
}



#' @title Solve Tukey \eqn{g}-&-\eqn{h} Parameters from Moments
#' 
#' @description
#' Solve Tukey \eqn{g}-, \eqn{h}- and \eqn{g}-&-\eqn{h} distribution parameters 
#' from mean, standard deviation, skewness and kurtosis.
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
#' Function [moment2GH()] solves the 
#' location \eqn{A}, scale \eqn{B}, skewness \eqn{g} 
#' and elongation \eqn{h} parameters of Tukey \eqn{g}-&-\eqn{h} distribution,
#' from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1), 
#' skewness and kurtosis.  
#' 
#' @returns
#' Function [moment2GH()] returns a \link[base]{length}-4 
#' \link[base]{numeric} \link[base]{vector} \eqn{(A, B, g, h)}.
#' 
#' @examples
#' moment2GH(skewness = .2, kurtosis = .3)
#' 
#' @name moment2GH
#' @importFrom stats optim
#' @export
moment2GH <- function(mean = 0, sd = 1, skewness, kurtosis) {
  optim(par = c(A = 0, B = 1, g = 0, h = 0), fn = function(x) {
    mm <- moment_GH_(g = x[3L], h = x[4L])
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm), kurtosis_moment_(mm)) - c(mean, sd, skewness, kurtosis))
  })$par
}


#' @rdname moment2GH
#' 
#' @details
#' An educational and demonstration function [moment2GH_h_demo()] solves 
#' \eqn{(B, h)} parameters of Tukey \eqn{h}-distribution,
#' from user-specified \eqn{\sigma} and kurtosis.
#' This is a non-skewed distribution, thus 
#' the location parameter \eqn{A=\mu=0}, and the skewness parameter \eqn{g=0}.
#' 
#' @returns
#' Function [moment2GH_h_demo()] returns a \link[base]{length}-2 
#' \link[base]{numeric} \link[base]{vector} \eqn{(B, h)}.
#' 
#' @examples
#' moment2GH_h_demo(kurtosis = .3)
#' 
#' @importFrom stats optim
#' @export
moment2GH_h_demo <- function(sd = 1, kurtosis) {
  optim(par = c(B = 1, h = 0), fn = function(x) {
    mm <- moment_GH_(g = 0, h = x[2L])
    crossprod(c(sd_moment_(mm, scale = x[1L]), kurtosis_moment_(mm)) - c(sd, kurtosis))
  })$par
}

#' @rdname moment2GH
#' 
#' @details
#' An educational and demonstration function [moment2GH_g_demo()] solves  
#' \eqn{(A, B, g)} parameters of Tukey \eqn{g}-distribution,
#' from user-specified \eqn{\mu}, \eqn{\sigma} and skewness.
#' For this distribution, the elongation parameter \eqn{h=0}.
#' 
#' @returns
#' Function [moment2GH_g_demo()] returns a \link[base]{length}-3 
#' \link[base]{numeric} \link[base]{vector} \eqn{(A, B, g)}.
#' 
#' @examples
#' moment2GH_g_demo(skewness = .2)
#' 
#' @importFrom stats optim
#' @export
moment2GH_g_demo <- function(mean = 0, sd = 1, skewness) {
  optim(par = c(A = 0, B = 1, g = 0), fn = function(x) {
    mm <- moment_GH_(g = x[3L], h = 0)
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm)) - c(mean, sd, skewness))
  })$par
}



