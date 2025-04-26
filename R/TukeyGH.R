



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
#' @keywords internal
#' @export
moment_GH <- function(A = 0, B = 1, g = 0, h = 0) {
  
  list(g = g, h = h) |>
    .mapply(FUN = moment_GH_, MoreArgs = NULL) |>
    c(list(FUN = c, SIMPLIFY = FALSE)) |>
    do.call(what = mapply) |>
    c(list(Class = 'moment', location = A, scale = B)) |>
    do.call(what = 'new')
  
}


moment_GH_ <- function(g, h) {
  
  if (g == 0) {
    return(list(
      raw1 = 0, 
      raw2 = 1 / (1-2*h) ^ (3/2), # (45a), page 502
      raw3 = 0,
      raw4 = 3 / (1-4*h) ^ (5/2) # (45b), page 502
    ))
  }
  
  r_g <- \(n) { # equation (47), page 503
    if (h > 1/n) return(NaN)
    tmp <- (0:n) |>
      vapply(FUN = \(i) {
        (-1)^i * choose(n,i) * exp((n-i)^2 * g^2 / 2 / (1-n*h))
      }, FUN.VALUE = NA_real_) |>
      sum()
    tmp / g^n / sqrt(1-n*h)
  }
  
  return(list(
    raw1 = r_g(1L),
    raw2 = r_g(2L),
    raw3 = r_g(3L),
    raw4 = r_g(4L)
  ))
  
}





moment_GH_OLD <- function(g = 0, h = 0) {
  tmp <- data.frame(g, h) # recycling
  g <- tmp[[1L]]
  h <- tmp[[2L]]
  
  g0 <- (g == 0)
  
  # 1st-4th raw moment E(Y^n), when `g = 0`
  r1 <- r3 <- rep(0, times = length(g))
  r2 <- 1 / (1-2*h) ^ (3/2) # (45a), page 502
  r4 <- 3 / (1-4*h) ^ (5/2) # (45b), page 502
  
  if (any(!g0)) {
    # {r}aw moment E(Y^n), when `g != 0`
    r_g <- \(n) { # equation (47), page 503
      tmp <- lapply(0:n, FUN = \(i) {
        (-1)^i * choose(n,i) * exp((n-i)^2 * g^2 / 2 / (1-n*h))
      })
      suppressWarnings(Reduce(f = `+`, tmp) / g^n / sqrt(1-n*h)) # warnings for `h > 1/n`
    }
    
    r1[!g0] <- r_g(1L)[!g0]
    r2[!g0] <- r_g(2L)[!g0]
    r3[!g0] <- r_g(3L)[!g0]
    r4[!g0] <- r_g(4L)[!g0]
  }
  
  moment_init(distname = 'GH', raw1 = r1, raw2 = r2, raw3 = r3, raw4 = r4)
  
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
#' @keywords internal
#' @name moment2GH
#' @importFrom stats optim
#' @export
moment2GH <- function(
    mean = 0, sd = 1, skewness, kurtosis,
    g,
    h,
    ...
) {
  
  g0 <- !missing(g) && identical(g, 0)
  h0 <- !missing(h) && identical(h, 0)
  
  if (g0 & h0) stop('just normal distribution')
  
  opt <- if (g0) {
    if (mean != 0) {
      target <- c(mean, sd, kurtosis)
      optim(par = c(A = 0, B = 1, h = 0), fn = function(x) {
        m <- moment_GH_(g = 0, h = x[3L]) |>
          do.call(what = moment_init)
        z <- c(
          mean = mean_moment_(m, location = x[1L], scale = x[2L]),
          sd = sd_moment_(m, scale = x[2L]),
          kurtosis = kurtosis_moment_(m)
        )
        crossprod(z - target)
      })
    } else {
      target <- c(sd, kurtosis)
      optim(par = c(B = 1, h = 0), fn = function(x) {
        m <- moment_GH_(g = 0, h = x[2L]) |>
          do.call(what = moment_init)
        z <- c(
          sd = sd_moment_(m, scale = x[1L]),
          kurtosis = kurtosis_moment_(m)
        )
        crossprod(z - target)
      })
    }
  } else if (h0) {
    target <- c(mean, sd, skewness)
    optim(par = c(A = 0, B = 1, g = 0), fn = function(x) {
      m <- moment_GH_(g = x[3L], h = 0) |>
        do.call(what = moment_init)
      z <- c(
        mean = mean_moment_(m, location = x[1L], scale = x[2L]),
        sd = sd_moment_(m, scale = x[2L]),
        skewness = skewness_moment_(m)
      )
      crossprod(z - target)
    })
  } else {
    target <- c(mean, sd, skewness, kurtosis)
    optim(par = c(A = 0, B = 1, g = 0, h = 0), fn = function(x) {
      m <- moment_GH_(g = x[3L], h = x[4L]) |>
        do.call(what = moment_init)
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






