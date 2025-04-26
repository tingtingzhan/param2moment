

#' @title Moments of Gamma Distribution
#' 
#' @param scale,shape \link[base]{numeric} scalars or \link[base]{vector}s, 
#' scale \eqn{\theta} and shape \eqn{\alpha},
#' see function \link[stats]{dgamma}.
#' 
#' @returns
#' Function [moment_gamma()] returns a \linkS4class{moment} object.
#' 
#' @keywords internal
#' @export
moment_gamma <- function(scale = 1, shape) {
  new(Class = 'moment', 
      distname = 'sn', 
      location = 0, scale = scale,
      raw1 = gamma(shape + 1) / gamma(shape),
      raw2 = gamma(shape + 2) / gamma(shape),
      raw3 = gamma(shape + 3) / gamma(shape),
      raw4 = gamma(shape + 4) / gamma(shape))
}







#' @title Solve Skew-Normal Parameters from Moments
#' 
#' @description
#' Solve skew-normal parameters from mean, standard deviation and skewness.
#' 
#' @param ... \link[base]{numeric} scalars, ***two*** of 
#' `mean` \eqn{\mu},
#' standard deviation `sd` \eqn{\sigma},
#' `skewness`, and
#' `kurtosis`
#' 
#' @details
#' Function [moment2gamma()] solves the 
#' location \eqn{\xi}, scale \eqn{\omega} and slant \eqn{\alpha} parameters 
#' of skew-normal distribution,
#' from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1) and 
#' skewness.  
#' 
#' @returns
#' Function [moment2gamma()] returns a \link[base]{length}-2
#' \link[base]{numeric} \link[base]{vector} \eqn{(\alpha, \theta)}.
#' 
#' @keywords internal
#' @importFrom stats optim
#' @export
moment2gamma <- function(...) {
  
  target <- c(...)
  id <- names(target) %in% c('mean', 'sd', 'skewness', 'kurtosis')
  if (sum(id) != 2L) stop('')
   
  opt <- optim(par = c(scale = 1, shape = 1), fn = \(x) {
    shape <- x[2L] |> unname()
    m <- moment_init(
      raw1 = gamma(shape + 1) / gamma(shape), 
      raw2 = gamma(shape + 2) / gamma(shape), 
      raw3 = gamma(shape + 3) / gamma(shape), 
      raw4 = gamma(shape + 4) / gamma(shape))
    z <- c(
      mean = mean_moment_(m, location = 0L, scale = unname(x[1L])),
      sd = sd_moment_(m, scale = unname(x[1L])),
      skewness = skewness_moment_(m),
      kurtosis = kurtosis_moment_(m)
    )
    crossprod(z[names(target)] - target)
  })
  
  return(opt$par)
  
}






