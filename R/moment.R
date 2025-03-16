

#' @title Raw, Central and Standardized Moments, and other Distribution Characteristics
#' 
#' @description
#' Up to 4th raw \eqn{\text{E}(Y^n)}, \href{https://en.wikipedia.org/wiki/Central_moment}{central} \eqn{\text{E}[(Y-\mu)^n]} and 
#' \href{https://en.wikipedia.org/wiki/Standardized_moment}{standardized moments} \eqn{\text{E}[(Y-\mu)^n/\sigma^n]} of the random variable
#' \deqn{Y = (X - \text{location})/\text{scale}}
#' 
#' Also, the mean, standard deviation, skewness and excess kurtosis of the random variable \eqn{X}.
#' 
#' @slot distname \link[base]{character} scalar, name of distribution,
#' e.g., `'norm'` for normal, `'sn'` for skew-normal, `'st'` for skew-\eqn{t}, 
#' and `'GH'` for Tukey \eqn{g}-&-\eqn{h} distribution,
#' following the nomenclature of \link[stats]{dnorm}, \link[sn]{dsn}, \link[sn]{dst} and `QuantileGH::dGH`
#' 
#' @slot location,scale \link[base]{numeric} scalars or \link[base]{vector}s, 
#' location and scale parameters
#' 
#' @slot mu \link[base]{numeric} scalar or \link[base]{vector}, 
#' 1st *raw* moment \eqn{\mu = \text{E}(Y)}. 
#' Note that the 1st central moment \eqn{\text{E}(Y-\mu)} and
#' standardized moment \eqn{\text{E}(Y-\mu)/\sigma} are both 0.
#' 
#' @slot raw2,raw3,raw4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 2nd or higher *raw* moments \eqn{\text{E}(Y^n)}, \eqn{n\geq 2}
#' 
#' @slot central2,central3,central4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 2nd or higher \href{https://en.wikipedia.org/wiki/Central_moment}{*central* moments}, \eqn{\sigma^2 = \text{E}[(Y-\mu)^2]} and 
#' \eqn{\text{E}[(Y-\mu)^n]}, \eqn{n\geq 3}
#' 
#' @slot standardized3,standardized4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 3rd or higher \href{https://en.wikipedia.org/wiki/Standardized_moment}{*standardized* moments}, 
#' \href{https://en.wikipedia.org/wiki/Skewness}{skewness} \eqn{\text{E}[(Y-\mu)^3]/\sigma^3} and
#' \href{https://en.wikipedia.org/wiki/Kurtosis}{kurtosis} \eqn{\text{E}[(Y-\mu)^4]/\sigma^4}.
# \eqn{\text{E}[(Y-\mu)^n]/\sigma^n}, \eqn{n=5,\cdots}. 
#' Note that the 2nd standardized moment is 1
#' 
#' @details
#' 
#' For \eqn{Y = (X - \text{location})/\text{scale}}, 
#' let \eqn{\mu = \text{E}(Y)}, then, according to  
#' \href{https://en.wikipedia.org/wiki/Binomial_theorem}{Binomial theorem},
#' the 2nd to 4th central moments of \eqn{Y} are,
#' \deqn{\text{E}[(Y-\mu)^2] = \text{E}(Y^2) - 2\mu \text{E}(Y) + \mu^2 = \text{E}(Y^2) - \mu^2}
#' \deqn{\text{E}[(Y-\mu)^3] = \text{E}(Y^3) - 3\mu \text{E}(Y^2) + 3\mu^2 \text{E}(Y) - \mu^3 = \text{E}(Y^3) - 3\mu \text{E}(Y^2) + 2\mu^3}
#' \deqn{\text{E}[(Y-\mu)^4] = \text{E}(Y^4) - 4\mu \text{E}(Y^3) + 6\mu^2 \text{E}(Y^2) - 4\mu^3 \text{E}(Y) + \mu^4 = \text{E}(Y^4) - 4\mu \text{E}(Y^3) + 6\mu^2 \text{E}(Y^2) - 3\mu^4}
#' 
#' The distribution characteristics of \eqn{Y} are,
#' \deqn{\mu_Y = \mu}
#' \deqn{\sigma_Y = \sqrt{\text{E}[(Y-\mu)^2]}}
#' \deqn{\text{skewness}_Y = \text{E}[(Y-\mu)^3] / \sigma^3_Y}
#' \deqn{\text{kurtosis}_Y = \text{E}[(Y-\mu)^4] / \sigma^4_Y - 3}
#' 
#' The distribution characteristics of \eqn{X} are
#' \eqn{\mu_X = \text{location} + \text{scale}\cdot \mu_Y},
#' \eqn{\sigma_X = \text{scale}\cdot \sigma_Y},
#' \eqn{\text{skewness}_X = \text{skewness}_Y}, and
#' \eqn{\text{kurtosis}_X = \text{kurtosis}_Y}.
#' 
#' 
#' @note
#' Potential name clash with function `e1071::moment`.
#' 
#' 
#' @importFrom methods setClass
#' @export
setClass(Class = 'moment', slots = c(
  distname = 'character',
  location = 'numeric', scale = 'numeric',
  mu = 'numeric',
  raw2 = 'numeric', raw3 = 'numeric', raw4 = 'numeric',
  central2 = 'numeric', central3 = 'numeric', central4 = 'numeric',
  standardized3 = 'numeric', standardized4 = 'numeric'
))







# converts raw moments to central (and standardized) moments, using binomial theorem
moment_int <- function(
    distname, 
    mu, raw2, raw3, raw4,
    ...
) {
  central2 <- raw2 - mu^2
  sigma <- sqrt(central2)
  central3 <- raw3 - 3*raw2*mu + 2*mu^3
  central4 <- raw4 - 4*raw3*mu + 6*raw2*mu^2 - 3*mu^4
  
  standardized3 <- central3 / sigma^3
  standardized4 <- central4 / sigma^4
  
  return(list(
    distname = distname, 
    mu = mu, raw2 = raw2, raw3 = raw3, raw4 = raw4,
    central2 = central2, central3 = central3, central4 = central4,
    standardized3 = standardized3, standardized4 = standardized4
  )) # do NOT return S4; saves time!
  
}


mean_moment_ <- function(x, location, scale) location + scale * x$mu
mean_moment <- function(x) x@location + x@scale * x@mu
sd_moment_ <- function(x, scale) scale * sqrt(x$central2)
sd_moment <- function(x) x@scale * sqrt(x@central2)
skewness_moment_ <- function(x) x$standardized3
skewness_moment <- function(x) x@standardized3
kurtosis_moment_ <- function(x) x$standardized4 - 3
kurtosis_moment <- function(x) x@standardized4 - 3




# stats::dt
# Raw moments of \eqn{t}-distribution: \url{https://en.wikipedia.org/wiki/Student%27s_t-distribution}
# Raw moments of non-central \eqn{t}-distribution: \url{https://en.wikipedia.org/wiki/Noncentral_t-distribution}





#' @title Show \linkS4class{moment}
#' 
#' @description
#' Print S4 object \linkS4class{moment} in a pretty manner.
#' 
#' @param object a \linkS4class{moment} object
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{moment} object 
#' does not have a returned value.
#' 
#' @keywords internal
#' @importFrom methods show signature
#' @export
setMethod(f = show, signature = signature(object = 'moment'), definition = function(object) {
  cat(sprintf('Distribution: %s\n', object@distname))
  cat(sprintf('Mean: %s\n', paste0(sprintf(fmt = '%.3f', mean_moment(object)), collapse = ', ')))
  cat(sprintf('Standard Deviation: %s\n', paste0(sprintf(fmt = '%.3f', sd_moment(object)), collapse = ', ')))
  cat(sprintf('Skewness: %s\n', paste0(sprintf(fmt = '%.3f', skewness_moment(object)), collapse = ', ')))
  cat(sprintf('Kurtosis: %s\n', paste0(sprintf(fmt = '%.3f', kurtosis_moment(object)), collapse = ', ')))
})

