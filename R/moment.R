

#' @title Raw, Central and Standardized Moments
#' 
#' @description ..
#' 
#' @slot distname \link[base]{character} scalar, name of distribution,
#' e.g., `'norm'` for normal, `'sn'` for skew-normal, `'st'` for skew-\eqn{t}, 
#' and `'GH'` for Tukey \eqn{g}-&-\eqn{h} distribution,
#' following the nomenclature of \link[stats]{dnorm}, \link[sn]{dsn}, \link[sn]{dst} and `QuantileGH::dGH`
#' 
#' @slot location,scale \link[base]{numeric} scalars or \link[base]{vector}s, 
#' location and scale parameters
#' 
#' @slot raw1,raw2,raw3,raw4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 1st or higher order *raw* moments
#' 
#' @slot central2,central3,central4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 2nd or higher order *central* moments
#' 
#' @slot standardized3,standardized4 \link[base]{numeric} scalars or \link[base]{vector}s, 
#' 3rd or higher *standardized* moments
# 
#' 
#' @note
#' Potential name clash with function `e1071::moment`.
#' 
#' 
#' @export
setClass(Class = 'moment', slots = c(
  distname = 'character',
  location = 'numeric', scale = 'numeric',
  raw1 = 'numeric', raw2 = 'numeric', raw3 = 'numeric', raw4 = 'numeric',
  central2 = 'numeric', central3 = 'numeric', central4 = 'numeric',
  standardized3 = 'numeric', standardized4 = 'numeric'
))



setMethod(f = initialize, signature = 'moment', definition = function(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  x@central2 <- x@raw2 - x@raw1^2
  sigma <- sqrt(x@central2)
  x@central3 <- x@raw3 - 3*x@raw2*x@raw1 + 2*x@raw1^3
  x@central4 <- x@raw4 - 4*x@raw3*x@raw1 + 6*x@raw2*x@raw1^2 - 3*x@raw1^4
  
  x@standardized3 <- x@central3 / sigma^3
  x@standardized4 <- x@central4 / sigma^4
  
  return(x)
  
})



# barebone initialize, in compute extensive jobs
moment_init <- function(raw1, raw2, raw3, raw4, ...) {
  
  central2 <- raw2 - raw1^2
  sigma <- sqrt(central2)
  central3 <- raw3 - 3*raw2*raw1 + 2*raw1^3
  central4 <- raw4 - 4*raw3*raw1 + 6*raw2*raw1^2 - 3*raw1^4
  
  return(list(
    raw1 = raw1, #raw2 = raw2, raw3 = raw3, raw4 = raw4,
    central2 = central2, central3 = central3, central4 = central4,
    standardized3 = central3 / sigma^3, 
    standardized4 = central4 / sigma^4
  ))
  
}


mean_moment_ <- function(x, location, scale) location + scale * x$raw1
sd_moment_ <- function(x, scale) scale * sqrt(x$central2)
skewness_moment_ <- function(x) x$standardized3
kurtosis_moment_ <- function(x) x$standardized4 - 3





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
#' @export
setMethod(f = show, signature = signature(object = 'moment'), definition = function(object) {
  
  object@distname |> sprintf(fmt = 'Distribution: %s\n') |> cat()
  
  mean_moment <- \(x) x@location + x@scale * x@raw1
  object |> mean_moment() |> sprintf(fmt = '%.3f') |> paste0(collapse = ', ') |> sprintf(fmt = 'Mean: %s\n') |> cat()
  
  sd_moment <- \(x) x@scale * sqrt(x@central2)
  object |> sd_moment() |> sprintf(fmt = '%.3f') |> paste0(collapse = ', ') |> sprintf(fmt = 'Standard Deviation: %s\n') |> cat()
  
  skewness_moment <- \(x) x@standardized3
  object |> skewness_moment() |> sprintf(fmt = '%.3f') |> paste0(collapse = ', ') |> sprintf(fmt = 'Skewness: %s\n') |> cat()
  
  kurtosis_moment <- \(x) x@standardized4 - 3
  object |> kurtosis_moment() |> sprintf(fmt = '%.3f') |> paste0(collapse = ', ') |> sprintf(fmt = '(Excess) Kurtosis: %s\n') |> cat()
})

