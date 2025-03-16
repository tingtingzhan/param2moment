
#' @title Moment to Parameters: A Batch Process
#' 
#' @description
#' Converts multiple sets of moments to multiple sets of distribution parameters.
#' 
#' @param distname \link[base]{character} scalar, distribution name.
#' Currently supported are `'GH'` for Tukey \eqn{g}-&-\eqn{h} distribution,
#' `'sn'` for skew-normal distribution and `'st'` for skew-\eqn{t} distribution
#' 
#' @param FUN \link[base]{name} or \link[base]{character} scalar,
#' (name of) \link[base]{function} used to solve the distribution parameters from moments.
#' Default is `paste0('moment2', distname)`, e.g., [moment2GH] will be used for `distname = 'GH'`.
#' To use one of the educational functions, specify
#' `FUN = moment2GH_g_demo` or `FUN = 'moment2GH_g_demo'`.
#' 
#' @param ... \link[base]{numeric} scalars, 
#' some or all of `mean`, `sd`, `skewness` and `kurtosis`
#' (length will be recycled).
#' 
#' @returns
#' Function [moment2param()] returns a \link[base]{list} of \link[base]{numeric} \link[base]{vector}s.
#' 
#' @examples
#' skw = c(.2, .5, .8)
#' krt = c(.5, 1, 1.5)
#' moment2param(distname = 'GH', skewness = skw, kurtosis = krt)
#' moment2param(distname = 'st', skewness = skw, kurtosis = krt)
#' 
#' @export
moment2param <- function(distname, FUN = paste0('moment2', distname), ...) {
  .mapply(
    FUN = FUN, 
    dots = as.list.data.frame(data.frame(..., check.names = FALSE, fix.empty.names = FALSE)), 
    # use ?base::data.frame to recycle the length
    MoreArgs = NULL
  )
}






