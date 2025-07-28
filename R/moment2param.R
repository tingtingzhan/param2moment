
#' @title Moment to Parameters: A Batch Process
#' 
#' @description
#' Converts multiple sets of moments to multiple sets of distribution parameters.
#' 
#' @param distname \link[base]{character} scalar, distribution name.
#' Currently supported are `'GH'` for Tukey \eqn{g}-&-\eqn{h} distribution,
#' `'sn'` for skew-normal distribution and `'st'` for skew-\eqn{t} distribution
#' 
#' @param MoreArgs ..
#' 
#' @param ... \link[base]{numeric} scalars, 
#' some or all of `mean`, `sd`, `skewness` and `kurtosis`
#' (length will be recycled).
#' 
#' @returns
#' Function [moment2param()] returns a \link[base]{list} of \link[base]{numeric} \link[base]{vector}s.
#' 
#' @keywords internal
#' @export
moment2param <- function(distname, ..., MoreArgs = NULL) {
  
  # length recycled!
  ret <- .mapply(FUN = paste0('moment2', distname), dots = list(...), MoreArgs = MoreArgs)
  
  attr(ret, which = 'distname') <- distname
  class(ret) <- 'moment2param'
  return(ret)
  
}

#' @title autoplot.moment2param
#' 
#' @param object ..
#' 
#' @param hjust ..
#' 
#' @param ... ..
#' 
#' @keywords internal
#' @importFrom ggplot2 autoplot ggplot aes stat_function labs
#' @importFrom geomtextpath geom_textpath 
#' @export autoplot.moment2param
#' @export
autoplot.moment2param <- function(object, hjust = .5, ...) {
  
  distname <- attr(object, which = 'distname', exact = TRUE)
  
  label <- object |> 
    lapply(FUN = \(x) { # tzh's ?gg.tzh:::getval_ function
      sprintf(fmt = '%s = %.3g', names(x), x) |>
        sub(pattern = '([-]?)0[.]', replacement = '\\1.') |> # remove leading zero
        paste0(collapse = '; ') |> 
        gsub(pattern = 'Inf', replacement = '\u221e') |>
        gsub(pattern = 'alpha', replacement = '\u03b1') |>
        gsub(pattern = 'lambda', replacement = '\u03bb') |>
        gsub(pattern = 'nu', replacement = '\u03bd') |>
        gsub(pattern = 'xi', replacement = '\u03be') |>
        gsub(pattern = 'omega', replacement = '\u03c9')
    })
  
  # tzh's ?gg.tzh::paths_function function
  
  mp <- object |>
    seq_along() |>
    sprintf(fmt = '%02d') |> # so that '10' is after '09'
    lapply(FUN = \(i) aes(color = i))
  
  lyr <- .mapply(
    FUN = stat_function, 
    dots = list(
      mapping = mp, 
      args = object,
      label = label, # recycled
      hjust = hjust # recycled
    ), 
    MoreArgs = list(
      fun = paste0('d', distname) |> get(), 
      geom = 'textpath',
      size = 2,
      show.legend = FALSE
    )
  )
  
  ggplot() + 
    lyr +
    labs(y = distname)
  
}


