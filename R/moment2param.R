
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
#' @param hjust ..
#' 
#' @param ... \link[base]{numeric} scalars, 
#' some or all of `mean`, `sd`, `skewness` and `kurtosis`
#' (length will be recycled).
#' 
#' @returns
#' Function [moment2param()] returns a \link[base]{list} of \link[base]{numeric} \link[base]{vector}s.
#' 
#' @keywords internal
#' @importFrom ggplot2 ggplot aes stat_function labs
#' @importFrom geomtextpath geom_textpath 
#' @export
moment2param <- function(
    distname, 
    MoreArgs = NULL,
    hjust = .5,
    ...
) {
  
  # length recycled!
  ret <- .mapply(FUN = paste0('moment2', distname), dots = list(...), MoreArgs = MoreArgs)
  
  label <- ret |> 
    lapply(FUN = \(x) { # tzh's ?gg.tzh:::getval_ function
      z <- sprintf(fmt = '%s = %.3g', names(x), x) |>
        sub(pattern = '([-]?)0[.]', replacement = '\\1.') |> # remove leading zero
        paste0(collapse = '; ')
      if (getOption('use_unicode')) {
        z <- z |> 
          gsub(pattern = 'Inf', replacement = '\u221e') |>
          gsub(pattern = 'alpha', replacement = '\u03b1') |>
          gsub(pattern = 'lambda', replacement = '\u03bb') |>
          gsub(pattern = 'nu', replacement = '\u03bd') |>
          gsub(pattern = 'xi', replacement = '\u03be') |>
          gsub(pattern = 'omega', replacement = '\u03c9')
      }
      return(z)
    })
  
  # tzh's ?gg.tzh::paths_function function
  
  mp <- ret |>
    seq_along() |>
    sprintf(fmt = '%02d') |> # so that '10' is after '09'
    lapply(FUN = \(i) aes(color = i))
  
  lyr <- .mapply(
    FUN = stat_function, 
    dots = list(
      mapping = mp, 
      args = ret,
      label = label, # recycled
      hjust = hjust # recycled
    ), 
    MoreArgs = list(
      fun = paste0('d', distname) |> get(), 
      geom = 'textpath', # must Depends geomtextpath; Imports does not work. why??
      size = 2,
      show.legend = FALSE
    )
  )
  
  p <- ggplot() + 
    lyr +
    labs(y = distname)
  
  attr(p, which = 'value') <- ret
  return(p)
  
}


