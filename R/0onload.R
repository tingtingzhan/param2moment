
.onLoad <- function(libname, pkgname = 'groupedHyperframe') {
  
  Sys.setenv('_R_CHECK_LIMIT_CORES_' = 'false') 
  # otherwise ?parallel:::.check_ncores causes error when ?devtools::check
  
  options(use_unicode = TRUE)
  
}

