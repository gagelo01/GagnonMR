#' Calculate Fstat from harmonised data with rsquared
#'
#' @param dat harmonised data
#'
#' @return a value that is the instrument Fstat
#' @export
fstat_fromdat <- function(dat){
  dat_split <- split(dat, dat$id.exposure)
  list_rep <- lapply(dat_split, function(dat) {
    k <- nrow(dat)
    n <- mean(dat$samplesize.exposure)
    r2sum <- sum(dat$rsq.exposure)
    Fstat <- ((n - k - 1)/k) * (r2sum/(1 - r2sum))
    return(rep(Fstat, k))
  })
 return( unlist(list_rep))
}

