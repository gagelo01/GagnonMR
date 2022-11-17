#' Calculate Fstat from harmonised data with rsquared
#'
#' @param dat_split list of harmonised data
#'
#' @return a value that is the instrument Fstat
#' @export
fstat_fromdat <- function(dat_split){
  list_rep <- lapply(dat_split, function(dat) {
    k <- nrow(dat)
    if(all(is.na(dat$ncase.exposure) & is.na(dat$ncontrol.exposure))) {
      n <- mean(dat$samplesize.exposure)
    } else {
      n <- mean( TwoSampleMR::effective_n(ncase = dat$ncontrol.exposure, ncontrol = dat$ncase.exposure))
    }
    r2sum <- sum(dat$rsq.exposure)
    Fstat <- ((n - k - 1)/k) * (r2sum/(1 - r2sum))
    return(Fstat)
  })
  return(unlist(list_rep))
}

