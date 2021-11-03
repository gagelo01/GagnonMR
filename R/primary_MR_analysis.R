#' Perform Primary MR analysis
#'
#' @param dat a harmonised data.frame same as in TwoSampleMR
#'
#' @return return a data.frame of either wald ratio or IVW.
#' @export

primary_MR_analysis <- function(dat) {

  if(nrow(dat) == 1) {
    res <- mr(dat, method_list=c("mr_wald_ratio"))
    res$nsnp <-dat$SNP
  } else {
    res <- mr(dat, method_list=c("mr_ivw"))
  }

  return(res)

}
