#' test_object
#'
#' @return a list containing useful test object
#'
#' @examples
test_object <- function() {
  region <- c(total = "20:15000000-17200000", subset= "20:15500000-16000000")
fn_exp1 <- system.file("extdata","bmi.vcf.gz", package="GagnonMR")
fn_exp2 <- system.file("extdata","whr.vcf.gz", package="GagnonMR")
fn_out <- system.file("extdata","cad.vcf.gz", package="GagnonMR")
fn_out2 <- system.file("extdata","t2d.vcf.gz", package="GagnonMR")
parameters <- GagnonMR::default_param()
parameters$path <- paste0(system.file("extdata", package = "GagnonMR"), "/")
return(list(region = region, fn_exp1 = fn_exp1, fn_exp2 = fn_exp2, fn_out = fn_out, fn_out2 = fn_out2, parameters = parameters))
}
