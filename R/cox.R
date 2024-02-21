#' Tidy coxph output
#'
#' @param output output from coxph
#'
#' @return copph output clean in a data.table
tidy_cox_output <- function(models) {
  k<-summary(models)
  k1 <- data.table(HR = k$coefficients[,c("exp(coef)")], pval = k$coefficients[,c("Pr(>|z|)")], exposure = rownames(k$coefficients))
  k2 <- data.table(lci = k$conf.int[, c("lower .95")], uci = k$conf.int[, c("upper .95")], exposure = rownames(k$conf.int))
  k3 <- merge(k1,k2,by="exposure") %>% as.data.table(.)
  return(k3)
}

#' proportional hazard assumption of coxph
#'
#' @param models coxph output
#'
#' @return coxph output proportionl hazard clean in a data.table

tidy_ph_assumption <- function(models) {
  test.ph <- survival::cox.zph(models)
  dt <- test.ph$table %>% as.data.frame()
  colnames(dt)<-paste0("ph_", colnames(dt))
  dt$exposure <- rownames(dt)
  return(dt)
}


#' Perform a coxph analysis and return in a tidy data.table
#'
#' @param dat the data.table to compute the model
#' @param IV The name of the variable you are interested in as a predictor
#' @param DV The variable name need to have *_time and *_censored in the data
#' @param cov_inc the character of the covariate such as " + age_enrollment + sex
#'
#' @return a tdy coxph ouput with coxph assumptions
#' @export

run_coxph_wrapper<- function(dat, #The data.frame
                             IV, #The name of the variable you are interested in as a predictor
                             DV, #The variable name need to have *_time and *_censored in the data
                             cov_inc #the character of the covariate such as " + age_enrollment + sex"
) {
  message(paste0("running on ", DV, "~", IV, cov_inc))
  formulas <- as.formula(paste0("survival::Surv(",DV,"_time,", DV, "_censored) ~ ",IV, cov_inc))
  models <-survival::coxph(formulas,data=dat)
  res <- GagnonMR:::tidy_cox_output(models)
  ph <- GagnonMR:::tidy_ph_assumption(models)
  res <- merge(res, ph, by = "exposure", all = TRUE)
  res<-res[!grepl("GLOBAL", res$exposure), ]
  res$IV<-IV
  res$outcome <- DV
  res$cov_inc<-cov_inc
  return(res)}
