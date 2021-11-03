#' Two Sample MR power analysis
#'
#' This method rely on the paper by Calculating statistical power in Mendelian
#' randomization studies Marie-Jo A Brion, Konstantin Shakhbazov, Peter M Visscher
#' International Journal of Epidemiology 2013 42: 1497-1501
#' I found the code on https://github.com/kn3in/mRnd/blob/master/functions.R
#'
#' @param N Sample size
#' @param alpha Type-I error rate (significance level) tipycally 0.05
#' @param R2xz Proportion of variance explained for the association between the SNP or allele score (Z) and the exposure variable (X)
#' @param varx Variance of the exposure variable (X)
#' @param vary Variance of the outcome variable (Y)
#' @param byx The regression coefficient βyx for the true underlying causal association between the exposure (X) and outcome (Y) variables
#' @param epower if epower is at NA it will return the power. if it is numeric between 0 and 1 it will return the samplesize.
#'
#' @return power for continuous outcome
#' @export

results_beta_based <- function(N, alpha, R2xz, varx, vary, byx, epower) {

  threschi <- qchisq(1 - alpha, 1)
  R2yz <- byx^2 * (varx / vary) * R2xz

  if(is.na(epower)) {

    NCP <- N * R2yz / (1 - R2yz)
    power <- 1 - pchisq(threschi, 1, NCP)
    dfpow<-data.frame(Parameter = c("Power", "NCP"), Value = c(power, NCP), Description = c("", "Non-Centrality-Parameter"))

  } else {

    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(epower)
    Z  <- (z1 + z2)^2
    N2 <- Z * (1 - R2yz) / R2yz
    dfpow<-data.frame(Parameter = "Sample Size", Value = N2)

  }
  return(dfpow[1,]$Value)
}


#' Calulate power for binary outcome
#'
#' @param N Sample size
#' @param alpha Type-I error rate (significance level) tipycally 0.05
#' @param R2xz Proportion of variance explained for the association between the SNP or allele score (Z) and the exposure variable (X)
#' @param K Proportion of cases in the study
#' @param OR True odds ratio of the outcome variable per standard deviation of the exposure variable
#' @param epower if epower is at NA it will return the power. if it is numeric between 0 and 1 it will return the samplesize.
#'
#' @return Power for binary outcome
#' @export

results_binary <- function(N, alpha, R2xz, K, OR, epower) {
  threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
  f.value <- 1 + N * R2xz / (1 - R2xz)

  if (is.na(epower)) {

    b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)

    v_MR <- (K * (1-K) - b_MR^2) / (N*R2xz)
    NCP <- b_MR^2 / v_MR

    # 2-sided test
    power <- 1 - pchisq(threschi, 1, NCP)
    dfpow<-data.frame(Parameter = c("Power", "NCP", "F-statistic"), Value = c(power, NCP, f.value), Description = c("", "Non-Centrality-Parameter", "The strength of the instrument"))


  } else {

    # Calculation of sample size given power
    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(epower)
    Z  <- (z1 + z2)^2

    b_01 <- K * ( OR/ (1 + K * (OR - 1)) -1)
    f <- K * (1-K) - b_01^2
    N1 <- Z * f / (b_01^2 * R2xz)
    N1 <- ceiling(N1)
    dfpow<-data.frame(Parameter = "Sample Size", Value = N1)

  }
  return(dfpow[1,]$Value)
}

#' Title
#'
#' @param dat an harmonised data sets
#' @param alpha the significance level you wish to use. tipycally 0.05
#' @param effect the true causal effect of the exposure on the outcome hypothetised.
#'
#' @return the power
#' @export

power.calculator_fromdat_vishner <- function(dat, alpha = 0.05, effect=0.1){

  N <- dat[1,]$samplesize.outcome
  R2xz <- sum(dat$rsq.exposure)
  epower <- NA

  #binary or continuous
  if(dat[1,]$units.outcome == "log odds") {

    K <- dat[1,]$ncase.outcome/N#proportion of case in the study
    OR <- exp(effect)
    results_binary(N, alpha, R2xz, K, OR, epower)

  }else {

    varx <- 1
    vary <- 1
    byx <- (effect*sqrt(varx)) / sqrt(vary)
    results_beta_based(N, alpha, R2xz, varx, vary, byx, epower)
  }
}
