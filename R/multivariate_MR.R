#' Analogous of MendelianRandomization::mr_mvinput, but on a list of harmonised data from TwoSampleMR
#'
#' @param list_harm3 the list of harmonised data. can be optain by running on your harmonised data
#' "harm". split(harm, harm$exposure)
#'
#' @return An MRMVInput object identical to those obtained from MendelianRandomization::mr_mvinput.
#' @export

mr_mvinput_wrapper <- function(list_harm3){
  lapply(list_harm3, setDT)
  rsid_include <- Reduce(intersect, lapply(list_harm3, function(x) x$SNP))

  list_harm3 <- lapply(list_harm3, function(x) x[SNP %in% rsid_include,])
  stopifnot(all(list_harm3[[1]]$SNP == list_harm3[[2]]$SNP & list_harm3[[1]]$SNP == list_harm3[[length(list_harm3)]]$SNP))
  bon <- map(list_harm3, `[`, "beta.exposure")
  bx <- matrix(unlist(bon), ncol = length(list_harm3), nrow = nrow(bon[[1]]) , byrow = FALSE)
  bon <- map(list_harm3, `[`,"se.exposure")
  bxse <- matrix(unlist(bon), ncol = length(list_harm3), nrow = nrow(bon[[1]]) , byrow = FALSE)

  stopifnot(all(lapply(list_harm3, function(x) x$beta.outcome == list_harm3[[1]]$beta.outcome ) %>% unlist(.)))
  by <-list_harm3[[1]]$beta.outcome
  byse<-list_harm3[[1]]$se.outcome


  stopifnot( all(list_harm3[[1]]$effect_allele == list_harm3[[2]]$effect_allele &
                   list_harm3[[1]]$effect_allele == list_harm3[[length(list_harm3)]]$effect_allele) &
               all(list_harm3[[1]]$other_allele == list_harm3[[2]]$other_allele &
                     list_harm3[[1]]$other_allele == list_harm3[[length(list_harm3)]]$other_allele)&
               all(list_harm3[[1]]$eaf.exposure == list_harm3[[2]]$eaf.exposure &
                     list_harm3[[1]]$eaf.exposuree == list_harm3[[length(list_harm3)]]$eaf.exposure))

  inpoot  <- mr_mvinput(bx = bx,
                        bxse = bxse,
                        by = by,
                        byse = byse,
                        exposure = map_chr(list_harm3, `[`, 1, "exposure"),
                        outcome = unique(list_harm3[[1]]$outcome),
                        snps = list_harm3[[1]]$SNP,
                        effect_allele = list_harm3[[1]]$effect_allele.exposure,
                        other_allele = list_harm3[[1]]$other_allele.exposure,
                        eaf = list_harm3[[1]]$eaf.exposure)

  return(inpoot)

}


#' perform MendelianRandomization::mr_mvivw, but automatically format the data in a data.frame
#'
#' @param inpoot MendelianRandomization::mr_mvinput object
#'
#' @return a data.frame of a multivariate MR with IVW.
#' @export

mr_mvivw_wrapper <- function(inpoot) {
result <- MendelianRandomization::mr_mvivw(inpoot)
result_formated <- data.frame(Exposure = result@Exposure, Outcome = result@Outcome,
                     Estimate = result@Estimate, StdError = result@StdError, CILower = result@CILower,
                     CIUpper = result@CIUpper, Pvalue = result@Pvalue)

return(result_formated)
}

#' Perform MendelianRandomization::mr_mvegger, but a return a data.frame
#'
#' @param inpoot MendelianRandomization::mr_mvinput object
#' @param distribution same as mr_mvegger
#'
#' @return data.frame of a multivariate MR with multivariable weighted linear regression.
#' @export

mr_mvegger_wrapper <- function(inpoot, distribution = "normal") {
  result<-mr_mvegger(inpoot, distribution = distribution)

  result_formated <- data.frame(Exposure = c(result@Exposure, "(intercept)"), Outcome = result@Outcome,
                                Estimate = c(result@Estimate, result@Intercept),
                                StdError = c(result@StdError.Est, result@StdError.Int),
                                CILower = c(result@CILower.Est, result@CILower.Int),
                                CIUpper = c(result@CIUpper.Est, result@CIUpper.Int),
                                Pvalue = c(result@Pvalue.Est, result@Pvalue.Int))

  return(result_formated)
}

#' Align a list of exposure with the same rsids on the same allele (this function has not been validated and is outdated)
#'
#' @param list_inst_mvmr a list of different exposure
#' @param action action of harmonise_data
#'
#' @return the same list, but align on the same allele
align_list_exposure_mvmr <- function(list_inst_mvmr, action = 2) {
  inst_mvmr_outcome <- lapply(list_inst_mvmr, function(x) {colnames(x) <- gsub("exposure", "outcome", colnames(x))
  return(x)})
  inst_mvmr_harm <- lapply(inst_mvmr_outcome, function(x) harmonise_data(list_inst_mvmr[[1]], x, action = action))
  outcome_mvmr_same_allele <- lapply( inst_mvmr_harm, function(x) x[,c("SNP", colnames(x)[grepl("outcome", names(x))])])
  inst_mvmr_same_allele    <-  lapply(outcome_mvmr_same_allele, function(x)
    setnames(x, old = colnames(x), new = gsub("outcome", "exposure", colnames(x))))
  lapply(inst_mvmr_same_allele, setDT)
  all_inst_dt <- data.table::rbindlist(inst_mvmr_same_allele, fill = TRUE)
  return(all_inst_dt)
}



#' First function to run on a data.frame of mvmr_inst
#'
#' @param exposure_dat A data.frame of significant instrument for exemple the output of TwoSampleMR::extract_instruments(c("ukb-b-19953", "ukb-b-9405"))
#' @param d1 A data.frame of the effect of all SNP present in exposure_dat,
#' @param clump_r2 The default is 0.001
#' @param clump_kb The default is 10000
#' @param harmonise_strictness See the action option of harmonise_data. The default is 2.
#' @param pval_threshold Instrument detection p-value threshold. Default = 5e-8
#' @param clump_exp if NULL clump the SNPs with respect to the lowest P-value corresponding to any of the exposures. Otherwise,
#' provide the name of the exposure as character. Will clump by selecting the SNP with the lowest p-value for that exposure.
#' This procedure should be implemented to select a maximum of strong genetic instruments for that exposure,
#' @param should_clump if TRUE clump if FALSE don't
#' @return an object very similar to exposure_dat (theinput), but ready to go further in the mvmr pipeline
#' @export
prepare_for_mvmr<- function (exposure_dat, d1, clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2,
                             pval_threshold = 1, clump_exp = NULL, should_clump = TRUE)  {
  if(should_clump) {
    if(is.null(clump_exp)) {
      temp <- exposure_dat
      temp$id.exposure <- 1
      temp <- temp[order(temp$pval.exposure, decreasing = FALSE),]
      temp <- subset(temp, !duplicated(SNP))
    }else{
      temp <- d1
      temp <- temp[exposure == clump_exp,]
      temp <- subset(temp, !duplicated(SNP))
    }
    inst_all_clump = ieugwasr::ld_clump(data.frame(rsid=temp$SNP,pval=temp$pval.exposure, id = temp$id.exposure),
                                        clump_kb=clump_kb,clump_r2=clump_r2,
                                        plink_bin=genetics.binaRies::get_plink_binary(),
                                        bfile="/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
    exposure_dat <- subset(exposure_dat, SNP %in%     inst_all_clump$rsid)
  }

  id_exposure <- unique(exposure_dat$id.exposure)
  colnames(d1) <- gsub("exposure", "outcome", colnames(d1))
  d1<-subset(d1, d1$SNP %in% exposure_dat$SNP)
  stopifnot(length(unique(d1$id)) == length(id_exposure))
  d1 <- subset(d1, mr_keep.outcome)
  d2 <- subset(d1, id.outcome != id_exposure[1])
  d1 <- subset(d1, id.outcome == id_exposure[1])
  colnames(d1) <- gsub("outcome", "exposure", colnames(d1))
  d <- TwoSampleMR::harmonise_data(d1, d2, action = harmonise_strictness)
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure) - 1]
  d <- subset(d, SNP %in% keepsnps)
  dh1 <- subset(d, id.outcome == id.outcome[1], select = c("SNP", colnames(d)[grepl("exposure", colnames(d))]))
  dh2 <- subset(d, select = c("SNP", colnames(d)[grepl("outcome", colnames(d))]))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  dh <- rbindlist(list(dh1, dh2), fill = TRUE)
  return(dh)
}


#' Perform mmvmr ivw and all mvmr robust analyses in MendelianRandomization
#'
#' @param exposure_outcome_harmonized the output from TwoSampleMR::mv_harmonise_data
#' @param pheno_cov_exp the phenotype covariance of the exposure
#' @param only_IVW whether to perform only to perform IVW
#'
#' @return a data.frame with the results of all the models
#' @export
mv_multiple_MendelianRandomization <- function(exposure_outcome_harmonized, pheno_cov_exp = 0, only_IVW = FALSE) {
  x <- exposure_outcome_harmonized
  mriobj <- MendelianRandomization::mr_mvinput(bx = x$exposure_beta,
                                               bxse = x$exposure_se, by = x$outcome_beta, byse = x$outcome_se,
                                               exposure = x$expname$exposure, outcome = x$outname$outcome)
  mvmr_res <- vector(mode = "list", length = 2)
  x <- MendelianRandomization::mr_mvivw(object = mriobj, model = "default", robust = FALSE)
  ivwestimate <- data.frame(exposure = x@Exposure,
                            outcome = x@Outcome, b = x@Estimate,
                            se = x@StdError, lci = x@CILower,
                            uci = x@CIUpper, pval = x@Pvalue,
                            cochranQ = x@Heter.Stat[1], cochranQpval = x@Heter.Stat[2],
                            nsnp = x@SNPs, method = "Multivariable IVW") %>% as.data.table

  sres2  <- obtain_conditionalFstatistics(exposure_outcome_harmonized, pheno_cov_exp)
  if(only_IVW == TRUE) { return(merge(ivwestimate, sres2, by = "exposure", sort = FALSE))}

  x <- MendelianRandomization::mr_mvegger(object = mriobj)
  eggerestimate <- data.frame(exposure = x@Exposure, outcome = x@Outcome, b = x@Estimate,
                              se = x@StdError.Est, lci = x@CILower.Est, uci = x@CIUpper.Est, pval = x@Pvalue.Est,
                              cochranQ = x@Heter.Stat[1], cochranQpval = x@Heter.Stat[2],
                              nsnp = x@SNPs, method = "Multivariable Egger") %>% as.data.table

  eggerintercept <- data.frame(exposure = x@Exposure, outcome = x@Outcome, b = x@Intercept,
                               se = x@StdError.Int, lci = x@CILower.Int, uci = x@CIUpper.Int, pval = x@Pvalue.Int,
                               nsnp = x@SNPs, method = "Multivariable Egger Intercept") %>% as.data.table

  mvmr_res[[1]] <- MendelianRandomization::mr_mvmedian(object = mriobj)
  mvmr_res[[2]] <- tryCatch(
    expr = {MendelianRandomization::mr_mvlasso(object = mriobj)   },
    error = function(e){return(NULL) })

  resmvmr  <- lapply(mvmr_res, function(x) {
    res_mvmr <- data.frame(exposure = x@Exposure,
                           outcome = x@Outcome, b = x@Estimate,
                           se = x@StdError, lci = x@CILower,
                           uci = x@CIUpper, pval = x@Pvalue,
                           nsnp = x@SNPs) %>% as.data.table
    return(res_mvmr)})


  methodvec <- c("Multivariable Median", "Multivariable Lasso")
  for(i in 1:length(resmvmr)){
    resmvmr[[i]][,method := methodvec[i]]
  }
  resmvmr <- rbindlist(resmvmr)

  resmvmr <-rbindlist(list(ivwestimate, resmvmr,eggerestimate, eggerintercept), fill = TRUE)
  resmvmr <- merge(resmvmr, sres2, by = "exposure", sort = FALSE)
  return(resmvmr)
}


#' Obtain conditional F statistcis using MVMR package from paper PMID: 34338327
#'
#' @param exposure_outcome_harmonized the output from TwoSampleMR::mv_harmonise_data
#' @param pheno_cov_exp phenotypic correlation between the exposures. It is a possbility to assume that these covariances are zero.
#'
#' @return A data.frame of F statistics
#' @export
obtain_conditionalFstatistics <- function(exposure_outcome_harmonized, pheno_cov_exp = 0) {
  mv_harm <- exposure_outcome_harmonized
  F.data <- MVMR::format_mvmr(BXGs = mv_harm$exposure_beta,
                              BYG = mv_harm$outcome_se,
                              seBXGs = mv_harm$exposure_se,
                              seBYG = mv_harm$outcome_se,
                              RSID = mv_harm$exposure_se %>% rownames) %>% setDT
  mvmrcovmatrix <- matrix(pheno_cov_exp, nrow = ncol(exposure_outcome_harmonized$exposure_beta), ncol = ncol(exposure_outcome_harmonized$exposure_beta)) #matrix of R
  diag(mvmrcovmatrix) <-  1
  Xcovmat<- MVMR::phenocov_mvmr(mvmrcovmatrix,F.data[, .SD,.SDcols = which(grepl("sebetaX", colnames(F.data)))])
  sres2 <- MVMR::strength_mvmr(r_input = F.data, gencov = Xcovmat)
  colnames(sres2) <- mv_harm$expname$exposure
  rownames(sres2) <- NULL
  sres2 <- data.frame(exposure = colnames(sres2), F_stastistics = as.numeric(sres2[1,]))
  return(sres2)
}
