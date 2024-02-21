#' A wrapper for contamination mixture sensitivity analysis
#'
#' A wrapper for contamination mixture sensitivity analysis that works with harmonised data with Two Sample MR.
#'
#' @param dat a harmonised data.frame same as TwoSampleMR.
#'
#' @return A contamination mixture model

mr_conmix_wrapper <- function (dat) {
  input <- MendelianRandomization::mr_input(bx = dat$beta.exposure,
                                            bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome,
                                            exposure = dat$exposure, outcome = dat$outcome, snps = dat$SNP,
                                            effect_allele = dat$effect_allele.exposure, other_allele = dat$other_allele.exposure,
                                            eaf = dat$eaf.exposure)
  conmixobj <- MendelianRandomization::mr_conmix(input)
  conmixres <- data.frame(id.exposure = unique(conmixobj@Exposure),
                          id.outcome = unique(conmixobj@Outcome),
                          outcome = unique(conmixobj@Outcome),
                          exposure = unique(conmixobj@Exposure), method = "Contamination mixture",
                          nsnp = conmixobj@SNPs, b = conmixobj@Estimate, se = NA,
                          pval = conmixobj@Pvalue, lci = min(conmixobj@CILower), uci = max(conmixobj@CIUpper))
  return(conmixres)
}

#' All sensitivity analysis
#'
#' Perform all sensitivity analysis and return in data.frame.
#'
#' @param dat An harmonised data.frame like in TwoSampleMR
#' @param SignifThreshold the significance threshold you want to use for MRpresso.
#' @param Primary_mr What type of primary_mr model should be used: "default", "random", "random_underdispersion" or "fixed".
#' The random-effects model ("random") is a multiplicative random-effects model, allowing overdispersion in the weighted linear regression (the residual standard error is not fixed to be 1, but is not allowed to take values below 1).
#' The fixed-effect model ("fixed") sets the residual standard error to be 1. The "default" setting is to use a fixed-effect model with 3 genetic variants or fewer, and otherwise to use a random-effects model.
#' The "random_underdispersion" is the default method in TwoSampleMR.
#' @param short if TRUE only runs egger, weighted mode median and egger if FALSE runs the rest
#' @param skip_presso default FALSE. will run MR presso. If True won't run MR presso. MR presso takes load of memory.
#' @param tsmr_robust will be passed on the TwoSampleMR::mr function
#' @return a data.frame with results of sensitivity analysis.
#'
#' @import data.table
#'
#' @export
all_mr_methods <- function(dat, SignifThreshold = 0.05, Primary_mr = "default", short = FALSE, skip_presso = FALSE,
                           tsmr_robust = c("mr_weighted_mode", "mr_weighted_median", "mr_egger_regression")) {
  dat<-dat[dat$mr_keep == TRUE,]
  if (nrow(dat) == 1) {
    res <- TwoSampleMR::mr(dat, method_list = c("mr_wald_ratio"))
    res$nsnp <- dat$SNP
    res$type_of_test <- "Primary analysis"
    setDT(res)
    res$lci<-(res$b-res$se)*1.96
    res$uci <- (res$b+res$se)*1.96
    return(res)
  }else{
    dat <- data.table::as.data.table(dat)
    if(Primary_mr == "random"){ primary_analysis <- "mr_ivw_mre"}
    if(Primary_mr == "fixed"){ primary_analysis <- "mr_ivw_fe"}
    if(Primary_mr == "random_underdispersion"){ primary_analysis <- "mr_ivw"}

    if(Primary_mr == "default") {
      if(nrow(dat) > 2){
        primary_analysis <- "mr_ivw"
      } else {
        primary_analysis <- "mr_ivw_fe"
      }
    }
    dat<- dat[dat$mr_keep,]
    res <- TwoSampleMR::mr(dat, method_list=c(primary_analysis, tsmr_robust))

    if(short == TRUE) {return(res)}
    #mr raps
    out <-  mr.raps::mr.raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure,
                             se_out = dat$se.outcome, over.dispersion = TRUE, loss.function = "tukey")

    out_raps <- data.frame(id.exposure = dat$id.exposure[1], id.outcome = dat$id.outcome[1], outcome = dat$outcome[1],
                           exposure = dat$exposure[1], method = "Robust adjusted profile score (RAPS)", nsnp = nrow(dat),
                           b = out$beta.hat, se = out$beta.se, pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2)
    res <- rbind(res, out_raps)
    #MR Presso

    if(nrow(dat) > 3 && !skip_presso){
      tryCatch({
        ok<-TwoSampleMR::run_mr_presso(as.data.frame(dat), NbDistribution = ifelse(nrow(dat)<900, 1000, nrow(dat)+100))


        if(!is.na(ok[[1]]$`Main MR results`[2, "Sd"])) {
          nsnpremoved <-length(ok[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
          numsnp <- ok[[1]]$`MR-PRESSO results`$`Outlier Test` %>%
            nrow(.)
          index <- 2
        } else{          index <- 1}

        globaltest_p <- ok[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
        if( globaltest_p != "<0.001" | globaltest_p > SignifThreshold) {
          print("MR-PRESSO global test indicates no significant pleiotropy")
        }

        dtf<- ok[[1]]$`Main MR results`[index,]
        data.table::setDT(dtf)
        data.table::setnames(dtf, c("MR Analysis", "Causal Estimate", "Sd", "P-value"), c("Analysis", "b", "se", "pval" ))
        dtf$id.exposure <- dat$id.exposure[1]
        dtf$id.outcome <- dat$id.outcome[1]
        dtf$outcome <- dat$outcome[1]
        dtf$exposure <- dat$exposure[1]
        dtf$method <- "MR-PRESSO (outlier-corrected)"

        if(dtf$Analysis == "Outlier-corrected") {
          dtf$nsnp <- numsnp - nsnpremoved
        } else {
          dtf$nsnp <- res[res$method == "Inverse variance weighted","nsnp"]
        }

        col <- colnames(res)
        dtf <- dtf[, ..col]
      },
      error=function(cond) {
        message(conditionMessage(cond))
      }
      )
      if(!exists("dtf")){
        dtf <- matrix(nrow=0,ncol=ncol(res)) %>% as.data.frame(.) %>% data.table::as.data.table(.)
        colnames(dtf)<-colnames(res)
      }
    } else{
      dtf <- matrix(nrow=0,ncol=ncol(res)) %>% as.data.frame(.) %>% data.table::as.data.table(.)
      colnames(dtf)<-colnames(res)
    }

    dt <- rbind(res, dtf[,c("id.exposure", "id.outcome", "outcome",  "exposure", "method", "nsnp", "b", "se", "pval")])
    data.table::setDT(dt)
    ####MR-radial
    if(!any(grepl("presso", tolower(dt$method)))){
      tryCatch({
        res <- TwoSampleMR::mr(dat, method_list = c("mr_ivw_radial"))
        dt <- rbindlist(list(dt, res), fill = TRUE)
      }, error = function(cond) {
        message(conditionMessage(cond))
      })
    }
    dt$lci <- dt$b-(dt$se*1.96)
    dt$uci <- dt$b+(dt$se*1.96)
    #contamination mixture
    tryCatch({ #this condition is not very good, but I do not have anything better for the moment
      conmixres <- GagnonMR:::mr_conmix_wrapper(dat)
      conmixres[, "id.exposure"] <- dt[1,"id.exposure"]
      conmixres[, "id.outcome"] <- dt[1,"id.outcome"]
    },
    error=function(cond) {
      message(conditionMessage(cond))
    }
    )
    if(!exists("conmixres")){
      conmixres <- matrix(nrow=0,ncol=ncol(dt)) %>% as.data.frame(.) %>% as.data.table(.)
      colnames(conmixres)<-colnames(dt)
    }
    #rbind everything
    result <- rbind(dt, conmixres)
    type_of_test <- result$method %>%
      ifelse(grepl("Inverse variance weighted", .), "Primary analysis", .) %>%
      ifelse(. %in% c("Weighted mode", "Weighted median"), "Consensus methods", .) %>%
      ifelse(. %in% c("Robust adjusted profile score (RAPS)", "MR Egger", "Contamination mixture"), "Pleiotropy test", .) %>%
      ifelse(grepl("MR-P|radial", .), "Outlier-robust methods", .)
    result$type_of_test <- type_of_test
    return(result)
  }
}


#' Perform a wald ratio
#'
#' @param inst_all_sign_clump_arg the instrument object
#' @param all_outcome_arg the outcome object
#' @param id_exposure_name the id.exposure of the exposure
#' @param id_outcome_name the id.outcome of the outcome
#'
#' @return
#' @export

wald_quick <- function(inst_all_sign_clump_arg, all_outcome_arg, id_exposure_name, id_outcome_name = ".*") {

  message(paste0("###### Analysing ", id_exposure_name, " #########"))

  harm_univariate <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure),],
                                                 outcome_dat = all_outcome_arg[ grep(paste0("^" ,id_outcome_name, "$"), id.outcome),],
                                                 action = 1)
  harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate)
  data.table::setDT(harm_univariate)
  harm_univariate <- harm_univariate[mr_keep==TRUE,]
  if (dim(harm_univariate)[1] == 0) {
    return(cbind(inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure), .(id.exposure,exposure)],
                 dplyr::distinct(all_outcome_arg[grep(paste0("^" ,id_outcome_name, "$"), id.outcome), .(id.outcome,outcome)])))
  }
  res <- TwoSampleMR::mr(harm_univariate, method_list = c("mr_wald_ratio"))
  res$nsnp <- harm_univariate$SNP
  res$type_of_test <- "Primary analysis"
  data.table::setDT(res)
  res$lci<-(res$b)-(res$se*1.96)
  res$uci <- (res$b)+(res$se*1.96)
  harm_univariate <- TwoSampleMR::add_rsq(harm_univariate) %>% data.table::as.data.table(.)
  res$fstat<-GagnonMR::fstat_fromdat(split(harm_univariate, 1:nrow(harm_univariate)))
  res$rsq <- harm_univariate$rsq.exposure
  res <- merge(res, harm_univariate[,.(id.exposure,id.outcome,steiger_dir,steiger_pval)], by = c("id.exposure", "id.outcome"))
  data.table::setDT(res)
  return(res)
}


#' Perform a wald ratio analyses and format the data
#'
#' @param inst_all_sign_clump_arg The instrument data.table
#' @param all_outcome_arg The all_outcome data.table
#' @param id_exposure_name A regular expression to select the id.exposure
#' @param id_outcome_name A regular expression to select the id.outcome
#'
#' @return
#' @export

wald_quick <- function(inst_all_sign_clump_arg, all_outcome_arg, id_exposure_name, id_outcome_name = ".*") {

  message(paste0("###### Analysing ", id_exposure_name, " #########"))
  stopifnot(inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure),.N]==1)
  harm_univariate <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure),],
                                                 outcome_dat = all_outcome_arg[ grep(paste0("^" ,id_outcome_name, "$"), id.outcome),],
                                                 action = 1)

  if (dim(harm_univariate)[1] == 0) {
    return(cbind(inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure), .(id.exposure,exposure)],
                 dplyr::distinct(all_outcome_arg[grep(paste0("^" ,id_outcome_name, "$"), id.outcome), .(id.outcome,outcome)])))
  }

  harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate)
  data.table::setDT(harm_univariate)
  harm_univariate <- harm_univariate[mr_keep==TRUE,]
  res <- TwoSampleMR::mr(harm_univariate, method_list = c("mr_wald_ratio"))
  res$nsnp <- harm_univariate$SNP
  res$type_of_test <- "Primary analysis"
  data.table::setDT(res)
  res$lci<-(res$b)-(res$se*1.96)
  res$uci <- (res$b)+(res$se*1.96)
  harm_univariate <- TwoSampleMR::add_rsq(harm_univariate) %>% data.table::as.data.table(.)
  res$fstat<-GagnonMR::fstat_fromdat(split(harm_univariate, 1:nrow(harm_univariate)))
  res$rsq <- harm_univariate$rsq.exposure
  res <- merge(res, harm_univariate[,.(id.exposure,id.outcome,steiger_dir,steiger_pval)], by = c("id.exposure", "id.outcome"))
  data.table::setDT(res)
  return(res)
}
