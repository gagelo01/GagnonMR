#' harmonise your instrument and your outcome
#'
#'   inst_IBD <- fread("/mnt/sde/gagelo01/Projects/Dysbiose_project/Data/Modified/Instruments/inst_IBD")
#'   harmonize_gwassumstat_with_instrument( inst = "inst_IBD",
#'   action = 1,
#' gwas_file = "/home/couchr02/Mendel_Commun/Christian/GWAS/40_diseases/CAD_META.gz",
#' snp_col = "oldID",
#' outcome_name = "van_der_Harst_CAD",
#' beta_col = "Effect",
#' se_col = "StdErr",
#' pval_col = "P-value",
#' eaf_col = "Freq1",
#' effect_allele_col = "Allele1",
#' other_allele_col = "Allele2",
#' ncase_col = 122733,
#' ncontrol_col = 424528,
#' samplesize_col = 122733+424528,
#' snp_format = "rsids",
#' units_col = "log odds",
#' prevalence_col = 0.11)
#'
#' @param inst the gwas sumstat of your instrument in character. Your instrument must exist in the global environnemnt.
#' @param action same as harmonise_data
#' @param gwas_file the directory where the gwas sumstat of the outcome lay.
#' @param snp_col the name of the Snp column.
#' @param outcome_name the name you want to give your outcome
#' @param beta_col the name of the beta column
#' @param se_col the name of the se coumn
#' @param pval_col the name of the pval column
#' @param eaf_col the name of the eaf column
#' @param effect_allele_col the name of the efffect_allele column
#' @param other_allele_col the name of the other_allele_col
#' @param ncase_col how many case? numeric. if continuous supply "ncase"
#' @param ncontrol_col how many control? numeric. if continuous supply "ncontrol".
#' @param samplesize_col total samplesize (ncase+nconrol). numeric. if the data has a samplesize column. simply provide the name of that column.
#' @param snp_format can take two values either "rsids" or "chr:bp". if "chr:bp" will attempt to convert to rsids.
#' @param units_col if the outcome is continuous units are in standard deviation supply "SD". if the outcome is dichotomous supply "log odds". if the outcome is continuous and units are not is SD (years, Cm , Mm ,etc.) supply a character of your units.
#' @param prevalence_col if outcome is binary provide the prevalence of the condition in the general population. if it is continuous supply as.numeric("1").
#'
#' @return a data.frame of harmonised data.
#' @export

harmonize_gwassumstat_with_instrument<- function(inst,
                                                 action = 1,
                                                 gwas_file,
                                                 snp_col,
                                                 outcome_name,
                                                 beta_col,
                                                 se_col,
                                                 pval_col,
                                                 eaf_col,
                                                 effect_allele_col,
                                                 other_allele_col,
                                                 ncase_col = "ncase",
                                                 ncontrol_col = "ncontrol",
                                                 samplesize_col = "samplesize",
                                                 snp_format = "rsids",
                                                 units_col,  #"SD", "log odds", ou autre.
                                                 prevalence_col = 1#prevalence of the disease in the population. 1 if continuous.
) {

  instrument <- get(inst)
  df_no_suitable_instrument <- data.frame(exposure = unique(instrument$exposure),
                                          outcome = outcome_name,
                                          steiger_dir = TRUE,
                                          SNP = 0)

  if(gwas_file == "ieugwasr"){
    outcome <- extract_outcome_data(instrument$SNP,  outcomes = outcome_name )
    if(is.null(outcome)) {
      return(df_no_suitable_instrument)
    }
  }else {
    #create the temp files with all SNPs used as instrument
    if(snp_format == "rsids") {
      snps <- instrument[,"SNP"]
    }

    if(snp_format == "chr:bp") {
      dfsnp <- ieugwasr::variants_rsid(instrument$SNP)
      snps <-data.frame(snps = paste0(dfsnp$chr, ":", dfsnp$pos))
      if(!length(instrument$SNP)==length(snps)) {
        print("there was a problem while converting to chr:pos")
      }
    }

    colnames(snps) <- snp_col		# Car on veut que zgrep trouve tous les items de la liste DONT l'entête qui contient colname_snps...
    data.table::fwrite(snps, file = paste0("snps_list_",inst,"_", outcome_name,".txt"))		# writing de la liste en format .txt in the current directory

    all_out <-  data.table::fread(cmd = paste0("zgrep -w -f snps_list_",inst,"_", outcome_name,".txt ", gwas_file))

    #remove file
    if (file.exists(paste0("snps_list_",inst,"_", outcome_name,".txt"))) {
      file.remove(paste0("snps_list_",inst,"_", outcome_name,".txt"))
    }

    if(colnames(all_out)[1] == "V1" | nrow(all_out) == 0) { #no matching SNP
      return(df_no_suitable_instrument)
    }

    if(snp_format == "chr:bp") {
      data.table::setnames(all_out,snp_col, "SNP")
      all_out[,SNP := dfsnp$name]
      data.table::setnames(all_out,"SNP", snp_col)
    }

    all_out[,outcome := outcome_name]
    all_out[,units := units_col]

    if(varhandle::check.numeric(ncase_col)) {
      all_out[,ncase := as.numeric(ncase_col)]
      ncase_col <- "ncase"
    }

    if(varhandle::check.numeric(ncontrol_col)) {
      all_out[,ncontrol := as.numeric(ncontrol_col)]
      ncontrol_col <- "ncontrol"
    }

    if(varhandle::check.numeric(samplesize_col)) {
      all_out[,samplesize := as.numeric(samplesize_col)]
      samplesize_col <- "samplesize"
    }

    outcome <- format_data(all_out,
                           type="outcome",
                           snp_col= snp_col,
                           beta_col= beta_col,
                           se_col= se_col,
                           effect_allele_col=effect_allele_col,
                           other_allele_col=other_allele_col,
                           pval_col=pval_col,
                           eaf_col = eaf_col,
                           phenotype_col = "outcome",
                           ncase_col = ncase_col,
                           ncontrol_col = ncontrol_col,
                           samplesize_col = samplesize_col,
                           units_col = "units")

    #add prevalence
    if(as.numeric(prevalence_col) != 1) {
      outcome$prevalence.outcome <- as.numeric(prevalence_col)
    }

  }

  dat <- TwoSampleMR::harmonise_data(exposure = instrument, outcome = outcome, action=action) %>%
    TwoSampleMR::clump_data(., clump_r2 = 0.01)


  if(any(is.na(dat$eaf.outcome))){
    dat[is.na(dat$eaf.outcome),"eaf.outcome"]<-dat[is.na(dat$eaf.outcome),"eaf.exposure"]
  }

  if(any(is.na(dat$eaf.exposure))){
    dat[is.na(dat$eaf.exposure),"eaf.exposure"]<-dat[is.na(dat$eaf.exposure),"eaf.outcome"]
  }

  if(!(all(is.na(dat$eaf.outcome)) | all(is.na(dat$eaf.outcome)))){
    dat <- TwoSampleMR::steiger_filtering(dat)
    dat <- dat[dat$steiger_dir,]
  }

  dat<-dat[dat$mr_keep,]

  if(nrow(dat) == 0){
    return(df_no_suitable_instrument)
  }

  return(dat)
}
