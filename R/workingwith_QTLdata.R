#' Calculate multicismr and correcting for ld matrix
#'
#' @param inst_pqtl one exposure clumped at r2<0.6.
#' @param outcome one outcome
#' @param ldmat the ldmatrix obtained by TwoSampleMR::ld_matrix()
#'
#' @return A data.frame object with the results of the multicis mendelian randomization investigation corrected for ld.
#' @export
calculate_multicis_mr <- function(inst_pqtl, outcome, ldmat ) {

  res_multicis <- data.frame(exposure = inst_pqtl[1,]$exposure,
                             outcome = outcome$outcome[1],
                             beta.multi_cis = NA,
                             se.multi_cis = NA,
                             pval.multi_cis = NA,
                             cochranQ.multi_cis = NA,
                             cochranQpval.multi_cis = NA, # rejection of the null is an indication that one or more variants may be pleiotropic.
                             nsnp.multi_cis = 1)
  ###harmonise
  if(nrow(inst_pqtl[SNP %in% outcome$SNP,]) == 0) {
    return(res_multicis)
  } else {
    harm_all <- harmonise_data(inst_pqtl, outcome, action = 2)
    setDT(harm_all)
    goodsnporder<-inst_pqtl[SNP %in% harm_all$SNP]$SNP
    harm_all<-harm_all[order(match(SNP, goodsnporder))]
    if(nrow(harm_all) == 1) {
      return(res_multicis)
    }
    ##harmonise same order as ldmat
    list_harmonised <- TwoSampleMR::harmonise_ld_dat(x = harm_all, ld = ldmat)
    harm_all <- list_harmonised$x
    ldmat <- list_harmonised$ld

    if(nrow(harm_all) == 0) {
      return(res_multicis)
    }
    if(!is.matrix(ldmat)) {
      return(res_multicis)
    } else {

      ####wrapper harm_all to mr_input
      MRInputobject.cor <- mr_input(bx = harm_all$beta.exposure,
                                    bxse = harm_all$se.exposure,
                                    by = harm_all$beta.outcome,
                                    byse = harm_all$se.outcome,
                                    corr = ldmat)

      IVWobject <- mr_ivw(MRInputobject.cor,
                          model = "default",
                          correl = TRUE,
                          distribution = "normal",
                          alpha = 0.05)


      res_multicis <- data.frame(exposure = inst_pqtl[1,]$exposure,
                                 outcome =  outcome$outcome[1],
                                 beta.multi_cis = IVWobject@Estimate,
                                 se.multi_cis = IVWobject@StdError,
                                 pval.multi_cis = IVWobject@Pvalue,
                                 cochranQ.multi_cis = IVWobject@Heter.Stat[1],
                                 cochranQpval.multi_cis = IVWobject@Heter.Stat[2], # rejection of the null is an indication that one or more variants may be pleiotropic.
                                 nsnp.multi_cis = nrow(harm_all))
      return(res_multicis)
    }
  }
}




#' Get eqtl data from GTEX
#'
#' @param tissue From which tissue do you want to extract the data
#' @param gene From which gene do you want to extract the data exemple SULT1A1
#' @param mywindow Which window to use
#' @param mythreshold Which p value threshold to use
#' @param gencode the GTEX v8 dictionnary of ensembl gene used.
#' @param vcf_file the path to the vcf file
#' @param expr_mat_file expression matrices. column are subjects, and row is their expression for each gene
#' @param cov_file the covariate file
#'
#' @return a data.frame
#' @export

get_eQTL<- function(tissue,
                    gene ,
                    mywindow = 1e6,
                    mythreshold = 1,
                    gencode = fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt"),
                    vcf_file = "/mnt/sda/pernic01/eQTL/GTEx8/genotypes/geno_693Indiv_maf_0.01_chr1_22.vcf.gz",
                    expr_mat_file = paste0("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/GTEx_EUR_Analysis_v8_eQTL_expression_matrices/",
                                           tissue,".v8.normalized_expression_EUR_chr1_22.bed"), #expression matrices. column are subjects, and row is their expression for each gene
                    cov_file = paste0("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/Covariates_GTEx8_EUR/",
                                      tissue,"_covariates_GTEx8_EUR.txt")) {



  gencode_small <- gencode[gene_name == gene, .(chr, start = min(start), end =max(end), gene_id, gene_name) ] %>% dplyr::distinct(.)
  gencode_small <- gencode_small[chr %in% 1:22, ]

  if(nrow(gencode_small)==0){
    return(data.table(tissue = tissue,
                      hgncgene = gene))
  }

  gene_id = unique(gencode_small$gene_id)
  chr = unique(gencode_small$chr)
  mystart = min(gencode_small$start)
  myend = max(gencode_small$end)

  ### treatment of the phenos file

  matr = as.data.frame(fread(expr_mat_file))
  num_suj = length(colnames(matr)[grep("GTEX", colnames(matr))])
  matr$"#chr" = sub("chr","",matr$"#chr")
  matr = cbind(matr[,c("#chr", "start", "end", "gene_id")], data.frame(gid = rep(".", dim(matr)[1]),
                                                                       strand = rep("+", dim(matr)[1])), matr[,grep("GTEX-", colnames(matr))])
  colnames(matr)[grep("GTEX-", colnames(matr))] = paste(colnames(matr)[grep("GTEX-", colnames(matr))],
                                                        colnames(matr)[grep("GTEX-", colnames(matr))], sep = "_")

  dir.create(tempdir())
  temp_file <- tempfile()
  fwrite(matr,paste0(temp_file,".bed"),sep="\t",quote=F,row.names=F)

  # indexing new file (needed for QTLTools)

  system2(command = "bgzip",
          args = c("-f",
                   paste0(temp_file,".bed"),
                   "&&",
                   "tabix",
                   "-p",
                   "bed",
                   paste0(temp_file,".bed.gz")))

  mymatr <-  paste0(temp_file,".bed.gz")

  ### traitement du fichier covar

  covar <- as.data.frame(fread(cov_file))
  colnames(covar)[grep("GTEX-", colnames(covar))] = paste(colnames(covar)[grep("GTEX-", colnames(covar))],
                                                          colnames(covar)[grep("GTEX-", colnames(covar))], sep = "_")
  colnames(covar)[colnames(covar) == "ID"] = "id"
  covar$id[grep("InferredCov", covar$id)] = gsub(covar$id[grep("InferredCov", covar$id)], pattern = " ", replacement = "_")

  fwrite(covar,paste0(temp_file, ".txt"),sep="\t",quote=F,row.names=F)

  mycovar = paste0(temp_file, ".txt")


  ###  QTLtools analysis

  system2(command = "QTLtools",
          args = c("cis",
                   "--vcf",
                   vcf_file,
                   "--bed",
                   mymatr,
                   "--cov",
                   mycovar,
                   "--nominal",
                   mythreshold,
                   "--region",
                   paste0(chr, ":", mystart, "-", myend),
                   "--window",
                   format(mywindow/2, scientific = F),
                   "--out",
                   paste0(temp_file, "out.txt")))

  path_qtl <- paste0(temp_file, "out.txt")
  if(!file.exists(path_qtl)) { return(data.table(tissue = tissue,
                                                 hgncgene = gene))}
  QTL <- fread(path_qtl)
  QTL <- QTL[sub("\\..*", "", V1) %in% gencode_small[, sub("\\..*", "", unique(gene_id))]] #make sure that even from ensembl it fits

  col_QTLtools_nominal = c(
    "probe",
    "chr",
    "start",
    "end",
    "strand",
    "n_cis_variant_tested",
    "distance_probe_variant",
    "variant",
    "chr_top_variant",
    "start_top_variant",
    "end_top_variant",
    "p_val",
    "beta",
    "top_variant"
  )

  colnames(QTL) = col_QTLtools_nominal

  covar = covar[,-which(colnames(covar) == "id")]
  QTL[, zscore := abs(qnorm(p_val/2))*(beta/abs(beta))]

  QTL[, se := beta / zscore]

  QTL[, tissue := tissue]
  QTL[,sample_size := num_suj]

  QTL <- QTL[!is.na(QTL$p_val),]

  unlink(mymatr)
  unlink(paste0(mymatr,".tbi"))
  unlink(mycovar)
  unlink(path_qtl)
  unlink(temp_file)
  unlink(tempdir(), recursive = TRUE)

  return(QTL)
}

#' to run after extracting data from gtex
#'
#' @param exposures the GTEx data.frame
#' @param translation the translation data.frame to get rsids.
#' @param gencode a data.frame to transform Ensemble into genecard name
#'
#' @return
#' @export
format_gtex_data <- function(exposures,
                             translation = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt"),
                             gencode = fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")) {
  setnames(gencode, "gene_id", "ensembl_gene_id", skip_absent = TRUE)
  data_eQTL <- separate(exposures, col = "variant", into = c("chr_todump", "pos_todump", "Allele1", "Allele2", "buildtodump"))
  data_eQTL <- merge(data_eQTL, translation, by.x = c("chr_top_variant", "start_top_variant"),
                     by.y = c("chr", "pos_b38"))

  #eaf
  data_eQTL <- data_eQTL[(Allele1 == a0 | Allele1 == a1) & (Allele2 == a0 | Allele2 == a1) & a0 != a1 & Allele1 != Allele2, ] #because low number removed, coded on the forward strand
  data_eQTL <- data_eQTL[chr_top_variant %in% 1:22, ]
  data_eQTL[, eaf := ifelse(Allele1 == a0, EUR, 1-EUR)]

  data_eQTL[, `:=`(probe, gsub(".", "_", probe, fixed = TRUE))]
  data_eQTL <- separate(data_eQTL, col = "probe", into = c("ensembl_gene_id", "to_dump"))
  data_eQTL[, to_dump := NULL]

  gencode[, `:=`(probe, gsub(".", "_", ensembl_gene_id, fixed = TRUE))]
  gencode <- separate(gencode, col = "ensembl_gene_id", into = c("ensembl_gene_id", "to_dump"))
  gencode[,to_dump := NULL]

  data_eQTL <- merge(data_eQTL, distinct(gencode[,.(ensembl_gene_id, gene_name)]), by = "ensembl_gene_id")

  data_eQTL[,Phenotype := paste0(tissue, "-", gene_name)]
  data_eQTL[,id := Phenotype]
  data_eQTL <- distinct(data_eQTL)
  data_eqtl <- TwoSampleMR::format_data(data_eQTL, type = "exposure",
                                        phenotype_col = "Phenotype", snp_col = "rsid", beta_col = "beta",
                                        se_col = "se", effect_allele_col = "Allele1", id_col = "id",
                                        other_allele_col = "Allele2", samplesize_col = "sample_size",
                                        pval_col = "p_val", gene_col = "ensembl_gene_id", chr_col = "chr",
                                        pos_col = "start_top_variant")
  setDT(data_eqtl)

  #return
  return(data_eqtl)
}



#' Format from TwoSampleMR to coloc
#'
#' @param dat harmonised dataset with one single exposure and outcome.
#'
#' @return the posterior probability of colocalisation
#' @export
format_to_coloc <- function(dat) {

  setDT(dat)
  mypossnp <- dat$pos.exposure
  dat_window <- dat
  toselect <- dat_window[,complete.cases(.SD), .SDcols = c("beta.exposure", "se.exposure", "samplesize.exposure", "beta.outcome", "se.outcome", "samplesize.outcome")]
  dat_window<-dat_window[toselect,]
  dat_window <- distinct(dat_window)
  if(dat_window[1,]$units.outcome != "log odds") {
    type <- "quant"
    s<-0.99999
  } else {
    type <- "cc"
    s<-dat_window[1, ncase.outcome / ncontrol.outcome]
  }


  if(!any(is.na(dat_window$eaf.exposure))) {
    dat_window <- dat_window[( eaf.exposure > 0) & eaf.outcome > 0, ]
    maf_exposure <- dat_window$eaf.exposure %>% ifelse(. < 0.5, ., 1-.)
  } else {
    dat_window <- dat_window[(maf.exposure > 0 ) & eaf.outcome > 0, ]
    maf_exposure <- dat_window$maf.exposure

  }


  D1<- list(exposure = dat_window$exposure[1],
            beta=dat_window$beta.exposure,
            varbeta = dat_window$se.exposure^2,
            N = dat_window$samplesize.exposure,
            MAF = maf_exposure,
            type="quant",
            snp = dat_window$SNP,
            position = dat_window$pos.exposure)
  D2 <- list(outcome = dat_window$outcome[1],
             beta=dat_window$beta.outcome,
             varbeta=dat_window$se.outcome^2,
             N = dat_window$samplesize.outcome,
             type=type,
             s = s,
             MAF = dat_window$eaf.outcome %>% ifelse(. < 0.5, ., 1-.),
             snp = dat_window$SNP,
             position = dat_window$pos.exposure)

  if(D2$type == "quant") {
    D2$list <- NULL
  }
  return(list(D1 = D1, D2 = D2))

}


#' Calculate coloc
#'
#' @param list_dataset the list of dataset
#' @param exp the name of the exposure in character
#' @param name_outcome The name that will go in the column "outcome". If outcome_id == "NA" then the name is used for the name of the outcome object
#' @return the posterior probability of colocalisation
#' @export
#'
calculate_coloc_from_list_dataset <- function(list_dataset, exp, name_outcome) {

  res_coloc <- coloc::coloc.abf(dataset1=list_dataset$D1,
                                dataset2=list_dataset$D2)$summary[6]

  df_res <-data.frame(exposure_outcome = paste0(list_dataset$D1$exposure, "_", list_dataset$D2$outcome),
                      posprob_coloc.mr = res_coloc)

}

#' a wrapper GagnonMR::calculate_multicis_mr for pmap
#'
#' @param Exposure as.character the name of your exposure data.frame
#' @param Outcome as.character the name of your outcome data.frame object
#' @param exposure_name the name of the precise exposure you want to do multi_cis mr on
#' @param outcome_name the name of the precise outcome you want to do multi_cis mr o
#'
#' @return a multicis mt results data.frame
#' @export

multicis_mr_forpmap <- function(Exposure, Outcome, exposure_name, outcome_name) {
  exp <- get(Exposure)
  out <- get(Outcome)
  inst_pqtl <- exp[exposure == exposure_name, ][order(SNP)]
  out <- out[outcome == outcome_name, ][SNP %in% exp$SNP][order(SNP)]
  ldmat <- list_ldmat[[exposure_name]]
  snpnames <- do.call(rbind, strsplit(rownames(ldmat), split = "_")) %>% as.data.frame(.)
  colnames(snpnames)<-c("SNP", "A1", "A2")
  setDT(snpnames)
  goodsnporder<-snpnames[SNP %in% inst_pqtl$SNP]$SNP
  inst_pqtl<-inst_pqtl[order(match(SNP, goodsnporder))]
  out<-out[order(match(SNP, goodsnporder))]
  snpnames[A1 == TRUE, A1 := TRUE %>% ifelse(A2 == inst_pqtl$effect_allele.exposure, inst_pqtl$other_allele.exposure, .) %>%
             ifelse(A2 == inst_pqtl$other_allele.exposure, inst_pqtl$effect_allele.exposure, .)]
  snpnames[A2 == TRUE, A2 := TRUE %>% ifelse(A1 == inst_pqtl$effect_allele.exposure, inst_pqtl$other_allele.exposure, .) %>%
             ifelse(A1 == inst_pqtl$other_allele.exposure, inst_pqtl$effect_allele.exposure, .)]

  rownames(ldmat) <- snpnames[,paste0(SNP, "_", A1, "_", A2)]
  colnames(ldmat) <- snpnames[,paste0(SNP, "_", A1, "_", A2)]
  setDT(inst_pqtl)
  setDT(out)
  res <- calculate_multicis_mr(inst_pqtl = inst_pqtl,
                               outcome =   out,
                               ldmat = ldmat)
  return(res)
}


#' Title run hypr on aligned
#'
#' @param df_aligned a data.frame with everything
#' @param ... arguments passed to te function hyprcoloc::hyprcoloc
#'
#' @return
#' @export

run_hypr_on_aligned <- function(df_aligned, ...) {
  df_reshaped <- reshape(df_aligned, idvar = c("SNP", "chr.exposure", "pos.exposure"), timevar = "id.exposure", direction = "wide")
  index <- df_reshaped[, .SD , .SDcols = colnames(df_reshaped)[grepl("^beta.exposure|^se.exposure", colnames(df_reshaped))]] %>% complete.cases(.)
  df_reshaped <- df_reshaped[index,]
  if(nrow(df_reshaped) > 2) {
    effect_est <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("^beta.exposure", names(df_reshaped))]] %>% as.matrix
    effect_se <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("^se.exposure", names(df_reshaped))]] %>% as.matrix

    res<-hyprcoloc::hyprcoloc(effect.est = effect_est,
                              effect.se = effect_se,
                              trait.names = unique(df_aligned$exposure),
                              snp.id = df_reshaped$SNP,
                              ...)
    toreturn <- res$results %>% as.data.table(.)
    toreturn[, nsnp := nrow(df_reshaped)]
    toreturn[, allexposures := paste(df_aligned[order(id.exposure),]$id.exposure %>% unique , collapse = ":")]
    toreturn[,gene_investigated := df_aligned$gene.exposure[1]]
    return(toreturn)
  } else {
    return(data.frame(traits = "", posterior_prob = NA, regional_prob = NA, candidate_snp = NA, posterior_explained_by_snp = NA, nsnp = nrow(df_reshaped)))
  }
}
