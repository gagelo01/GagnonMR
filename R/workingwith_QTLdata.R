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
      break()
    }
    ##harmonise same order as ldmat
    list_harmonised <- TwoSampleMR:::harmonise_ld_dat(x = harm_all, ld = ldmat)
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
#' @param gencode the path to the gencode file
#' @param vcf_file the path to the vcf file
#' @param temp_dir the path to the temporary directory whare all files will at the end of the process get erased.
#'
#' @return a data.frame
#' @export
#'
#' @examples

get_eQTL<- function(tissue,
                    gene ,
                    mywindow = 1000000,
                    mythreshold = 1,
                    gencode = as.data.frame(fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")),
                    vcf_file = "/mnt/sde/pernic01/eQTL/GTEx8/genotypes/geno_693Indiv_maf_0.01_chr1_22.vcf.gz",
                    temp_dir = "/mnt/sde/gagelo01/temp/" ) {
  expr_mat_file = paste0("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/GTEx_EUR_Analysis_v8_eQTL_expression_matrices/",
                         tissue,".v8.normalized_expression_EUR_chr1_22.bed")
  cov_file = paste0("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/Covariates_GTEx8_EUR/",
                    tissue,"_covariates_GTEx8_EUR.txt")


  gencode_small = gencode[which(gencode$gene_name %in% gene),]
  if(nrow(gencode_small) == 0) {
    return(data.table(tissue = tissue, probe= gencode_small[1,"gene_id"]))
  } else {
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
    fwrite(matr,paste0(temp_dir,"mymatr.bed"),sep="\t",quote=F,row.names=F)

    # indexing new file (needed for QTLTools)

    system2(command = "bgzip",
            args = c("-f",
                     paste0(temp_dir, "mymatr.bed"),
                     "&&",
                     "tabix",
                     "-p",
                     "bed",
                     paste0(temp_dir, "mymatr.bed.gz")))

    mymatr =  paste0(temp_dir, "mymatr.bed.gz")

    ### traitement du fichier covar

    covar = as.data.frame(fread(cov_file))
    colnames(covar)[grep("GTEX-", colnames(covar))] = paste(colnames(covar)[grep("GTEX-", colnames(covar))],
                                                            colnames(covar)[grep("GTEX-", colnames(covar))], sep = "_")
    colnames(covar)[colnames(covar) == "ID"] = "id"
    covar$id[grep("InferredCov", covar$id)] = gsub(covar$id[grep("InferredCov", covar$id)], pattern = " ", replacement = "_")
    fwrite(covar,paste0(temp_dir, "mycovar.txt"),sep="\t",quote=F,row.names=F)

    mycovar = paste0(temp_dir, "mycovar.txt")


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
                     paste0(temp_dir, tissue,"_", gene,"_QTLtools_nominal_",mythreshold,"_window_",mywindow / 1000,"kb",".txt")))

    path_qtl <- paste0(temp_dir, tissue,"_", gene,"_QTLtools_nominal_",mythreshold,"_window_",mywindow / 1000,"kb",".txt")
    if(!file.exists(path_qtl)) { return(data.table(tissue = tissue, probe= gencode_small[1,"gene_id"]))}
    QTL = fread(path_qtl)

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

    file.remove(mymatr)
    file.remove(paste0(mymatr,".tbi"))
    file.remove(mycovar)
    file.remove(paste0(temp_dir, tissue,"_", gene,"_QTLtools_nominal_",mythreshold,"_window_",mywindow / 1000,"kb",".txt"))

    return(QTL)
  }
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
                             translation = fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt"),
                             gencode = fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")) {
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
#' @param exp the name of one exposure in the inst_pQTL object to analyse.
#' @param inst_pQTL the inst_pqtl with all exposures
#' @param outcome_id if "NA" then will use the argument name outcome to determine the outcome object
#' @param name_outcome the name as character of the outcome object
#' @param outcome_column The name that will go in the column "outcome". If outcome_id == "NA" then the name is used for the name of the outcome object
#'
#'
#' @return the posterior probability of colocalisation
#' @export
format_to_coloc <- function(exp, inst_pQTL, outcome_id, name_outcome, outcome_column = "NA") {

  inst_pqtl <- get(inst_pQTL)
  dt_exposure<-inst_pqtl[exposure == exp,]

  #outcome
  if(outcome_id == "NA"){
    dt_outcome <- get(name_outcome)
    dt_outcome <- dt_outcome[SNP %in% dt_exposure$SNP,]
    dt_outcome <- dt_outcome[outcome %in% outcome_column, ]
    name_outcome<-dt_outcome[1,]$outcome
  } else {
    dt_outcome<- TwoSampleMR::extract_outcome_data(outcomes = outcome_id, snps = dt_exposure$SNP, proxies = FALSE)
  }

  dat <- harmonise_data(dt_exposure, dt_outcome, action=1)
  setDT(dat)
  mypossnp <- dat[exposure == exp]$pos.exposure

  if(length(mypossnp)== 0) {
    if(dat[, max(pos.exposure) - min(pos.exposure)] > 2*10^6) {
      return("not sure how to select the region")
    }

    dat_window <- dat
  } else {
    window<-1*10^6
    dat_window <- dat[(pos.exposure >= mypossnp-window) & (pos.exposure <= mypossnp+window),]
  }

  toselect <- dat_window[,complete.cases(.SD), .SDcols = c("beta.exposure", "se.exposure", "samplesize.exposure", "beta.outcome", "se.outcome", "samplesize.outcome")]
  dat_window<-dat_window[toselect,]
  dat_window <- distinct(dat_window)
  if(dt_outcome[1,]$units.outcome != "log odds") {
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



#' Title
#'
#' @param SOMAmerID the SOMAmerID
#' @param genecard_name the name of the Gene. exemple GSTA1
#' @param window The window. note that 2e5 window means the function will return the SNP 1e5 KB upstream and 1e5KB downstream of the gene coding region
#' @param minimum_maf the minimum maf SNP you want to include (default 0.01)
#' @param gencode_chr the name in character of the gencode. load gencode prior to running the function gencode <- data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
#' @param traduction_chr the name in character of the traduction data.table. load traduction prior to running the function traduction <-  fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
#' @param ensembl to convert to genecard name
#' @param type if cis selection only cis
#' @param path_sun the path towards the sun data

#'
#' @return
#' @export
#'
#' @examples

extract_all_cis_instrument_sun <- function(genecard_name = NULL, SOMAmerID = NULL, window = 1e6,
                                           minimum_maf =0.001, gencode_chr="gencode", traduction_chr="traduction",
                                           ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl"),
                                           type = "cis", path_sun = "/home/couchr02/Mendel_Commun/dbs_Web/meta_filtered_final") {

  conversion <- readxl::read_excel("/home/gagelo01/workspace/Projects/small_MR_exploration/Plasma_proteome/Data/Raw/media-1.xlsx", sheet = 2, skip = 2) %>% as.data.table(.)

  if(!is.null(genecard_name)){
    SOMAmerID <- conversion[EntrezGeneSymbol == genecard_name, ID]
  }

  if(!is.null(SOMAmerID)){
    genecard_name<-conversion[ID == SOMAmerID, EntrezGeneSymbol]
    genecard_name <- sub(".", "_", genecard_name, fixed = TRUE) %>% strsplit(., "_") %>% unlist
  }
  if(is.character(gencode_chr)) {gencode  <- get(gencode_chr)}  else {gencode <- gencode_chr}
  if(is.character(traduction_chr))  {traduction <- get(traduction_chr)} else {traduction<-traduction_chr}
  gene_coords <-gencode[gene_name %in% genecard_name,]
  if(nrow(gene_coords) == 0) {
    gene_coords <- biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"),
                                  filters = "hgnc_symbol", values = genecard_name,
                                  mart = ensembl) %>% as.data.table(.)
    gene_coords <- gene_coords[chromosome_name %in% 1:22,]
    colnames(gene_coords)<-c("gene_name", "chr", "start", "end")
  } else {
    gene_coords[, start := min(start)]
    gene_coords[, end := max(end)]
    gene_coords[,gene_name := paste0(unique(gene_name, collapse = "-"))]
    gene_coords <- gene_coords %>% distinct(.)}

  name_dir_prot <- SOMAmerID

  sun <- fread(paste0(path_sun, "/",name_dir_prot, "/", name_dir_prot, "_chrom_", gene_coords$chr[1], "_meta_final_v1.tsv"))
  setnames(sun, "log(P)", "log_p")
  sun[, p:=exp(log_p)]

  if(type == "pan") {
    sun_list<-vector(mode="list", length=22)
    for(i in 1:22) {
      sun_list[[i]]<-fread(paste0(path_sun, "/",name_dir_prot, "/", name_dir_prot, "_chrom_", i, "_meta_final_v1.tsv"))
    }
    sun_full<-rbindlist(sun_list)
    setnames(sun_full, "log(P)", "log_p")
    sun_full[, p:=exp(log_p)]
    sun_trans <-  sun_full[chromosome != gene_coords$chr[1] & p<5e-8,]
    sun <- rbind(sun, sun_trans)
  }


  #getrsisd
  #in which imputed variants had the same genomic position (GRCh37) and alleles
  traduction[,chr:=as.integer(chr)]
  traduction[,pos:=as.integer(position)]
  sun <- merge(sun, traduction, by.x = c("chromosome", "position"), by.y = c("chr", "position"))
  sun[,a0 := tolower(a0)]
  sun[,a1 := tolower(a1)]
  #eaf
  sun<-sun[Allele1 == a0 | Allele1 == a1,]
  sun[, eaf := ifelse(Allele1 == a0, EUR, 1-EUR)]

  #transform  cis-trans
  gene_coords[, start_cis := (start - window/2) %>% ifelse(.<1, 1, .) ]
  gene_coords[, end_cis := end + window/2]
  sun <- cbind(sun, gene_coords[,.(start_cis, end_cis)])
  sun[, pqtl_type := ifelse( (position > start_cis) & (position < end_cis), "cis", "trans")]
  sun <- sun[, -c("a0", "a1", "EUR" ,"maf", "start_cis", "end_cis" )]
  sun[,Phenotype := paste0("sun_",   gene_coords[1,gene_name])]
  sun<-distinct(sun)
  sun[, maf := ifelse(eaf > 0.5, 1-eaf, eaf)]
  sun[maf> minimum_maf,]
  sun[rsid %in% sun[, .N , by ="rsid"][N > 1,]$rsid, ] #many snp with different maf
  sun[, N := 3301]
  if(type=="cis") {sun <- sun[pqtl_type == "cis",]}
  if(type=="pan") {sun <- sun[pqtl_type == "cis" | (pqtl_type == "trans" & p<5e-8)]}
  sun_format <-TwoSampleMR::format_data(sun,
                                        type="exposure",
                                        snp_col="rsid",
                                        phenotype_col = "Phenotype",
                                        beta_col="Effect",
                                        effect_allele_col="Allele1",
                                        other_allele_col="Allele2",
                                        se_col="StdErr",
                                        pval_col="p",
                                        chr_col = "chromosome",
                                        pos_col = "position",
                                        samplesize_col = "N")
  sun_format <- merge(sun_format, sun[, .(rsid, pqtl_type)], by.x = "SNP", by.y = "rsid", all.y = FALSE)
  return(sun_format)
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
