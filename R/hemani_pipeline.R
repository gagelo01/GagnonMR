#' Title
#'
#' @param genecard_name A vector the names of the gene genecard name
#' @param window the window 1e6 == 5e5 upward 5e5 downward
#' @param gencode the gencode data.frame
#' @param ensembl ensembl
#'
#' @return
#' @export

from_genecard_to_generegion <- function(genecard_name, window = 1e6,
                                        gencode = data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt"),
                                        ensembl = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) {

  genecard_name <- genecard_name %>% strsplit(., split = "\\.|,") %>% unlist
  gene_coords <-gencode[gene_name %in% genecard_name,]

  if(!all(genecard_name %in% gene_coords$gene_name) | nrow(gene_coords) == 0) {
    gene_coords_mart <- biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"),
                                       filters = "hgnc_symbol", values = genecard_name[!(genecard_name %in% gene_coords$gene_name)],
                                       mart = ensembl) %>% as.data.table(.)
    gene_coords_mart <- gene_coords_mart[chromosome_name %in% 1:22,]
    colnames(gene_coords_mart)<-c("gene_name", "chr", "start", "end")
    gene_coords <- rbindlist(list(gene_coords, gene_coords_mart), fill = TRUE)
  }
  gene_coords[, start := min(start), by = "gene_name"]
  gene_coords[, end := max(end), by = "gene_name"]
  gene_coords[,gene_name := paste0(unique(gene_name), collapse = "-"), by = "gene_name"]
  gene_coords <- gene_coords %>% distinct(.)
  gene_coords <- gene_coords[chr %in% 1:22,]

  setDT(gene_coords)
  if(gene_coords[,.N]==0) {  return(NA)} else {
    gene_coords[, start_cis := (start - window/2) %>% ifelse(.<1, 1, .) ]
    gene_coords[, end_cis := end + window/2]
    gene_to_keep <- gene_coords[,length(unique(chr)), by = "gene_name"][V1 == 1]$gene_name
    gene_region <- gene_coords[gene_name %in% gene_to_keep, paste0(chr, ":", start_cis, "-", end_cis)]
    names(gene_region) <- gene_coords[gene_name %in% gene_to_keep, gene_name]
  }

  k <- genecard_name[!(genecard_name %in% gene_coords$gene_name)]
  vecna <- rep(NA, times = length(k))
  names(vecna)<-k
  gene_region <- c(gene_region, vecna)
  gene_region <- gene_region[genecard_name] #order
  return(gene_region)
}




#' This functions format to vcf and create an index. THIS FUNCTION CANNOT BE RUN IN PARALLEL
#'
#' @param all_out the data.frame object
#' @param snp_col the name of the snp col, NULL if absent
#' @param outcome_name the name of the traits
#' @param beta_col the name of the beta col, Must be present
#' @param se_col the name of the se col, Must be presnt
#' @param pval_col the name of the pval col, Must be present
#' @param eaf_col the name of the eaf col, NULL if absent. If present will flip eaf according to reference panel.
#' @param maf_col the name of the maf col, NULL if absent. If present, will use maf to infer eaf.
#' @param effect_allele_col the name of the effect allele, must be present
#' @param other_allele_col the name of the other allele, must be present
#' @param ncase_col The number of ncase (in numeric) or the name of the column where the info can be retrieved (in character) otherwise NULL
#' @param ncontrol_col  The number of ncontrol (in numeric) or the name of the column where the info can be retrieved (in character) otherwise NULL
#' @param samplesize_col The samplesize (in numeric) or the name of the column where the info can be retrieved (in character)
#' @param chr_col the name of the chr col, NULL if absent
#' @param pos_col the name of the pos col, NULL if absent
#' @param units ideally "SD", "log odds". If this is not the case, can put any character such as "year" or "Cm"
#' @param traduction the traduction file to map to the reference panel. IF NULL, all harmonisation steps will be skipped. This take less time, and must be used only on already harmonised data.
#' @param out_wd where you want your vcf.gz and vcf.gz.tbi to be stored
#' @param df_index the data.frame with all id.
#' @param group_name will be paste in df_index
#' @param year will be paste in df_index, must be numeric
#' @param author will be paste in df_index
#' @param consortium will be paste in df_index
#' @param sex Either "Males and Females", "Males", "Females"
#' @param population either "European", "African", "South Asian", "East Asian", "Mix"
#' @param initial_build the initial build ideally "HG19/GRCh37"
#' @param category either protein", "eqtl", "Metabolites","Trait","Disease"
#' @param pmid the pmid as.numeric
#' @param note any aditional note in character.
#' @param should_create_id IF TRUE will automatiquely create an ID. If false will not, but will search for ID in the df_index as specified in the argument ID.
#' @param ID Must be NA when create ID == TRUE. Otherwise, the ID name will be used to name the files.
#' @return a vcf.gz and vcf.gz.tbi with all necessary column, in a standard format and mapped to reference panel
#' @export
formattovcf_createindex <- function(all_out,
                                    snp_col,
                                    outcome_name,
                                    beta_col,
                                    se_col,
                                    pval_col,
                                    eaf_col,
                                    maf_col = NULL,
                                    effect_allele_col,
                                    other_allele_col,
                                    ncase_col = "ncase",
                                    ncontrol_col = "ncontrol",
                                    samplesize_col = "samplesize",
                                    chr_col,
                                    pos_col,
                                    units,  #"SD", "log odds", ou autre.
                                    traduction = fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt"),
                                    out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf",
                                    df_index = fread( "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt"),
                                    group_name ,
                                    year ,
                                    author ,
                                    consortium ,
                                    sex ,
                                    population ,
                                    initial_build ,
                                    category ,
                                    pmid ,
                                    note ,
                                    should_create_id = TRUE,
                                    ID = NA) {


  if(should_create_id == TRUE) {message("When creating ID, this function CANNOT BE RUN IN PARALLEL")}
  stopifnot(is.numeric(year))
  stopifnot(sex %in% c("Males and Females", "Males", "Females"))
  stopifnot(population %in% c("European", "African", "South Asian", "East Asian", "Mix"))
  if(initial_build != "HG19/GRCh37") {message("input must be in HG19/GRCh37. Please map the summary statistic to this build")}
  stopifnot(initial_build %in% c("HG19/GRCh37", "HG38/GRCh38", "HG18/Build36"))
  stopifnot(category %in% c("protein", "eqtl", "Metabolites","Trait","Disease"))
  stopifnot(is.numeric(pmid))
  if(!is.null(traduction)) {
    if(0 %in% traduction$EUR) {stop("Error: traduction cannot contain 0. You should replace by 0.001") }
    if( 1 %in% traduction$EUR) {stop("Error: traduction cannot contain 1. You should replace by 0.999")}
  }
  if(grepl(" ", outcome_name)) {stop("Error: outcome_name cannot contain space")}
  if(should_create_id == FALSE & df_index[,!(ID %in% id)]) {stop("Error: when should_create_id == FALSE, ID must exist in df_index$id")}
  if(grepl( "\\W", outcome_name) & !grepl("-", outcome_name)) {stop("Error: outcome_name can only include alphanumeric (upppercase or lowercase), underscore and hyphen")}
  if(!is.null(eaf_col) & !is.null(maf_col)) {stop("Error: if eaf_col is not NULL, then maf_col must be NULL")}

  #including all_out
  all_out[,outcome := outcome_name]
  oldname <- c(snp_col, beta_col, se_col, pval_col, effect_allele_col, other_allele_col, chr_col, pos_col, eaf_col, maf_col)
  newname <- c(unlist(ifelse(is.null(snp_col),list(NULL),"snp")), "beta", "se", "pval", "effect_allele", "other_allele", unlist(ifelse(is.null(chr_col),list(NULL),"chrom")), unlist(ifelse(is.null(pos_col),list(NULL),"pos")), unlist(ifelse(is.null(eaf_col),list(NULL),"eaf")), unlist(ifelse(is.null(maf_col),list(NULL),"maf")))
  setnames(all_out, oldname, newname)
  all_out[, c("beta", "se", "pval") := lapply(.SD, as.numeric), .SDcols = c("beta", "se", "pval")]

  if(is.null(ncase_col)){
    ncase_col<-ncase_col} else if(varhandle::check.numeric(ncase_col)) {
      all_out[,ncase := as.numeric(ncase_col)]
      ncase_col <- "ncase"
    }

  if(is.null(ncontrol_col)){
    ncontrol_col<-ncontrol_col}else if(varhandle::check.numeric(ncontrol_col)) {
      all_out[,ncontrol := as.numeric(ncontrol_col)]
      ncontrol_col <- "ncontrol"
    }

  if(varhandle::check.numeric(samplesize_col)) {
    all_out[,samplesize := as.numeric(samplesize_col)]
    samplesize_col <- "samplesize"
  }

  nrow_init <- all_out[,.N]
  all_out[, effect_allele := toupper(effect_allele)]
  all_out[, other_allele := toupper(other_allele)]

  if(!is.null(traduction)) {
    if(is.null(chr_col)) {all_out <- merge(all_out, traduction, by.x = c("snp"), by.y = c("rsid"), all = FALSE)} else {
      all_out <- merge(all_out, traduction, by.x = c("chrom", "pos"), by.y = c("chr", "position"), all = FALSE) }

    arguments <- data.table::data.table( ao_col = c( "snp_col", "chr_col", "pos_col"),
                             trad_col = c("rsid", "chr",  "position"))

    for(i in 1:nrow(arguments)) {
      if(is.null(get(arguments$ao_col[i]))) {
        setnames(all_out, arguments[i, trad_col], arguments$ao_col[i] %>% ifelse(. == "snp_col", "snp", . ) %>%
                   ifelse(. == "chr_col", "chrom", . ) %>% ifelse(. == "pos_col", "pos", . ))
      } }

    all_out <- all_out[(effect_allele == a0 | effect_allele == a1) & (other_allele == a0 | other_allele == a1) & a0 != a1  & effect_allele != other_allele, ] #because low number removed, coded on the forward strand
    all_out <- all_out[chrom %in% 1:22, ]
    nrow_harm <- all_out[, .N]
    all_out[, chrom := as.integer(chrom)]
    SwitchedAlleles <- all_out[effect_allele !=  a1, .N]

    all_out[effect_allele == a0, beta := beta*-1]
    if(!is.null(eaf_col)) { all_out[effect_allele == a0, eaf := 1-eaf] }
    all_out[effect_allele == a0, effect_allele := a1]
    all_out[other_allele == a1, other_allele := a0] #less than mrbase, possibly because I do not use the same traduction file
    if(!is.null(maf_col)) {all_out[, eaf := ifelse(EUR < 0.5, maf, 1-maf)]}
    if(is.null(maf_col) & is.null(eaf_col)) {all_out[, eaf := EUR]}


  } else {
    nrow_harm <- all_out[, .N]
    SwitchedAlleles <- 0
  }
  all_out_vcf <- gwasvcf::create_vcf(chrom= all_out[,chrom], pos= all_out[, pos], nea= all_out[,other_allele],
                                     ea= all_out[,effect_allele], snp= all_out[, snp], ea_af= all_out[, eaf],
                                     effect= all_out[, beta], se= all_out[,se], pval= all_out[,pval],
                                     n= all_out[, get(samplesize_col)], ncase = if(is.null(ncase_col)){NULL}else{all_out[, get(ncase_col)]},
                                     name= outcome_name)

  StudyType <- ifelse(units == "log odds", "CaseControl", "Continuous")


  df <- S4Vectors::DataFrame(TotalVariants = as.character(nrow_init), VariantsNotRead = as.character(0), HarmonisedVariants = as.character(nrow_harm),
                             VariantsNotHarmonised = as.character(nrow_init-nrow_harm), SwitchedAlleles = as.character(SwitchedAlleles),
                             StudyType = StudyType, row.names = outcome_name)

  df$TotalControls <- if(is.null(ncontrol_col)){NULL}else{all_out[, max(get(ncontrol_col))]}
  df$TotalCases <-if(is.null(ncase_col)){NULL}else{all_out[, max(get(ncase_col))]}

  VariantAnnotation::meta(VariantAnnotation::header(all_out_vcf))[["SAMPLE"]]  <- df

  VariantAnnotation::info(VariantAnnotation::header(all_out_vcf)) <- S4Vectors::DataFrame(Number = c("A", "A", "1"),
                                                                                          Type = c("Float", "Integer", "Integer"),
                                                                                          Description = c("Allele Frequency", "Allele count in genotypes", "Total number of alleles in called genotypes"),
                                                                                          row.names = c("AF", "AC", "AN"))
  if(should_create_id == TRUE) {
    fwrite(df_index,paste0(out_wd, "/server_gwas_id_copy.txt" )) #create a copy for safety purpose
    newrow <- data.frame(id = NA, trait = outcome_name, group_name = group_name, year = year, author = author, consortium = consortium,
                         sex = sex, population = population, unit = units, nsnp = nrow_harm, sample_size = all_out[, max(get(samplesize_col))],
                         initial_build = initial_build, category = category, pmid = pmid, ncase = 1,
                         sd = ifelse(units == "SD", 1, NA), note = note, ncontrol = 1)
    newrow$ncase <- if(is.null(ncase_col)){NA}else{all_out[, max(get(ncase_col))]}
    newrow$ncontrol <- if(is.null(ncontrol_col)){NA}else{all_out[, max(get(ncontrol_col))]}

    df_index <- rbind(df_index, newrow)
    df_index<- create_id(df_index)

    ID <- df_index[.N, id]
  }

  system2(command = "mkdir", args = paste0(out_wd, "/", ID))
  VariantAnnotation::writeVcf(all_out_vcf, file=paste0(out_wd, "/", ID, "/",ID, ".vcf"))
  current_wd<-getwd()
  setwd(paste0(out_wd, "/", ID))
  system2(command = "bgzip", args = c("-f", paste0(ID, ".vcf")))
  system2(command = "tabix", args = c("-f" ,"-p", "vcf", paste0(ID, ".vcf.gz")))
  setwd(current_wd)
  if(should_create_id == TRUE){ fwrite(df_index, "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt") }
  message("This script finished without errors")
}
#' Attribute a new id to your new phenotype
#'
#' @param df the df_index, with the new row
#'
#' @return the df_index with a new row and a new index
#' @export

create_id <- function(df) {

  id1 <- df[is.na(id), ]$category %>% ifelse(. == "protein", "prot", .) %>%
    ifelse(. == "eqtl", "eqtl", .) %>%
    ifelse(. == "Metabolites", "met", .) %>%
    ifelse(. == "Trait", "trait", .) %>%
    ifelse(. == "Disease", "dis", .)

  test <- df[group_name %in% df[is.na(id), ]$group_name & year %in% df[is.na(id), ]$year & author %in% df[is.na(id), ]$author
             & consortium %in% df[is.na(id), ]$consortium & category %in% df[is.na(id), ]$category, ][1,]
  if(!is.na(test$id)) {
    id2 <- test[, stringr::str_match(id, paste0(id1,"-(.*?)-"))[,2] ]
  } else {
    id2 <- df[category %in% df[is.na(id)]$category, .N, by = c("group_name", "year", "author", "consortium")][,.N]  #combien de différente combinaison existante incluant la row NA
  }

  id3 <- (sum(grepl(paste0("^",id1,"-",id2), df[,id])) + 1):(sum(grepl(paste0("^",id1,"-",id2), df[,id])) + length(id1))

  df[,id:=as.character(id)]
  df[is.na(id), id := paste0(id1, "-", id2, "-", id3)]
  return(df)
}


#' Get vcf file and its index directly from MRBAse
#'
#' @param id the MRBase id
#' @param out_wd where to store the data
#' @param nthreads the number of threads
#'
#' @return
#' @export

get_data_from_mrbase <- function(id, out_wd = "/home/couchr02/Mendel_Commun/Vcffile/MRBase_vcf", nthreads =3) {
  plan(multisession, workers = nthreads)
  map(as.list(id), function(x) system2(command = "mkdir", args = c("-p", paste0(out_wd, "/", x))))

  furrr::future_map(as.list(id), function(x)
    system2(command = "wget", args = c( "-O", paste0(out_wd, "/", x, "/", x, ".vcf.gz"), paste0("https://gwas.mrcieu.ac.uk/files/", x, "/", x, ".vcf.gz"))))

  furrr::future_map(as.list(id), function(x) system2(command = "tabix", args = c("-p", "vcf", paste0(out_wd, "/",x, "/", x, ".vcf.gz"))))

  message("script finished without errors")
}


#' Perform Coloc analysis
#'
#' @param vcffile_exp the path to the exposure
#' @param vcffile_out the path to the outcome
#' @param chrompos the genetic region ex : "22:25115489-26127836"
#' @param ldref the ld reference panel
#'
#' @return
#' @export

get_coloc <- function(vcffile_exp, vcffile_out, chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  message("initializing coloc")
  vout <- GagnonMR::gwasvcf_to_coloc_eloi(vcffile_exp, vcffile_out, chrompos)
  if(is.null(vout)) {return( data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples,
                                        outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples))}
  vout[[1]]$MAF <- vout[[1]]$MAF %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)
  vout[[2]]$MAF <- vout[[2]]$MAF %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)
  vres <- coloc::coloc.abf(vout[[1]], vout[[2]])
  vresres <- vres$results %>% as.data.table(.)
  df_res <- data.frame(exposure = vout$dataset1$id[1], outcome = vout$dataset2$id[1], posprob_coloc.mr = vres$summary[6],
                       posprob_coloc.SNP = vresres[which.max(SNP.PP.H4), snp] ,posprob_coloc.SNPexplained_var = vresres[which.max(SNP.PP.H4), SNP.PP.H4] )

  return(df_res)
}

#' Perform lead variant cis analysis
#'
#' @param vcffile_exp the path to the exposure .vcf.gz file
#' @param vcffile_out the path to the outcome .vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param ldref the path to the ldreference panel
#'
#' @return a data.table of results
#' @export

get_uni_cis <-  function(vcffile_exp, vcffile_out, chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  message("initializing uni-cis MR")
  dat_vcf <- gwasvcf::query_gwas(vcffile_exp, chrompos = chrompos)
  dat_tsmr <- gwasglue::gwasvcf_to_TwoSampleMR(dat_vcf)
  data.table::setDT(dat_tsmr)
  out_vcf <- tryCatch(
    expr = {  gwasvcf::query_gwas(vcf = vcffile_out, rsid = dat_tsmr[which.min(pval.exposure), ]$SNP, proxies = "yes", bfile = ldref)   },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]==0) {return(data.frame(exposure= dat_tsmr$exposure[1], outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, b.wald=NA,se.wald=NA,pval.wald=NA))}
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)

  harm <- TwoSampleMR::harmonise_data(dat_tsmr, out_tsmr, action = 1)
  harm <- TwoSampleMR::add_rsq(harm)
  harm$fstat.exposure <- fstat_fromdat(harm)
  harm <- TwoSampleMR::steiger_filtering( harm )
  res_all <- TwoSampleMR::mr(harm, method_list = "mr_wald_ratio") %>% data.table::as.data.table(.)
  if(dim(res_all)[1]==0) {return(data.frame(exposure= dat_tsmr$exposure[1], outcome = out_tsmr$outcome[1], b.wald=NA,se.wald=NA,pval.wald=NA))}
  res_all <- res_all[, .(outcome,  exposure, b, se ,pval)]
  data.table::setnames(res_all, old = c("b", "se", "pval"), new = paste0(c("b", "se", "pval"), ".wald"))
  res_wald <- cbind(res_all, data.frame(lead_snp.wald = harm$SNP))
  res_wald <- cbind(res_wald, pval_exposure.wald = harm$pval.exposure)
  res_wald <- cbind(res_wald, steiger_pval.wald = harm$steiger_pval)
  res_wald <- cbind(res_wald, steiger_dir.wald = harm$steiger_dir)
  col_to_select<-c("chr.exposure", "pos.exposure", "rsq.exposure", "fstat.exposure")
  res_wald <- cbind(res_wald, harm[,col_to_select])
  setnames(res_wald, col_to_select, paste0(col_to_select, ".wald"))
  return(res_wald)
}

#' Perform pan (cis + trans) analysis
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param ldref the path to the ldreference panel
#' @param chrompos the genetic region (only a dummy variable to make consitent form with get_uni_cis)
#'
#' @return a data.table of pan anaysis
#' @export

get_pan <- function(vcffile_exp, vcffile_out, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", chrompos = NA) {
  message("initializing pan mr")
  chrompos<-NA
  df_null <- data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples,
                        outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, b.ivw=NA,se.ivw=NA,pval.ivw=NA)

  if(file.exists(gsub(".vcf.gz", "_inst.txt", vcffile_exp))) {
    dat_sign <- fread(gsub(".vcf.gz", "_inst.txt", vcffile_exp))
    retain_snp <- dat_sign$SNP
  } else {

    dat_vcf <- gwasvcf::query_gwas(vcffile_exp, pval = 5e-8)
    if(dim(dat_vcf)[1]==0) {return(df_null)}
    dat_sign <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(. , "exposure") %>% data.table::as.data.table(.)

    retain_snp  <- tryCatch(
      expr = {dat_sign %>% dplyr::select(rsid=SNP, pval=pval.exposure, id = id.exposure) %>%
          ieugwasr::ld_clump(., plink_bin=genetics.binaRies::get_plink_binary(), bfile=ldref) %>%
          {.$rsid} },
      error = function(e){return(matrix(nrow = 0, ncol = 1)) })
  }

  out_vcf <- tryCatch(
    expr = {gwasvcf::query_gwas(vcf = vcffile_out, rsid = retain_snp, proxies = "yes", bfile = ldref) },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]==0) {return(df_null)}

  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)

  harm <- TwoSampleMR::harmonise_data(dat_sign, out_tsmr, action = 1)
  if(nrow(harm)>1){
    res_all <- TwoSampleMR::mr(harm, method_list = "mr_ivw")  %>% as.data.table
    res_all <- res_all[, .(outcome,  exposure, b, se ,pval, nsnp)]
    data.table::setnames(res_all, old = c("b", "se", "pval", "nsnp"), new = paste0(c("b", "se", "pval", "nsnp"), ".ivw"))
  } else{ res_all <- data.frame(outcome = harm$outcome[1], exposure = harm$exposure[1], b.ivw = NA, se.ivw  = NA, pval.ivw = NA)}
  return(res_all)
}


#' Get reverse MR
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region (dummy variable) only there for consistency with get_uni_cis
#'
#' @return a data.table of results
#' @export

get_reverseMR <- function(vcffile_exp, vcffile_out, chrompos = NA) {
  message("initializin reverseMR")
  res <- GagnonMR::get_pan(vcffile_exp = vcffile_out, vcffile_out = vcffile_exp)
  data.table::setnames(res, c("outcome", "exposure","b.ivw", "se.ivw", "pval.ivw"), c("exposure", "outcome", "b.reverseMR", "se.reverseMR", "pval.reverseMR"))
  return(res)
}



#' pipeline clump at 0.6 and correct for ld matrix
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param clumping_treshold the clumping threshold to do the analysis
#' @param ldref the path to the ldreferene panel
#'
#' @return a data.table of the results
#' @export

get_multicis <- function(vcffile_exp, vcffile_out,  chrompos, clumping_treshold = 0.6, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  message("initializing multicis mr")
  dat_tsmr <- gwasvcf::query_gwas(vcffile_exp, chrompos = chrompos) %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% data.table::as.data.table(.)
  df_null <- data.frame(exposure= dat_tsmr$exposure[1], outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, b.multi_cis=NA,se.multi_cis=NA,pval.multi_cis=NA, nsnip.multi_cis = 0)
  if(dat_tsmr[pval.exposure < 5e-8,.N] == 0) {return(df_null)}

  retain_snp  <- dat_tsmr[pval.exposure < 5e-8,] %>% dplyr::select(rsid=SNP, pval=pval.exposure, id = id.exposure) %>%
    ieugwasr::ld_clump(., clump_r2 = 0.6, clump_p = 5e-8, plink_bin=genetics.binaRies::get_plink_binary(), bfile=ldref) %>%
    {.$rsid}

  out_vcf <- tryCatch(
    expr = { gwasvcf::query_gwas(vcffile_out, rsid = retain_snp, proxies = "yes", bfile = ldref )  },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]<2) {return(df_null)}
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  harm <- TwoSampleMR::harmonise_data(dat_tsmr, out_tsmr, 1)
  ldmat <- ieugwasr::ld_matrix_local(harm$SNP, plink_bin = genetics.binaRies::get_plink_binary(), bfile = ldref)

  if (nrow(harm) < 2) {
    return(df_null)
  }
  else {
    x <- harm
    mriobj <- MendelianRandomization::mr_input(bx = x$beta.exposure,
                                               bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome,
                                               exposure = x$exposure[1], outcome = x$outcome[1],
                                               snps = x$SNP, effect_allele = x$effect_allele.exposure,
                                               other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)
    mriobj@correlation<-ldmat

    IVWobject <- MendelianRandomization::mr_ivw(mriobj, model = "default",
                                                correl = TRUE, distribution = "normal", alpha = 0.05)

    res_multicis <- data.frame(exposure = dat_tsmr[1, ]$exposure, outcome = out_tsmr$outcome[1], beta.multi_cis = IVWobject@Estimate,
                               se.multi_cis = IVWobject@StdError, pval.multi_cis = IVWobject@Pvalue,
                               cochranQ.multi_cis = IVWobject@Heter.Stat[1], cochranQpval.multi_cis = IVWobject@Heter.Stat[2],
                               nsnp.multi_cis = nrow(harm))
    return(res_multicis)
  }
}

#' Run all pqtl analyses with control on which analysis specifically to run
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param method_list the list of method to include "get_pan", "get_reverseMR", "get_uni_cis", "get_coloc", "get_multicis", "get_susie_coloc", "get_multicis_susie"
#'
#' @return a data.table with all results
#' @export

run_all_pqtl_analyses <- function(vcffile_exp, vcffile_out,  chrompos,
                                  method_list = list("get_pan", "get_reverseMR", "get_uni_cis", "get_coloc", "get_multicis")) {
  message(paste0("***********initilalizing all qtl analysis for ", vcffile_exp, " and ", vcffile_out, "************"))
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  if(is.na(chrompos)) {
    method_list <- method_list[sapply(method_list, function(x) (x %in% c("get_pan", "get_reverseMR")))]
  }

  if(length(method_list) == 0) {
    df_null <- data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples,
                          outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples) %>% as.data.table(.)
    message(paste0("no suitable method for ", df_null$exposure))
    return(df_null)
  }

  res <- lapply(method_list, function(meth) {
    get(meth)(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, chrompos = chrompos) })

  res_all <- res %>% purrr::reduce(merge, by = c("exposure", "outcome"))
  return(res_all)
}



#' Same as TwoSampleMR::clump_data, but using local clumping
#'
#' @param dat a TSMR format data.frame
#' @param ldref the path to the ldref panel
#' @param clump_kb same as in TwoSampleMR::clump_data
#' @param clump_r2 same as in TwoSampleMR::clump_data
#' @param clump_p same as in TwoSampleMR::clump_data
#'
#' @return a clumped data.frame
#' @export

clump_data_local <- function(dat, ldref ="/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", clump_kb = 1e4, clump_r2 = 0.001,
                             clump_p =1){
  d <- dat %>% dplyr::select(rsid=SNP, pval=pval.exposure, id = id.exposure)
  out <- ieugwasr::ld_clump(d, plink_bin=genetics.binaRies::get_plink_binary(), bfile=ldref, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = clump_p)
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
  return(dat[keep, ])
}


#' Obtain susie object from vcf
#'
#' @param vcffile a vcf-file path
#' @param chrompos chrompos
#' @param ldref th ld reference panel
#' @param L the maximum number of causal variants
#'
#' @return an object of class susie
#' @export
get_susie <- function(vcffile, chrompos,
                      ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", L = 10) {
  dat <- GagnonMR::gwasvcf_to_finemapr_gagnon(region = chrompos, vcf = vcffile, bfile = ldref)


  fitted_rss <- vector(mode = "list", length = L)
  for(i in 1:L) {
    tryCatch(expr = {
      fitted_rss[[i]] <- susieR::susie_rss(
        z = dat[[1]]$z$zscore,
        R = dat[[1]]$ld,
        L=i)
    }, error = function(e) {
      return(NULL)
    })
  }
  index <- sapply(fitted_rss, function(x) !is.null(x))
  if(sum(index)==0){return(NULL)}
  fitted_rss <- fitted_rss[[which(index)[sum(index) ]]]
  return(fitted_rss)
}



#' Obtain susie coloc results from wallace 2021.
#'
#' @param vcffile_exp dah
#' @param vcffile_out dah
#' @param chrompos dah
#' @param ldref dah
#'
#' @return
#' @export

get_susie_coloc <-function(vcffile_exp, vcffile_out,  chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  message("initializing susie coloc")

  df_null <- data.table(  exposure = gwasvcf::query_gwas(vcffile_exp,  chrompos = "1:40000-80000")@metadata$header@samples,
               outcome = gwasvcf::query_gwas(vcffile_out,  chrompos = "1:40000-80000")@metadata$header@samples)
  sus1 <- GagnonMR::get_susie(vcffile = vcffile_exp, chrompos = chrompos)
  sus2 <- GagnonMR::get_susie(vcffile = vcffile_out, chrompos = chrompos)
  if(is.null(sus1)|is.null(sus2)){return(df_null)}
  susie.res <- coloc::coloc.susie(sus1,sus2)

  susie_res_summary <-susie.res$summary
  data.table(setDT(susie_res_summary))
  if(nrow(susie_res_summary) == 0) { return(df_null)}

  k <- susie_res_summary[, .(hit1, PP.H4.abf, PP.H3.abf)]
  k<-k[,.SD[which.max(PP.H4.abf)],by="hit1"]
  k$index<-1:nrow(k)
  k$ID<-1
  data.table::setnames(k, c("hit1","PP.H4.abf", "PP.H3.abf"), c("susiecoloc.hit", "susiecoloc.PP.H4.abf", "susiecoloc.PP.H3.abf"))
  susie_res_summary <- data.table::dcast(k,ID  ~  index, value.var = c("susiecoloc.hit", "susiecoloc.PP.H4.abf", "susiecoloc.PP.H3.abf"))
  susie_res_summary[,ID := NULL]
  susie_res_summary$exposure <- gwasvcf::query_gwas(vcffile_exp,  chrompos = "1:40000-80000")@metadata$header@samples
  susie_res_summary$outcome <- gwasvcf::query_gwas(vcffile_out,  chrompos = "1:40000-80000")@metadata$header@samples

  return(susie_res_summary)

}


#' Same as get_multicis, but instead of selecting instruments based on LD clumping select with finemapping using SuSiE.
#'
#' @param vcffile_exp dah
#' @param vcffile_out dah
#' @param chrompos dah
#' @param ldref dah
#'
#' @return select as instrument the SNP with the highest PIP in each credible set. If no SNP has a p-value under 5e-8 or number of credible sets < 2 return empty data.frame.
#' @export

get_multicis_susie <- function(vcffile_exp, vcffile_out,  chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  message("initializing multicis mr with susie finemapping")
  df_null <- data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples, outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, nsnp.multi_cis_susie = 0)
  fitted_rss <- GagnonMR::get_susie(vcffile = vcffile_exp, chrompos = chrompos)
  if(is.null(fitted_rss)) {return(df_null)}
  dt_sus <- summary(fitted_rss)$vars %>% data.table::as.data.table(.)
  dt_sus <- dt_sus[cs != -1, ]
  retain_snp <- names(fitted_rss$pip)[dt_sus[,.SD[which.max(variable_prob)] , by = "cs"]$variable]
  df_null$nsnp.multi_cis_susie <- length(retain_snp)

  dat_tsmr <- gwasvcf::query_gwas(vcffile_exp, chrompos = chrompos) %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% data.table::as.data.table(.)
  # if(dat_tsmr[pval.exposure < 5e-8,.N] == 0) {return(df_null)}

  out_vcf <- tryCatch(
    expr = { gwasvcf::query_gwas(vcffile_out, rsid = retain_snp, proxies = "yes", bfile = ldref )  },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]<2) {return(df_null)}
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  harm <- TwoSampleMR::harmonise_data(dat_tsmr, out_tsmr, 1)
  harm <- TwoSampleMR::steiger_filtering(harm) %>% as.data.table()
ndirectionalyinconsistent <- harm[, sum(steiger_dir == FALSE)]
harm <- harm[steiger_dir == TRUE,]

  ldmat <- ieugwasr::ld_matrix_local(harm$SNP, plink_bin = genetics.binaRies::get_plink_binary(), bfile = ldref)
  ldmat_test<-ldmat ; diag(ldmat_test) <- 0
  if (nrow(harm) < 2 | any(ldmat_test == 1) ) {
    return(df_null)
  }
  else {
    x <- harm
    mriobj <- MendelianRandomization::mr_input(bx = x$beta.exposure,
                                               bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome,
                                               exposure = x$exposure[1], outcome = x$outcome[1],
                                               snps = x$SNP, effect_allele = x$effect_allele.exposure,
                                               other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)
    mriobj@correlation<-ldmat

    IVWobject <- MendelianRandomization::mr_ivw(mriobj, model = "default",
                                                correl = TRUE, distribution = "normal", alpha = 0.05)

    res_multicis <- data.frame(exposure = dat_tsmr[1, ]$exposure, outcome = out_tsmr$outcome[1], beta.multi_cis_susie = IVWobject@Estimate,
                               se.multi_cis_susie = IVWobject@StdError, pval.multi_cis_susie = IVWobject@Pvalue,
                               cochranQ.multi_cis_susie = IVWobject@Heter.Stat[1], cochranQpval.multi_cis_susie = IVWobject@Heter.Stat[2],
                               nsnp.multi_cis_susie = nrow(harm), minpvalexposure.multi_cis_susie = dat_tsmr[,min(pval.exposure)],
                               nsteigerfalse.multi_cis_susie = ndirectionalyinconsistent)

    return(res_multicis)
  }
}




#' Get instrument ld clump = .001 and pvalue < 5e-8
#'
#' @param vcffile the path to the vcffile
#' @param should_write if FALSE, simply return the value if TRUE write in the directory
#'
#' @return
#' @export

get_inst <- function(vcffile, pval = 5e-08, clump = TRUE, r2 = 0.001,
                     kb = 10000, should_write = FALSE) {
  dat_vcf <- gwasvcf::query_gwas(vcffile, pval = pval)
  dat_sign <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "exposure") %>% data.table::as.data.table(.)

  if(clump==TRUE) {
    dat_sign <- GagnonMR::clump_data_local(dat_sign, clump_kb = kb, clump_r2 = r2,)
  }
  if(should_write == FALSE)  {
    fwrite(dat_sign, gsub(".vcf.gz", "_inst.txt", vcffile))
  }
  return(dat_sign)
}
