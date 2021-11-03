#' Title
#'
#' @param genecard_name the name of the gene genecard name
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

if(nrow(gene_coords) == 0) {
  gene_coords <- biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"),
                                filters = "hgnc_symbol", values = genecard_name,
                                mart = ensembl) %>% as.data.table(.)
  gene_coords <- gene_coords[chromosome_name %in% 1:22,]
  colnames(gene_coords)<-c("gene_name", "chr", "start", "end")
} else {
  gene_coords[, start := min(start)]
  gene_coords[, end := max(end)]
  gene_coords[,gene_name := paste0(unique(gene_name), collapse = "-")]
  gene_coords <- gene_coords %>% distinct(.)
  gene_coords <- gene_coords[chr %in% 1:22,]
}
setDT(gene_coords)
if(gene_coords[, length(unique(chr)) > 1] | gene_coords[,.N]==0) { gene_region <- NULL} else {
  gene_coords[, start_cis := (start - window/2) %>% ifelse(.<1, 1, .) ]
  gene_coords[, end_cis := end + window/2]
  gene_region <- gene_coords[, paste0(chr, ":", start_cis, "-", end_cis)]
}

return(gene_region)
}




#' Title
#'
#' @param all_out the data.frame object
#' @param snp_col the name of the snp col, NULL if absent
#' @param outcome_name the name of the traits
#' @param beta_col the name of the beta col, Must be present
#' @param se_col the name of the se col, Must be presnt
#' @param pval_col the name of the pval col, Must be present
#' @param eaf_col the name of the eaf col, NULL if absent
#' @param effect_allele_col the name of the effect allele, must be present
#' @param other_allele_col the name of the other allele, must be present
#' @param ncase_col The number of ncase (in numeric) or the name of the column where the info can be retrieved (in character) otherwise NULL
#' @param ncontrol_col  The number of ncontrol (in numeric) or the name of the column where the info can be retrieved (in character) otherwise NULL
#' @param samplesize_col The samplesize (in numeric) or the name of the column where the info can be retrieved (in character)
#' @param chr_col the name of the chr col, NULL if absent
#' @param pos_col the name of the pos col, NULL if absent
#' @param units ideally "SD", "log odds". If this is not the case, can put any character such as "year" or "Cm"
#' @param traduction the traduction file to map to the reference panel
#' @param out_wd where you want your vcf.gz and vcf.gz.tbi to be stored
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
#'
#' @return a vcf.gz and vcf.gz.tbi with all necessary column, in a standard format and mapped to reference panel
#' @export

formattovcf_createindex <- function(all_out,
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
                                     chr_col,
                                     pos_col,
                                     units,  #"SD", "log odds", ou autre.
                                     traduction = fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt"),
                                     out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf/",
                                     group_name ,
                                     year ,
                                     author ,
                                     consortium ,
                                     sex ,
                                     population ,
                                     initial_build ,
                                     category ,
                                     pmid ,
                                     note ) {
  #test
  stopifnot(is.numeric(year))
  stopifnot(sex %in% c("Males and Females", "Males", "Females"))
  stopifnot(population %in% c("European", "African", "South Asian", "East Asian", "Mix"))
  if(initial_build != "HG19/GRCh37") {message("input must be in HG19/GRCh37. Please map the summary statistic to this build")}
  stopifnot(initial_build %in% c("HG19/GRCh37", "HG38/GRCh38"))
  stopifnot(category %in% c("protein", "eqtl", "Metabolites","Trait","Disease"))
  stopifnot(is.numeric(pmid))

  #putting traduction in right format
  traduction[,chr:=as.integer(chr)]
  traduction[,pos:=as.integer(position)]
  setorder(traduction, chr, pos)
  traduction <- distinct(traduction)
  traduction <- traduction[, .(rsid, chr, position, a0, a1, EUR)]
  traduction[ , EUR := EUR %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)]

  #including all_out
  all_out <- fread(gwas_file)
  all_out[,outcome := outcome_name]

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

  nrow_init <- all_out[,.N]
  all_out[, c(effect_allele_col) := toupper(get(effect_allele_col))]
  all_out[, c(other_allele_col) := toupper(get(other_allele_col))]

  if(is.null(chr_col)) {all_out <- merge(all_out, traduction, by.x = c(snp_col), by.y = c("rsid"), all = FALSE)} else {
    all_out <- merge(all_out, traduction, by.x = c(chr_col, pos_col), by.y = c("chr", "position"), all = FALSE) }

  arguments <- data.frame( ao_col = c( "snp_col", "chr_col", "pos_col", "other_allele_col", "effect_allele_col" ,"eaf_col"),
                           trad_col = c("rsid", "chr",  "position", "a0", "a1", "EUR"))

  for(i in 1:ncol(arguments)) {
    if(is.null(get(arguments$ao_col[i]))) {
      assign(arguments$ao_col[i], ao$trad_col[i])
    } }

  all_out <- all_out[(get(effect_allele_col) == a0 | get(effect_allele_col) == a1) & (get(other_allele_col) == a0 | get(other_allele_col) == a1) & a0 != a1 & get(effect_allele_col) != get(other_allele_col), ] #because low number removed, coded on the forward strand
  all_out <- all_out[get(chr_col) %in% 1:22, ]
  nrow_harm <- all_out[, .N]
  all_out[, c(chr_col) := as.integer(get(chr_col))]
  SwitchedAlleles <- all_out[get(effect_allele_col) !=  a1, .N]

  all_out[get(effect_allele_col) == a0, c(beta_col) := get(beta_col)*-1]
  all_out[get(effect_allele_col) == a0, c(effect_allele_col) := a1]
  all_out[get(other_allele_col) == a1, c(other_allele_col) := a0] #less than mrbase, possibly because I do not use the same traduction file

  debugonce(create_vcf)
  all_out_vcf <- gwasvcf::create_vcf(chrom= all_out[,get(chr_col)], pos= all_out[, get(pos_col)], nea= all_out[,get(other_allele_col)],
                                     ea= all_out[,get(effect_allele_col)], snp= all_out[, get(snp_col)], ea_af= all_out[, get(eaf_col)],
                                     effect= all_out[, get(beta_col)], se= all_out[,get(se_col)], pval= all_out[, get(pval_col)],
                                     n= all_out[, get(samplesize_col)], ncase = if(is.null(ncase_col)){NULL}else{all_out$ncase},
                                     name= outcome_name)

  StudyType <- ifelse(units == "log odds", "CaseControl", "Continuous")

  df <- DataFrame(TotalVariants = as.character(nrow_init), VariantsNotRead = as.character(0), HarmonisedVariants = as.character(nrow_harm),
                  VariantsNotHarmonised = as.character(nrow_init-nrow_harm), SwitchedAlleles = as.character(SwitchedAlleles),
                  TotalControls = if(is.null(ncontrol_col)){NULL}else{all_out[, max(get(ncontrol_col))]},
                  TotalCases = if(is.null(ncase_col)){NULL}else{all_out[, max(get(ncase_col))]},
                  StudyType = StudyType, row.names = outcome_name)


  VariantAnnotation::meta(VariantAnnotation::header(all_out_vcf))[["SAMPLE"]]  <- df

  info(VariantAnnotation::header(all_out_vcf)) <- DataFrame(Number = c("A", "A", "1"),
                                                            Type = c("Float", "Integer", "Integer"),
                                                            Description = c("Allele Frequency", "Allele count in genotypes", "Total number of alleles in called genotypes"),
                                                            row.names = c("AF", "AC", "AN"))

  df_index <- fread( "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")

  newrow <- data.frame(id = NA, trait = outcome_name, group_name = group_name, year = year, author = author, consortium = consortium,
                       sex = sex, population = population, unit = units, nsnp = nrow_harm, sample_size = all_out[, max(get(samplesize_col))],
                       initial_build = initial_build, category = category, pmid = pmid, ncase = all_out[, max(get(ncase_col))],
                       sd = ifelse(units == "SD", 1, NA), note = note, ncontrol = all_out[, max(get(ncontrol_col))])

  df_index <- rbind(df_index, newrow)
  df_index<- create_id(df_index)
  fwrite(df_index, "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")

  ID <- df_index[trait == SOMAmerID, id]
  system2(command = "mkdir", args = c("-p", paste0(out_wd, ID)))
  writeVcf(sun_full_vcf, file=paste0(out_wd, ID, "/",ID, ".vcf"))
  current_wd<-getwd()
  setwd(paste0(out_wd, ID))
  system2(command = "bgzip", args = c("-f", paste0(ID, ".vcf")))
  system2(command = "tabix", args = c("-f" ,"-p", "vcf", paste0(ID, ".vcf.gz")))
  setwd(current_wd)
  message("This script finished without errors")
}

#' Title attribute a new id to your new phenotype
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

  id2 <- df[!is.na(id) & category %in% df[is.na(id)]$category, .N, by = c("group_name", "year", "author", "consortium")][,.N + 1]  #combien de combinaison existante déjà

  id3 <- 1:df[is.na(id), .SD, by = c("category", "group_name", "year", "author", "consortium"), .SDcols = "trait"][,.N]

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

get_data_from_mrbase <- function(id, out_wd = "/home/couchr02/Mendel_Commun/Vcffile/MRBase_vcf/", nthreads =3) {
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
  dat_tsmr <- gwasglue::gwasvcf_to_TwoSampleMR(dat_vcf) %>% as.data.table(.)
  out_vcf <- tryCatch(
    expr = {  gwasvcf::query_gwas(vcf = vcffile_out, rsid = dat_tsmr[which.min(pval.exposure), ]$SNP, proxies = "yes", bfile = ldref)   },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]==0) {return(data.frame(exposure= dat_tsmr$exposure[1], outcome = query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, b.wald=NA,se.wald=NA,pval.wald=NA))}
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)

  harm <- TwoSampleMR::harmonise_data(dat_tsmr, out_tsmr, action = 1)
  res_all <- TwoSampleMR::mr(harm, method_list = "mr_wald_ratio") %>% data.table::as.data.table(.)
  res_all <- res_all[, .(outcome,  exposure, b, se ,pval)]
  data.table::setnames(res_all, old = c("b", "se", "pval"), new = paste0(c("b", "se", "pval"), ".wald"))
  res_wald <- cbind(res_all, data.frame(lead.snp.wald = harm$SNP))
  res_wald <- cbind(res_wald, pval.exposure.wald = harm$pval.exposure)
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

get_pan <- function(vcffile_exp, vcffile_out, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", chrompos = NULL) {
  message("initializing pan mr")
  df_null <- data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples,
                        outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, b.ivw=NA,se.ivw=NA,pval.ivw=NA)
  dat_vcf <- gwasvcf::query_gwas(vcffile_exp, pval = 5e-8)
  if(dim(dat_vcf)[1]==0) {return(df_null)}
  dat_sign <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(. , "exposure") %>% data.table::as.data.table(.)

  retain_snp  <- tryCatch(
    expr = {dat_sign %>% dplyr::select(rsid=SNP, pval=pval.exposure, id = id.exposure) %>%
        ieugwasr::ld_clump(., plink_bin=genetics.binaRies::get_plink_binary(), bfile=ldref) %>%
        {.$rsid} },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  out_vcf <- tryCatch(
    expr = {gwasvcf::query_gwas(vcf = vcffile_out, rsid = retain_snp, proxies = "yes", bfile = ldref) },
    error = function(e){return(matrix(nrow = 0, ncol = 1)) })

  if(dim(out_vcf)[1]==0) {return(df_null)}

  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)

  harm <- TwoSampleMR::harmonise_data(dat_sign, out_tsmr, action = 1)
  if(nrow(harm)>1){
    res_all <- TwoSampleMR::mr(harm, method_list = "mr_ivw")  %>% as.data.table
    res_all <- res_all[, .(outcome,  exposure, b, se ,pval)]
   data.table::setnames(res_all, old = c("b", "se", "pval"), new = paste0(c("b", "se", "pval"), ".ivw"))
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

get_reverseMR <- function(vcffile_exp, vcffile_out, chrompos = NULL) {
  message("initializin reverseMR")
  res <- GagnonMR::get_pan(vcffile_exp = vcffile_out, vcffile_out = vcffile_exp, chompos = chrompos)
  data.table::setnames(res, colnames(res), c("exposure", "outcome", "b.reverseMR", "se.reverseMR", "pval.reverseMR"))
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
#' @param method_list the list of method to include
#'
#' @return a data.table with all results
#' @export

run_all_pqtl_analyses <- function(vcffile_exp, vcffile_out,  chrompos,
                                  method_list = list("get_pan", "get_reverseMR", "get_uni_cis", "get_coloc", "get_multicis")) {
  message(paste0("***********initilalizing all qtl analysis for ", vcffile_exp, " and ", vcffile_out, "************"))

  if(is.null(chrompos)) {
    method_list <- method_list[sapply(method_list, function(x) (x %in% c("get_pan", "get_reverseMR")))]
  }

  if(length(method_list) == 0) {
    df_null <- data.frame(exposure= query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples,
                          outcome = query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples)
    message(paste0("no suitable method for ", df_null$exposure))
    return(df_null)
  }

  res <- lapply(method_list, function(meth) {
    get(meth)(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, chrompos = chrompos) })

  res_all <- res %>% purrr::reduce(merge, by = c("exposure", "outcome"))
  return(res_all)
}

