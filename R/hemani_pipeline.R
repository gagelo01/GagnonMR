#' get hgnc nomenclature from uniprot id assessing multiple nomenclature c("uniprotswissprot", "uniprot_gn_id", "uniprot_gn_symbol", "uniprot_isoform",  "uniprotsptrembl") and NCBI database
#' @param uniprot_vec A vector of Uniprot ID
#' @param ensembl ensembl database
#' @param uniprot_dictionnary uniprot dictionnary downloaded from NCBI
#'
#' @return a data.table of uniprot and hgnc_id
#' @export
from_uniprot_to_hgnc <- function(uniprot_vec,
                                 ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl"),
                                 uniprot_dictionnary = fread("/mnt/sda/boujer01/DATA/Drug_Targets/Uniprot_ID/uniprot_final.txt")) {

  colnames(uniprot_dictionnary)<- c("uniprot", "gene_name", "n_aa")
  uniprot_filters <- c("uniprotswissprot", "uniprot_gn_id", "uniprot_gn_symbol", "uniprot_isoform",  "uniprotsptrembl")
  uniprot_biomart <- data.table(hgnc_symbol = "hgnc_symbol", uniprot_gn_id = "uniprot_gn_id")
  for(i in 1:length(uniprot_filters)) {
    values <- uniprot_vec[!(uniprot_vec %in% uniprot_biomart$uniprot_gn_id)]
    if(length(values)>0){
    k <- biomaRt::getBM(attributes=c(uniprot_filters[i], "hgnc_symbol"),
                        filters =  uniprot_filters[i] , values=values, mart=ensembl)
    setDT(k)
    setnames(k, uniprot_filters[i], "uniprot_gn_id")
    uniprot_biomart <- rbind(uniprot_biomart, k)
    }
  }
  values <- uniprot_vec[!(uniprot_vec %in% uniprot_biomart$uniprot_gn_id)]
if(length(values)>0) {
  vec_res <- vector(mode = "list", length = length(values))
  for(i in 1:length(values)) {
    x  <- uniprot_dictionnary[grepl(gsub(" ", "|", values[i],),uniprot),] #gsub("-.*", "",values[i])
    x[,uniprot_gn_id := values[i]]
    setnames(x, "gene_name",c("hgnc_symbol"))
    x <- distinct(x[, .(hgnc_symbol, uniprot_gn_id)]) %>% as.data.table(.)
    x[, hgnc_symbol := paste(hgnc_symbol, collapse = " ")]
    vec_res[[i]] <- distinct(x[, .(hgnc_symbol, uniprot_gn_id)]) %>% as.data.table(.)
  }
  dt_res <- rbindlist(vec_res)
  uniprot_biomart <- rbind(uniprot_biomart, dt_res)
}
  uniprot_biomart <- uniprot_biomart[!(hgnc_symbol == "hgnc_symbol"),]
  uniprot_biomart <- uniprot_biomart[hgnc_symbol != "",]

  double_hgnc <-uniprot_biomart[, .N, by = "uniprot_gn_id"][N>1]$uniprot_gn_id
  uniprot_biomart[ , hgnc_symbol := unique(paste(hgnc_symbol, collapse = " ")), by = "uniprot_gn_id"]
  uniprot_biomart <- distinct(uniprot_biomart) %>% as.data.table(.)

  return(uniprot_biomart)
}

#' Title
#'
#' @param genecard_name A vector the names of the gene genecard name
#' @param window the window 1e6 == 5e5 upward 5e5 downward
#' @param gencode the gencode data.frame
#' @param ensembl ensembl
#'
#' @return
#' @export

from_genecard_to_generegion <- function(genecard_name, window = 1e+06, gencode = data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt"),
                                        ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                                                   host = "grch37.ensembl.org", path = "/biomart/martservice",
                                                                   dataset = "hsapiens_gene_ensembl"))
{
  genevec <- genecard_name %>% strsplit(., split = " ") %>% unlist(.)
  gene_coords <- gencode[gene_name %in% genevec, ]
  if (!all(genevec %in% gene_coords$gene_name) | nrow(gene_coords) ==  0) {
    gene_coords_mart <- biomaRt::getBM(attributes = c("hgnc_symbol",
                                                      "chromosome_name",
                                                      "start_position",
                                                      "end_position"),
                                       filters = "hgnc_symbol",
                                       values = genevec[!(genevec %in% gene_coords$gene_name)], mart = ensembl) %>%
      as.data.table(.)
    # gene_coords_mart <- gene_coords_mart[chromosome_name %in%  1:22, ]
    colnames(gene_coords_mart) <- c("gene_name", "chr", "start","end")
    gene_coords <- rbindlist(list(gene_coords, gene_coords_mart),
                             fill = TRUE)
  }
  gene_coords[, `:=`(start, min(start)), by = "gene_name"]
  gene_coords[, `:=`(end, max(end)), by = "gene_name"]
  gene_coords[, `:=`(gene_name, paste0(unique(gene_name), collapse = "-")),
              by = "gene_name"]
  gene_coords <- gene_coords %>% distinct(.)
  # gene_coords <- gene_coords[chr %in% 1:22, ]
  setDT(gene_coords)
  if (gene_coords[, .N] == 0) {
    return(NA)
  }

  gene_coords[, start_cis := (start - window/2) %>% ifelse(.<1, 1, .) ]
  gene_coords[, end_cis := (end + window/2)]
  gene_region <- sapply(genecard_name, function(x) {
    k <- gene_coords[gene_name %in% strsplit(x, split = " ")[[1]],]
    if(nrow(k)==0){return(NA)}
    gene_region <- k[,paste0(unique(chr), ":", min(start_cis), "-", max(end_cis))]
    return(gene_region)})

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
                                    traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt"),
                                    out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                    df_index = fread( "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt"),
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
  stopifnot(population %in% c("European", "African", "South Asian", "Asian unspecified",
                              "East Asian","Mixed", "Mix", "Hispanic or Latin American"))
  if(initial_build != "HG19/GRCh37") {message("input must be in HG19/GRCh37. Please map the summary statistic to this build")}
  stopifnot(initial_build %in% c("HG19/GRCh37", "HG38/GRCh38", "HG18/Build36"))
  stopifnot(category %in% c("protein", "eqtl", "Metabolites","Trait","Disease"))
  stopifnot(is.numeric(pmid))
  if(!is.null(traduction)) {
    if(0 %in% traduction[1:1000,]$EUR) {stop("Error: traduction cannot contain 0. You should replace by 0.001") }
    if( 1 %in% traduction[1:1000,]$EUR) {stop("Error: traduction cannot contain 1. You should replace by 0.999")}
  }
  if(grepl(" ", outcome_name)) {stop("Error: outcome_name cannot contain space")}
  if(should_create_id == FALSE & df_index[, !any(id %in% ID)]) {stop("Error: when should_create_id == FALSE, ID must exist in df_index$id")}
  if(grepl( "\\W", outcome_name) & !grepl("-", outcome_name) & !grepl("\\.", outcome_name)) {stop("Error: outcome_name can only include alphanumeric (upppercase or lowercase), underscore, hyphen and dot")}
  if(!is.null(eaf_col) & !is.null(maf_col)) {stop("Error: if eaf_col is not NULL, then maf_col must be NULL")}
  k <- all_out[, sapply(.SD, function(x) any(is.na(x)))]
  if(any(k)){stop(paste0("Error: column --", paste(names(k)[k==TRUE], collapse = "-- and --"), "-- contains NA. Please verify and potentially remove these rows."))}
  if(!is.null(snp_col)){
    if(all_out[,any(!grepl("^rs", get(snp_col)))]) {stop("Error: SNP column contains rsid that do not start by 'rs'")}}
  if(!(is.null(ncase_col)&is.null(ncontrol_col)) & units != "log odds") {stop("Error: When trait is continuous, ncase_col and ncontrol_col should be NULL")}

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
      all_out[,c("chrom", "pos") := lapply(.SD, as.integer), .SDcols = c("chrom", "pos")]
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
                                     n= all_out[, get(samplesize_col)], ncase = if(is.null(ncase_col)){NA}else{all_out[, get(ncase_col)]},
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
  system2(command = "bgzip", args = c("-f", paste0(out_wd, "/", ID, "/", ID, ".vcf")))
  system2(command = "tabix", args = c("-f" ,"-p", "vcf", paste0(out_wd, "/", ID, "/", ID, ".vcf.gz")))
  if(should_create_id == TRUE){ fwrite(df_index, "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt") }
  message("This script finished without errors")
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
#' @param pos_b37_col the name of the pos col in build 37, NULL if absent
#' @param pos_b38_col the name of the pos col in build 38, NULL if absent
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
#' @param ID Must be NA when create ID == TRUE. Otherwise, the ID name will be used to name the files.
#' @return a vcf.gz and vcf.gz.tbi with all necessary column, in a standard format and mapped to reference panel
#' @export

formattovcf_createindex2 <- function(all_out,
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
                                     pos_b37_col,
                                     pos_b38_col,
                                     units,  #"SD", "log odds", ou autre.
                                     traduction,
                                     out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                     df_index = fread( "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt"),
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
                                     ID = NA) {


  stopifnot(is.numeric(year))
  stopifnot(sex %in% c("Males and Females", "Males", "Females"))
  stopifnot(population %in% c("European", "African", "South Asian", "Asian unspecified",
                              "East Asian","Mixed", "Mix", "Hispanic or Latin American"))
  if(initial_build != "HG19/GRCh37") {message("input must be in HG19/GRCh37. Please map the summary statistic to this build")}
  stopifnot(initial_build %in% c("HG19/GRCh37", "HG38/GRCh38", "HG18/Build36"))
  stopifnot(category %in% c("protein", "eqtl", "Metabolites","Trait","Disease"))
  stopifnot(is.numeric(pmid))
  if(!is.null(traduction)) {
    if(0 %in% traduction[1:1000,]$eur_eaf) {stop("Error: traduction cannot contain 0. You should replace by 0.001") }
    if( 1 %in% traduction[1:1000,]$eur_eaf) {stop("Error: traduction cannot contain 1. You should replace by 0.999")}
  }
  if(grepl(" ", outcome_name)) {stop("Error: outcome_name cannot contain space")}
  if(df_index[, !any(id %in% ID)]) {stop("Error:  ID must exist in df_index$id")}
  if(grepl( "\\W", outcome_name) & !grepl("-", outcome_name) & !grepl("\\.", outcome_name)) {stop("Error: outcome_name can only include alphanumeric (upppercase or lowercase), underscore, hyphen and dot")}
  if(!is.null(eaf_col) & !is.null(maf_col)) {stop("Error: if eaf_col is not NULL, then maf_col must be NULL")}
  k <- all_out[, sapply(.SD, function(x) any(is.na(x)))]
  if(any(k)){stop(paste0("Error: column --", paste(names(k)[k==TRUE], collapse = "-- and --"), "-- contains NA. Please verify and potentially remove these rows."))}
  if(!is.null(snp_col)){
    if(all_out[,any(!grepl("^rs", get(snp_col)))]) {stop("Error: SNP column contains rsid that do not start by 'rs'")}}
  if(!(is.null(ncase_col)&is.null(ncontrol_col)) & units != "log odds") {stop("Error: When trait is continuous, ncase_col and ncontrol_col should be NULL")}
  if(!is.null(pos_b37_col)&!is.null(pos_b38_col)) {stop("Error: you can only provide one of pos_b37_col and pos_b38_col. The other must be NULL")}
  stopifnot(all(c("chr","pos_b38","other_all","effect_all", "rsid", "pos_b37", "eur_eaf") %in% colnames(traduction)))
  #including all_out
  all_out[,outcome := outcome_name]
  oldname <- c(snp_col, beta_col, se_col, pval_col, effect_allele_col, other_allele_col, chr_col, pos_b37_col, pos_b38_col, eaf_col, maf_col)
  newname <- c(unlist(ifelse(is.null(snp_col),list(NULL),"snp")), "beta", "se", "pval", "effect_allele", "other_allele", unlist(ifelse(is.null(chr_col),list(NULL),"chrom")), unlist(ifelse(is.null(pos_b37_col),list(NULL),"pos")), unlist(ifelse(is.null(pos_b38_col),list(NULL),"pos")), unlist(ifelse(is.null(eaf_col),list(NULL),"eaf")), unlist(ifelse(is.null(maf_col),list(NULL),"maf")))
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
    if(is.null(chr_col)) {all_out <- merge(all_out, traduction, by.x = c("snp"), by.y = c("rsid"), all = FALSE)}
    else if (!is.null(pos_b37_col)) {
      all_out[,c("chrom", "pos") := lapply(.SD, as.integer), .SDcols = c("chrom", "pos")]
      all_out <- merge(all_out, traduction, by.x = c("chrom", "pos"), by.y = c("chr", "pos_b37"), all = FALSE) }
    else if(!is.null(pos_b38_col)) {
      all_out <- merge(all_out, traduction, by.x = c("chrom", "pos"), by.y = c("chr", "pos_b38"), all = FALSE)
      all_out[,pos := pos_b37]
    }

    all_out <- all_out[(effect_allele == other_all | effect_allele == effect_all) & (other_allele == other_all | other_allele == effect_all) & other_all != effect_all  & effect_allele != other_allele, ] #because low number removed, coded on the forward strand
    all_out <- all_out[chrom %in% 1:22, ]
    nrow_harm <- all_out[, .N]
    all_out[, chrom := as.integer(chrom)]
    SwitchedAlleles <- all_out[effect_allele !=  effect_all, .N]

    all_out[effect_allele == other_all, beta := beta*-1]
    if(!is.null(eaf_col)) { all_out[effect_allele == other_all, eaf := 1-eaf] }
    all_out[effect_allele == other_all, effect_allele := effect_all]
    all_out[other_allele == effect_all, other_allele := other_all]
    if(!is.null(maf_col)) {all_out[, eaf := ifelse(eur_eaf < 0.5, maf, 1-maf)]}
    if(is.null(maf_col) & is.null(eaf_col)) {all_out[, eaf := eur_eaf]}


  } else {
    nrow_harm <- all_out[, .N]
    SwitchedAlleles <- 0
  }
  all_out_vcf <- gwasvcf::create_vcf(chrom= all_out[,chrom], pos= all_out[, pos], nea= all_out[,other_allele],
                                     ea= all_out[,effect_allele], snp= all_out[, snp], ea_af= all_out[, eaf],
                                     effect= all_out[, beta], se= all_out[,se], pval= all_out[,pval],
                                     n= all_out[, get(samplesize_col)], ncase = if(is.null(ncase_col)){NA}else{all_out[, get(ncase_col)]},
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

  system2(command = "mkdir", args = paste0(out_wd, "/", ID))
  VariantAnnotation::writeVcf(all_out_vcf, file=paste0(out_wd, "/", ID, "/",ID, ".vcf"))
  system2(command = "bgzip", args = c("-f", paste0(out_wd, "/", ID, "/", ID, ".vcf")))
  system2(command = "tabix", args = c("-f" ,"-p", "vcf", paste0(out_wd, "/", ID, "/", ID, ".vcf.gz")))
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
#'
#' @return
#' @export

get_data_from_mrbase <- function(id, out_wd = "/mnt/sda/gagelo01/Vcffile/MRBase_vcf") {
  purrr::map(as.list(id), function(x) {
    message(paste0("Initializing ", x))
    system2(command = "mkdir", args = c("-p", paste0(out_wd, "/", x)))
    system2(command = "wget", args = c( "-O", paste0(out_wd, "/", x, "/", x, ".vcf.gz"), paste0("https://gwas.mrcieu.ac.uk/files/", x, "/", x, ".vcf.gz")))
    system2(command = "tabix", args = c("-p", "vcf", paste0(out_wd, "/",x, "/", x, ".vcf.gz")))
    })
  message("script finished without errors")
}


#' Perform Coloc analysis
#'
#' @param vcffile_exp the path to the exposure
#' @param vcffile_out the path to the outcome
#' @param chrompos the genetic region ex : "22:25115489-26127836"
#' @param ldref the ld reference panel
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return
#' @export

get_coloc <- function(vcffile_exp, vcffile_out, chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                      parameters = default_param()) {
  message("initializing coloc")
  stopifnot(length(vcffile_exp) == 1)
  dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out,
                                       parameters = parameters)
  vout <- GagnonMR::gwasvcf_to_coloc_eloi(vcffile_exp, vcffile_out,
                                          chrompos)
  if (is.null(vout)) {
    return(dt_null)
  }
  dt_res <- get_coloc_intern(vout = vout, dt_null = dt_null)
  return(dt_res)
}


#' Perform lead variant cis analysis
#'
#' @param vcffile_exp the path to the exposure .vcf.gz file
#' @param vcffile_out the path to the outcome .vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param ldref the path to the ldreference panel
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return a data.table of results
#' @export

get_uni_cis <-  function(vcffile_exp, vcffile_out, chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                         parameters = default_param()) {
  message("initializing uni-cis MR")
  stopifnot(gwasvcf::check_bcftools()&gwasvcf::check_plink())
  dat_tsmr <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = vcffile_exp, chrompos = chrompos, parameters = parameters)
  dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, parameters = parameters)
  rsid <- dat_tsmr[dat_tsmr$pval.exposure<parameters$uni_cis_minpval, ]$SNP
  if(!is.null(parameters$snp_bim)){rsid<-rsid[rsid%in%parameters$snp_bim]}
  dat_ld<-dat_tsmr[dat_tsmr$SNP%in%rsid,]
  if (dim(dat_ld)[1] == 0) {return(dt_null)}
  dat_ld <- dat_ld[,.SD[which.min(pval.exposure)], by = "id.exposure"]
  out_tsmr <- GagnonMR::extract_outcome_variant(snps =  dat_ld$SNP, outcomes = vcffile_out,
                                                ldref = ldref, parameters = parameters)

  if(ncol(out_tsmr)==2){
    k <- GagnonMR::convert_exposure_to_outcome(dat_ld)
    setDT(k)
    k[, c("id.outcome", "outcome")]<-NULL
    k[,]<-NA
    out_tsmr<-cbind(out_tsmr, k)
  }
  arguments<-tidyr::crossing(id.exposure = unique(dat_ld$id.exposure), id.outcome = unique(out_tsmr$id.outcome))

  res <- map(split(arguments, 1:nrow(arguments)), function(x)
    get_uni_cis_intern(dat_exposure = dat_ld[dat_ld$id.exposure==x$id.exposure],
                       dat_outcome = out_tsmr[out_tsmr$id.outcome==x$id.outcome,],
                       dt_null = dt_null[dt_null$id.exposure==x$id.exposure & dt_null$id.outcome==x$id.outcome])) %>%
    rbindlist(., fill = TRUE)
  res <- merge(res, dt_null, by = c("id.exposure", "exposure", "id.outcome", "outcome"), all = TRUE)

  return(res)
}

#' Perform pan (cis + trans) analysis
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param ldref the path to the ldreference panel
#' @param chrompos the genetic region (only a dummy variable to make consitent form with get_uni_cis)
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return a data.table of pan anaysis
#' @export

get_pan <- function(vcffile_exp, vcffile_out, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", chrompos = NA,
                    parameters = default_param()) {
  message("initializing pan mr")
  stopifnot(gwasvcf::check_bcftools()&gwasvcf::check_plink())
  chrompos<-NA
  dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, parameters = parameters)

  rsid <- map(vcffile_exp, function(x) {
    dat_vcf <- gwasvcf::query_gwas(vcf = x, pval = parameters$pan_clumping["pval"])
    if(dim(dat_vcf)[1]==0){return(NULL)}
    return(dat_vcf@assays@data@listData$ID %>% unlist)}) %>%
    unlist(.)
  if(is.null(rsid)){return(dt_null)}
  dat_tsmr <- map(vcffile_exp, function(x) {
    dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = x, vcffile_out = NULL, parameters = parameters)
    dat_vcf <- gwasvcf::query_gwas(vcf = x, rsid = rsid,  proxies = "no" )
    if(dim(dat_vcf)[1]==0){return(dt_null)}
    dat_tsmr <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "exposure") %>% data.table::as.data.table(.)
    dat_tsmr$id.exposure <- gsub(paste(parameters$path, collapse = "|"), "", x) %>% gsub("/.*$|.vcf.gz$", "", .)
    dat_tsmr <- GagnonMR::clump_data_local(dat_tsmr, ldref = ldref,
                                           clump_kb = parameters$pan_clumping["kb"],
                                           clump_r2 = parameters$pan_clumping["r2"])
    return(dat_tsmr)}) %>% data.table::rbindlist(., fill = TRUE)

  out_tsmr <- GagnonMR::extract_outcome_variant(snps = rsid, outcomes = vcffile_out, ldref = ldref, parameters = parameters)

  harm <- TwoSampleMR::harmonise_data(dat_tsmr, out_tsmr, action = 1)
  harm <- TwoSampleMR::add_rsq(harm)
  harm$id.exposureoutcome <- paste0(harm$id.exposure,harm$id.outcome)
  harm$fstat.exposure <- GagnonMR::fstat_fromdat(split(harm, harm$id.exposureoutcome))
  harm <- TwoSampleMR::steiger_filtering( harm )
  res_all <- map(split(harm, f = harm$id.exposureoutcome), function(x) {
    GagnonMR::all_mr_methods(x) %>%
      data.table::as.data.table(.)}) %>% data.table::rbindlist(., fill = TRUE)
  res_all[,c("lci", "uci", "type_of_test"):=NULL]
  k<-c("method", "nsnp", "b", "se", "pval")
  setnames(res_all, k, paste0("pan_", k))
  return(res_all)
}

#' Get reverse MR
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region (dummy variable) only there for consistency with get_uni_cis
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return a data.table of results
#' @export

get_reverseMR <- function(vcffile_exp, vcffile_out, chrompos = NA, parameters = default_param()) {
  message("initializin reverseMR")
  res <- GagnonMR::get_pan(vcffile_exp = vcffile_out, vcffile_out = vcffile_exp, parameters = parameters)
  k<-c("pan_method", "pan_nsnp", "pan_b", "pan_se", "pan_pval")
  setnames(res, k, paste0("rev", k))
  return(res)
}

#' pipeline clump at 0.6 and correct for ld matrix
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param ldref the path to the ldreferene panel
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return a data.table of the results
#' @export

get_multicis <- function(vcffile_exp, vcffile_out,  chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                         parameters = default_param()) {
  message("initializing multicis mr")
  stopifnot(gwasvcf::check_bcftools()&gwasvcf::check_plink())
  dat_tsmr <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = vcffile_exp, chrompos = chrompos, parameters = parameters)
  dt_null <- intern_dt_null(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, parameters = parameters)
  rsid <- dat_tsmr[dat_tsmr$pval.exposure<5e-8, ]$SNP
  dat_tsmr<-dat_tsmr[dat_tsmr$SNP%in%rsid,]
  if (dim(dat_tsmr)[1] == 0) {return(dt_null)}
  dat_tsmr <- GagnonMR::clump_data_local(dat_tsmr, ldref =ldref, clump_r2 = parameters$multicis_clumping["clump_r2"],
                                         clump_kb = parameters$multicis_clumping["clump_kb"],
                                         clump_p = parameters$multicis_clumping["clump_p"])

  out_tsmr <- extract_outcome_variant(snps =  dat_tsmr$SNP, outcomes = vcffile_out,
                                      ldref = ldref, parameters = parameters)

  ldmat <- ieugwasr::ld_matrix_local(unique(dat_tsmr$SNP), plink_bin = genetics.binaRies::get_plink_binary(), bfile = ldref)

  arguments<-  dt_null
  res_multicis <- lapply(split(arguments, 1:arguments[,.N]), function(x) {
    harm <- TwoSampleMR::harmonise_data(dat_tsmr[exposure==x$exposure, ], out_tsmr[outcome==x$outcome, ], 1)
    harm <- TwoSampleMR::steiger_filtering(harm) %>% as.data.table()
    if (nrow(harm) == 0) { return(x)}
    ndirectionalyinconsistent <- harm[, sum(steiger_dir == FALSE & steiger_pval < 0.05)]
    harm <- harm[!(steiger_dir == FALSE & steiger_pval < 0.05),]
    harm$id.exposure_outcome <- paste0(harm$id.exposure, "_", harm$id.outcome)

    if (nrow(harm) < 2) {return(x)}
    x<-harm
    mriobj <- MendelianRandomization::mr_input(bx = x$beta.exposure,
                                               bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome,
                                               exposure = x$exposure[1], outcome = x$outcome[1],
                                               snps = x$SNP, effect_allele = x$effect_allele.exposure,
                                               other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)

    k<-x[,c(paste(SNP, effect_allele.exposure, other_allele.exposure, sep = "_"),
            paste(SNP, other_allele.exposure, effect_allele.exposure, sep = "_"))]
    mriobj@correlation<-ldmat[k[k%in%rownames(ldmat)],k[k%in%colnames(ldmat)]]

    IVWobject <- MendelianRandomization::mr_ivw(mriobj, model = "default",
                                                correl = TRUE, distribution = "normal", alpha = 0.05)

    res_multicis <- cbind(x[1,c("id.exposure", "id.outcome", "exposure", "outcome")],
                          data.frame(beta.multi_cis = IVWobject@Estimate,
                                     se.multi_cis = IVWobject@StdError, pval.multi_cis = IVWobject@Pvalue,
                                     cochranQ.multi_cis = IVWobject@Heter.Stat[1], cochranQpval.multi_cis = IVWobject@Heter.Stat[2],
                                     nsnp.multi_cis = nrow(x), nsteigerfalse.multi_cis = ndirectionalyinconsistent))
    return(res_multicis)
  }) %>% rbindlist(., fill = TRUE)
  res_multicis <- merge(res_multicis, dt_null, by = c("id.exposure", "exposure", "id.outcome", "outcome"), all = TRUE)
  return(res_multicis)
}


#' Run all pqtl analyses with control on which analysis specifically to run
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param method_list the list of method to include "get_pan", "get_reverseMR", "get_uni_cis", "get_coloc", "get_multicis", "get_susie_coloc", "get_multicis_susie"
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return a data.table with all results
#' @export

run_all_pqtl_analyses <- function(vcffile_exp, vcffile_out,  chrompos,
                                  method_list = c("get_coloc", "get_uni_cis"),
                                  parameters = default_param()) {
  message(paste0("***********initilalizing all qtl analysis for ", vcffile_exp, " and ", vcffile_out, "************"))
  stopifnot(gwasvcf::check_bcftools()&gwasvcf::check_plink())
  if(is.na(chrompos)) {
    method_list <- method_list[sapply(method_list, function(x) (x %in% c("get_pan", "get_reverseMR")))]
  }

  if(length(method_list) == 0) {
    dt_null <- intern_dt_null(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, parameters = parameters)
    message(paste0("no suitable method for ", dt_null$id.exposure))
    return(dt_null)
  }

  vecexp<-c("get_uni_cis", "get_multicis", "get_pan")
  arguments<-tidyr::crossing(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out) %>% split(., 1:nrow(.))
  res <- lapply(method_list[!(method_list%in%vecexp)], function(meth) {
    map(arguments, function(x)
      get(meth)(vcffile_exp = x$vcffile_exp, vcffile_out = x$vcffile_out, chrompos = chrompos, parameters = parameters)) %>%
      rbindlist(., fill = TRUE)})

  res2 <- lapply(method_list[(method_list%in%vecexp)], function(meth) {
    get(meth)(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, chrompos = chrompos, parameters = parameters) })

  res_all <- c(res, res2) %>% purrr::reduce(merge, by = c("id.exposure", "id.outcome", "exposure", "outcome"), all = TRUE)
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
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return an object of class susie
#' @export
get_susie <- function(vcffile, chrompos,
                      ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                      parameters = default_param()) {
  dat <- GagnonMR::gwasvcf_to_finemapr_gagnon(region = chrompos, vcf = vcffile, bfile = ldref)
  L<- parameters$L

  fitted_rss <- vector(mode = "list", length = length(L))
  for(i in L) {
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
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return
#' @export

get_susie_coloc <-function(vcffile_exp, vcffile_out,  chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                           parameters=default_param()) {
  message("initializing susie coloc")
  sus1 <- GagnonMR::get_susie(vcffile = vcffile_exp, chrompos = chrompos, parameters = parameters)
  sus2 <- GagnonMR::get_susie(vcffile = vcffile_out, chrompos = chrompos, parameters = parameters)
  df_null <- data.table(  exposure = gwasvcf::query_gwas(vcffile_exp,  chrompos = "1:40000-80000")@metadata$header@samples,
                          outcome = gwasvcf::query_gwas(vcffile_out,  chrompos = "1:40000-80000")@metadata$header@samples,
                          susiecoloc.n.credibleset.exposure = length(sus1$sets$cs), susiecoloc.n.credibleset.outcome = length(sus2$sets$cs))
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
  susie_res_summary$susiecoloc.n.credibleset.exposure <- length(sus1$sets$cs)
  susie_res_summary$susiecoloc.n.credibleset.outcome <- length(sus2$sets$cs)
  return(susie_res_summary)

}

#' Same as get_multicis, but instead of selecting instruments based on LD clumping select with finemapping using SuSiE.
#'
#' @param vcffile_exp dah
#' @param vcffile_out dah
#' @param chrompos dah
#' @param ldref dah
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return select as instrument the SNP with the highest PIP in each credible set. If no SNP has a p-value under 5e-8 or number of credible sets < 2 return empty data.frame.
#' @export

get_multicis_susie <- function(vcffile_exp, vcffile_out,  chrompos, ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                               parameters = default_param()) {
  message("initializing multicis mr with susie finemapping")
  df_null <- data.frame(exposure= gwasvcf::query_gwas(vcffile_exp, chrompos = "1:40000-80000")@metadata$header@samples, outcome = gwasvcf::query_gwas(vcffile_out, chrompos = "1:40000-80000")@metadata$header@samples, nsnp.multi_cis_susie = 0)
  fitted_rss <- GagnonMR::get_susie(vcffile = vcffile_exp, chrompos = chrompos, parameters = parameters)
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
ndirectionalyinconsistent <- harm[, sum(steiger_dir == FALSE & steiger_pval < 0.05)]
harm <- harm[!(steiger_dir == FALSE & steiger_pval < 0.05),]

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


#' The default param for pQTL analyses
#'
#' @return the list of default param
#' @export
default_param <- function() {
  return(list(multicis_clumping = c(clump_kb = 1000, clump_r2 = 0.6, clump_p = 1),
              pan_clumping = c(pval = 5e-08, r2 = 0.01, kb = 1000),
              L = 1:10, minpval1 = TRUE, uni_cis_minpval = 5e-8, snp_bim = NULL,
              path = paste0("/mnt/sda/gagelo01/Vcffile/", c("MRBase", "Server"), "_vcf/")))
}


#' Get instrument
#'
#' @param vcffile the path to the vcffile
#' @param pval the minimum pvalue of the instruments
#' @param clump should clump TRUE or FALSE
#' @param r2 the maximum ld r2
#' @param kb the clumping kilobases
#' @param should_write should write (TRUE or FALSE) the instrument in the vcf directory.
#' @param parameters default_param(), that it is possible to change.
#'
#' @return
#' @export

get_inst <- function(vcffile, pval = 5e-08, clump = TRUE, r2 = 0.001,
                     kb = 10000, should_write = FALSE, parameters = default_param()) {
  stopifnot(gwasvcf::check_plink()&gwasvcf::check_bcftools())
  instfile <- gsub(".vcf.gz", paste0("_inst_rsquare", r2, "_kb", kb,".txt"), vcffile)
  if(file.exists(instfile)&clump == TRUE){
    dat_sign <- read.table(instfile, sep = "\t")
    data.table::setDT(dat_sign)
  } else {
    dat_vcf <- gwasvcf::query_gwas(vcf = vcffile, pval = pval)
    dat_sign <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "exposure") %>% data.table::as.data.table(.)
    dat_sign[,id.exposure := gsub(paste(parameters$path, collapse = "|"), "", vcffile) %>% gsub("/.*$|.vcf.gz$", "", .)]
    if(clump==TRUE) {
      dat_sign <- GagnonMR::clump_data_local(dat_sign, clump_kb = kb, clump_r2 = r2)
    }
  }
  if(should_write == TRUE)  {
    write.table(dat_sign, file = instfile, sep = "\t")
  }
  return(dat_sign)
}


#' extract outcome
#'
#' @param snps a vector of SNPs to extract
#' @param outcomes the path to the vcf file
#' @param proxies passed to the function gwasvcf::query_gwas
#' @param rsq passed to the function gwasvcf::query_gwas
#' @param ldref path to plink bed/bim/fam ld reference panel
#' @param parameters default_param()
#'
#' @return
#' @export

extract_outcome_variant <- function(snps,
                                    outcomes,
                                    proxies = "yes",
                                    rsq = 0.8,
                                    ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                                    parameters = default_param()) {

  stopifnot(gwasvcf::check_bcftools()&gwasvcf::check_plink())
  dt_null<-GagnonMR:::intern_dt_null(vcffile_exp = NULL, vcffile_out = outcomes, parameters = parameters)
  res <- map(outcomes, function(x) {
    out_vcf <- tryCatch(
      expr = {gwasvcf::query_gwas(vcf = x, rsid = unique(snps), proxies = proxies, bfile = ldref, tag_r2 = rsq) },
      error = function(e){return(matrix(nrow = 0, ncol = 1))})
    if (dim(out_vcf)[1] == 0) {return(dt_null)}
    out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.,
                                                             "outcome") %>% data.table::as.data.table(.)
    out_tsmr$id.outcome <-  gsub(paste(parameters$path, collapse = "|"), "", x) %>% gsub("/.*$|.vcf.gz$", "", .)
    return(out_tsmr)
  }) %>% data.table::rbindlist(., fill = TRUE)

  return(res)
}


#' Function that exits if you code takes 95 of total memory
#'
#' @return
#' @export
#'
#' @examples protect_memory()
protect_memory <- function() {
  test <- as.data.frame(system2("free", " --mega", stdout = TRUE))
  mem.col <- strsplit(test[1,], " +")[[1]]
  mem.total <- as.numeric(strsplit(test[2,], " +")[[1]][which(mem.col=='total')])
  mem.curr <- as.numeric(strsplit(test[2,], " +")[[1]][which(mem.col=='used')])
  mem.used <- as.numeric(sapply(strsplit(as.character(pryr::mem_used()), ' '), `[[`, 1))/1e06
  if(((mem.curr / mem.total) >= 0.9)){
    message(paste0("Warning : memory usage above ", 0.9*100 ," % . This script represents ", (mem.used/mem.curr)*100, "% of current memory usage."))
    if(((mem.curr / mem.total) >= 0.95) | ((mem.used/mem.curr) > opt$check_memory_lim)){
      stop("Memory usage too high. Stopping analysis.")
      quit(save = "no", status=1)
    }
  }
}


#' intern_dt_null
#'
#' @param vcffile_exp the path to the exposure vcf.gz file, can be a vector
#' @param vcffile_out the path to the outcome vcf.gz file, can be a vector
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return

intern_dt_null <- function(vcffile_exp, vcffile_out, parameters = default_param()) {
  exp <- data.table(id.exposure = if(is.null(vcffile_exp)){vcffile_exp}else{gsub(paste(parameters$path,  collapse = "|"), "", vcffile_exp) %>% gsub("/.*$|.vcf.gz","", .)},
                    exposure = if(is.null(vcffile_exp)){vcffile_exp}else{sapply(as.list(vcffile_exp), function(x) gwasvcf::query_gwas(x, chrompos = "1:40000-40000")@metadata$header@samples)})

  out <- data.table(id.outcome = if(is.null(vcffile_out)){vcffile_out}else{gsub(paste(parameters$path, collapse = "|"), "", vcffile_out) %>% gsub("/.*$|.vcf.gz", "", .)},
                    outcome = if(is.null(vcffile_out)){vcffile_out}else{sapply(as.list(vcffile_out), function(x) gwasvcf::query_gwas(x, chrompos = "1:40000-40000")@metadata$header@samples)})

  res <- if (is.null(vcffile_exp)) {
    out
  } else if (is.null(vcffile_out)) {
    exp
  } else {
    tidyr::crossing(exp, out)
  }

  return(data.table::as.data.table(res))

}


#' intern_vcfpath_to_TwoSampelMR_region
#'
#' @param vcf the path to the vcf.gz file, can be a vector
#' @param chrompos  the gene region e.g. 1:30000-40000
#' @param type either "exposure" or "outcome"
#' @param parameters parameters as obtained by GagnonMR::default_param()
#'
#' @return
#' @export
#'
#' @examples
intern_vcfpath_to_TwoSampelMR_region <- function(vcf, chrompos, type = "exposure", parameters = default_param()) {
  dat_tsmr <- map(as.list(vcf), function(x) {
    dt_null<-GagnonMR:::intern_dt_null(vcffile_exp = if(type=="exposure"){x}else{NULL},
                                       vcffile_out = if(type=="outcome"){x}else{NULL},
                                       parameters = parameters)
    dat_vcf <- tryCatch(
      expr = {gwasvcf::query_gwas(vcf = x, chrompos = chrompos)},
      error = function(e){return(matrix(nrow = 0, ncol = 1))})
    if (dim(dat_vcf)[1] == 0) {return(dt_null)}

    dat_tsmr <- gwasglue::gwasvcf_to_TwoSampleMR(dat_vcf, type = type)
    data.table::setDT(dat_tsmr)
    colnom<-paste0("id.", type)
    dat_tsmr[,(colnom) := gsub(paste(parameters$path, collapse = "|"), "", x) %>% gsub("/.*$|.vcf.gz$", "", .)]
    return(dat_tsmr)
  }) %>% rbindlist(., fill = TRUE)

  return(dat_tsmr)
}


#' Title
#'
#' @param dat_exposure TwoSampleMR exposure
#' @param dat_outcome TwoSampleMR outcome
#'
#' @return
#' @export

from_tsmr_to_coloc <- function(dat_exposure, dat_outcome) {
  type1 <- ifelse(all(is.na(dat_exposure$ncase.exposure)&is.na(dat_exposure$ncontrol.exposure)), "quant", "cc")
  type2 <- ifelse(all(is.na(dat_exposure$ncase.outcome)&is.na(dat_exposure$ncontrol.exposure)), "quant", "cc")

  harm <- TwoSampleMR::harmonise_data(dat_exposure, dat_outcome, action = 1)

  out1 <- data.table(pvalues = harm$pval.exposure,
                     N = harm$samplesize.exposure,
                     MAF = harm$eaf.exposure,
                     beta = harm$beta.exposure,
                     varbeta = (harm$se.exposure)^2,
                     type = type1,
                     snp = harm$SNP,
                     id = harm$exposure[1])

  out2 <- data.table(pvalues = harm$pval.outcome,
                     N = harm$samplesize.outcome,
                     MAF = harm$eaf.outcome,
                     beta = harm$beta.outcome,
                     varbeta = (harm$se.outcome)^2,
                     type = type2,
                     snp = harm$SNP,
                     id = harm$outcome[1])

  if (type1 == "cc") {
    out1$s <- mean(tab1$NC/tab1$SS, na.rm = TRUE)
  }
  if (type2 == "cc") {
    out2$s <- mean(tab2$NC/tab2$SS, na.rm = TRUE)
  }

  vout<-list(dataset1 = out1, dataset2 = out2)
  return(vout)
}


#' Title
#'
#' @param vout the object that from_tsmr_to_coloc returns
#' @param dt_null the object that intern_dt_null returns
#'
#' @return
#' @export
#' @example

get_coloc_intern <- function(vout, dt_null) {

  vout[[1]]$MAF <- vout[[1]]$MAF %>% ifelse(. == 0, 0.001,
                                            .) %>% ifelse(. == 1, 0.999, .)
  vout[[2]]$MAF <- vout[[2]]$MAF %>% ifelse(. == 0, 0.001,
                                            .) %>% ifelse(. == 1, 0.999, .)
  vres <- coloc::coloc.abf(vout[[1]], vout[[2]])
  pp <- as.data.table(as.list(vres$summary[2:6]))
  colnames(pp) <- paste0("posprob_coloc_PPH", 0:4)
  vresres <- vres$results %>% as.data.table(.)
  df_res1 <- data.frame(exposure = vout$dataset1$id[1], outcome = vout$dataset2$id[1])
  df_res2 <- data.frame(posprob_colocH4.SNP = vresres[which.max(SNP.PP.H4),
                                                      snp], posprob_colocH4.SNPexplained_var = vresres[which.max(SNP.PP.H4),
                                                                                                       SNP.PP.H4])
  dt_res <- do.call(cbind, list(df_res1, pp, df_res2)) %>%
    data.table::as.data.table(.)
  dt_res <- merge(dt_null, dt_res, by = c("exposure", "outcome"))
  return(dt_res)
}


#' get_uni_cis_intern, this function cannot receive several exposure and outcome at once
#'
#' @param dat_exposure only top SNP one exposure at one
#' @param dat_outcome must contain all the SNP in exposure out outcome at once
#' @param dt_null the object that intern_dt_null returns
#'
#' @return
#' @export
#'
#' @examples
get_uni_cis_intern <- function(dat_exposure, dat_outcome, dt_null) {
  dt_null[, lead_snp.wald := dat_exposure$SNP[1]]
  #harmonisation
  harm <- TwoSampleMR::harmonise_data(dat_exposure, dat_outcome, action = 1)
  harm <- TwoSampleMR::add_rsq(harm)
  harm$id.exposureoutcome <- paste0(harm$id.exposure,harm$id.outcome)
  harm$fstat.exposure <- GagnonMR::fstat_fromdat(split(harm, harm$id.exposureoutcome))
  harm <- TwoSampleMR::steiger_filtering( harm )
  harm <- harm[harm$mr_keep == TRUE,]
  if(nrow(harm)==0){return(dt_null)}
  res_all <- GagnonMR::all_mr_methods(harm) %>% data.table::as.data.table(.)
  if (dim(res_all)[1] == 0) {return(dt_null)}
  res_all[, c("method", "type_of_test", "lci", "uci") := NULL]
  data.table::setnames(res_all, old = c("nsnp", "b",  "se", "pval"), new = paste0(c("lead_snp", "b",  "se", "pval"), ".wald"))

  col_to_select<-c("pval.exposure", "steiger_pval", "steiger_dir", "chr.exposure", "pos.exposure", "rsq.exposure", "fstat.exposure")
  res_all <- cbind(res_all, harm[,col_to_select])
  data.table::setnames(res_all, col_to_select, paste0(col_to_select, ".wald"))
  dt_null$lead_snp.wald<-NULL
  res_all <- merge(res_all, dt_null, by = c("id.exposure", "exposure", "id.outcome", "outcome"), all = TRUE)
  return(res_all)
}


