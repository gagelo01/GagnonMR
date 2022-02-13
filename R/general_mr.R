#' get ld from ldlink
#'
#' @param SNP the list of rsids you want to have the proxies from
#'
#' @return the data.frame of proxies
#' @export

LD_proxy_wrapper <- function(SNP) {

  if (file.exists("combined_query_snp_list.txt")) {
    file.remove("combined_query_snp_list.txt")
  }


  LDlinkR::LDproxy_batch(snp = unique(SNP), token = Sys.getenv("LDLINK_TOKEN"),
                         append = TRUE)
  combined_query <- fread("combined_query_snp_list.txt")

return(combined_query)
}

#' Title
#'
#' @param gwas_file dah
#' @param snp_col dah
#' @param outcome_name dag
#' @param beta_col dah
#' @param se_col dah
#' @param pval_col dah
#' @param eaf_col dah
#' @param effect_allele_col dah
#' @param other_allele_col dah
#' @param ncase_col dah
#' @param ncontrol_col dah
#' @param samplesize_col dah
#' @param units_col dah
#' @param prevalence_col dag
#'
#' @return
#' @export

fromdisease_tooutcome <- function (gwas_file, snp_col, outcome_name,
                                   beta_col, se_col, pval_col, eaf_col, effect_allele_col, other_allele_col,
                                   ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                                   units_col, prevalence_col)
{

  disease <- fread(gwas_file)
  disease[, Phenotype := outcome_name]
  disease[, ncase := as.numeric(ncase_col)]
  disease[, ncontrol := as.numeric(ncontrol_col) ]
  disease[, samplesize := as.numeric(samplesize_col)]
  disease[, units := "log odds"]

  disease_outcome <- TwoSampleMR::format_data(disease,
                                 type = "outcome",
                                 snp_col = snp_col,
                                 beta_col = beta_col,
                                 se_col = se_col,
                                 pval_col = pval_col,
                                 eaf_col = eaf_col,
                                 effect_allele_col = effect_allele_col,
                                 other_allele_col = other_allele_col)
  setDT(disease_outcome)
  disease_outcome[ , prevalence.outcome := prevalence_col]
  return(disease_outcome)
}



#' Title remove SNPs close to pleiotropic gene region
#'
#' @param inst you instrument data.frame
#' @param to_exclude the name of gene to exclude
#' @param window your window you want to remove around the gene
#'
#' @return
#' @export
#'
#' @examples
#'
remove_gene_region <- function(inst, to_exclude = c("APOE", "ABO", "HLA-A"), window = 2e6) {

gencode <- as.data.table(fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt"))

list <- vector(mode = "list", length = length(to_exclude))

for(i in 1:length(to_exclude)) {
  bon <- gencode[gene_name == to_exclude[i],]
  list[[i]] <- data.frame(chr = bon[1,]$chr,
                          start = min(bon$start) - window/2,
                          end = max(bon$end) + window/2,
                          gene_name = bon[1,]$gene_name)

}

region_df <- rbindlist(list)

for(i in 1:nrow(region_df)) {
  inst <- inst[!((chr.exposure == region_df[i,]$chr) &
                             (pos.exposure >= region_df[i,]$start) &
                             (pos.exposure <= region_df[i,]$end)), ]

}
return(inst)
}

#multi-cis

