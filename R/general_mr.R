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

#' In TwoSampleMR format convert exposure to outcome
#'
#' @param exposure_dat A TwoSampleMR data.frame
#'
#' @return
#' @export
#'
#' @examples
convert_exposure_to_outcome <- function (exposure_dat) {
  id <- subset(exposure_dat, !duplicated(exposure), select = c(exposure,
                                                               id.exposure))
  outcome_dat <- TwoSampleMR::format_data(exposure_dat, beta_col = "beta.exposure", type = "outcome",
                                          se_col = "se.exposure", pval_col = "pval.exposure", phenotype_col = "exposure",
                                          effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",
                                          eaf_col = "eaf.exposure", units_col = "units.exposure")
  outcome_dat <- base::merge(outcome_dat, id, by.x = "outcome", by.y = "exposure")
  outcome_dat <- subset(outcome_dat, select = -c(id.outcome))
  names(outcome_dat)[names(outcome_dat) == "id.exposure"] <- "id.outcome"
  return(outcome_dat)
}

