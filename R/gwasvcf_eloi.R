#' same as the gwasvcf function except fix bug
#'
#' @param vcf1 dah
#' @param vcf2 dah
#' @param chrompos dah
#'
#' @return dah
#' @export

gwasvcf_to_coloc_eloi <- function (vcf1, vcf2, chrompos) {
  o <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
  vcf1 <- o[[1]]
  vcf2 <- o[[2]]
  if (length(vcf1) == 0) {
    message("No overlaps in vcf1")
    return(NULL)
  }
  if (length(vcf2) == 0) {
    message("No overlaps in vcf2")
    return(NULL)
  }
  stopifnot(length(vcf1) == length(vcf2))
  tab1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  tab2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  index <- as.character(tab1$REF) == as.character(tab2$REF) &
    as.character(tab1$ALT) == as.character(tab2$ALT) & as.character(tab1$seqnames) ==
    as.character(tab2$seqnames) & tab1$start == tab2$start
  stopifnot(sum(index) > 0)
  type1 <- ifelse(VariantAnnotation::header(vcf1) %>% VariantAnnotation::meta() %>%
                    {
                      .[["SAMPLE"]][["StudyType"]]
                    } == "Continuous", "quant", "cc")
  type2 <- ifelse(VariantAnnotation::header(vcf2) %>% VariantAnnotation::meta() %>%
                    {
                      .[["SAMPLE"]][["StudyType"]]
                    } == "Continuous", "quant", "cc")
  tab1$AF[is.na(tab1$AF)] <- 0.5
  tab2$AF[is.na(tab2$AF)] <- 0.5
  out1 <- tab1[index, ] %>% {
    list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES,
         varbeta = .$SE^2, type = type1, snp = names(vcf2)[index],
         z = .$ES/.$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf1))[1])
  }
  out2 <- tab2[index, ] %>% {
    list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES,
         varbeta = .$SE^2, type = type2, snp = names(vcf2)[index],
         z = .$ES/.$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf2))[1])
  }
  if (type1 == "cc") {
    out1$s <- mean(tab1$NC/tab1$SS, na.rm = TRUE)
  }
  if (type2 == "cc") {
    out2$s <- mean(tab2$NC/tab2$SS, na.rm = TRUE)
  }
  return(list(dataset1 = out1, dataset2 = out2))
}



#' Title dah
#'
#' @param region dah
#' @param vcf dah
#' @param bfile dah
#' @param plink_bin dah
#'
#' @return
#' @export

gwasvcf_to_finemapr_gagnon <- function (region, vcf, bfile, plink_bin = genetics.binaRies::get_plink_binary())
{
  message("Extracting data from vcf")
  ext <- gwasvcf::query_gwas(vcf = vcf, chrompos = region)
  out <- lapply(unique(region), function(i) {
    message(i)
    m <- list()
    temp <- gwasvcf::query_gwas(vcf = ext, chrompos = i)
    m[["ld"]] <- ieugwasr::ld_matrix(names(temp), bfile = bfile,
                                     plink_bin = plink_bin, with_alleles = FALSE) %>%
      gwasglue:::greedy_remove()
    tib <- gwasvcf::vcf_to_tibble(temp)
    m[["z"]] <- tib %>% subset(rsid %in% rownames(m[["ld"]])) %>%
      dplyr::mutate(z = ES/SE) %>% dplyr::select(snp = rsid,
                                                 zscore = z)
    m[["n"]] <- tib[["SS"]]
    return(m)
  })
  class(out) <- "FinemaprList"
  return(out)
}

