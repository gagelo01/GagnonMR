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


#' provide vcf path and get coloc input
#'
#' @param vcf1 dah
#' @param vcf2 if vcf2 is NULL will simply format vcf1 in the coloc format which is usefull for CoPheScan
#' @param chrompos  dah
#' @param parameters dah
#'
#' @return
#' @export
#'

gwasvcf_to_coloc_input <- function (vcf1, vcf2 = NULL, chrompos, parameters = default_param()) {
  if(is.null(vcf2)) {
    o <-  list(gwasvcf::query_gwas(vcf = vcf1, chrompos = chrompos))
  } else {
    o <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
    stopifnot(length(o[[1]]) == length(o[[2]]))}

  if (any(map_lgl(o, function(x) length(x)==0))) {
    message("No SNPs")
    return(NULL)
  }

  dtnull <- lapply(list(vcf1,vcf2), function(x) {
    GagnonMR:::intern_dt_null(vcffile_exp = x, vcffile_out = NULL,
                              parameters = parameters)
  }) %>% data.table::rbindlist(., fill = TRUE)


  dt <- lapply(o, gwasglue::gwasvcf_to_TwoSampleMR) %>% data.table::rbindlist(.,
                                                                              fill = TRUE)
  data.table::setDT(dt)
  dt[,chr.exposure:=chr.exposure%>%as.character(.)%>%as.integer(.)]
  dtnull[, todump := dt$id.exposure %>% unique]
  dt <- merge(dtnull[, .(id.exposure, todump)], dt, by.x = "todump",
              by.y = "id.exposure", sort = FALSE)
  dt[, todump := NULL]

  if(!is.null(vcf2)) {
    k <- GagnonMR::prepare_for_mvmr(dt, dt, harmonise_strictness = 1,
                                    should_clump = FALSE)} else{k<-dt}
  k <- split(dt, by = "id.exposure")


  type <-  map(o, function(x) ifelse(VariantAnnotation::header(x) %>% VariantAnnotation::meta() %>%
                                       {
                                         .[["SAMPLE"]][["StudyType"]]
                                       } == "Continuous", "quant", "cc"))

  out<-   map(1:length(o), function(i) {
    out <- k[[i]] %>% {
      list(pvalues = .$pval.exposure, N = .$samplesize.exposure,
           MAF = .$eaf.exposure %>% ifelse(.<0.5, ., 1-.),
           beta = .$beta.exposure,
           varbeta = .$se.exposure^2, type = type[[i]], snp = .$SNP,
           z = .$beta.exposure/.$se.exposure, chr = .$chr.exposure,
           pos = .$pos.exposure, id = .$id.exposure)
    }

    if (type[[i]] == "cc") {
      out$s <- mean(k[[i]]$ncase.exposure/k[[i]]$samplesize.exposure, na.rm = TRUE)
    }
    return(out)
  })
  # names(out)<-c("dataset1", "dataset2")
  return(out)
}


