#' Extract outcome and proxy from an output of format_data("outcome") object
#'
#' @param snps The snps (RSID) to extract
#' @param phenotype A complete summarystatistics object that was obtained from format_data("outcome") or format_data("outcome").
#' @param type whether the summarystatistics object is an outcome or an exposure
#' @param proxies if TRUE will attempt to find proxy if FALSE will not.
#' @param rsq Minimum LD rsq value (if proxies = TRUE). Default = 0.8.
#' @param combined_query A combined_querry data.frame from LDproxy_batch(snp = snps,  token = Sys.getenv("LDLINK_TOKEN"), append = TRUE)
#' if NULL it will automatically estimate it. But this step is long. So if you already have it or you plan using the function often it is
#' more efficient to first to save the combined_querry object and use it each time. I personnaly always use a combined_query that I already saved.
#'
#' @return Dataframe of summary statistics for all available outcomes in the format recognised by TwoSampleMR
#' @export

extract_and_proxy_fromsumstat <- function (snps, phenotype, type = "outcome", proxies = TRUE,
                                           rsq = 0.8, combined_query = NULL)
{
  snps_toproxy <- snps[!(snps %in% phenotype$SNP)]
  proxies <- ifelse(length(snps_toproxy) == 0, FALSE, proxies)
  if (proxies) {
    if (is.null(combined_query)) {
      if (file.exists("combined_query_snp_list.txt")) {
        file.remove("combined_query_snp_list.txt")
      }
      LDlinkR::LDproxy_batch(snp = snps_toproxy, token = Sys.getenv("LDLINK_TOKEN"),
                             append = TRUE)
      combined_query <- fread("combined_query_snp_list.txt")
    }
    combined_query <- combined_query[query_snp %in% snps_toproxy,
    ]
    best_proxies <- combined_query[, .SD[min(which(RS_Number %in%
                                                     phenotype$SNP))], by = query_snp, .SDcols = c("RS_Number",
                                                                                                   "R2", "Distance", "Alleles", "Correlated_Alleles")][R2 >=
                                                                                                                                                         rsq, ]
    best_proxies <- merge(best_proxies, phenotype[, .SD,
                                                  .SDcols = c("SNP", paste0(c("effect_allele.", "other_allele."),
                                                                            type))], by.x = "RS_Number", by.y = "SNP", all = FALSE)
    best_proxies <- tidyr::separate(best_proxies, col = "Correlated_Alleles",
                                    into = c("targetV1", "proxyV1", "targetV2", "proxyV2"))
    if (nrow(best_proxies) == 0) {
      message("no available proxy in sumstat")
      outcome_to_return <- phenotype[SNP %in% snps, ]
      return(outcome_to_return)
    }
    else {
      part_phenotype_proxy <- data.frame(V1 = TRUE, V2 = best_proxies$query_snp,
                                         V3 = best_proxies$RS_Number, V4 = best_proxies[,
                                                                                        get(paste0("effect_allele.", type))], V5 = best_proxies[,
                                                                                                                                                get(paste0("other_allele.", type))], V6 = best_proxies[,
                                                                                                                                                                                                       ifelse(targetV1 == get(paste0("effect_allele.",
                                                                                                                                                                                                                                     type)), proxyV1, proxyV2)], V7 = best_proxies[,
                                                                                                                                                                                                                                                                                   ifelse(targetV1 == get(paste0("other_allele.",
                                                                                                                                                                                                                                                                                                                 type)), proxyV1, proxyV2)], V8 = best_proxies$R2,
                                         V9 = best_proxies$Distance)
      colnames(part_phenotype_proxy) <- paste0(c("proxy.",
                                                 "target_snp.", "proxy_snp.", "target_a1.", "target_a2.",
                                                 "proxy_a1.", "proxy_a2.", "proxy_R2.", "proxy_distance."),
                                               type)
      setDT(part_phenotype_proxy)
      mirge_proxy <- merge(phenotype[SNP %in% part_phenotype_proxy[,
                                                                   get(paste0("proxy_snp.", type))], ], part_phenotype_proxy,
                           by.x = "SNP", by.y = paste0("proxy_snp.", type),
                           all = FALSE)
      mirge_proxy[, `:=`(SNP, get(paste0("target_snp.",
                                         type)))]
      mirge_proxy[, `:=`(paste0("effect_allele.", type),
                         get(paste0("target_a1.", type)))]
      mirge_proxy[, `:=`(paste0("other_allele.", type),
                         get(paste0("target_a2.", type)))]
      outcome_to_return <- plyr::rbind.fill(phenotype[SNP %in%
                                                        snps, ], mirge_proxy)
    }
  }
  else {
    outcome_to_return <- phenotype[SNP %in% snps, ]
  }
  return(outcome_to_return)
}

