#' Title
#'
#' @param vcffile_exp the path to the exposure vcf.gz file
#' @param vcffile_out the path to the outcome vcf.gz file
#' @param chrompos the gene region e.g. 1:30000-40000
#' @param parameters Parameters to be used for various cis-MR methods. Default is output from default_param.
#'
#' @return
#' @export

get_hyprcoloc <- function(vcffile_exp, vcffile_out = NULL, chrompos,
                          parameters = default_param()) {
  list_vcf <- as.list(c(vcffile_exp, vcffile_out))
  o <- gwasvcf::vcflist_overlaps(vcflist = list_vcf, chrompos = chrompos)
  dtnull <- lapply(list_vcf, function(x) {
    GagnonMR:::intern_dt_null(vcffile_exp = x, vcffile_out = NULL,
                              parameters = parameters)
  }) %>% data.table::rbindlist(., fill = TRUE)
  dt <- lapply(o, gwasglue::gwasvcf_to_TwoSampleMR) %>% data.table::rbindlist(.,
                                                                              fill = TRUE)
  data.table::setDT(dt)
  dt[,chr.exposure:=chr.exposure%>%as.character(.)%>%as.integer(.)]
  dtnull[, todump := dt$id.exposure %>% unique]
  dt <- merge(dtnull[,.(id.exposure, todump)], dt, by.x = "todump", by.y = "id.exposure")
  dt[, todump := NULL]
  k <- GagnonMR::prepare_for_mvmr(dt, dt, harmonise_strictness = 1,
                                  should_clump = FALSE)
  res <- GagnonMR::run_hypr_on_aligned(k)
  return(list(res = res, align = k))
}

#' Title
#'
#' @param df_aligned A data.table of exposure harmonised as outputted by prepare_for_mvmr or get_hyprcoloc
#' @param res_hypr1 The results of hyprcoloc as outputted by get_hyprcoloc
#' @param ldref the path to the ld reference panel either in grch37 or grch38
#' @param traits_inorder a vector of the name of the name of the exposure to stipulate the order of the exposure for the plot
#' @param build must be either 37 or 38
#'
#' @return
#' @export

stack_assoc_plot_wrapper <- function(df_aligned,
                                     res_hypr1,
                                     ldref = default_param()$ldref,
                                     traits_inorder = unique(df_aligned$exposure),
                                     build = 37) {
  stopifnot(res_hypr1[,.N==1])
  stopifnot(build %in% c(37))
  stopifnot("gassocplot" %in% (.packages()))
  stopifnot(is.integer(df_aligned$chr.exposure)&is.integer(df_aligned$chr.exposure))

  # stopifnot(is.factor(df_aligned$exposure))

  df_order <- data.frame(exposure=rev(traits_inorder), var_order = c(1:length(traits_inorder)))
  df_aligned <- merge(df_aligned, df_order, by = "exposure")
  df_aligned <- df_aligned[order(var_order), ]
  df_aligned[,z:=beta.exposure/se.exposure]
  df_reshaped <- reshape(df_aligned, idvar = c("SNP", "chr.exposure", "pos.exposure"), timevar = "exposure", direction = "wide")

  ldmat <- ieugwasr::ld_matrix_local(
    df_reshaped$SNP,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = ldref,
    with_alleles = FALSE)

  top_snp <- res_hypr1[, candidate_snp]
  if(!(top_snp %in% rownames(ldmat))) {
    ldmat<-NULL} else {
      df_reshaped <- df_reshaped[SNP %in% rownames(ldmat), ]
      ldmat <- ldmat[df_reshaped$SNP, df_reshaped$SNP] #order in the same way as the marker data.frame
    }
  markers <- df_reshaped[, .(SNP, chr.exposure, pos.exposure)]
  data.table::setnames(markers, colnames(markers), c("marker", "chr", "pos"))

  colnom <- colnames(df_reshaped)[grepl("^z\\.", colnames(df_reshaped))]
  zscores<-as.matrix(df_reshaped[, .SD, .SDcols = colnom])


  res <- gassocplot::stack_assoc_plot(markers = markers,
                                      z = zscores,
                                      corr = ldmat,
                                      traits=rev(traits_inorder),
                                      top.marker= top_snp)

  return(res)

}

#' Title
#'
#' @param df_aligned the data.frame aligned as outputed by get_hyprcoloc()$align
#' @param traits_inorder a vector of the name of the name of the exposure to stipulate the order of the exposure for the plot
#'
#' @return
#' @export

sensitivity.plot_wrapper <- function(df_aligned, traits_inorder =  levels(df_aligned$exposure)) {
  df_aligned[, exposure := factor(exposure, levels = traits_inorder)]
  df_aligned <- df_aligned[order(exposure), ]
  df_reshaped <- reshape(df_aligned, idvar = c("SNP", "chr.exposure", "pos.exposure"), timevar = "exposure", direction = "wide")
  effect_est <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("beta.exposure", names(df_reshaped))]] %>% as.matrix
  effect_se <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("^se.exposure", names(df_reshaped))]] %>% as.matrix
  res<- hyprcoloc::sensitivity.plot(effect.est = effect_est,
                                    effect.se = effect_se,
                                    trait.names = traits_inorder,
                                    snp.id = df_reshaped$SNP,
                                    similarity.matrix = TRUE)
  return(res)
}


#' Title
#'
#' @param mat_cor a matrix of correlation as outputted by sensitivity.plot_wrapper
#' @param should_hclust if true will perform a dendrogramm. If false will not.
#'
#' @return
#' @export

drawheatmap<-function(mat_cor, should_hclust = TRUE) {

  levels <- colnames(mat_cor)
  heat <- as.data.frame(mat_cor)
  heat$row <- rownames(heat)
  rownames(heat) <- NULL
  data.table::setDT(heat)
  heat <- data.table::melt(heat, id.vars = "row")
  heat[, `:=`(row, factor(row, levels = rev(levels)))]
  heat[, `:=`(variable, factor(variable, levels = rev(levels)))]

  if(should_hclust == TRUE) {
    otter_dendro <- as.dendrogram(hclust(d = dist(x = mat_cor)))
    otter_order <- order.dendrogram(otter_dendro)
    heat[,row:=factor(row, levels = row[otter_order], ordered = TRUE)]
    heat[,variable:=factor(variable, levels = row[otter_order], ordered = TRUE)]
    heat<- heat[order(row,variable)]
  }

  g <- ggplot2::ggplot(heat, aes(x = variable, tissue, y = row, fill = value))  +
    geom_tile() +
    # scale_fill_gradient(low = "lightblue", high = "blue3",limits=c(0,1)) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF",limits=c(0,1),name = "Colocalisation rate") +
    labs(fill = "") +
    theme(
      panel.background = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
      legend.position = "top",
      legend.title = element_text(
        color = "gray20",
        size = 12
      ),
      legend.text = element_text(
        color = "gray20",
        size = 10
      ),
      legend.title.align = 0.5,
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      axis.title = element_blank(),
      axis.line = element_line(size = 1, colour = "gray20"),
      axis.ticks = element_line(size = 1, colour = "gray20"),
      # axis.text.y = element_text(
      #   size = 10,
      #   colour = "gray20"
      # ),
      axis.text.y=element_blank(),
      axis.text.x = element_text(
        angle = 60,
        size = 10,
        hjust = 1,
        colour = "gray20"
      ),
      axis.ticks.length = unit(.25, "cm"))

  g
}

