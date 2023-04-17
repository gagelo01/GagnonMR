#' This function take as input a data.frame of primary MR analysis and return a balloon plot. credit Nicholas.
#'
#' @param dat_plot a data.frame containing the column "exposure" containing all the name of your exposure (or outcomes) as factor. the names will appear in the order you specified in the leevels of your factor from the bottom to up.
#' containing the column "outcome"  containing all the name of your outcomes (or exposures) as factor. the names will appear in the order you specified in the leevels of your factorfrom left to right.
#' containing the column "pval" a numeric column containing the p-value of the primary mR causal estimates.
#' containing the column "z_score" a numeric column containing the z-score of the primary MR causal estimates
#' containing the column "Category"  a factor column the name of the category. the names will appear in the order you specified in the levels of your factor from left to right on top of the graph
#' If the function throws an error it is because you did not specified dat_plot accordingly so read this section carefully.
#' @param bonferroni_threshold the pvalue threshold of significance you want to use.
#' @param subdivide_by_category if TRUE will divide by category. If FALSE won't
#' @return A Ballon plot of the primary MR analysis.
#' @export
plot_balloon <- function(dat_plot, bonferroni_threshold = 0.05, subdivide_by_category = TRUE) {

  dat_plot$pval = as.numeric(dat_plot$pval)
  dat_plot$log10_pval = -log10(dat_plot$pval)
  dat_plot$shape_point = sapply(dat_plot$pval, FUN = function(x) {ifelse( x < bonferroni_threshold, "rond", "Non-significant")})

  dat_plot_rond = dat_plot
  dat_plot_rond$z_score[which(dat_plot_rond$pval >= bonferroni_threshold)] = NA
  dat_plot_rond$log10_pval[which(dat_plot_rond$pval >= bonferroni_threshold)] = NA
  dat_plot_rond$shape_point = sapply(dat_plot_rond$pval, FUN = function(x) {ifelse( x < bonferroni_threshold, "rond", NA)})


  dat_plot_croix = dat_plot
  dat_plot_croix$z_score[which(dat_plot_croix$pval < bonferroni_threshold)] = NA
  dat_plot_croix$log10_pval[which(dat_plot_croix$pval < bonferroni_threshold)] = NA
  dat_plot_croix$shape_point = sapply(dat_plot_croix$pval, FUN = function(x) {ifelse( x >= bonferroni_threshold, "Non-significant", NA)})

  ordered_category = dat_plot %>% dplyr::arrange(Category) %>% dplyr::select(Category) %>% dplyr::distinct()
  breaks_category = sapply(unique(dat_plot$Category), FUN = function(x){dim(dat_plot %>% dplyr::select(outcome , Category) %>% dplyr::distinct() %>% dplyr::filter(Category ==  !!x))[1]})

  for (i in 1:length(ordered_category$Category[!is.na(ordered_category$Category)])) {
    if (i == 1)
    {
      list_breaks = as.numeric(breaks_category[i])
      label_breaks = as.numeric(breaks_category[i])/2
    } else {
      list_breaks = c(list_breaks, (data.table::last(list_breaks) + as.numeric(breaks_category[i])))
      label_breaks = c(label_breaks, (list_breaks[i - 1] + (as.numeric(breaks_category[i])/2 + 0.5)))
    }
  }

  balloon_plot = ggplot2::ggplot() +
    ggplot2::geom_point(data = dat_plot_croix, ggplot2::aes(x = outcome, y = exposure, shape = factor(shape_point)), size = 2, color = "gray20")

  when_annotate = seq(from = 1 , to  = length(ordered_category$Category[!is.na(ordered_category$Category)]), by = 2)

  cpt = 0
  for (i in c(1:length(ordered_category$Category[!is.na(ordered_category$Category)]))) {
    cpt = cpt + 1
    if (i == 1)
    {
      list_breaks_background = as.numeric(breaks_category[i]) + 0.5
      balloon_plot = balloon_plot +
        ggplot2::annotate(
          "rect",
          xmin = 0,
          xmax = list_breaks_background[i],
          ymin = 0,
          ymax = length(unique(dat_plot$exposure)) + 0.5,
          fill = "gray5",
          alpha = 0.1
        )
    } else {
      list_breaks_background = c(list_breaks_background, (data.table::last(list_breaks_background) + as.numeric(breaks_category[i])))
      if (cpt %in% when_annotate)
      {
        balloon_plot = balloon_plot +
          annotate(
            "rect",
            xmin = list_breaks_background[i - 1],
            xmax = list_breaks_background[i],
            ymin = 0,
            ymax = length(unique(dat_plot$exposure)) + 0.5,
            fill = "gray5",
            alpha = 0.1
          )
      }
    }
  }


  balloon_plot = balloon_plot +
    ggplot2::geom_point(data = dat_plot_rond, ggplot2::aes(x = outcome, y = exposure, size = log10_pval, color = z_score)) +
    ggplot2::scale_color_gradient2(name = "Z-score",
                          low = scales::muted("#5884E5"),
                          mid = "white",
                          high = scales::muted("#9E131E"),
                          midpoint = 0
    ) +
    ggplot2::scale_shape_manual(name = "", values = c(4,1)) +
    ggplot2::scale_size(name = expression(-Log[10](P)), range = c(4,7)) +
    ggplot2::coord_fixed(clip = "off", ratio = 1) +
    ggplot2::guides(size = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    ggplot2::theme(
      panel.grid.major.y = element_line(size = 0.25, colour = "gray60"),
      panel.grid.major.x = element_line(size = 0.25, colour = "gray60"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 2, r = 0.5, b = 0.5, l = 0.5, "cm"),
      legend.position = "right",
      legend.text = element_text(
        color = "gray20",
        size = 10,
        margin = margin(l = 0.2, r = 0.2)
      ),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      axis.title = element_blank(),
      axis.line = element_line(size = 0.5, colour = "gray20"),
      axis.ticks = element_line(size = 0.5, colour = "gray20"),
      axis.text.y = element_text(
        size = 10,
        colour = "gray20"
      ),
      axis.text.x = element_text(
        angle = 60,
        size = 8,
        hjust = 1,
        face = "plain",
        colour = "gray20"
      ))

  if(subdivide_by_category){
    for (i in 1:length(ordered_category$Category[!is.na(ordered_category$Category)])) {

      balloon_plot = balloon_plot + ggplot2::annotation_custom(grob =  grid::textGrob(label = ordered_category$Category[!is.na(ordered_category$Category)][i],
                                                                       hjust = 0,
                                                                       vjust = 0.5,
                                                                       rot = 18,
                                                                       gp = grid::gpar(cex = 0.75, fontface = "bold")),
                                                      ymin = length(unique(dat_plot$exposure)) + 1,
                                                      ymax = length(unique(dat_plot$exposure)) + 1,
                                                      xmin = label_breaks[i],
                                                      xmax = label_breaks[i])
    }
  }

  print(balloon_plot)

}


#' plot a scatter plot that is publication ready
#'
#' @param mr_results results from GagnonMR::all_mr_methods
#' @param dat results from TwoSampleMR::harmonise_data
#' @param equation_facet_grid an equation to facet
#' @param legend.position where to put the legend
#'
#' @return
#' @export

my_mr_scatter_plot <- function (mr_results, dat,
                                equation_facet_grid = "", legend.position = c(.32,.89) ) { #equation to parse into facet grid example includes "outcome ~  align", " ~  align", and "", if we do not want to have a grid

  d <- dat
  if(!("align" %in% colnames(d))) { d$align<-""}
  if(!("align" %in% colnames(mr_results))) { mr_results$align<-""}
  if(!("Locus" %in% colnames(d))) { d$Locus<-""}


  d <- plyr::mutate(d)
  if (nrow(d) < 1 | sum(d$mr_keep) == 0) {
    return(blank_plot("Insufficient number of SNPs"))
  }
  d <- subset(d, mr_keep)
  index <- d$beta.exposure < 0
  d$beta.exposure[index] <- d$beta.exposure[index] *
    -1
  d$beta.outcome[index] <- d$beta.outcome[index] *
    -1
  mrres <- mr_results
  mrres$a <- 0
  if ("MR Egger" %in% mrres$method) {
    d_split <- split(d, d$align)
    k_intercept <- lapply(d_split, function(data){
      temp <- TwoSampleMR::mr_egger_regression(data$beta.exposure,
                                               data$beta.outcome, data$se.exposure, data$se.outcome,
                                               default_parameters())
      return(data.table(align = unique(data$align), intercept = temp$b_i))}) %>%
      rbindlist(., fill = TRUE)
    mrres <- merge(mrres, k_intercept, by = "align")
    mrres[method == "MR Egger",a:=intercept];mrres[,intercept := NULL]
  }
  if ("MR Egger (bootstrap)" %in% mrres$method) {
    temp <- TwoSampleMR::mr_egger_regression_bootstrap(d$beta.exposure,
                                                       d$beta.outcome, d$se.exposure, d$se.outcome,
                                                       default_parameters())
    mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
  }


  ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure,
                                         y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome -
                                                                                                    se.outcome, ymax = beta.outcome + se.outcome),
                                                                                     colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure -
                                                                                                                                                          se.exposure, xmax = beta.exposure + se.exposure),
                                                                                                                                           colour = "grey", height = 0) + ggplot2::geom_point() +
    ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a,
                                                    slope = b, colour = method), show.legend = TRUE) +
    # ggplot2::scale_colour_manual(values = c("#a6cee3",
    #                                         "#1f78b4", "#b2df8a", "#6a3d9a","#33a02c", "#fb9a99",
    #                                         "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
    #                 "#b15928")) +
    ggplot2::scale_colour_manual(values = c(rgb(112, 54, 153, maxColorValue = 255),
                                            rgb(66, 94, 191, maxColorValue = 255),
                                            rgb(84,201,237, maxColorValue = 255),
                                            rgb(59,166,18, maxColorValue = 255),
                                            rgb(255,110,26, maxColorValue = 255),
                                            rgb(149,199,71, maxColorValue = 255),
                                            rgb(161,15,125, maxColorValue = 255),
                                            rgb(249, 106, 27, maxColorValue = 255),
                                            rgb(214,15,102, maxColorValue = 255),
                                            rgb(8, 161, 217, maxColorValue = 255),
                                            rgb(255,186,51, maxColorValue = 255),
                                            rgb(54, 150, 214, maxColorValue = 255))) +
    ggplot2::labs(colour = "MR Test",
                  x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on",
                                                                       d$outcome[1])) + ggplot2::theme(legend.position = "top",
                                                                                                       legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
    expand_limits(x = 0, y = 0) +
    ggrepel::geom_text_repel(label=d$Locus, check_overlap = T) +
    theme(
      legend.background = element_rect(fill = "white", size = 4, colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.2),
      # panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.key=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      # panel.border = element_rect(colour = "grey70", fill=NA, size=1),
      axis.line.x.bottom = element_line(color = "black", size = 0.2),
      axis.line.y.left   = element_line(color = "black", size = 0.2),
      axis.line.y.right  = element_blank(),
      axis.text.y.right  = element_blank(),
      panel.border       = element_blank(),
      legend.position=legend.position) +
    facet_grid(eval(parse(text = equation_facet_grid))) +
    expand_limits(x = 0, y = 0)
}
