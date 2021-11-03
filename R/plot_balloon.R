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
#' @examples
#' dat_plot <- readRDS("/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Primary/dat_balloon_plot")
#' plot_balloon(dat_plot, bonferroni_threshold = 0.05)
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
    ggplot2::scale_size(name = expression(-Log[10](P)), range = c(4,10)) +
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
