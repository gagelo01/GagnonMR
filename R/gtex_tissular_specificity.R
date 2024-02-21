#' obtain the tpm for all Gex tissues
#'
#' @param gene_list a vector of hgnc genes
#' @param ensembl_id_list a vector of ensembl_id can be given of gene_list.
#' @param gencode the file used by GTEX
#' @param tissue_gtex_wd the path to all tissue expression from GTEX
#' @param gtex_tpm_file_name the name of te file with all TPM
#' @param individus_EUR the path to the file with the id of European participants
#'
#' @return
#' @export

get_tpm_for_genes_on_all_tissues <- function(gene_list,
                                             ensembl_id_list = NULL,
                                             gencode = data.table::fread("/home/couchr02/workspace/GTEx_v8/gencode.v26.GRCh38.genes.txt", stringsAsFactors = F),
                                             tissue_gtex_wd =  "/home/couchr02/workspace/GTEx_v8/list/",
                                             gtex_tpm_file_name = "/home/couchr02/Mendel_Commun/GTEx_v8_TPM/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
                                             individus_EUR = data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/EUR_GTEx8.txt", stringsAsFactors = F, data.table = F, header = F)
){

  if(!("ensembl_id"%in%colnames(gencode))){ gencode[,c("ensembl_id", "suffix") :=tstrsplit(gene_id, split= ".", fixed = TRUE)]}
  if(is.null(ensembl_id_list)){
    gene_loop = gencode[gene_name %in% gene_list[2:length(gene_list)], gene_id] %>% unique
    if(gene_list[1]!="Name"){gene_list <- c("Name",gene_list)}
  }  else{
    gene_list <- ensembl_id_list
    if(gene_list[1]!="Name"){gene_list <- c("Name",gene_list)}
    gene_loop = gencode[ensembl_id %in% gene_list[2:length(gene_list)], gene_id] %>% unique
  }
  file_temp <- paste0(tempfile(), ".txt")

  data.table::fwrite(list(gene_list), file_temp, quote = F, col.names = F, row.names = F, sep = "\t")

  tpm_count = data.table::fread(cmd = paste0("zgrep -w -f ", file_temp, " ", gtex_tpm_file_name),
                                data.table = F,
                                stringsAsFactors = F)
  # tpm_count = tpm_count[which(tpm_count$Description %in% gene_list[-1]),]   # gene_list[-1] pour enlever "Name" de la liste

  for (gene in gene_loop) {
    gene_name = unique(gencode$gene_name[which(gencode$gene_id == gene)])
    message("\n#### Gene in progress: ", gene_name," (", which(gene_loop == gene), "/", length(gene_loop), ") ####")

    if (length(which(tpm_count$Description == gene_name)) > 0){
      tpm_count_gene = subset(tpm_count, Description == gene_name)
    } else {
      message("   Missing gene")
      next()
    }

    gene_id = tpm_count_gene$Name

    if (gene == gene_loop[1]){
      data_tau = data.frame(Ensembl.Gene.ID = gene, stringsAsFactors = F)
    } else {
      data_tau = rbind(data_tau, NA)
      data_tau$Ensembl.Gene.ID[which(is.na(data_tau$Ensembl.Gene.ID) == T)] = gene
    }



    tissues_loop = list.files(tissue_gtex_wd, pattern="_liste.txt")
    # tissues_loop = tissues_loop[-which(tissues_loop %in% c("Vagina_liste.txt", "Uterus_liste.txt", "Testis_liste.txt", "Prostate_liste.txt", "Ovary_liste.txt", "Breast_Mammary_Tissue_liste.txt", "Fallopian_Tube_liste.txt", "Cervix_Ectocervix_liste.txt", "Cervix_Endocervix_liste.txt", "Bladder_liste.txt", "Kidney_Cortex_liste.txt", "Kidney_Medulla_liste.txt"))]

    #tissue_list = tissues_loop[1]
    for (tissue_list in tissues_loop) {
      tissue = gsub(tissue_list, pattern = "_liste.txt", replacement = "")

      # message("     ## Tissue in progress: ", tissue," (", which(tissues_loop == tissue_list), "/", length(tissues_loop), ") ####")

      list_tissue = data.table::fread(paste0(tissue_gtex_wd, tissue_list), data.table = F, stringsAsFactors = F, nThread = 6, header = F)

      list_tissue$ID = sapply(list_tissue$V1, FUN = function(x) {paste(unlist(strsplit(x, "-"))[1], unlist(strsplit(x, "-"))[2], sep = "-")})

      list_tissue_EUR = subset(list_tissue, ID %in% individus_EUR$V1)

      tpm_count_tissue_EUR = tpm_count_gene[,c("Name", "Description", colnames(tpm_count_gene)[which(colnames(tpm_count_gene) %in% list_tissue_EUR$V1)])]

      tpm_count_tissue_EUR = t(tpm_count_tissue_EUR)
      tpm_count_tissue_EUR = tpm_count_tissue_EUR[grep("GTEX-", rownames(tpm_count_tissue_EUR)),, drop = F]

      if (gene == gene_loop[1]) {
        data_tau = cbind(data_tau, data.frame(tissue = mean(as.numeric(tpm_count_tissue_EUR[,1])), stringsAsFactors = F))
        colnames(data_tau)[which(colnames(data_tau) == "tissue")] = paste0("Averaged_TPM.", tissue)
      } else {
        data_tau[data_tau$Ensembl.Gene.ID == gene_id, paste0("Averaged_TPM.", tissue)] = mean(as.numeric(tpm_count_tissue_EUR[,1]))
      }

    } # end tissues loop
  }# end genes loop


  fTau_TPM <- function(average_TPM)
  {
    if (all(!is.na(average_TPM)))
    {
      if (min(average_TPM, na.rm = TRUE) >= 0)
      {
        if (max(average_TPM) != 0)
        {
          x =  1-(average_TPM/max(average_TPM))
          res <- sum(x, na.rm = TRUE)
          res <- res/(length(x) - 1)
        } else {
          res <- 0
        }
      } else {
        res <- NA
        print("Expression values have to be positive!")
      }
    } else {
      res <- NA
      print("No data for this gene available.")
    }
    return(res)
  }

  tpm = 1
  x = data_tau[,c(-1)]
  x[x < tpm] = 1
  x = log2(x)

  data_tau$fTau = apply(x, 1, fTau_TPM)
  data_tau = merge(data_tau, tpm_count[,c("Description", "Name")], by.x = "Ensembl.Gene.ID", by.y = "Name")
  setDT(data_tau)
  dt_united <- data.table::melt(data_tau, id.vars = c("Ensembl.Gene.ID", "fTau", "Description"))
  dt_united[, variable := gsub("Averaged_TPM.", "", variable, fixed = TRUE)]
  data.table::setnames(dt_united, c("variable", "value", "Description"), c("tissue", "TPM", "hgnc"))
  dt_united[, tissue := gsub(tissue, pattern = "_", replacement = " ")]
  dt_united[, tissue := gsub(tissue, pattern = "([a-z])([A-Z])", replacement = "\\1 \\2", perl = T)]

  return(dt_united)
}


#' Plot a publication ready heatmap of TPM
#'
#' @param dt_united the object from get_tpm_for_genes_on_all_tissues
#'
#' @return a publication ready heatmap of TPM accross tissues
#' @export

heatmap_tissular_specificity <- function(dt_united ) {
  # sorted_gene = unique(dplyr::arrange(data_united, chromosome_name, band) %>% dplyr::select(Description, chr_band))
  # x = rev(sorted_gene$Description)[1]
  label_names = sapply(dt_united$hgnc, FUN = function(x){bquote(paste(.(x), " (",tau,": ", .(format(unique(dt_united$fTau[which(dt_united$hgnc == x)]),  digits = 2, nsmall = 2)), ")"))})



  g <- ggplot(dt_united, aes(x = tissue, y = hgnc, fill = TPM))  +
    geom_tile() +
    labs(fill = "Averaged TPM\n") +
    viridis::scale_fill_viridis(discrete = FALSE) +
    scale_y_discrete(labels = label_names) +
    # guides(colour = guide_colourbar(title.vjust = 1)) +
    theme(
      # legend.position = "none",
      # panel.border = element_rect(colour = "gray20", fill = "transparent", size = 1),
      # panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor.y = element_blank(),
      # panel.grid.minor.x = element_blank(),
      panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
      #legend.position = "right",
      legend.position = "top",
      legend.title = element_text(
        color = "gray20",
        size = 12
        # margin = margin(l = 0.2, r = 0.2, t = 0.2, b = 0.8)
      ),
      legend.text = element_text(
        color = "gray20",
        size = 10
        # margin = margin(l = 0.2, r = 0.2)
      ),
      legend.title.align = 0.5,
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      # panel.grid.major.x = element_blank(),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 1, colour = "gray20"),
      axis.ticks = element_line(linewidth = 1, colour = "gray20"),
      axis.text.y = element_text(
        # angle = 60,
        size = 10,
        # vjust = 0.5,
        colour = "gray20"
      ),
      axis.text.x = element_text(
        angle = 60,
        size = 10,
        # vjust = 0.55,
        hjust = 1,
        colour = "gray20"
      ),
      # axis.text = element_text(
      #   # angle = 60,
      #   size = 20,
      #   # vjust = 0.5,
      #   colour = "gray20"
      # ),
      axis.ticks.length = unit(.25, "cm"))

  print(g)
}
