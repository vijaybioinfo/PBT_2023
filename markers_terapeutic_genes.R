
{ cat(redb("### Target genes heatmap Expanded cells - CD4 & CD8 ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  setwd('/home/kmlanderos/tmp_large/pbtumor-all/results/figures/CD4_CD8_subclustering')

  sc_cd3p_cd4 <- readRDS("data/sc_cd3p_cd4_seurat_object.rds")
  sc_cd3p_cd8 <- readRDS("data/sc_cd3p_cd8_seurat_object.rds")

  source("/home/ciro/scripts/clustering/R/utilities.R")
  source("/home/ciro/scripts/dgea/R/group_specific.R")
  source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R")

  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
    "/home/ciro/scripts/figease/source.R", # fig_global_objects
    "ggplot2", "cowplot", "patchwork", "Seurat", "stringr", "tidyverse", "data.table", "pheatmap", "viridis"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())

  result_id = "markers_clusters/"; dir.create("markers_clusters")
  fconfigs = list(
    list(object = "sc_cd3p_cd4", name = "target_genes", hname = "orig.donor",
      genes = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "ICOS", "AMICA1", "TNFRSF9"), #  TIM3 == HAVCR2, 41BB == TNFSF9, OX40 == CD134 == TNFRSF4
      filters = list(c("expDegree", "Expanded"), c("cell_classification", "-TREG"))
    ),
    list(object = "sc_cd3p_cd8", name = "target_genes", hname = "orig.donor",
      genes = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "ICOS", "AMICA1", "TNFRSF9"), #
      filters = list(c("expDegree", "Expanded"))
    )
  )

  # fconfig = fconfigs[[2]]; fconfig = fconfigs[[1]]
  for (fconfig in fconfigs) {
    cat("===================", fconfig$object, fconfig$name, "\n")
    command = paste0("Idents(", fconfig$object, ") <- '", fconfig$hname, "'")
    cat(" -", command, "\n"); eval(parse(text = command))
    scells = filters_subset_df(
      fconfig$filters, eval(parse(text = fconfig$object))@meta.data, v = TRUE)
    path_i <- paste0(result_id, fconfig$object, "_", fconfig$name, "/")
    if(!dir.exists(path_i)) dir.create(path_i)

    object = eval(parse(text = fconfig$object))[, scells]
    expr_info <- object@assays$RNA[fconfig$genes,]

    mdata <- object@meta.data
    donors <- c("BT1", "BT2", "BT3", "BT4", "BT5", "BT7", "BT8", "BT9", "BT10", "BT11", "BT12", "BT13", "BT15", "BT17_brain", "BT18", "BT19", "BT20", "BT21", "BT22", "BT23", "BT24", "BT25", "BT26", "BT27")
    grade <- c("LG","HG","LG","LG","LG","LG","LG","LG","HG","HG","HG","HG","HG","HG","HG","LG","HG","HG","HG","HG","LG","LG","HG","LG")
    if(fconfig$object == "sc_cd3p_cd4") {
      grade <- grade[which(!donors %in% c("BT2", "BT9", "BT10"))]
      donors <- donors[which(!donors %in% c("BT2", "BT9", "BT10"))]
    }
    if(fconfig$object == "sc_cd3p_cd8"){
      grade <- grade[which(!donors %in% c("BT2", "BT8", "BT10"))]
      donors <- donors[which(!donors %in% c("BT2", "BT8", "BT10"))]
    }

    # donor = "BT1"
    rm(result_matrix)
    for (donor in donors){
      cells_per_donor <- rownames(mdata[mdata$orig.donor == donor & !is.na(mdata$orig.donor),])

      if(length(cells_per_donor)>1){
        sgenes_mean <- rowMeans(as.matrix(expr_info[, cells_per_donor]), na.rm = FALSE)
        prct_cell <- apply(as.matrix(expr_info[, cells_per_donor]) > 0, MARGIN= 1, function(x){100*sum(x)/length(x)})
        sgenes_mean_positive <- apply(as.matrix(expr_info[, cells_per_donor]), MARGIN= 1, function(x){pos_cells <- x>0; ifelse(sum(pos_cells)>=1, mean(x[pos_cells]), 0)  })
      } else{
        sgenes_mean <- rep(NA, length(fconfig$genes))
        prct_cell <- rep(NA, length(fconfig$genes))
        sgenes_mean_positive <- rep(NA, length(fconfig$genes))
      }

      if(donor == donors[1]){
      result_matrix <- matrix(c(rbind(sgenes_mean, prct_cell, sgenes_mean_positive)), nrow=1, dimnames = list(c(donor), c(rbind(paste0(names(sgenes_mean),"_meanSeuratNormalized"), paste0(names(prct_cell),"_percentageSeuratNormalized"), paste0(names(sgenes_mean_positive),"_meanOfPositiveSeuratNormalized") ))))
      } else{result_matrix <- rbind(result_matrix, c(rbind(sgenes_mean, prct_cell, sgenes_mean_positive))); pos <- which(donors == donor); rownames(result_matrix) <- donors[1:pos]}
    }
    write.csv(as.data.frame(result_matrix), file=paste0(path_i, "gene_expression.csv"))

    # -- Plot

    pallettes <- c(pal1 = colorRampPalette(c("yellow", "red")), pal2 = colorRampPalette(c("purple", "yellow")), pal3 = colorRampPalette(c("blue", "yellow")), pal4 = colorRampPalette(c("blue", "red")), pal5 = colorRampPalette(c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')),
    pal6 = colorRampPalette(c("blue", "yellow", "red")), pal7 = colorRampPalette(c("blue", "white", "red")), pal8 = colorRampPalette(viridis(3)), pal9 = colorRampPalette(c("#5F9DF7", "white", "red")), pal10 = colorRampPalette(c("#87A2FB", "white", "red")),
    pal11 = colorRampPalette(c("#277BC0", "white", "red")), pal12 = colorRampPalette(c("#1F4690", "white", "red")), pal13 = colorRampPalette(c("#1363DF", "white", "red"))
    )

    for(paletteVersion in names(pallettes)){ # paletteVersion = "pal11"
      # ------------- Heatmaps -------------- #
      dir.create(paste0(path_i, paletteVersion))
      palette <- pallettes[[paletteVersion]]

      # Mean
      rm_mean <- result_matrix[,grepl("_meanSeuratNormalized", colnames(result_matrix))]
      colnames(rm_mean) <- sub("_meanSeuratNormalized", "", colnames(rm_mean))
      donor_order <- names(sort(rm_mean[,"PDCD1"], decreasing = TRUE)) # Arrange based on PDCD1
      grade_df = data.frame("Grade" = grade); rownames(grade_df) = rownames(rm_mean) # name matching
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_scaled_rows.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_scaled_cols.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      # showing numbers
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = round(t(rm_mean[donor_order,]),2), fontsize_number = 2.5)
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_scaled_rows_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)
      pheatmap(t(rm_mean[donor_order,]), fontsize_col = 7, main = "Mean Expression", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_scaled_cols_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)

      # Percentage
      rm_mean_perct <- result_matrix[,grepl("_percentageSeuratNormalized", colnames(result_matrix))]
      colnames(rm_mean_perct) <- sub("_percentageSeuratNormalized", "", colnames(rm_mean_perct))
      donor_order <- names(sort(rm_mean_perct[,"PDCD1"], decreasing = TRUE)) # Arrange based on PDCD1
      grade_df = data.frame("Grade" = grade); rownames(grade_df) = rownames(rm_mean_perct) # name matching

      paletteLength <- 100
      capped_mat <- t(rm_mean_perct[donor_order,]); capped_mat[which(capped_mat>50)] = 50
      myBreaks <- c(seq(min(capped_mat), 20, length.out=ceiling(paletteLength/2) + 1),
              seq(20, max(capped_mat), length.out=floor(paletteLength/2))[-1])
      pheatmap(capped_mat, fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA",
        breaks = myBreaks, legend_breaks = c(0,10,20,30,40,50))
      pheatmap(t(rm_mean_perct[donor_order,]), fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage_scaled_rows.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      pheatmap(t(rm_mean_perct[donor_order,]), fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage_scaled_cols.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      # showing numbers
      pheatmap(t(rm_mean_perct[donor_order,]), fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = round(t(rm_mean_perct[donor_order,]),2), fontsize_number = 2.5)
      pheatmap(t(rm_mean_perct[donor_order,]), fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage_scaled_rows_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)
      pheatmap(t(rm_mean_perct[donor_order,]), fontsize_col = 7, main = "Percentage of Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_cell_percentage_scaled_cols_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)

      # Mean Positive
      rm_mean_pos <- result_matrix[,grepl("_meanOfPositiveSeuratNormalized", colnames(result_matrix))]
      colnames(rm_mean_pos) <- sub("_meanOfPositiveSeuratNormalized", "", colnames(rm_mean_pos))
      donor_order <- names(sort(rm_mean_pos[,"PDCD1"], decreasing = TRUE)) # Arrange based on PDCD1
      grade_df = data.frame("Grade" = grade); rownames(grade_df) = rownames(rm_mean_perct) # name matching
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive_scaled_rows.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive_scaled_cols.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
      # showing numbers
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = round(t(rm_mean_pos[donor_order,]),2), fontsize_number = 2.5)
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "row", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive_scaled_rows_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)
      pheatmap(t(rm_mean_pos[donor_order,]), fontsize_col = 7, main = "Mean Expression of Positive Cells", annotation_col = grade_df, cluster_rows=F, cluster_cols=F, scale = "column", filename = paste0(path_i, paletteVersion, "/heatmap_mean_expression_positive_scaled_cols_wNumbers.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA", display_numbers = TRUE, fontsize_number = 2.5)
    }

  }

}
