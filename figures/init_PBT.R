#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2023-04-20
# ---

### ==== CD4/CD8 clustering ==== ####
# CD4 - pct25_pc15_res0.6
# CD8 - pct20_pc15_res0.6

### ================== Figures ================== ###
setwd('/home/kmlanderos/tmp_large/pbtumor-all-Batch2'); here::here()

# Read CD4 and CD8 Seurat Objects
sc_cd3p_cd8 <- readRDS(file = "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_cd3p_cd8_Batch1_Batch2_seurat_object_subset.rds")
sc_cd3p_cd4 <- readRDS(file = "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_cd3p_cd4_Batch1_Batch2_seurat_object_subset.rds")
gzmk_trm_cd8_20pct_30pc <- readRDS(file = "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_cd3p_cd8_gzmk_trm_Batch1_Batch2_seurat_object_subset.rds")

figures_path <- here::here("results", "figures")
setwd(figures_path);

{ cat(redb("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # As you progress they appear in many tasks and so become more or less global

  redu = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2"))
  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/clustering/R/utilities.R",
    "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
    "/home/ciro/scripts/figease/source.R", # fig_global_objects
    "ggplot2", "cowplot", "patchwork", "Seurat", "stringr", "tidyverse", "data.table", "pheatmap", "fmsb", "grid",
    "ggpubr", "ComplexHeatmap", "RColorBrewer", "ggsci", "scales", "UpSetR", "igraph"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())

  # CD8 cluster colours
  cd8_cluster_color <- c( "#a3e3d4", "#72e5d2", "#3bb8ae", "#007e74",   "#ffa7d7",   "#fcee62",   "#c5b8ff", "#8c50cf",  "#ff971f",   "#4da7fd",   "#9c0006")
  names(cd8_cluster_color) <- c("0", "1", "2", "8",   "7",   "6",   "5", "9",   "4",   "3",   "10")
  # CD4 cluster colours
  cd4_cluster_color <- c("#c5b8ff", "#d26dec", "#a084ff", "#8c50cf", "#5d1d88",   "#52edf4", "#97ccf5", "#379ee4",   "#ffa7d7", "#ff717c", "#9c0006")
  names(cd4_cluster_color) <- c("1", "2", "4", "5", "8",   "0", "3", "7",   "6", "10", "11")

  # CD8 paper cluster colours
  cd8_paper_cluster_color <- c("#72e5d2", "#ffa7d7", "#4da7fd", "#ff971f", "#c5b8ff", "#fcee62", "#8c50cf", "#9c0006")
  names(cd8_paper_cluster_color) <- c("0", "1", "2", "3", "4", "5", "6", "7")

  # CD8 cell_classification colours
  cd8_cell_classification_color <- c("Cell_Cycle" = "#9c0006", "GZMK_HI" = "#72e5d2", "TRM" = "#ffa7d7", "CD16p_Effector" = "#fcee62", "Effector" = "#ff971f", "MAIT" = "#9e89e8", "TCF7_HI" = "#4da7fd")
  # CD4 cell_classification colours
  # cd4_cell_classification_color <- c("Cell_Cycle" = "#9c0006", "CTL" = "#a084ff", "TFH" = "#ff717c", "TREG" = "#ffa7d7", "TCM_TN" = "#97ccf5")
  cd4_cell_classification_color <- c("Cell_Cycle" = "#9c0006", "CTL" = "#a084ff", "TFH" = "#ff717c", "TREG" = "#ffa7d7", "TN" = "#97ccf5", "TCM" = "#379ee4")


}

{ cat(redb("### Signature Analyses ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # cp /Users/kmlanderos/Documents/PBT/Data/signatures_cd3p*csv /Volumes/kmlanderos/pbtumor-all-Batch2/info/
  # Rscript /home/ciro/scripts/functions/csvCorrect.R /home/kmlanderos/pbtumor-all-Batch2/info/signatures_cd3p.csv

  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R") # stats_summary_table
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
  source("/home/ciro/scripts/handy_functions/devel/utilities.R")  # vlist2df
  source("/home/ciro/scripts/figease/figease.R") # fig_module_score
  source("/home/ciro/scripts/handy_functions/devel/filters.R")

  dir.create("gsea")

  ## ---- GSEA ---- ###

  fname <- "/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/signatures_cd3p_human_v2.csv"
  slists <- readfile(fname, stringsAsFactors = FALSE)

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  gzmk_trm_cd8_20pct_30pc@meta.data <- gzmk_trm_cd8_20pct_30pc@meta.data %>% mutate(TRM_GZMK_subcluster_0.6 = case_when(
    RNA_snn_res.0.6 %in% c("3", "4", "5") ~ "TRM",
    RNA_snn_res.0.6 %in% c("0", "1", "2", "6", "7") ~ "GZMK_HI",
    TRUE ~ "Other"
  ), TRM_GZMK_subcluster_0.8 = case_when(
    RNA_snn_res.0.8 %in% c("2", "5", "6", "9") ~ "TRM",
    RNA_snn_res.0.8 %in% c("0", "1", "3", "4", "7", "8", "10") ~ "GZMK_HI",
    TRUE ~ "Other"
  ))

  fconfigs = list(
    ## ------ CD8
    ## CD8 subsets
    list(result_id = "gsea/CD8_cell_classification/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cell_classification")
    ),
    ## CD8 clusters
    list(result_id = "gsea/CD8_cell_cluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cluster")
    ),
    ## CD8 PDCD1p vs PDCD1n
    list(result_id = "gsea/CD8_PDCD1_tag/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("PDCD1_tag", "PDCD1n", "PDCD1p"),
      comparisons = c("PDCD1_tag")
    ),
    ## CD8 Expansion comparison
    list(result_id = "gsea/CD8_expansion/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("NeoTCR8..ref.Lowery.et..Science.2022.", "CD8_GZMK_Guo2018")], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
    ),
    ## -------- CD4
    ## CD4 subsets
    list(result_id = "gsea/CD4_cell_classification/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Treg_Schmiedel2018", "TFH_Locci2013", "Hu.CD4.Naive.CM", "CellCycle_Best2013", "CD4_NAIVE_outofTcells_R24", "CD4_NAIVE_outofotherImmunecells_R24", "CD4_TCM_UP202_Broad_DS", "Oja_Braga_TCM", "GSE11057_TN_Down", "GSE11057_TN_Up", "GSE22886_TN_Down")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cell_classification_deconv")
    ),
    ## CD4 clusters
    list(result_id = "gsea/CD4_cluster/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Treg_Schmiedel2018", "TFH_Locci2013", "Hu.CD4.Naive.CM", "CellCycle_Best2013", "CD4_NAIVE_outofTcells_R24", "CD4_NAIVE_outofotherImmunecells_R24", "CD4_TCM_UP202_Broad_DS", "Oja_Braga_TCM", "GSE11057_TN_Down", "GSE11057_TN_Up", "GSE22886_TN_Down")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cluster")
    ),
    ## nonTREG CD4 Expansion comparison
    list(result_id = "gsea/nonTREG_CD4_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("NeoTCR4..ref.Lowery.et..Science.2022.", "CD4CTLvsTCM_Patil2018", "Treg_Schmiedel2018", "TFH_Locci2013", "Hu.CD4.Naive.CM", "CellCycle_Best2013")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("expDegree")
    ),
    ## CD4 Expansion comparison
    list(result_id = "gsea/CD4_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("NeoTCR4..ref.Lowery.et..Science.2022.", "CD4CTLvsTCM_Patil2018", "Treg_Schmiedel2018", "TFH_Locci2013", "Hu.CD4.Naive.CM", "CellCycle_Best2013")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("expDegree", "Expanded", "Non_expanded"),
      comparisons = c("expDegree")
    ),
    ## Identify TCM and TN (clusters)
    list(result_id = "gsea/CD4_TCM_TN/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4_NAIVE_outofTcells_R24", "CD4_NAIVE_outofotherImmunecells_R24", "CD4_TCM_UP202_Broad_DS", "Oja_Braga_TCM", "GSE11057_TN_Down", "GSE11057_TN_Up", "GSE22886_TN_Down")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification", "TCM_TN"),
      comparisons = c("cluster")
    ),
    ## Identify TCM and TN (grouped clusters)
    list(result_id = "gsea/CD4_TCM_TN_grouped_clusters/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4_NAIVE_outofTcells_R24", "CD4_NAIVE_outofotherImmunecells_R24", "CD4_TCM_UP202_Broad_DS", "Oja_Braga_TCM", "GSE11057_TN_Down", "GSE11057_TN_Up", "GSE22886_TN_Down")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification_deconv", "TCM", "TN"),
      comparisons = c("cell_classification_deconv")
    ),
    ## -------- CD8 subclustering
    ## TRM/GZMK signatures cluster
    list(result_id = "gsea/CD8_GZMK_TRM_0.6_clusters/",
      edata = "expm1(gzmk_trm_cd8_20pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "TRM_CLARKE_small", "DICE_C1.Ag.spec.TRM_logFC0.3", "DICE_C1.Ag.spec.TRM_logFC0.4", "DICE_C2.TRM_logFC0.3", "DICE_C2.TRM_logFC0.4")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("RNA_snn_res.0.6")
    ),
    list(result_id = "gsea/CD8_GZMK_TRM_0.6_grouped_clusters/",
      edata = "expm1(gzmk_trm_cd8_20pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "TRM_CLARKE_small", "DICE_C1.Ag.spec.TRM_logFC0.3", "DICE_C1.Ag.spec.TRM_logFC0.4", "DICE_C2.TRM_logFC0.3", "DICE_C2.TRM_logFC0.4")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("TRM_GZMK_subcluster_0.6")
    ),
    ## TRM/GZMK signatures grouped clusters
    list(result_id = "gsea/CD8_GZMK_TRM_0.8_clusters/",
      edata = "expm1(gzmk_trm_cd8_20pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "TRM_CLARKE_small", "DICE_C1.Ag.spec.TRM_logFC0.3", "DICE_C1.Ag.spec.TRM_logFC0.4", "DICE_C2.TRM_logFC0.3", "DICE_C2.TRM_logFC0.4")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("RNA_snn_res.0.8")
    ),
    list(result_id = "gsea/CD8_GZMK_TRM_0.8_grouped_clusters/",
      edata = "expm1(gzmk_trm_cd8_20pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "TRM_CLARKE_small", "DICE_C1.Ag.spec.TRM_logFC0.3", "DICE_C1.Ag.spec.TRM_logFC0.4", "DICE_C2.TRM_logFC0.3", "DICE_C2.TRM_logFC0.4")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("TRM_GZMK_subcluster_0.8")
    ),
    # Main UMAP clusters
    list(result_id = "gsea/CD8_GZMK_TRM_MainUMAP_cell_classification/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cell_classification_deconv")
    ),
    # Main UMAP cell_classification
    list(result_id = "gsea/CD8_GZMK_TRM_MainUMAP_cell_cluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cluster_deconv")
    )
  )

  # CD8
  pp_gsea = fig_gsea(fconfigs[c(1:4)], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # CD4
  pp_gsea = fig_gsea(fconfigs[c(5:10)], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  # CD8 TRM GZMK subclustering
  pp_gsea = fig_gsea(fconfigs[c(11:14)], features = rownames(gzmk_trm_cd8_20pct_30pc), return_plot = TRUE, verbose = T)
  # CD8 TRM GZMK subclustering
  pp_gsea = fig_gsea(fconfigs[c(15:16)], features = rownames(mdata_sc_cd3p_cd8), return_plot = TRUE, verbose = T)

  ## ---- Module scoring analysis ---- ###

  source("/home/ciro/scripts/handy_functions/R/gsea_signature.R") # clean_feature_list, signature_scoring
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  source("/home/kmlanderos/scripts/handy_functions/R/moduleScore_functions.R") # signature_scoring (updated)
  dir.create("modulescores");

  # sc_cd3p_lists <- readfile("/home/ciro/pbtumor/large/results/figures/gsea/sc_cd3p_mod_signatures.csv", stringsAsFactors = FALSE)
  sc_cd3p_lists <- readfile("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/signatures_cd3p_human_v2.csv", stringsAsFactors = FALSE)

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

  fconfigs = list(
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4", object = "sc_cd3p_cd4",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2),c("CD4CTLvsTCM_Patil2018", "TFH_Locci2013", "Treg_Schmiedel2018", "TH17_Seumois2020", "tcell_activation_GO.0042110.1", "tcr_signaling_RHSA202403.1", "CellCycle_Best2013.1", "TypeIandII_IFNsignaling_Seumois2020.1", "th1_signature1_arlehamn", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "CD4_NAIVE_outofTcells_R24", "CD4_NAIVE_outofotherImmunecells_R24", "CD4_TCM_UP202_Broad_DS", "Oja_Braga_TCM", "GSE11057_TN_Down", "GSE11057_TN_Up", "GSE22886_TN_Down", "GSE22886_TN_Up")]),
      axis_x = list(col = "cell_classification"),# order=c("TCM", "CTL", "1", "10", "11", "5")),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    ),
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "sc_cd3p_cd8@meta.data", object = "sc_cd3p_cd8",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2), c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "CYTOTOXICITY.SIGNATURE_CD161.paper_matthewson..Cell", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like", "Corgnac2020_CD103.TIL_vs_KLRG1.TIL", "Kumar2017_Lung.CD69..vs.CD69.", "Milner2017_Core.ms.TRM")]),
      axis_x = list(col = "cell_classification"),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    ),
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(gzmk_trm_cd8_20pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data", object = "gzmk_trm_cd8_20pct_30pc",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2), c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "CYTOTOXICITY.SIGNATURE_CD161.paper_matthewson..Cell", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like", "Corgnac2020_CD103.TIL_vs_KLRG1.TIL", "Kumar2017_Lung.CD69..vs.CD69.", "Milner2017_Core.ms.TRM")]),
      axis_x = list(col = "cluster"),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    ),
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(gzmk_trm_cd8_15pct_30pc@assays$RNA@data)",
      metadata = "gzmk_trm_cd8_15pct_30pc@meta.data", object = "gzmk_trm_cd8_15pct_30pc",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2), c("CD8_GZMK_Guo2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "MAIT_Yao2020", "CD8.CX3CR1_Guo2018", "Sade.feldman_Hu.CD8.TSCM..repsonders", "CellCycle_Best2013", "CYTOTOXICITY.SIGNATURE_CD161.paper_matthewson..Cell", "Krishna_Hu.CD8.TSCM..responders", "HU_CD8.TSCM.vs.TN", "HU_CD8.TSCM.vs.TCM", "Yao2019_Progenitor.like", "Corgnac2020_CD103.TIL_vs_KLRG1.TIL", "Kumar2017_Lung.CD69..vs.CD69.", "Milner2017_Core.ms.TRM")]),
      axis_x = list(col = "cluster"),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    )
  )
  pp_ms = fig_module_score(fconfigs[1], reductions = redu[1], violins_color = 'mean') # calls signature_scoring
  pp_ms = fig_module_score(fconfigs[2], reductions = redu[1], violins_color = 'mean') # calls signature_scoring
  pp_ms = fig_module_score(fconfigs[3:4], reductions = redu[1], violins_color = 'mean') # calls signature_scoring

  # Adjust scale of the CD4-CTL signature (0-0.2)
  signature_scores <- read.csv("modulescores/sc_cd3p_cd4_test/signatures.csv", row.names=1)
  if( all(rownames(mdata_sc_cd3p_cd4) == rownames(signature_scores)) ) {
    mdata_sc_cd3p_cd4$CD4CTL_signature <- signature_scores[,"CD4CTLvsTCM_Patil2018"]
  }

  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(CD4CTL_signature_adjusted = case_when(
      CD4CTL_signature > 0.15 ~ 0.15, # 0.20, Good 0.15
      TRUE ~ CD4CTL_signature)
    )

  # Using Vicente's color scale 0.15
  pdf(paste0("modulescores/sc_cd3p_cd4_test/", "cd4ctlvstcm_patil2018_0.15_umap_v3.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()


}

{ cat(redb("### V and J genes usage ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  # ==================== CD8

  # ---> Beta Chain
  # trb.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.v)) %>% select(cell_classification_deconv, trb.v) %>% group_by(cell_classification_deconv, trb.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_beta_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.j)) %>% select(cell_classification_deconv, trb.j) %>% group_by(cell_classification_deconv, trb.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_beta_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.v)) %>% select(cell_classification_deconv, tra.v) %>% group_by(cell_classification_deconv, tra.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_alpha_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.j)) %>% select(cell_classification_deconv, tra.j) %>% group_by(cell_classification_deconv, tra.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_alpha_usage_perCelltype.tsv", sep = "\t", row.names = F)

  # trb.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.v)) %>% select(cluster_deconv, trb.v) %>% group_by(cluster_deconv, trb.v) %>% summarize(n_cells = n()) %>% arrange(cluster_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.j)) %>% select(cluster_deconv, trb.j) %>% group_by(cluster_deconv, trb.j) %>% summarize(n_cells = n()) %>% arrange(cluster_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.v)) %>% select(cluster_deconv, tra.v) %>% group_by(cluster_deconv, tra.v) %>% summarize(n_cells = n()) %>% arrange(cluster_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.j)) %>% select(cluster_deconv, tra.j) %>% group_by(cluster_deconv, tra.j) %>% summarize(n_cells = n()) %>% arrange(cluster_deconv, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)

  # ==================== CD4

  # ---> Beta Chain
  # trb.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.v)) %>% select(cell_classification, trb.v) %>% group_by(cell_classification, trb.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_beta_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.j)) %>% select(cell_classification, trb.j) %>% group_by(cell_classification, trb.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRJ_beta_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.v)) %>% select(cell_classification, tra.v) %>% group_by(cell_classification, tra.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_alpha_usage_perCelltype.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.j)) %>% select(cell_classification, tra.j) %>% group_by(cell_classification, tra.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRJ_alpha_usage_perCelltype.tsv", sep = "\t", row.names = F)

  # trb.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.v)) %>% select(cluster, trb.v) %>% group_by(cluster, trb.v) %>% summarize(n_cells = n()) %>% arrange(cluster, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.j)) %>% select(cluster, trb.j) %>% group_by(cluster, trb.j) %>% summarize(n_cells = n()) %>% arrange(cluster, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRJ_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.v)) %>% select(cluster, tra.v) %>% group_by(cluster, tra.v) %>% summarize(n_cells = n()) %>% arrange(cluster, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.j)) %>% select(cluster, tra.j) %>% group_by(cluster, tra.j) %>% summarize(n_cells = n()) %>% arrange(cluster, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRJ_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)

}

{ cat(redb("### Markers: Dot/curtain ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dotplots");

  fconfigs = list(
      list(result_id = "dotplots/CD8/supplementary_2b",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "3", "2", "5", "4", "9")),
       features = c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "CCL4", "CCL5", "EOMES", "SH2D1A", "ITM2C", "CD74", "CXCR4", "ZNF683", "ITGAE", "ITGA1",
       "FCGR3A", "FGFBP2", "GZMB", "GZMA", "CX3CR1", "KLRG1", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "TCF7", "SELL", "CCR7", "STMN1", "TOP2A", "TUBA1B")
     ),
    list(result_id = "dotplots/CD4/supplementary_",
       edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd4", axis_x = list(col = "RNA_snn_res.0.4",
         order = c("0", "2", "9", "4", "1", "10", "5", "3", "6", "7", "8", "11")),
       features = c("CD3D", "CD4", "CD8A", "CD8B", "GZMA", "GZMK", "PRF1", "TNF", "NFKB1", "CXCL13", "BTLA", "SH2D1A", "PDCD1", "FOXP3", "CTLA4", "IL2RA", "ENTPD1", "S1PR1", "SELL", "CCR7", "TCF7", "STMN1", "TOP2A", "TUBA1B")
     )
   )
  pp_curtains = fig_plot_curtain(fconfigs, dot.scale = 12, col.min = -1.5, col.max = 1.5)
}

{ cat(redb("### Markers: violins ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  dir.create("violins")

  fconfigs = list(
      list(result_id = "violins/sc_cd8/cluster_deconv/violin_",
        edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
        axis_x = list(col = "cluster_deconv"),
        # sample_filter = c("celltype", "CD8"),
        features = c("KLRB1", "TCF7", "GZMK", "IFNG", "TNF", "PDCD1", "PRF1", "ZNF683", "S1PR1", "S1PR5", "ITGAE", "ITGA1", "CCL4", "CCL5", "EOMES", "CXCR6", "LAG3")
      ),
      list(result_id = "violins/sc_cd8/cell_classification_deconv/violin_",
        edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
        axis_x = list(col = "cell_classification_deconv", order = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI")),
        # sample_filter = c("celltype", "CD8"),
        features = c("KLRB1", "TCF7", "GZMK", "IFNG", "TNF", "PDCD1", "PRF1", "ZNF683", "S1PR1", "S1PR5", "ITGAE", "ITGA1", "CCL4", "CCL5", "EOMES", "CXCR6", "LAG3")
      )
  )

  pp_violins = fig_plot_violins(
    fconfigs,
    theme_extra = theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10)),
    colour_by = "pct", couls = couls_opt$red_gradient$white
  )

{ cat(redb("### Markers: dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dim_reduction_features", showWarnings = FALSE)

  fconfigs = list(
    list(result_id = "dim_reduction_features/sc_cd3p_cd4_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "PRF1", "GZMA", "IFNG", "TBX21", "EOMES", "GZMB", "SLAMF7", "CRTAM", "TNF", "GZMH", "GNLY", "CCL4", "CCL5",
      "CTSW", "ZNF683", "FCRL6", "SPON2", "CX3CR1", "S1PR5", "NKG7", "CD244", "CXCR3", "S1PR1", "SELL", "TCF7", "FOXP3", "RORA", "RORC",
      "CCR6", "KLRB1", "IRF7", "CXCL13", "BTLA", "CD40LG", "ITGAE", "CD69", "S1PR5", "S1PR1", "KLRB1"),
      col = c("lightgrey", "blue")
    ),
    list(result_id = "dim_reduction_features/sc_cd3p_cd8_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd8@assays$RNA@data",
      metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "GZMK", "GZMA", "PRF1", "LAG3", "IFNG", "TNF", "CCL4", "CCL3", "FASLG", "HLA-DRB1", "TOX", "S1PR5", "S1PR1", "KLRB1",
      "ITGAE", "ZNF683", "CXCR6", "ITGA6", "ITGA1", "AMICA1", "GPR25", "RBPJ", "RGS1", "KLF2", "S1PR5", "S1PR1", "TCF7"),
      col = c("lightgrey", "blue")
    ),
    list(result_id = "dim_reduction_features/gzmk_trm_cd8_20pct_30pc_", sufix = "exm_markers/",
      edata = "gzmk_trm_cd8_20pct_30pc@assays$RNA@data",
      metadata = "gzmk_trm_cd8_20pct_30pc@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "GZMK", "GZMA", "PRF1", "LAG3", "IFNG", "TNF", "CCL4", "CCL3", "FASLG", "HLA-DRB1", "TOX", "S1PR5", "S1PR1", "KLRB1",
      "ITGAE", "ZNF683", "CXCR6", "ITGA6", "ITGA1", "AMICA1", "GPR25", "RBPJ", "RGS1", "KLF2", "S1PR5", "S1PR1"),
      col = c("lightgrey", "blue")
    ),
    list(result_id = "dim_reduction_features/gzmk_trm_cd8_15pct_30pc_", sufix = "exm_markers/",
      edata = "gzmk_trm_cd8_15pct_30pc@assays$RNA@data",
      metadata = "gzmk_trm_cd8_15pct_30pc@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "GZMK", "GZMA", "PRF1", "LAG3", "IFNG", "TNF", "CCL4", "CCL3", "FASLG", "HLA-DRB1", "TOX", "S1PR5", "S1PR1", "KLRB1",
      "ITGAE", "ZNF683", "CXCR6", "ITGA6", "ITGA1", "AMICA1", "GPR25", "RBPJ", "RGS1", "KLF2", "S1PR5", "S1PR1"),
      col = c("lightgrey", "blue")
    ))

  pp_markers = fig_plot_scatters(fconfigs[2])

  # Change scale to plots
  plot_blank_fun = function(x){
      plot_rm_layer(x) + theme(legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank())
  }
  # ---> CD4
  # PDCD1
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/PDCD1_new_colorScale.pdf")
  FeaturePlot(sc_cd3p_cd4, features = c("PDCD1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/PDCD1_new_colorScale_blank.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd4, features = c("PDCD1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) ))
  dev.off()
  # LAG3
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/LAG3_new_colorScale.pdf")
  FeaturePlot(sc_cd3p_cd4, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/LAG3_new_colorScale_blank.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd4, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) ))
  dev.off()

  # ---> CD8
  # PDCD1
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/PDCD1_new_colorScale.pdf")
  FeaturePlot(sc_cd3p_cd8, features = c("PDCD1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/PDCD1_new_colorScale_blank.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd8, features = c("PDCD1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) ))
  dev.off()
  # LAG3
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/LAG3_new_colorScale_v2.pdf")
  FeaturePlot(sc_cd3p_cd8, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "gray90", "#7933ff", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/LAG3_new_colorScale_blank_v2.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd8, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "gray90", "#7933ff", "#5800FF"))(100) ))
  dev.off()
  # TCF7
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/TCF7_new_colorScale.pdf")
  FeaturePlot(sc_cd3p_cd8, features = c("TCF7"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/TCF7_new_colorScale_blank.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd8, features = c("TCF7"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) ))
  dev.off()
  # KLRB1
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/KLRB1_new_colorScale_v2.pdf")
  FeaturePlot(sc_cd3p_cd8, features = c("KLRB1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "gray90", "#7933ff", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/KLRB1_new_colorScale_blank_v2.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd8, features = c("KLRB1"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "gray90", "#7933ff", "#5800FF"))(100) ))
  dev.off()

  # ------------------- Inset plot
  # ---> CD4
  # PDCD1
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/PDCD1_inset_plot.pdf")
  mdata <- sc_cd3p_cd4@meta.data; edata <- sc_cd3p_cd4@assays$RNA
  tmp <- as.vector(as.matrix(edata["PDCD1",rownames(mdata)]))
  mdata <- cbind(mdata, PDCD1_expr = tmp)
  mdata$tag_PDCD1 <- add_gene_tag(c("PDCD1"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(PDCD1_expr),
    pct.expr = 100*sum(tag_PDCD1 == "PDCD1p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM", "TN"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30, 50)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()
  # PDCD1 pct_expr capped
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/PDCD1_inset_plot_capped_pct_expr.pdf")
  mdata <- sc_cd3p_cd4@meta.data; edata <- sc_cd3p_cd4@assays$RNA
  tmp <- as.vector(as.matrix(edata["PDCD1",rownames(mdata)]))
  mdata <- cbind(mdata, PDCD1_expr = tmp)
  mdata$tag_PDCD1 <- add_gene_tag(c("PDCD1"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(PDCD1_expr),
    pct.expr = 100*sum(tag_PDCD1 == "PDCD1p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM", "TN"))
  a <- a %>% mutate(pct.expr = case_when(
    pct.expr > 30 ~ 30,
    TRUE ~ pct.expr
  ))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30, 50)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()
  # LAG3
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/LAG3_inset_plot.pdf")
  mdata <- sc_cd3p_cd4@meta.data; edata <- sc_cd3p_cd4@assays$RNA
  tmp <- as.vector(as.matrix(edata["LAG3",rownames(mdata)]))
  mdata <- cbind(mdata, LAG3_expr = tmp)
  mdata$tag_LAG3 <- add_gene_tag(c("LAG3"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(LAG3_expr),
    pct.expr = 100*sum(tag_LAG3 == "LAG3p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM", "TN"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(5, 10, 15, 20, 25)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()


  # ---> CD8
  # PDCD1
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/PDCD1_inset_plot.pdf")
  mdata <- sc_cd3p_cd8@meta.data; edata <- sc_cd3p_cd8@assays$RNA
  tmp <- as.vector(as.matrix(edata["PDCD1",rownames(mdata)]))
  mdata <- cbind(mdata, PDCD1_expr = tmp)
  mdata$tag_PDCD1 <- add_gene_tag(c("PDCD1"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(PDCD1_expr),
    pct.expr = 100*sum(tag_PDCD1 == "PDCD1p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()
  # LAG3
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/LAG3_inset_plot.pdf")
  mdata <- sc_cd3p_cd8@meta.data; edata <- sc_cd3p_cd8@assays$RNA
  tmp <- as.vector(as.matrix(edata["LAG3",rownames(mdata)]))
  mdata <- cbind(mdata, LAG3_expr = tmp)
  mdata$tag_LAG3 <- add_gene_tag(c("LAG3"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(LAG3_expr),
    pct.expr = 100*sum(tag_LAG3 == "LAG3p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30, 40, 50)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()
  # TCF7
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/TCF7_inset_plot.pdf")
  mdata <- sc_cd3p_cd8@meta.data; edata <- sc_cd3p_cd8@assays$RNA
  tmp <- as.vector(as.matrix(edata["TCF7",rownames(mdata)]))
  mdata <- cbind(mdata, TCF7_expr = tmp)
  mdata$tag_TCF7 <- add_gene_tag(c("TCF7"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(TCF7_expr),
    pct.expr = 100*sum(tag_TCF7 == "TCF7p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(15, 30, 45, 60, 70)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()
  # KLRB1
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/KLRB1_inset_plot.pdf")
  mdata <- sc_cd3p_cd8@meta.data; edata <- sc_cd3p_cd8@assays$RNA
  tmp <- as.vector(as.matrix(edata["KLRB1",rownames(mdata)]))
  mdata <- cbind(mdata, KLRB1_expr = tmp)
  mdata$tag_KLRB1 <- add_gene_tag(c("KLRB1"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification_deconv) %>% summarize(mean_expr = mean(KLRB1_expr),
    pct.expr = 100*sum(tag_KLRB1 == "KLRB1p")/n())
  a$category <- "Expansion"; a$cell_classification_deconv <- factor(a$cell_classification_deconv, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification_deconv, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(10, 30, 50, 70, 90)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  print(p)
  dev.off()



}

{ cat(redb("### Markers: contour ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("scatter_contour")

  # ---> Define Function
  scatter_contour <- function(
    data, axis_x, axis_y, axis_z = NULL, dp = TRUE) {
    if(is.null(data$Density) && is.null(axis_z)){
      data$Density <- 0
      dps = if(isTRUE(dp)){
        rowSums(data[, c(axis_x, axis_y)] > 0) > 1
      }else{ rownames(data) }
      tmp = try(MASS_kde2d(
        x = data[dps, axis_x],
        y = data[dps, axis_y]
      ), silent = TRUE)
      if(class(tmp) != "try-error"){
        data[dps, ]$Density <- tmp; axis_z = "Density"
      }else{ data$Density = NULL }
    }
    aesy = if(!is.null(axis_z)){
      aes_string(x = axis_x, y = axis_y, color = axis_z)
      aes(x = !!rlang::sym(axis_x), y = !!rlang::sym(axis_y), color = !!rlang::sym(axis_z))
    }else{ aes(x = !!rlang::sym(axis_x), y = !!rlang::sym(axis_y)) }
    p <- ggplot(data = data, mapping = aesy) + geom_point(size = 0.5)
    if(!is.null(axis_z))
      p <- p + geom_density2d(data = data[dps, ], colour="white") # Making the contour line disappear
    return(p)
  }

  # ----> CD8

  # Adjust the metadata
  sc_cd3p_cd8_mdata <- sc_cd3p_cd8@meta.data

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  sc_cd3p_cd8_mdata[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  sc_cd3p_cd8_mdata <- sc_cd3p_cd8_mdata %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))


  # Expanded
  fconfigs = list(
    list(result_id = "scatter_contour/sc_cd8_Expanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("GZMA"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("expDegree", "Expanded"))
    ),
    list(result_id = "scatter_contour/sc_cd8_Expanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("TNF"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("expDegree", "Expanded"))
    )
  )
  pp_contour = fig_plot_contour(fconfigs[1],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4.5)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )
  pp_contour = fig_plot_contour(fconfigs[2],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.5)) + scale_x_continuous(limits = c(0, 6)) + scale_y_continuous(limits = c(0, 6)) + geom_point(size=3), limits = list(0.5, 0.5), type = "percent")
  )

  # Non_expanded
  fconfigs = list(
    list(result_id = "scatter_contour/sc_cd8_Non_expanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("GZMA"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("expDegree", "Non_expanded"))
    ),
    list(result_id = "scatter_contour/sc_cd8_Non_expanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("TNF"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("expDegree", "Non_expanded"))
    )
  )
  pp_contour = fig_plot_contour(fconfigs[1],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4.5)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )
  pp_contour = fig_plot_contour(fconfigs[2],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.5)) + scale_x_continuous(limits = c(0, 6)) + scale_y_continuous(limits = c(0, 6)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # PDCD1p
  fconfigs = list(
    list(result_id = "scatter_contour/sc_cd8_PDCD1p/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("GZMA"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("PDCD1_tag", "PDCD1p"))
    ),
    list(result_id = "scatter_contour/sc_cd8_PDCD1p/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("TNF"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("PDCD1_tag", "PDCD1p"))
    )
  )
  pp_contour = fig_plot_contour(fconfigs[1],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4.5)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )
  pp_contour = fig_plot_contour(fconfigs[2],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.4)) + scale_x_continuous(limits = c(0, 6)) + scale_y_continuous(limits = c(0, 6)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # PDCD1n
  fconfigs = list(
    list(result_id = "scatter_contour/sc_cd8_PDCD1n/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("GZMA"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("PDCD1_tag", "PDCD1n"))
    ),
    list(result_id = "scatter_contour/sc_cd8_PDCD1n/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata",
      features = list(x = c("TNF"), y = c("PRF1", "IFNG")),
      sample_filter = list(c("PDCD1_tag", "PDCD1n"))
    )
  )
  pp_contour = fig_plot_contour(fconfigs[1],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4.5)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )
  pp_contour = fig_plot_contour(fconfigs[2],
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.4)) + scale_x_continuous(limits = c(0, 6)) + scale_y_continuous(limits = c(0, 6)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # ----> CD4

  # Adjust the metadata
  sc_cd3p_cd4_mdata <- sc_cd3p_cd4@meta.data
  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd4@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd4@assays$RNA@data)
  sc_cd3p_cd4_mdata[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  sc_cd3p_cd4_mdata <- sc_cd3p_cd4_mdata %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  fconfigs = list(
      list(result_id = "scatter_contour/sc_cd4_Expanded/",
        edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster!='5',]", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
        features = list(x = "GZMA", y = c("PRF1", "IFNG")),
        sample_filter = list(c("expDegree", "Expanded"))
      )
  )
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # Non_expanded
  fconfigs = list(
      list(result_id = "scatter_contour/sc_cd4_Non_expanded/",
        edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster!='5',]", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
        features = list(x = "GZMA", y = c("PRF1", "IFNG")),
        sample_filter = list(c("expDegree", "Non_expanded"))
      )
  )
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 4)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # PDCD1p
  fconfigs = list(
      list(result_id = "scatter_contour/sc_cd4_PDCD1p/",
        edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster!='5',]",
        features = list(x = "GZMA", y = c("PRF1", "IFNG")),
        sample_filter = list(c("PDCD1_tag", "PDCD1p"))
      )
  )
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

  # PDCD1n
  fconfigs = list(
      list(result_id = "scatter_contour/sc_cd4_PDCD1n/",
        edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster!='5',]",
        features = list(x = "GZMA", y = c("PRF1", "IFNG")),
        sample_filter = list(c("PDCD1_tag", "PDCD1n"))
      )
  )
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent")
  )

}

{ cat(redb("### Blanks volcano ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R") # volplot
  source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
  source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
  source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
  source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes
  source("/home/ciro/scripts/handy_functions/devel/volcano.R")

  fconfigs = list(
    list(file = here::here("results/dgea/CD8/CD8_Expansion/ExpandedvsNon_expanded/mastlog2cpm_results.csv"),
      group1 = "Expanded", group2 = "Non_expanded"
    ),
    list(file = here::here("results/dgea/CD8/CD8_PDCD1_tag/PDCD1pvsPDCD1n/mastlog2cpm_results.csv"),
      group1 = "PDCD1p", group2 = "PDCD1n",
      showgenes = c("GZMA", "GZMK", "GZMH", "CCL4", "CCL5", "CD74", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "RGS1", "CST7", "NKG7", "IFNG", "CXCR6", "TIGIT", "PRDM1", "ICOS", "LAG3", "CCL3", "CXCR3",
        "TNFSF9", "EOMES", "ITGAE", "XCL2", "TNF", "ITGA1", "CCR5", "FASLG", "NFATC2", "RBPJ", "CD44", "PRF1", "BATF", "SLAMF7", "FAS", "CXCR4", "IL7R", "KLF2", "TCF7", "CD55", "CD7", "TMEM123", "FLT3LG", "LEF1",
        "SELL", "FAM65B", "CCR7")
    ),
    list(file = here::here("results/dgea/CD4/CD4_nonTREG_Expansion/ExpandedvsNon_expanded/mastlog2cpm_results.csv"),
      group1 = "Expanded", group2 = "Non_expanded",
      showgenes = c("CCL4","CCL5","GZMA","GZMH","GZMK","IFNG","TBX21","PDCD1","TOX","LAG3","RGS1","CXCR6","RBPJ","ITGAE","PRDM1","RUNX2","HOPX","BATF","RUNX3","HLA-DRB1","NKG7","THEMIS","TCF7","LEF1","IL7R","SELL","KLF2","CD7","CD55")
      # showgenes = c("GZMA","GZMH","GZMK","PRDM1","HLA-DRB1","IFNG","CCL4","CCL5","IFNG","CXCR6","PDCD1","RGS1","CST7","CD2","CXCR3","CD74","KLRB1","IL2RB","RBPJ","CTLA4","CD40LG","CXCL13","BATF","KLF2","SELL","FAM65B","CCR7","TCF7","MAL","LEF1","CD55","CD7","SOCS3","IL7R","LTB","FOS","JUNB","IER2")
    ),
    list(file = here::here("results/dgea/CD4/CD4_nonTREG_PDCD1_tag/PDCD1pvsPDCD1n/mastlog2cpm_results.csv"),
      group1 = "PDCD1p", group2 = "PDCD1n",
      showgenes = c("GZMA", "CCL4", "CCL5", "GZMH", "HLA-DRB1", "TOX", "CXCR6", "PRDM1", "IFNG", "ITGAE", "RBPJ", "PRF1", "ZEB2", "ICOS")
    ),
    list(file = here::here("results/dgea/CD8/CD8_LAG3_tag/LAG3pvsLAG3n/mastlog2cpm_results.csv"),
      group1 = "LAG3p", group2 = "LAG3n",
      showgenes = c("GZMA", "CCL4", "CCL5", "GZMH", "HLA-DRB1", "TOX", "CXCR6", "PRDM1", "IFNG", "ITGAE", "RBPJ", "PRF1", "ZEB2", "ICOS", "GZMB")
    ),
    list(file = here::here("results/dgea/CD8/CD8_LAG3_tag_PDCD1low_donors/LAG3pvsLAG3n/mastlog2cpm_results.csv"),
      group1 = "LAG3p", group2 = "LAG3n",
      showgenes = c("GZMA", "CCL4", "CCL5", "GZMH", "HLA-DRB1", "TOX", "CXCR6", "PRDM1", "IFNG", "ITGAE", "RBPJ", "PRF1", "ZEB2", "ICOS", "GZMB")
    )
  )

  # fconfig=fconfigs[[6]]
  for (fconfig in fconfigs){
    file <- fconfig$file
    showgenes <- fconfig$showgenes
    group1 <- fconfig$group1
    group2 <- fconfig$group2
    results <- read.csv(file)

    fcthr = 0.35
    padjthr = 0.05
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    # heatmap parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    # Remove Unwanted genes
    mygenes = grep("XIST|^RP", rownames(resdf), value = TRUE, invert = TRUE)
    resdf <- resdf[mygenes,]
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=10),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 50),
        head(resdf[genes2plot, "gene"], 50)))
    }

    cat("---------------------- Volcano -------------------------\n")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    source("/home/kmlanderos/scripts/handy_functions/devel/volplot.R")

    # pdf("CD8_covsex_volcano_blankq.pdf") #Change name for CD8
    # pdf(paste0(dirname(file),"/mastlog2cpm_volcano_blank.pdf"), 10, 10)
    pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots_FC", as.character(fcthr), "_blank.pdf"), 10, 10)
      print(volplot( # volplot, volplot_wnames
        resdf,
        pvalth = padjthr,
        lfcth = fcthr,
        pvaltype = "padj",
        lfctype = "log2FoldChange",
        col_feature = means,
        size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
        gene_name = "gene",
        group = "degs",
        legends = FALSE, #TRUE, #
        blank = TRUE, #FALSE, #
        check_genes = list(text = features_parse_ensembl(showgenes)),
        return_plot = TRUE,
        clipp = FALSE,
        verbose = verbose,
        dot_size_low = 1, # 0
        dot_size_high = 4 # 3
      ) + xlim(-4,4)
      )
      # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
      #   title = paste(group2, "(-) vs ", group1, "(+)"))
    dev.off()

    # With names (big dots)
    pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots_FC", as.character(fcthr), ".pdf"), 10, 10)
      print(volplot_wnames(
        resdf,
        pvalth = padjthr,
        lfcth = fcthr,
        pvaltype = "padj",
        lfctype = "log2FoldChange",
        col_feature = means,
        size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
        gene_name = "gene",
        group = "degs",
        # legends = FALSE, #TRUE, #
        # blank = TRUE, #FALSE, #
        check_genes = list(text = features_parse_ensembl(showgenes)),
        return_plot = TRUE,
        clipp = FALSE,
        verbose = verbose,
        dot_size_low = 1, # 0
        dot_size_high = 4 # 3
      ) + xlim(-4,4)
      )
      # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
      #   title = paste(group2, "(-) vs ", group1, "(+)"))
    dev.off()

  }

}

{ cat(redb("### Main dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  dir.create("dim_reduction")
  source("/home/ciro/scripts/clustering/R/plotting.R")

  sc_cd3p_cd8@meta.data$subcluster_tag.0.6 <- "Other"
  sc_cd3p_cd8@meta.data[rownames(gzmk_trm_cd8_20pct_30pc@meta.data), "subcluster_tag.0.6"] <- as.character(gzmk_trm_cd8_20pct_30pc@meta.data$RNA_snn_res.0.6)
  sc_cd3p_cd8@meta.data$subcluster_tag.0.8 <- "Other"
  sc_cd3p_cd8@meta.data[rownames(gzmk_trm_cd8_20pct_30pc@meta.data), "subcluster_tag.0.8"] <- as.character(gzmk_trm_cd8_20pct_30pc@meta.data$RNA_snn_res.0.8)
  sc_cd3p_cd8@meta.data <- sc_cd3p_cd8@meta.data %>% mutate(TRM_GZMK_subcluster_0.6 = case_when(
    subcluster_tag.0.6 == "3" ~ "TRM_1",
    subcluster_tag.0.6 == "4" ~ "TRM_2",
    subcluster_tag.0.6 == "5" ~ "TRM_3",
    subcluster_tag.0.6 %in% c("0", "1", "2", "6", "7") ~ "GZMK_HI",
    TRUE ~ "Other"
  ), TRM_GZMK_subcluster_0.8 = case_when(
    subcluster_tag.0.8 == "2" ~ "TRM_1",
    subcluster_tag.0.8 == "5" ~ "TRM_2",
    subcluster_tag.0.8 == "6" ~ "TRM_3",
    subcluster_tag.0.8 == "9" ~ "TRM_4",
    subcluster_tag.0.8 %in% c("0", "1", "3", "4", "7", "8", "10") ~ "GZMK_HI",
    TRUE ~ "Other"
  ))

  TRM_GZMK_subcluster_color_0.6 <- c(TRM_1 = "#ffa7d7", TRM_2 = "#ff6ede", TRM_3 = "#FFE6F7", GZMK_HI = "#72e5d2", Other = "gray85")
  TRM_GZMK_subcluster_color_0.8 <- c(TRM_1 = "#ffa7d7", TRM_2 = "#ff6ede", TRM_3 = "#FFE6F7", TRM_4 = "#FBC5C5", GZMK_HI = "#72e5d2", Other = "gray85")

  fconfigs = list(
    # CD8
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster",
      colors = cd8_cluster_color,
      redu = redu[[1]]
    ),
    # CD4
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster",
      colors = cd4_cluster_color,
      redu = redu[[1]]
    ),
    # CD8 TRM/GZMK subclusters in main UMAP
    list(result_id = "dim_reduction/sc_cd3p_cd8_gzmk_trm_0.6",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "TRM_GZMK_subcluster_0.6",
      colors = TRM_GZMK_subcluster_color_0.6,
      redu = redu[[1]]
    ),
    # CD8 TRM/GZMK subclusters in main UMAP
    list(result_id = "dim_reduction/sc_cd3p_cd8_gzmk_trm_0.8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "TRM_GZMK_subcluster_0.8",
      colors = TRM_GZMK_subcluster_color_0.8,
      redu = redu[[1]]
    ),
    # CD8 paper clusters
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "paper_cluster",
      colors = cd8_paper_cluster_color,
      redu = redu[[1]]
    )
  )

  pp_main_dr = fig_plot_base(
    fconfigs[5],
    theme_extra = function(x){
      p <- x + geom_point(size = 0.4, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2", color = NULL) +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left")) #+
      # xlim(-10, 10) + ylim(-10, 10)
      # Seurat::LabelClusters(plot = p, id = "cluster", color = "black") # , fontface = "bold"  , point.padding = 1, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc"))
    }
  )

}

{ cat(redb("### TCR plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source('/home/ciro/scripts/handy_functions/devel/overlap.R') # fig_plot_overlaps
  #source("/home/kmlanderos/scripts/figease/figease.R")
  dir.create("tcr")

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data # [!sc_cd3p_cd8@meta.data$cluster %in% c("6", "10", "16"),]
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data # [!sc_cd3p_cd4@meta.data$cluster %in% c("12", "14"),]

  pie_fun = function(x){
      axes <- c(rlang::as_name(x$mapping$y), rlang::as_name(x$mapping$x))
      y <- reshape2::melt(table(x$data[, axes]))
      y[, 2] <- factor(as.character(y[, 2]), levels(x$data[, axes[2]]))
      data.table::setorderv(y, cols = rev(axes), order = c(1, -1))
      yy <- y %>% group_by(.dots = axes[2]) %>%
        mutate(prop = value / sum(value) * 100) %>%
        mutate(ypos = cumsum(prop) - (0.5*prop) ) %>% as.data.frame() %>%
        ggplot(aes(x="", y=prop, fill=TCR.tag)) +
        geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) +
        facet_wrap(paste0("~", axes[2])) + theme_void() +
        geom_text(aes(y = ypos, label = round(prop, 1)), color = "white") +
        scale_fill_brewer(palette="Set1")
  }
  bar_fun = function(x){
    axes <- c(rlang::as_name(x$mapping$y), rlang::as_name(x$mapping$x))
    y <- reshape2::melt(table(x$data[, axes]))
    y[, 2] <- factor(as.character(y[, 2]), rev(levels(x$data[, axes[2]])))
    data.table::setorderv(y, cols = rev(axes), order = c(1, -1))
    y %>% group_by(.dots = axes[2]) %>%
      mutate(prop = value / sum(value) * 100, ypos = cumsum(prop) - (0.5*prop)) %>%
      ggplot(aes(x=!!rlang::sym(axes[2]), y=value, fill=!!rlang::sym(axes[1]))) +
      geom_bar(stat='identity', position='fill', width=0.8) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      coord_flip() +
      geom_text(aes(y = ypos/100, label = round(prop, 1)), color = "white") +
      scale_fill_brewer(palette = "Set1") + labs(x = NULL, y = NULL, fill = NULL)
  }

  fconfigs = list(
    list(result_id = "tcr/sc_cd3p_cd8/", edata = "sc_cd3p_cd8@assays$RNA@data",
      metadata = "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag), ]",
      axis_x = "cluster_deconv", axis_y = "clonotype.tag", # filters = c("celltype", "CD8"),
      # col = c("GZMK" = "#54D8FD", "TRM" = "#FFB94F", "INNATE-LIKE_NKR" = "#D655A7", "INNATE-LIKE_XCL1-2" = "#44A306", "TCM" = "#969696", "MAIT" = "#C994C7", "Cell_Cycle" = "#464E2E"),
      col = cd8_cell_classification_color,
      vars2col = c("cell_classification_deconv"),
      element_type = c("clones", "cells"),
      clust_ord = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI")
    ),
    list(result_id = "tcr/sc_cd3p_cd4/", edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag), ]",
      axis_x = "cluster_deconv", axis_y = "clonotype.tag", # filters = c("celltype", "CD4"),
      col = cd4_cell_classification_color,
      vars2col = c("cell_classification_deconv"),
      element_type = c("clones", "cells"),
      clust_ord = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM", "TN")
    )
  )
  # ## 6) Upsetplot: Clones shared  ## ---------------------------
  # pp_overlaps = fig_plot_overlaps(fconfigs)

  # ## 6) Upsetplot: Clones shared between cell types (custom) ## ---------------------------
  if(!exists("gcolours")) gcolours = NULL
  dfplot <- mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag), ]
  verbose=T

  sharing_data_list <- list()
  for(i in 1:length(fconfigs[1])){ # fconfig=fconfigs[[1]]; i=1
    # Setting data
    fconfig = fig_config_check(fconfigs[[i]])
    fig_data = fig_set_data(fconfig,
      return_pdata = TRUE, verbose = verbose)
    ddf = fig_data$pdata
    axes = c(fconfig$axis_x[[1]][[1]], fconfig$axis_y[[1]][[1]])
    # feats <- fconfig$features[fconfig$features %in% fig_data$features]

    vars2col_data_list <- list()
    for(i in unique(c(fconfig$vars2col, "cluster_deconv"))){ # i="cell_classification";  i="cluster_deconv"
      if(verbose) cat("-", i, "\n")
      fname0 <- paste0(fig_data$name, gsub("orig\\.", "", i))
      if(!is.numeric(ddf[, i])) ddf[, i] <- droplevels(factor(ddf[, i]))
      elist = base::split(x = ddf[[axes[2]]], f = ddf[[i]])
      cols_i = v2cols(names(elist), fconfig$col)
      if(is.null(fconfig$element_type)) fconfig$element_type = c("collapsed", "complete")

      element_type_data_list <- list()
      for(type in fconfig$element_type){ # type = fconfig$element_type[2]
        if(verbose) cat(" *", type, "\n")
        elist_i = if(grepl("cells|complete", type)){ # clones size as the overlap
          lapply(X = setNames(nm = names(elist)), FUN = function(x){
            cells_with_clones <- ddf[[axes[2]]] %in% elist[[x]]
            rownames(ddf)[which(cells_with_clones)]
          })
        }else{ elist }

        fname <- paste0(fconfig$result_id, gsub("orig|\\.", "", i), "_", "custom_", type, "_")
        cat("* Upset\n")
        uddf <- as.data.frame(ComplexHeatmap::list_to_matrix(elist_i))
        pdf(paste0(fname, "upset.pdf"), onefile = FALSE, width = 10)
        if(i == "cell_classification_deconv") print(UpSetR::upset(data = uddf, sets.bar.color = rev(fconfig$col), sets = rev(fconfig$clust_ord), keep.order = TRUE)) else print(UpSetR::upset(data = uddf, sets = names(elist_i), keep.order = TRUE))
        graphics.off()
        pdf(paste0(fname, "upset_decreasing.pdf"), onefile = FALSE, width = 10)
        if(i == "cell_classification_deconv") print(UpSetR::upset(data = uddf, sets = rev(fconfig$clust_ord), order.by = "freq")) else print(UpSetR::upset(data = uddf, sets = names(elist_i), order.by = "freq"))
        graphics.off()

        element_type_data_list <- append(element_type_data_list, list(uddf))
      }; names(element_type_data_list) <- fconfig$element_type
      vars2col_data_list <- append(vars2col_data_list, list(element_type_data_list))
    }; names(vars2col_data_list) <- unique(c(fconfig$vars2col, "cluster_deconv"))
    sharing_data_list <- append(sharing_data_list, list(vars2col_data_list))
  }

  ### Get the number of cells of the celltypes of the combinations
  df <- copy(sharing_data_list[[1]]$cell_classification_deconv$cells)
  df <- df[rowSums(df)>1,] # Keep sahring in at least 2 groups
  df <- setDT(df)[,list(Count=.N),names(df)]
  df <- df[order(Count,decreasing=TRUE),]
  df <- df[Count > 100,]

  tmp <- apply(df, MARGIN=1, function(x) names(df)[as.logical(x)])
  combinations <- lapply(tmp, function(x) x[x!="Count"])

  dir.create("tcr/sc_cd3p_cd8/sharing_contributions_cell_classification_deconv")
  df2 <- copy(sharing_data_list[[1]]$cell_classification_deconv$cells)
  for(combination in combinations){ # combination <- combinations[[5]]
    vec <- c("CD16p_Effector" = 0, "Cell_Cycle" = 0, "Effector" = 0, "GZMK_HI" = 0, "MAIT" = 0, "TCF7_HI" = 0, "TRM" = 0)
    vec[combination] <- 1
    # Get the cells with the desired combination
    cells <- apply(df2, MARGIN=1, FUN = function(x) all(x == vec)); cells <- cells[cells]
    contribution <- table(mdata_sc_cd3p_cd8[names(cells),"cell_classification_deconv"])
    cat(combination, "\n", contribution, "\n")
    write.csv(contribution, paste0("tcr/sc_cd3p_cd8/sharing_contributions_cell_classification_deconv/", paste(combination, collapse="_"), ".csv"), row.names = FALSE, col.names = FALSE)
  }

  ### Get the number of cells of the celltypes of the combinations
  df <- copy(sharing_data_list[[1]]$cluster_deconv$cells)
  df <- df[rowSums(df)>1,] # Keep sahring in at least 2 groups
  df <- setDT(df)[,list(Count=.N),names(df)]
  df <- df[order(Count,decreasing=TRUE),]
  df <- df[Count > 100,]

  tmp <- apply(df, MARGIN=1, function(x) names(df)[as.logical(x)])
  combinations <- lapply(tmp, function(x) x[x!="Count"])

  dir.create("tcr/sc_cd3p_cd8/sharing_contributions_cluster_deconv")
  df2 <- copy(sharing_data_list[[1]]$cluster_deconv$cells)
  for(combination in combinations){ # combination <- combinations[[5]]
    vec <- c("10" = 0, "3" = 0, "4" = 0, "5" = 0, "6" = 0, "9" = 0, "GZMK_HI.0" = 0, "GZMK_HI.1" = 0, "GZMK_HI.2" = 0, "GZMK_HI.6" = 0, "GZMK_HI.7" = 0, "TRM.3" = 0, "TRM.4" = 0, "TRM.5" = 0)
    vec[combination] <- 1
    # Get the cells with the desired combination
    cells <- apply(df2, MARGIN=1, FUN = function(x) all(x == vec)); cells <- cells[cells]
    contribution <- table(mdata_sc_cd3p_cd8[names(cells),"cluster_deconv"])
    cat(combination, "\n", contribution, "\n")
    write.csv(contribution, paste0("tcr/sc_cd3p_cd8/sharing_contributions_cluster_deconv/", paste(combination, collapse="_"), ".csv"), row.names = FALSE, col.names = FALSE)
  }


  ## Scatter plot (umap) for cluster clone size ## ----------------------------------
  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$sufix <- "cluster_clone_size_v2_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification_deconv";# fconfigs[[i]]$vars2col = "cluster"
  }

  # Create the "cluster_clones_size" column in the metadata of the 2 objects.
  clonotype_cd4 <- sc_cd3p_cd4@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd4) <- clonotype_cd4$clonotype.tag; clonotype_cd4$clonotype.tag <- NULL

  clonotype_cd8 <- sc_cd3p_cd8@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster_deconv) %>% group_by(clonotype.tag, cluster_deconv) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster_deconv) %>%
    pivot_wider(names_from = cluster_deconv, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd8) <- clonotype_cd8$clonotype.tag; clonotype_cd8$clonotype.tag <- NULL

  # Update CD4 mdata
  cluster_clone_size <- c()
  for (i in 1:dim(mdata_sc_cd3p_cd4)[1]){
    clonotype <- mdata_sc_cd3p_cd4[i,"clonotype.tag"]
    c <- as.character(mdata_sc_cd3p_cd4[i,"cluster"])
    cluster_clone_size <- c(cluster_clone_size, clonotype_cd4[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  mdata_sc_cd3p_cd4$cluster_clone_size <- as.numeric(cluster_clone_size) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(cluster_clone_size = case_when(
    cluster_clone_size >= 30 ~ 30,
    TRUE ~ cluster_clone_size
  ))

  # Update CD8 mdata
  cluster_clone_size <- c()
  for (i in 1:dim(mdata_sc_cd3p_cd8)[1]){
    clonotype <- mdata_sc_cd3p_cd8[i,"clonotype.tag"]
    c <- as.character(mdata_sc_cd3p_cd8[i,"cluster_deconv"])
    cluster_clone_size <- c(cluster_clone_size, clonotype_cd8[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  mdata_sc_cd3p_cd8$cluster_clone_size <- as.numeric(cluster_clone_size) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(cluster_clone_size = case_when(
    cluster_clone_size >= 40 ~ 40,
    TRUE ~ cluster_clone_size
  ))

  mdata_sc_cd3p_cd4$cluster <- factor(mdata_sc_cd3p_cd4$cluster, levels = unique(mdata_sc_cd3p_cd4$cluster))
  mdata_sc_cd3p_cd8$cluster_deconv <- factor(mdata_sc_cd3p_cd8$cluster_deconv, levels = unique(mdata_sc_cd3p_cd8$cluster_deconv))
  pp_clones = fig_plot_base(
    fconfigs[1], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        guides(size = "none") +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray35", alpha = 0.1, stroke = 0.4) + # gray47
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6)), size = guide_legend(override.aes = list(color = "black", shape = 1, fill = "black"))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  pp_clones = fig_plot_base(
    fconfigs[2], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray35", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 30), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6)), size = guide_legend(override.aes = list(color = "black", shape = 1, fill = "black"))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  # -- Only clonally-expanded
  fconfigs2 <- fconfigs
  fconfigs2[[1]]$sufix = "cluster_clone_size_v2_onlyExp_"
  fconfigs2[[1]]$metadata = "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag) & mdata_sc_cd3p_cd8$clon.size.tag != 1, ]"
  pp_clones = fig_plot_base(
    fconfigs2[1], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        guides(size = "none") +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray35", alpha = 0.1, stroke = 0.4) + # gray47
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6)), size = guide_legend(override.aes = list(color = "black", shape = 1, fill = "black"))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  fconfigs2[[2]]$sufix = "cluster_clone_size_v2_onlyExp_"
  fconfigs2[[2]]$metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag) & mdata_sc_cd3p_cd4$clon.size.tag != 1, ]"
  pp_clones = fig_plot_base(
    fconfigs2[2], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray35", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6)), size = guide_legend(override.aes = list(color = "black", shape = 1, fill = "black"))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  # Inset plot
  # -- CD8
  pdf(paste0(fconfigs[[1]]$result_id, "inset_plot_clone_size.pdf"))
  a <- sc_cd3p_cd8@meta.data %>% select(clonotype.tag, cell_classification, expDegree) %>% filter(!is.na(clonotype.tag)) %>% group_by(cell_classification) %>%
    summarize(pct.exp = 100*sum(expDegree == "Expanded")/length(expDegree)) %>% data.frame()
  a$category <- "Expansion"; a$cell_classification <- factor(a$cell_classification, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification, y = category, color = cell_classification)) +
    geom_point(aes(size = pct.exp), alpha = 0.8, sho1w.legend = TRUE) +
    scale_radius(breaks = c(10, 30, 50, 70)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_manual(values = c("Cell_Cycle" = "#9c0006", "GZMK_HI" = "#3bb8ae", "TRM" = "#ff6ede", "CD16p_Effector" = "#fcee62", "Effector" = "#ff971f", "MAIT" = "#9e89e8", "TCF7_HI" = "#4da7fd"))
  print(p)
  dev.off()
  # --- CD4
  pdf(paste0(fconfigs[[2]]$result_id, "inset_plot_clone_size.pdf"))
  a <- sc_cd3p_cd4@meta.data %>% select(clonotype.tag, cell_classification, expDegree) %>% filter(!is.na(clonotype.tag)) %>% group_by(cell_classification) %>%
    summarize(pct.exp = 100*sum(expDegree == "Expanded")/length(expDegree)) %>% data.frame()
  a$category <- "Expansion"; a$cell_classification <- factor(a$cell_classification, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM_TN"))
  p <- a %>% ggplot(aes(x = cell_classification, y = category, color = cell_classification)) +
    geom_point(aes(size = pct.exp), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 20, 40, 60, 80)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_manual(values = c("Cell_Cycle" = "#9c0006", "CTL" = "#a084ff", "TFH" = "#ff717c", "TREG" = "#ffa7d7", "TCM_TN" = "#97ccf5"))
  print(p)
  dev.off()

  #########################################
  # Stacked Barplots Clonotypes per donor #
  #########################################

  output.dir <- "tcr"
  srt.objs.list <- list(sc_cd3p_cd4, sc_cd3p_cd8); names(srt.objs.list) <- c("CD4", "CD8")
  table(sc_cd3p_cd8@meta.data$orig.donor)
  table(sc_cd3p_cd8@meta.data[!is.na(sc_cd3p_cd8@meta.data$clonotype.tag), "orig.donor"])
  # --------> Clone Size

  all.clonotypes <- lapply(X = srt.objs.list, FUN = function(obj){ # obj = sc_cd3p_cd4
    meta.data <- as.data.table(obj@meta.data)

    tmp.data.1 <- meta.data[!is.na(orig.donor) & !is.na(clonotype.tag), .(cell.count=.N), by=.(donor=orig.donor, clonotype=clonotype.tag)]
    setorderv(x=tmp.data.1, cols='cell.count', order=-1) # Maybe want to order by PDCD1_pct_byDonor as a 2nd variable
    uniq.donors <- tmp.data.1[, unique(donor)]
    # This function might not be necessary.
    tmp.all.clonotypes <- lapply(X=uniq.donors, FUN=function(tmp.donor){
      tmp.data <- tmp.data.1[donor==tmp.donor]
      tmp.data <- tmp.data[, .(donor, clonotype, cell.count)]
      tmp.data <- na.omit(tmp.data)
      # if(tmp.data[])
      return(tmp.data)
    })
    tmp.all.clonotypes <- rbindlist(l=tmp.all.clonotypes, use.names=TRUE, idcol=NULL)
    tmp.all.clonotypes[,donor:=factor(x = donor,
      levels = uniq.donors)]
  })
  # c("BT1", "BT26", "BT4", "BT24", "BT13", "BT20", "BT5", "BT21", "BT3", "BT12", "BT15", "BT19", "BT23", "BT44", "BT7", "BT8", "BT25", "BT46", "BT18", "BT11", "BT9", "BT27", "BT2", "BT17_brain", "BT22", "BT40", "BT32", "BT39", "BT10", "BT31", "BT37", "BT35", "BT30", "BT38", "BT45")

  # Save CD4 and CD8 XL tables
  write.csv(all.clonotypes[["CD4"]], file = "tcr/cellCount_allClones_PBT_CD4.csv", row.names = FALSE)
  write.csv(all.clonotypes[["CD8"]], file = "tcr/cellCount_allClones_PBT_CD8.csv", row.names = FALSE)

  ################ - Collapsing clonotypes with cs == 1

  # df <- all.clonotypes[["CD4"]]
  all.clonotypes.1 <- lapply(X = all.clonotypes, FUN = function(df){
    df.new <- df %>% filter(cell.count != 1)
    df.tmp <- df %>% filter(cell.count == 1) %>% group_by(donor) %>% summarize(
      cell.count = sum(cell.count)
    )
    df.tmp$clonotype <- "clonotypes_nonExp"; df.tmp <- df.tmp[,c("donor", "clonotype", "cell.count")]
    df.final <- rbind(df.new, df.tmp)
    df.final
  })

  # Plot

  # donors <- c("BT1", "BT26", "BT4", "BT24", "BT13", "BT20", "BT5", "BT21", "BT3", "BT12", "BT15", "BT19", "BT23", "BT44", "BT7", "BT8", "BT25", "BT46", "BT18", "BT11", "BT9", "BT27", "BT2", "BT17_brain",
  #   "BT22", "BT40", "BT32", "BT39", "BT10", "BT31", "BT37", "BT35", "BT30", "BT38", "BT45")

  # CD4
  tmp <- all.clonotypes.1[["CD4"]] %>% unite(col = donor_clonotype, donor, clonotype, sep = "_", remove = FALSE)
  tmp.1 <- tmp %>% filter(clonotype != "clonotypes_nonExp") %>% arrange(desc(cell.count))
  tmp.2 <- tmp %>% filter(clonotype == "clonotypes_nonExp") %>% arrange(-desc(cell.count))
  tmp.final <- rbind(tmp.1, tmp.2)
  tmp.final$donor_clonotype <- factor(tmp.final$donor_clonotype, levels = tmp.final$donor_clonotype)

  tmp.matrix <- apply(tmp.final, MARGIN = 1, function(x) { # x <- tmp.final[2,]
    nCells <- x["cell.count"]
    tmp <- Reduce("rbind", replicate(nCells, x[c("donor", "clonotype", "donor_clonotype")], simplify = FALSE))
  })
  CD4.matrix <- as.data.frame(Reduce("rbind", tmp.matrix))
  CD4.matrix$donor_clonotype <- factor(CD4.matrix$donor_clonotype, levels = tmp.final$donor_clonotype)
  donors <- unique(sc_cd3p_cd4$orig.donor)[!is.na(unique(sc_cd3p_cd4$orig.donor))]
  # CD4.matrix$donor <- factor(CD4.matrix$donor, levels = rev(donors))

  # CD8
  tmp <- all.clonotypes.1[["CD8"]] %>% unite(col = donor_clonotype, donor, clonotype, sep = "_", remove = FALSE)
  tmp.1 <- tmp %>% filter(clonotype != "clonotypes_nonExp") %>% arrange(desc(cell.count))
  tmp.2 <- tmp %>% filter(clonotype == "clonotypes_nonExp") %>% arrange(-desc(cell.count))
  tmp.final <- rbind(tmp.1, tmp.2)
  tmp.final$donor_clonotype <- factor(tmp.final$donor_clonotype, levels = tmp.final$donor_clonotype)

  tmp.matrix <- apply(tmp.final, MARGIN = 1, function(x) { # all.clonotypes.1[["CD8"]]
    nCells <- x["cell.count"]
    tmp <- Reduce("rbind", replicate(nCells, x[c("donor", "clonotype", "donor_clonotype")], simplify = FALSE))
  })
  CD8.matrix <- as.data.frame(Reduce("rbind", tmp.matrix))
  CD8.matrix$donor_clonotype <- factor(CD8.matrix$donor_clonotype, levels = tmp.final$donor_clonotype)
  donors <- unique(sc_cd3p_cd8$orig.donor)[!is.na(unique(sc_cd3p_cd4$orig.donor))]
  # CD8.matrix$donor <- factor(CD8.matrix$donor, levels = rev(donors))

  # donors <- c("BT1", "BT7", "BT32", "BT19", "BT3", "BT44", "BT27", "BT35", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT38", "BT47", "BT22", "BT40", "BT30", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT37", "BT39", "BT23", "BT13", "BT20", "BT17_brain", "BT31", "BT46", "BT2", "BT45", "BT48", "BT29", "BT43")
  donors <- c("BT1", "BT7", "BT32", "BT19", "BT3", "BT44", "BT27", "BT35", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT38", "BT47", "BT22", "BT40", "BT30", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT37", "BT39", "BT23", "BT13", "BT20", "BT17_brain", "BT31", "BT46", "BT2", "BT45", "BT48", "BT29", "BT43")
  CD4.matrix <- CD4.matrix %>% filter(donor %in% donors)
  CD8.matrix <- CD8.matrix %>% filter(donor %in% donors)
  CD4.matrix$donor <- factor(CD4.matrix$donor, levels = rev(donors))
  CD8.matrix$donor <- factor(CD8.matrix$donor, levels = rev(donors))
  donors_exclude <- c("BT43")


  # - Plot
  pdf(paste0("tcr/cellCount_allClones_PBT_CD4_CD8_nonExpCollapsed_v5.pdf"), 5, 6)
  tmp.y.axis <- "CD4+ cells"
  tmp.x.axis <- "Patient"
  # n_non_Exp <- sum(grepl("clonotypes_nonExp", cd4.matrix$donor_clonotype))
  p1 <- ggplot(CD4.matrix, aes(x = donor, fill = donor_clonotype)) + geom_bar(position = "stack", color = "cadetblue4", size = 0.0001) + # gray60
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_reverse(limits = c(2000,0)) + coord_flip() +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(CD4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(CD4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"), legend.position="none")

  tmp.y.axis <- "CD8+ cells"
  p2 <- ggplot(CD8.matrix, aes(x = donor, fill = donor_clonotype)) + geom_bar(position = "stack", color = "royalblue4", size = 0.0001) + # gray60
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_continuous(limits = c(0,2000)) + coord_flip() +
    scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(CD8.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(CD8.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y = element_blank(), axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"), legend.position="none")

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
  dev.off()


}

{ cat(redb("### GLIPH2 plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("tcr/gliph")

  donor_order <- c("BT1", "BT2", "BT3", "BT4", "BT5", "BT7", "BT8", "BT9", "BT10", "BT11", "BT12", "BT13", "BT15", "BT17_brain", "BT18", "BT19", "BT20", "BT21", "BT22", "BT23", "BT24", "BT25", "BT26", "BT27", "BT29",
    "BT30", "BT31", "BT32", "BT35", "BT37", "BT38", "BT39", "BT40", "BT43", "BT44", "BT45", "BT46", "BT47", "BT48")

  ############################################################
  # NUMBER of Specificity groups & clones per donor (tables)
  ############################################################

  # ---> All

  # CD4
  tmp <- sc_cd3p_cd4@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor)) %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
  tmp.1 <- tmp %>% filter(str_detect(pattern, "_", negate = TRUE))
  tmp.2 <- tmp %>% filter(str_detect(pattern, "_"))
  donor_ct_2 <- c()
  for (i in 1:nrow(tmp.2)){
    motifs  <- unlist(str_split(tmp.2[i, "pattern"], "_"))
    for(motif in motifs){
      donor_ct_2 <- rbind(donor_ct_2, c(as.character(tmp.2[i, "orig.donor"]), motif))
    }
  }; colnames(donor_ct_2) <- c("orig.donor", "pattern")
  df <- rbind(tmp.1, donor_ct_2) %>% mutate(orig.donor = factor(orig.donor, levels = donor_order))
  # SG per donor
  df_cd4 <- df %>% group_by(orig.donor) %>% summarize(n_SG = length(unique(pattern))) %>% as.data.frame()
  df_cd4

  # CD8
  tmp <- sc_cd3p_cd8@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor)) %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
  tmp.1 <- tmp %>% filter(str_detect(pattern, "_", negate = TRUE))
  tmp.2 <- tmp %>% filter(str_detect(pattern, "_"))
  donor_ct_2 <- c()
  for (i in 1:nrow(tmp.2)){
    motifs  <- unlist(str_split(tmp.2[i, "pattern"], "_"))
    for(motif in motifs){
      donor_ct_2 <- rbind(donor_ct_2, c(as.character(tmp.2[i, "orig.donor"]), motif))
    }
  }; colnames(donor_ct_2) <- c("orig.donor", "pattern")
  df <- rbind(tmp.1, donor_ct_2) %>% mutate(orig.donor = factor(orig.donor, levels = donor_order))
  # SG per donor
  df_cd8 <- df %>% group_by(orig.donor) %>% summarize(n_SG = length(unique(pattern))) %>% as.data.frame()
  df_cd8

  # Save tables
  write.table(df_cd4, "tcr/gliph/specificityGroups_per_Donor_CD4_simplified.csv", sep = ",", row.names = F, quote = F)
  write.table(df_cd8, "tcr/gliph/specificityGroups_per_Donor_CD8_simplified.csv", sep = ",", row.names = F, quote = F)

  # # ------------------------- Number of Expanded specificity groups -------------------------

  GLIPH2Output_CD4 <- read.csv(here::here("results", "tcr/gliph/2023-04-19/GLIPH2Output_TCR-Data_CD4_hto.csv")) # 2022-02-07
  GLIPH2Output_CD4 <- GLIPH2Output_CD4 %>% mutate(id = paste(GLIPH2Output_CD4$TcRb, GLIPH2Output_CD4$V, GLIPH2Output_CD4$J, sep = "_"))
  tmp <- GLIPH2Output_CD4 %>% group_by(pattern) %>% summarize(exp_tag = any(Freq >1))
  sum(tmp$exp_tag)

  GLIPH2Output_CD8 <- read.csv(here::here("results", "tcr/gliph/2023-04-19/GLIPH2Output_TCR-Data_CD8_hto.csv")) # 2022-02-07
  GLIPH2Output_CD8 <- GLIPH2Output_CD8 %>% mutate(id = paste(GLIPH2Output_CD8$TcRb, GLIPH2Output_CD8$V, GLIPH2Output_CD8$J, sep = "_"))
  tmp <- GLIPH2Output_CD8 %>% group_by(pattern) %>% summarize(exp_tag = any(Freq >1))
  sum(tmp$exp_tag)
}

{ cat(redb("### GG-paired plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("paired_plots")

  # ---> PDCD1n and PDCD1p Gene expression (CD8) - MANY GENES

  tmp <- as.matrix(sc_cd3p_cd8@assays$RNA["PDCD1",])
  cells <- colnames(sc_cd3p_cd8@assays$RNA["PDCD1",])
  PDCD1p_cells <- cells[tmp>0]
  PDCD1n_cells <- cells[tmp==0]

  genes <- c("CCL4", "GZMK", "GZMB", "GZMA", "HLA-DRB1", "HLA-DPA1", "HLA-DRA", "CCL5", "COTL1", "TOX", "ITM2A", "IFNG", "CD2", "TIGIT",
  "PRDM1", "ICOS", "LAG3", "CCL3", "EOMES", "XCL2", "TNF", "ITGAE", "CCR5", "FASLG", "ARAP2", "CD70", "IL32", "CD44", "TANK", "PRF1",
  "SLAMF7", "CRTAM", "TNFRSF9", "TNFSF9", "CCR2", "CCR4", "CTLA4", "MAF", "IL2")
  donors <- unique(sc_cd3p_cd8@meta.data$orig.donor)[!is.na(unique(sc_cd3p_cd8@meta.data$orig.donor))]

  df <- c()
  # donor <- "BT5"
  for (donor in donors){
    donor_cells <- rownames(sc_cd3p_cd8@meta.data[sc_cd3p_cd8@meta.data$orig.donor == donor & !is.na(sc_cd3p_cd8@meta.data$orig.donor),])
    # PDCD1p cells
    PDCD1p_exp <- sc_cd3p_cd8@assays$RNA[genes,intersect(PDCD1p_cells,donor_cells)]
    PDCD1p_exp <- apply(PDCD1p_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    # PDCD1n cells
    PDCD1n_exp <- sc_cd3p_cd8@assays$RNA[genes,intersect(PDCD1n_cells,donor_cells)]
    PDCD1n_exp <- apply(PDCD1n_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    tmp <- data.frame(donor,PDCD1p_exp,PDCD1n_exp, genes)
    df <- rbind(df,tmp)
  }
  colnames(df) <- c("donor","PDCD1p","PDCD1n","genes")

  # Plot
  pdf("paired_plots/paired_dots_CD8_PDCD1_gene_expression_pct.pdf", 15, 15)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", facet.by = "genes", palette = "jco", ylab = "Percentage")
  dev.off()
  write.csv(df, file = "paired_plots/paired_dots_CD8_PDCD1_gene_expression_pct.csv", row.names = FALSE)

  # ---> PDCD1n and PDCD1p clonal expansion percentage (%) (CD8)

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))
  mdata_sc_cd3p_cd8$expansion_tag <- ifelse(mdata_sc_cd3p_cd8$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(orig.donor)) %>% group_by(orig.donor, PDCD1_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, PDCD1_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = PDCD1_tag, values_from = pct_exp_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD8_PDCD1_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD8_PDCD1_expansion_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # ---> Expanded and Non-Expanded cells PDCD1 expression (CD8)

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  mdata_sc_cd3p_cd8$expansion_tag <- sc_cd3p_cd8$expDegree
  # mdata_sc_cd3p_cd8$expansion_tag <- ifelse(mdata_sc_cd3p_cd8$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor)) %>% group_by(orig.donor, expansion_tag, PDCD1_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, expansion_tag) %>% summarize(pct_PDCD1_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = expansion_tag, values_from = pct_PDCD1_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD8_expansion_PDCD1_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD8_expansion_PDCD1_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "Non_expanded", cond2 = "Expanded", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of cells expressing PDCD1")
  dev.off()

  # CD4

  # ---> PDCD1n and PDCD1p Gene expression (nonTREG CD4) - MANY GENES

  tmp <- as.matrix(sc_cd3p_cd4@assays$RNA["PDCD1",])
  cells <- colnames(sc_cd3p_cd4@assays$RNA["PDCD1",])
  PDCD1p_cells <- cells[tmp>0]
  PDCD1n_cells <- cells[tmp==0]

  genes <- c("CCL4", "GZMK", "GZMB", "GZMA", "HLA-DRB1", "HLA-DPA1", "HLA-DRA", "CCL5", "COTL1", "TOX", "ITM2A", "IFNG", "CD2", "TIGIT",
  "PRDM1", "ICOS", "LAG3", "CCL3", "EOMES", "XCL2", "TNF", "ITGAE", "CCR5", "FASLG", "ARAP2", "CD70", "IL32", "CD44", "TANK", "PRF1",
  "SLAMF7", "CRTAM", "TNFRSF9", "TNFSF9", "CCR2", "CCR4", "CTLA4", "MAF")
  donors <- unique(sc_cd3p_cd4@meta.data$orig.donor)[!is.na(unique(sc_cd3p_cd4@meta.data$orig.donor))]

  df <- c()
  # donor <- "BT5"
  for (donor in donors){
    donor_cells <- rownames(sc_cd3p_cd4@meta.data[sc_cd3p_cd4@meta.data$orig.donor == donor & !is.na(sc_cd3p_cd4@meta.data$orig.donor) & sc_cd3p_cd4@meta.data$cell_classification != "TREG",])
    # PDCD1p cells
    PDCD1p_exp <- sc_cd3p_cd4@assays$RNA[genes,intersect(PDCD1p_cells,donor_cells)]
    PDCD1p_exp <- apply(PDCD1p_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    # PDCD1n cells
    PDCD1n_exp <- sc_cd3p_cd4@assays$RNA[genes,intersect(PDCD1n_cells,donor_cells)]
    PDCD1n_exp <- apply(PDCD1n_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    tmp <- data.frame(donor,PDCD1p_exp,PDCD1n_exp, genes)
    df <- rbind(df,tmp)
  }
  colnames(df) <- c("donor","PDCD1p","PDCD1n","genes")

  # Plot
  pdf("paired_plots/paired_dots_CD4_PDCD1_gene_expression_pct.pdf", 15, 15)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", facet.by = "genes", palette = "jco", ylab = "Percentage")
  dev.off()
  write.csv(df, file = "paired_plots/paired_dots_CD4_PDCD1_gene_expression_pct.csv", row.names = FALSE)

  # ---> PDCD1n and PDCD1p clonal expansion percentage (%) (nonTREG CD4)

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd4@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd4@assays$RNA@data)
  mdata_sc_cd3p_cd4[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))
  mdata_sc_cd3p_cd4$expansion_tag <- ifelse(mdata_sc_cd3p_cd4$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd4 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd4$orig.donor) & cell_classification != "TREG") %>% group_by(orig.donor, PDCD1_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, PDCD1_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = PDCD1_tag, values_from = pct_exp_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD4_PDCD1_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD4_PDCD1_expansion_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # ---> Expanded and Non-Expanded cells PDCD1 expression (nonTREG CD4)

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd4@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd4@assays$RNA@data)
  mdata_sc_cd3p_cd4[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  mdata_sc_cd3p_cd4$expansion_tag <- sc_cd3p_cd4$expDegree
  # mdata_sc_cd3p_cd4$expansion_tag <- ifelse(mdata_sc_cd3p_cd4$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd4 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd4$orig.donor) & cell_classification != "TREG") %>% group_by(orig.donor, expansion_tag, PDCD1_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, expansion_tag) %>% summarize(pct_PDCD1_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = expansion_tag, values_from = pct_PDCD1_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD4_expansion_PDCD1_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD4_expansion_PDCD1_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "Non_expanded", cond2 = "Expanded", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of cells expressing PDCD1")
  dev.off()


  ######################################## LAG3 ################################################

  # ------------ CD8

  # ---> LAG3n and LAG3p Gene expression (CD8) - MANY GENES

  tmp <- as.matrix(sc_cd3p_cd8@assays$RNA["LAG3",])
  cells <- colnames(sc_cd3p_cd8@assays$RNA["LAG3",])
  LAG3p_cells <- cells[tmp>0]
  LAG3n_cells <- cells[tmp==0]

  genes <- c("PRF1", "GZMA", "GZMB", "GZMK", "GZMH", "IFNG", "TNF", "CCL4", "CCL5", "HLA-DRB1", "PDCD1", "EOMES", "PRDM1", "TOX", "ITGAE", "CCR5", "FASLG")
  donors <- unique(sc_cd3p_cd8@meta.data$orig.donor)[!is.na(unique(sc_cd3p_cd8@meta.data$orig.donor))]

  df <- c()
  # donor <- "BT5"
  for (donor in donors){
    donor_cells <- rownames(sc_cd3p_cd8@meta.data[sc_cd3p_cd8@meta.data$orig.donor == donor & !is.na(sc_cd3p_cd8@meta.data$orig.donor),])
    # LAG3p cells
    LAG3p_exp <- sc_cd3p_cd8@assays$RNA[genes,intersect(LAG3p_cells,donor_cells)]
    LAG3p_exp <- apply(LAG3p_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    # LAG3n cells
    LAG3n_exp <- sc_cd3p_cd8@assays$RNA[genes,intersect(LAG3n_cells,donor_cells)]
    LAG3n_exp <- apply(LAG3n_exp, MARGIN=1, FUN=function(x){
      pCells <- 100*sum(x>0)/length(x)
    })
    tmp <- data.frame(donor,LAG3p_exp,LAG3n_exp, genes)
    df <- rbind(df,tmp)
  }
  colnames(df) <- c("donor","LAG3p","LAG3n","genes")

  # # save table
  # write.csv(df, file = "paired_plots/paired_dots_CD8_LAG3_percentageExpr.csv", row.names = FALSE)
  # Plot
  pdf("paired_plots/paired_dots_CD8_LAG3_gene_expression_pct.pdf", 15, 15)
  ggpaired(df, cond1 = "LAG3n", cond2 = "LAG3p", line.color = "gray", #color = "condition",
    fill = "condition", facet.by = "genes", palette = "jco", ylab = "Percentage")
  dev.off()
  write.csv(df, file = "paired_plots/paired_dots_CD8_LAG3_gene_expression_pct.csv", row.names = FALSE)

  # ---> LAG3n and LAG3p clonal expansion percentage (%) (CD8)

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  LAG3_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["LAG3",]) > 0
  LAG3_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[LAG3_cells, "LAG3_tag_tmp"] <- LAG3_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(LAG3_tag = case_when(
    LAG3_tag_tmp > 0 ~ "LAG3p",
    LAG3_tag_tmp == 0 ~ "LAG3n"
  ))
  mdata_sc_cd3p_cd8$expansion_tag <- ifelse(mdata_sc_cd3p_cd8$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor)) %>% group_by(orig.donor, LAG3_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, LAG3_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = LAG3_tag, values_from = pct_exp_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD8_LAG3_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD8_LAG3_expansion_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "LAG3n", cond2 = "LAG3p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # ---> Expanded and Non-Expanded cells LAG3 expression (CD8)

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  LAG3_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["LAG3",]) > 0
  LAG3_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[LAG3_cells, "LAG3_tag_tmp"] <- LAG3_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(LAG3_tag = case_when(
    LAG3_tag_tmp > 0 ~ "LAG3p",
    LAG3_tag_tmp == 0 ~ "LAG3n"
  ))

  mdata_sc_cd3p_cd8$expansion_tag <- sc_cd3p_cd8$expDegree
  # mdata_sc_cd3p_cd8$expansion_tag <- ifelse(mdata_sc_cd3p_cd8$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor)) %>% group_by(orig.donor, expansion_tag, LAG3_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, expansion_tag) %>% summarize(pct_LAG3_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = expansion_tag, values_from = pct_LAG3_cells)

  write.csv(df, file = "paired_plots/paired_dots_CD8_expansion_LAG3_pct.csv", row.names = FALSE)

  # Plot
  pdf("paired_plots/paired_dots_CD8_expansion_LAG3_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "Non_expanded", cond2 = "Expanded", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of cells expressing LAG3")
  dev.off()


}

{ cat(redb("### Cluster Markers ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("cluster_markers")

  fnames <- paste0(
    "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/clustering/",
    c("CD45pCD3p_CD4_Batch1_Batch2_manuscript_v11pm_clean1_harmony/seurat_mean0.01_pct15_pc15_res0.6/dgea_MAST_fc0.25_padj0.05_summary_stats.csv",
      "CD45pCD3p_CD8_Batch1_Batch2_manuscript_clean1_harmony/seurat_mean0.01_pct20_pc30_res0.6/dgea_MAST_fc0.25_padj0.05_summary_stats.csv")
  );names(fnames) <- c("CD4", "CD8")

  ## ---> CD8

  celltype <- "CD8"
  n_features <- 50

  markers_cd8_paper_cluster <- dgea_seurat(
      object = sc_cd3p_cd8,
      results_prefix = paste0("cluster_markers", "/dgea_", "MAST"),
      config_markers = list(test="MAST", avg_logFC="0.25", "0.05"),
      cluster_column = "paper_cluster",
      return_table = TRUE,
      verbose = TRUE
    )
  cmarkers <- readRDS("cluster_markers/.dgea_MAST.rds")
  cmarkers <- markers_summary(
        marktab = cmarkers,
        annot = sc_cd3p_cd8@meta.data,
        datavis = GetAssayData(sc_cd3p_cd8),
        cluster_column = "paper_cluster",
        datatype = "SeuratNormalized",
        verbose = TRUE
      )
  write.table(cmarkers, "cluster_markers/paper_cluster_markers_CD8_summary_stats.csv", sep = ",")

  Idents(sc_cd3p_cd8) <- "paper_cluster"
  df <- FindAllMarkers(
          object = sc_cd3p_cd8,
          # assay = "RNA",
          # slot = "data",
          only.pos = TRUE,
          min.pct = 0.25,
          logfc.threshold = 0.25,
          min.diff.pct = 0.05,
          test.use = "MAST",
          return.thresh = 0.2,
          verbose = TRUE
        )
  # Idents(sc_cd3p_cd8) <- "cluster"
  # df <- FindAllMarkers(sc_cd3p_cd8, test.use = "MAST")

  df <- df %>% select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  top10 <- df %>%
    group_by(cluster) %>%
    top_n(n = n_features, wt = avg_logFC)
  write.table(top10, paste0("cluster_markers/paper_cluster_markers_", celltype, "_top", as.character(n_features), "features.csv"), sep = ",", row.names = F, quote = F)

  pdf(paste0("cluster_markers/paper_cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"))
  sc_cd3p_cd8$paper_cluster <- factor(sc_cd3p_cd8$paper_cluster, levels = names(cd8_paper_cluster_color))
  top10 <- top10 %>% arrange(factor(cluster, levels = names(cd8_paper_cluster_color)))
  DoHeatmap(sc_cd3p_cd8, features = top10$gene, slot = "scale.data", group.by = "paper_cluster", draw.lines = F, group.colors = cd8_paper_cluster_color) +
  scale_fill_gradientn(colors = c("blue", "black", "yellow")) #+ NoLegend()
  dev.off()


  ## ---> CD4

  celltype <- "CD4"
  n_features <- 50

  df <- read.csv(fnames[celltype], stringsAsFactors = FALSE, row.names = 1, check.names = FALSE) %>%
    select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  head(df); dim(df)

  top10 <- df %>%
    group_by(cluster) %>%
    top_n(n = n_features, wt = avg_logFC)
  write.table(top10, paste0("cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features.csv"), sep = ",", row.names = F, quote = F)

  pdf(paste0("cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"))
  sc_cd3p_cd4$cluster <- factor(sc_cd3p_cd4$cluster, levels = names(cd4_cluster_color))
  top10 <- top10 %>% arrange(factor(cluster, levels = names(cd4_cluster_color)))
  top10 <- top10 %>% filter(!cluster %in% c("9", "12"))
  DoHeatmap(sc_cd3p_cd4, features = top10$gene, slot = "scale.data", group.by = "cluster", draw.lines = F, group.colors = cd4_cluster_color) +
  scale_fill_gradientn(colors = c("blue", "black", "yellow")) #+ NoLegend()
  dev.off()

}

{ cat(redb("### Supplementary QC Plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("QC")

  # --- CD8
  df <- sc_cd3p_cd8@meta.data %>%
    mutate( paper_cluster=factor(paper_cluster,levels=names(cd8_paper_cluster_color)) )

  p1 <- ggplot(df, aes(x=paper_cluster, y=nCount_RNA, fill = paper_cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd8_paper_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  p2 <- ggplot(df, aes(x=paper_cluster, y=nFeature_RNA, fill = paper_cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd8_paper_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  p3 <- ggplot(df, aes(x=paper_cluster, y=percent.mt, fill = paper_cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd8_paper_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  df <- sc_cd3p_cd8@meta.data %>% select(paper_cluster, origlib) %>% group_by(paper_cluster, origlib) %>% summarize(n = n()) %>%
    mutate( paper_cluster=factor(paper_cluster,levels=names(cd8_paper_cluster_color)) )
  p4 <- ggplot(df, aes(fill=origlib, y=n, x=paper_cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Libraries", fill = "Library")

  df <- sc_cd3p_cd8@meta.data %>% select(paper_cluster, orig.Run) %>% group_by(paper_cluster, orig.Run) %>% summarize(n = n()) %>%
    mutate( paper_cluster=factor(paper_cluster,levels=names(cd8_paper_cluster_color)) )
  p5 <- ggplot(df, aes(fill=orig.Run, y=n, x=paper_cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Runs", fill = "Run")

  df <- sc_cd3p_cd8@meta.data %>% select(paper_cluster, orig.donor) %>% group_by(paper_cluster, orig.donor) %>% summarize(n = n()) %>%
    mutate( paper_cluster=factor(paper_cluster,levels=names(cd8_paper_cluster_color)) )
  p6 <- ggplot(df, aes(fill=orig.donor, y=n, x=paper_cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Donors", fill = "Donor")

  # Change violin plot colors by groups
  pdf("QC/qc_cd8.pdf", 10, 7)
  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
  dev.off()

  # --- CD4

  df <- sc_cd3p_cd4@meta.data %>%
    mutate( cluster=factor(cluster,levels=c("1", "2", "4", "5", "8",   "0", "3", "7",   "6", "10", "11")) )

  p1 <- ggplot(df, aes(x=cluster, y=nCount_RNA, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd4_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  p2 <- ggplot(df, aes(x=cluster, y=nFeature_RNA, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd4_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  p3 <- ggplot(df, aes(x=cluster, y=percent.mt, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = cd4_cluster_color) +
    theme_classic() +
    theme(legend.position="none")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster, origlib) %>% group_by(cluster, origlib) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("1", "2", "4", "5", "8",   "0", "3", "7",   "6", "10", "11")) )
  p4 <- ggplot(df, aes(fill=origlib, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Libraries", fill = "Library")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster, orig.Run) %>% group_by(cluster, orig.Run) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("1", "2", "4", "5", "8",   "0", "3", "7",   "6", "10", "11")) )
  p5 <- ggplot(df, aes(fill=orig.Run, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Runs", fill = "Run")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster, orig.donor) %>% group_by(cluster, orig.donor) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("1", "2", "4", "5", "8",   "0", "3", "7",   "6", "10", "11")) )
  p6 <- ggplot(df, aes(fill=orig.donor, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Donors", fill = "Donor")

  # Change violin plot colors by groups
  pdf("QC/qc_cd4.pdf", 10, 7)
  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
  dev.off()

}

{ cat(redb("### Clonality dot plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

  donor = "BT1" # donor = "BT5"
  donors <- unique(sc_cd3p_cd8@meta.data$orig.donor); donors <- donors[!is.na(donors)]

  for (donor in donors){
    cat("\nDonor: ", donor)

    df_cd8 <- sc_cd3p_cd8@meta.data %>% filter(orig.donor == donor & !is.na(orig.donor) & !is.na(clonotype.tag)) %>% select(clonotype.tag, clon.size.tag) %>%
      distinct() %>% arrange(clon.size.tag)
    df_cd4 <- sc_cd3p_cd4@meta.data %>% filter(orig.donor == donor & !is.na(orig.donor) & !is.na(clonotype.tag)) %>% select(clonotype.tag, clon.size.tag) %>%
      distinct() %>% arrange(clon.size.tag)
    min.val=1; max.val=max(df_cd8$clon.size.tag, df_cd4$clon.size.tag); min_dotsize=1.5; max_dotsize=25

    # -----------> CD8
    if(max(df_cd8$clon.size.tag) != 1){
      pdf(paste0("tcr/sc_cd3p_cd8/", donor, "_clones_dotplot.pdf"))
      net <- make_empty_graph(nrow(df_cd8))
      V(net)$size <- ifelse(df_cd8$clon.size.tag == 1, df_cd8$clon.size.tag, df_cd8$clon.size.tag+1.5) # df_cd8$clon.size.tag # Set node size based on clone size
      # min.val=min(df_cd8$clon.size.tag); max.val=max(df_cd8$clon.size.tag); min_dotsize=1.5; max_dotsize=25
      V(net)$size <- ( (V(net)$size - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # Scale so that it goes from 2 to 10.
      plot(net, vertex.label=NA, vertex.color=c("royalblue4", "gray80")[1+(V(net)$size== 1.5)], vertex.frame.color	= c("white", "white")[1+(V(net)$size== 1.5)]) # c("none", "gray")[1+(V(net)$size== 1.5)] # "gray"
      dev.off()
    }

    # -----------> CD4
    if(max(df_cd4$clon.size.tag) != 1){
      pdf(paste0("tcr/sc_cd3p_cd4/", donor, "_clones_dotplot.pdf"))
      net <- make_empty_graph(nrow(df_cd4))
      V(net)$size <- ifelse(df_cd4$clon.size.tag == 1, df_cd4$clon.size.tag, df_cd4$clon.size.tag+1.5) # df_cd4$clon.size.tag  # Set node size based on clone size
      # min.val=min(df_cd4$clon.size.tag); max.val=max(df_cd4$clon.size.tag); min_dotsize=1.5; max_dotsize=25
      V(net)$size <- ( (V(net)$size - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # Scale so that it goes from 2 to 10.
      plot(net, vertex.label=NA, vertex.color=c("cadetblue4", "gray80")[1+(V(net)$size== 1.5)], vertex.frame.color	= c("white", "white")[1+(V(net)$size== 1.5)])
      dev.off()
    }
  }

}

{ cat(redb("### PDCD1 responsiveness Signature Heatmap ### %%%%%%%%%%%%%%%%%%\n"))
  dir.create("ICB_response/")

  fname <- "/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/signatures_cd3p_human_v2.csv"
  sc_cd3p_lists <- readfile(fname, stringsAsFactors = FALSE)
  signatures <- c("Du_PD1_responders", "Guide_PD1_responders", "Sade_PD1_NON_responders")
  signatures.list <- as.list(sc_cd3p_lists[,signatures]) # , "Sade_PD1_responders"
  signatures.list <- lapply(signatures.list, function(x) x[x!=""])
  lapply(signatures.list, length)
  lapply(signatures.list, head)

  # ---> CD8
  # grep("ALB", rownames(sc_cd3p_cd8@assays$RNA@data), value=T, ignore.case=T)
  set.seed(42)
  sc_cd3p_cd8_mdata <- AddModuleScore(
    object = sc_cd3p_cd8,
    features = signatures.list,
    # ctrl = 150,
    name = names(signatures.list)
  )

  sig_col <- paste0(signatures, 1:3)
  table <- sc_cd3p_cd8_mdata@meta.data[!is.na(sc_cd3p_cd8_mdata@meta.data$orig.donor) & !is.na(sc_cd3p_cd8_mdata@meta.data$clonotype.tag) & sc_cd3p_cd8_mdata@meta.data$expDegree == "Expanded",c("orig.donor", "cell_classification", sig_col)] %>%
    pivot_longer(!c(orig.donor, cell_classification), names_to = "signature", values_to = "value") %>%
    group_by(orig.donor, signature) %>% summarize(value = mean(value)) %>% filter(signature == "Guide_PD1_responders2") %>%
    data.frame()
  rownames(table) <- table$orig.donor

  donor_order <- c("BT8", "BT21", "BT20", "BT44", "BT12", "BT1", "BT9", "BT25", "BT23", "BT38", "BT3", "BT5", "BT7", "BT24", "BT19", "BT46", "BT4", "BT13", "BT35", "BT31", "BT30", "BT32", "BT37", "BT26", "BT27", "BT40", "BT17_brain", "BT15", "BT11", "BT18", "BT22", "BT2")
  # donors2exclude <- c("BT10", "BT29", "BT39", "BT45", "BT47", "BT48")
  # donor_order <- donor_order[which(!donor_order %in% donors2exclude)]
  palette = colorRampPalette(c("#277BC0", "white", "red"))
  matrix <- matrix(table[donor_order,"value"], ncol = 1); rownames(matrix) <- donor_order
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
  write.csv(matrix, file = "ICB_response/CD8_donor_signature.csv", row.names = TRUE)

  # Move middle color position
  paletteLength <- 100
  myBreaks <- c(seq(min(matrix), 0.26, length.out=ceiling(paletteLength/2) + 1),
          seq(0.26, max(matrix), length.out=floor(paletteLength/2))[-1])
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap_scaleChanged.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA",
    breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))
  palette2 = colorRampPalette(c("gray80", "gray80", "#3E7C17"))
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap_scaleChanged_v2.pdf"), cellwidth = 10, cellheight = 10, color = palette2(100), border_color = "NA",
      breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))
  palette3 = colorRampPalette(c("gray80", "gray80", "#146C94"))
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap_scaleChanged_v3.pdf"), cellwidth = 10, cellheight = 10, color = palette3(100), border_color = "NA",
      breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))
  palette4 = colorRampPalette(c("gray80", "gray80", "#764AF1"))
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap_scaleChanged_v4.pdf"), cellwidth = 10, cellheight = 10, color = palette4(100), border_color = "NA",
      breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))
  palette5 = colorRampPalette(c("gray80", "gray80", "#5D3891"))
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD8_donor_signature_heatmap_scaleChanged_v5.pdf"), cellwidth = 10, cellheight = 10, color = palette5(100), border_color = "NA",
      breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))

  # ---> CD4
  # grep("ALB", rownames(sc_cd3p_cd4@assays$RNA@data), value=T, ignore.case=T)
  set.seed(42)
  sc_cd3p_cd4_mdata <- AddModuleScore(
    object = sc_cd3p_cd4,
    features = signatures.list,
    # ctrl = 150,
    name = names(signatures.list)
  )

  # sig_col <- paste0(signatures, 1:3)
  table <- sc_cd3p_cd4_mdata@meta.data[!is.na(sc_cd3p_cd4_mdata@meta.data$orig.donor) & !is.na(sc_cd3p_cd4_mdata@meta.data$clonotype.tag) & sc_cd3p_cd4_mdata@meta.data$expDegree == "Expanded" & sc_cd3p_cd4_mdata@meta.data$cell_classification != "TREG",c("orig.donor", "cell_classification", sig_col)] %>%
    pivot_longer(!c(orig.donor, cell_classification), names_to = "signature", values_to = "value") %>%
    group_by(orig.donor, signature) %>% summarize(value = mean(value)) %>% filter(signature == "Guide_PD1_responders2") %>%
    data.frame()
  rownames(table) <- table$orig.donor

  donor_order <- c("BT5", "BT12", "BT26", "BT8", "BT25", "BT1", "BT21", "BT17_brain", "BT40", "BT46", "BT9", "BT15", "BT23", "BT3", "BT20", "BT24", "BT19", "BT27", "BT11", "BT7", "BT4", "BT32", "BT13", "BT18", "BT2", "BT22")
  # donors2exclude <- c("BT10", "BT29", "BT30", "BT31", "BT35", "BT37", "BT38", "BT39", "BT44", "BT45", "BT47", "BT48")
  # donor_order <- donor_order[which(!donor_order %in% donors2exclude)]
  palette = colorRampPalette(c("#277BC0", "white", "red"))
  matrix <- matrix(table[donor_order,"value"], ncol = 1); rownames(matrix) <- donor_order
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD4_donor_signature_heatmap.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
  write.csv(matrix, file = "ICB_response/CD4_donor_signature.csv", row.names = TRUE)

  # Move middle color position
  paletteLength <- 100
  myBreaks <- c(seq(min(matrix), 0.18, length.out=ceiling(paletteLength/2) + 1),
          seq(0.18, max(matrix), length.out=floor(paletteLength/2))[-1])
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/CD4_donor_signature_heatmap_scaleChanged.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA",
    breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))


  # ---------------------- ABT patients

  abt_abt_combined <- readRDS(file = "/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/integration/seurat/hvg1500_pc30_2studies/figures/data/sc_combined_ABT_PBT_seurat_object.rds")
  DefaultAssay(abt_abt_combined) <- "RNA"
  abt_cd8_cells <- abt_abt_combined@meta.data %>% filter(orig.study == "mathewson_cell_2021" & orig.celltype %in% c("CD8")) %>% pull(cellname)
  abt_cd4_cells <- abt_abt_combined@meta.data %>% filter(orig.study == "mathewson_cell_2021" & orig.celltype %in% c("CD4", "Treg")) %>% pull(cellname)
  abt_cd8 <- subset(x = abt_abt_combined, subset = cellname %in% abt_cd8_cells)
  abt_cd4 <- subset(x = abt_abt_combined, subset = cellname %in% abt_cd4_cells)

  # ---> CD8
  # grep("ALB", rownames(sc_cd3p_cd8@assays$RNA@data), value=T, ignore.case=T)
  set.seed(42)
  abt_cd8_mdata <- AddModuleScore(
    object = abt_cd8,
    features = signatures.list,
    # ctrl = 150,
    name = names(signatures.list)
  )

  sig_col <- paste0(signatures, 1:3)
  table <- abt_cd8_mdata@meta.data[!is.na(abt_cd8_mdata@meta.data$orig.subject) & !is.na(abt_cd8_mdata@meta.data$clonotype.tag) & abt_cd8_mdata@meta.data$exp.Degree == "Expanded",c("orig.subject", "cell_classification", sig_col)] %>%
    pivot_longer(!c(orig.subject, cell_classification), names_to = "signature", values_to = "value") %>%
    group_by(orig.subject, signature) %>% summarize(value = mean(value)) %>% filter(signature == "Guide_PD1_responders2") %>%
    data.frame()
  rownames(table) <- table$orig.subject

  donor_order <- table %>% arrange(desc(value)) %>% pull(orig.subject)
  palette = colorRampPalette(c("#277BC0", "white", "red"))
  matrix <- matrix(table[donor_order,"value"], ncol = 1); rownames(matrix) <- donor_order
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/ABT_CD8_donor_signature_heatmap.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
  write.csv(matrix, file = "ICB_response/ABT_CD8_donor_signature.csv", row.names = TRUE)

  # Move middle color position
  paletteLength <- 100
  myBreaks <- c(seq(min(matrix), 0.25, length.out=ceiling(paletteLength/2) + 1),
          seq(0.25, max(matrix), length.out=floor(paletteLength/2))[-1])
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/ABT_CD8_donor_signature_heatmap_scaleChanged.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA",
    breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))

  # ---> CD4
  # grep("ALB", rownames(sc_cd3p_cd4@assays$RNA@data), value=T, ignore.case=T)
  set.seed(42)
  abt_cd4_mdata <- AddModuleScore(
    object = abt_cd4,
    features = signatures.list,
    # ctrl = 150,
    name = names(signatures.list)
  )

  sig_col <- paste0(signatures, 1:3)
  table <- abt_cd4_mdata@meta.data[!is.na(abt_cd4_mdata@meta.data$orig.subject) & !is.na(abt_cd4_mdata@meta.data$clonotype.tag) & abt_cd4_mdata@meta.data$exp.Degree == "Expanded",c("orig.subject", "cell_classification", sig_col)] %>%
    pivot_longer(!c(orig.subject, cell_classification), names_to = "signature", values_to = "value") %>%
    group_by(orig.subject, signature) %>% summarize(value = mean(value)) %>% filter(signature == "Guide_PD1_responders2") %>%
    data.frame()
  rownames(table) <- table$orig.subject

  donor_order <- table %>% arrange(desc(value)) %>% pull(orig.subject)
  palette = colorRampPalette(c("#277BC0", "white", "red"))
  matrix <- matrix(table[donor_order,"value"], ncol = 1); rownames(matrix) <- donor_order
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/ABT_CD4_donor_signature_heatmap.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA")
  write.csv(matrix, file = "ICB_response/ABT_CD4_donor_signature.csv", row.names = TRUE)

  # Move middle color position
  paletteLength <- 100
  myBreaks <- c(seq(min(matrix), 0.25, length.out=ceiling(paletteLength/2) + 1),
          seq(0.25, max(matrix), length.out=floor(paletteLength/2))[-1])
  pheatmap(matrix, fontsize_col = 7, main = "PDCD1 responsiveness", cluster_rows=F, cluster_cols=F, filename = paste0("ICB_response/ABT_CD4_donor_signature_heatmap_scaleChanged.pdf"), cellwidth = 10, cellheight = 10, color = palette(100), border_color = "NA",
    breaks = myBreaks, legend_breaks = c(0,0.2,0.2,0.3,0.4))

}

{ cat(redb("### PBT vs ABT cell type proportions ### %%%%%%%%%%%%%%%%%%%%%%%%\n"))

  #  ---> Functions
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  data_summary_collapsed <- function(x) {
    m <- mean(x)
    ymin <- m
    ymax <- m
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  # Read ABT separate clustering
  abt_ss2_cd4_mdata <- readRDS("~/tmp_large/tmp/abt/results/clustering/mathewson_Smartseq_CD4/.object_meta.data_seurat_mean0.01_pct25_pc25.rds")
  abt_ss2_cd4 <- readRDS("~/tmp_large/tmp/abt/results/clustering/mathewson_Smartseq_CD4/.object_stem_seurat_mean0.01_pct25.rds")

  abt_ss2_cd4@meta.data <- abt_ss2_cd4_mdata
  abt_ss2_cd4@meta.data$cluster <- abt_ss2_cd4@meta.data$RNA_snn_res.0.2

  # ---> CTL proportions

  # ABT
  donor_cells <- abt_ss2_cd4@meta.data %>% filter(!is.na(orig.subject)) %>% group_by(orig.subject) %>% summarize(cells = n())
  donor_CTLcells <- table(abt_ss2_cd4@meta.data$orig.subject, abt_ss2_cd4@meta.data$cluster) %>% melt() %>% as_tibble() %>% rename(orig.subject = Var1, cluster = Var2, CTL_cells = value) %>%
    filter(cluster == "1")

  table <- merge(donor_cells, donor_CTLcells, by = "orig.subject")
  table <- table %>% mutate(CTL_pct = 100*CTL_cells/cells)
  table$orig.study <- "mathewson_cell_2021"

  # PBT
  donor_cells.pbt <- sc_cd3p_cd4@meta.data %>% filter(!is.na(orig.donor)) %>% group_by(orig.donor) %>% summarize(cells = n()) %>% rename(orig.subject = orig.donor)
  donor_CTLcells.pbt <- table(sc_cd3p_cd4@meta.data$orig.donor, sc_cd3p_cd4@meta.data$cell_classification) %>% melt() %>% rename(orig.subject = Var1, cluster = Var2, CTL_cells = value) %>% filter(cluster == "CTL")

  table.pbt <- merge(donor_cells.pbt, donor_CTLcells.pbt, by = "orig.subject")
  table.pbt <- table.pbt %>% mutate(CTL_pct = 100*CTL_cells/cells)
  table.pbt$orig.study <- "upadhye_NNNN_202N"

  table_full_final <- rbind(table, table.pbt)

  # Dotplot
  pdf(paste0("ICB_response/CTL_proportions_study_full_cohort.pdf"))
  p <- ggplot(table_full_final, aes(x=orig.study, y=CTL_pct, fill=orig.study, color=orig.study)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + theme_classic() +
    scale_fill_manual(values=c("white", "white")) +
    scale_color_manual(values=c("gray50", "indianred3")) + # "cadetblue4", "royalblue4"
    # ylim(0,80) +
    stat_summary(
      geom='crossbar',
      fun.data=data_summary_collapsed,
      width=0.3,
      color=c("gray50", "indianred3")
    ) +
    stat_summary(
      geom='errorbar',
      fun.data=data_summary,
      width=0.2,
      color=c("gray50", "indianred3"),
      size=1.2
    ) +
    stat_compare_means(paired = F) +
    labs(title = paste0("CTL percentage"), x = "Study", y = "Percentage of CTLs among CD4+ cells (%)")
  print(p)
  dev.off()
  write.csv(table_full_final[,-3], file = "ICB_response/CTL_proportions_study_full_cohort.csv", row.names = FALSE)


  # ------------------------------------ CTL proportions (CD4 universe)

  # ABT (ss2)
  abt_ss2_cd4_mdata <- readRDS("~/tmp_large/tmp/abt/results/clustering/mathewson_Smartseq_CD4/.object_meta.data_seurat_mean0.01_pct25_pc25.rds")
  abt_ss2_cd4 <- readRDS("~/tmp_large/tmp/abt/results/clustering/mathewson_Smartseq_CD4/.object_stem_seurat_mean0.01_pct25.rds")

  abt_ss2_cd4@meta.data <- abt_ss2_cd4_mdata
  abt_ss2_cd4@meta.data$cluster <- abt_ss2_cd4@meta.data$RNA_snn_res.0.2

  donor_cells <- abt_ss2_cd4@meta.data %>% filter(!is.na(orig.subject)) %>% group_by(orig.subject) %>% summarize(cells = n())
  donor_TREGcells <- table(abt_ss2_cd4@meta.data$orig.subject, abt_ss2_cd4@meta.data$orig.celltype) %>% melt() %>% as_tibble() %>% rename(orig.subject = Var1, cluster = Var2, TREG_cells = value) %>%
    filter(cluster == "Treg")

  table <- merge(donor_cells, donor_TREGcells, by = "orig.subject")
  table <- table %>% mutate(TREG_pct = 100*TREG_cells/cells)
  table$orig.study <- "mathewson_cell_2021"

  # ABT (10x)
  mdata_10X <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/info/ABT/Mathewson_GSE163108_metadata_10x.csv", row.names=1) %>% filter(!is.na(annotate_Tcelltype) & annotate_Tcelltype %in% c("CD4+", "Treg")) %>%
    rename(orig.subject = sampleid)
  mdata_10X$cellname <- rownames(mdata_10X)

  donor_cells.abt_10x <- mdata_10X %>% filter(!is.na(orig.subject)) %>% group_by(orig.subject) %>% summarize(cells = n())
  donor_TREGcells <- table(mdata_10X$orig.subject, mdata_10X$annotate_Tcelltype) %>% melt() %>% as_tibble() %>% rename(orig.subject = Var1, cluster = Var2, TREG_cells = value) %>%
    filter(cluster == "Treg")

  table.abt_10x <- merge(donor_cells.abt_10x, donor_TREGcells, by = "orig.subject")
  table.abt_10x <- table.abt_10x %>% mutate(TREG_pct = 100*TREG_cells/cells)
  table.abt_10x$orig.study <- "mathewson_cell_2021"

  # PBT
  donor_cells.pbt <- sc_cd3p_cd4@meta.data %>% filter(!is.na(orig.donor)) %>% group_by(orig.donor) %>% summarize(cells = n()) %>% rename(orig.subject = orig.donor)
  donor_TREGcells.pbt <- table(sc_cd3p_cd4@meta.data$orig.donor, sc_cd3p_cd4@meta.data$cell_classification) %>% melt() %>% rename(orig.subject = Var1, cluster = Var2, TREG_cells = value) %>% filter(cluster == "TREG")

  table.pbt <- merge(donor_cells.pbt, donor_TREGcells.pbt, by = "orig.subject")
  table.pbt <- table.pbt %>% mutate(TREG_pct = 100*TREG_cells/cells)
  table.pbt$orig.study <- "upadhye_NNNN_202N"

  table_full_final <- rbind(table, table.abt_10x, table.pbt)

  # Dotplot
  pdf(paste0("ICB_response/TREG_proportions_study_full_cohort_including_ABT_10x.pdf"))
  p <- ggplot(table_full_final, aes(x=orig.study, y=TREG_pct, fill=orig.study, color=orig.study)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + theme_classic() +
    scale_fill_manual(values=c("white", "white")) +
    scale_color_manual(values=c("gray50", "indianred3")) + # "cadetblue4", "royalblue4"
    # ylim(0,80) +
    stat_summary(
      geom='crossbar',
      fun.data=data_summary_collapsed,
      width=0.3,
      color=c("gray50", "indianred3")
    ) +
    stat_summary(
      geom='errorbar',
      fun.data=data_summary,
      width=0.2,
      color=c("gray50", "indianred3"),
      size=1.2
    ) +
    stat_compare_means(paired = F) +
    labs(title = paste0("TREG percentage"), x = "Study", y = "Percentage of TREGs among CD4+ cells (%)")
  print(p)
  dev.off()
  write.csv(table_full_final[,-3], file = "ICB_response/TREG_proportions_study_full_cohort_including_ABT_10x.csv", row.names = FALSE)

  table_full_final <- rbind(table, table.pbt)

  # Dotplot
  pdf(paste0("ICB_response/TREG_proportions_study_full_cohort.pdf"))
  p <- ggplot(table_full_final, aes(x=orig.study, y=TREG_pct, fill=orig.study, color=orig.study)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + theme_classic() +
    scale_fill_manual(values=c("white", "white")) +
    scale_color_manual(values=c("gray50", "indianred3")) + # "cadetblue4", "royalblue4"
    # ylim(0,80) +
    stat_summary(
      geom='crossbar',
      fun.data=data_summary_collapsed,
      width=0.3,
      color=c("gray50", "indianred3")
    ) +
    stat_summary(
      geom='errorbar',
      fun.data=data_summary,
      width=0.2,
      color=c("gray50", "indianred3"),
      size=1.2
    ) +
    stat_compare_means(paired = F) +
    labs(title = paste0("TREG percentage"), x = "Study", y = "Percentage of TREGs among CD4+ cells (%)")
  print(p)
  dev.off()
  write.csv(table_full_final[,-3], file = "ICB_response/TREG_proportions_study_full_cohort.csv", row.names = FALSE)

}

{ cat(redb("### PBT and BCC  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  yost_bcc <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/integration/seurat/yost_GSE123813/yost_bcc_seurat_object.rds")
  yost_scc <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/integration/seurat/yost_GSE123813/yost_scc_seurat_object.rds")

  table(yost_bcc$orig.response)
  table(yost_bcc$expDegree)
  table(yost_bcc@meta.data[yost_bcc@meta.data$orig.response == "Responder","expDegree"])
  table(yost_bcc@meta.data[yost_bcc@meta.data$orig.response == "Non-Responder","expDegree"])

  donor_response_bcc_df <- yost_bcc@meta.data %>% select(orig.Patient_ID, orig.response) %>% distinct()
  expPct_perDonor_bcc <- yost_bcc@meta.data %>% filter(!is.na(clonotype.tag)) %>% group_by(orig.Patient_ID) %>% summarize(exp=sum(expDegree == "Expanded"), non_exp=sum(expDegree == "Non_expanded"),expanded_pct = 100*sum(expDegree == "Expanded")/(sum(expDegree == "Expanded")+sum(expDegree == "Non_expanded")) ) %>%
    arrange(desc(expanded_pct))
  expPct_perDonor_bcc <- merge(donor_response_bcc_df, expPct_perDonor_bcc, by = "orig.Patient_ID")
  p <- expPct_perDonor_bcc %>% ggplot(aes(x=orig.response, y=expanded_pct)) +
    geom_dotplot(binaxis='y', stackdir='center')
  pdf("expPct_perDonor_bcc.pdf")
  print(p)
  dev.off()

  donor_response_scc_df <- yost_scc@meta.data %>% select(orig.Patient_ID, orig.response) %>% distinct()
  expPct_perDonor_scc <- yost_scc@meta.data %>% filter(!is.na(clonotype.tag)) %>% group_by(orig.Patient_ID) %>% summarize(exp=sum(expDegree == "Expanded"), non_exp=sum(expDegree == "Non_expanded"),expanded_pct = 100*sum(expDegree == "Expanded")/(sum(expDegree == "Expanded")+sum(expDegree == "Non_expanded")) ) %>%
    arrange(desc(expanded_pct))
  expPct_perDonor_scc <- merge(donor_response_scc_df, expPct_perDonor_scc, by = "orig.Patient_ID")
  p <- expPct_perDonor_scc %>% ggplot(aes(x=orig.response, y=expanded_pct)) +
    geom_dotplot(binaxis='y', stackdir='center')
  pdf("expPct_perDonor_scc.pdf")
  print(p)
  dev.off()

  figures_path <- here::here("results/integration/seurat/yost_GSE123813/figures")
  if(!dir.exists(figures_path)) dir.create(figures_path)
  setwd(figures_path); cat("Working at:", getwd(), "\n")
  system("ls -loh")

  ### Yost create Seurat Objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  mdata_bcc <- get(load("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all-Batch2/raw/sc_integration/yost_GSE123813_cd8_bcc_preTreatment_10x_annotation.rdata"))
  edata_bcc <- get(load("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all-Batch2/raw/sc_integration/yost_GSE123813_cd8_bcc_preTreatment_10x_tpm.rdata"))
  yost_bcc <- CreateSeuratObject(counts=edata_bcc, meta.data=mdata_bcc)

  yost_bcc <- NormalizeData(
      object = yost_bcc,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = TRUE
    )
  yost_bcc <- FindVariableFeatures(
      object = yost_bcc,
      selection.method = "vst",
      nfeatures = 2000,
      mean.cutoff = c(0.01, 8),
      dispersion.cutoff = c(1, Inf),
      verbose = TRUE
    )
  yost_bcc <- ScaleData(
      object = yost_bcc,
      vars.to.regress = c("nCount_RNA", "percent.mt"),
      block.size = 2000,
      verbose = TRUE
    )
  yost_bcc <- RunPCA(
    object = yost_bcc,
    features = VariableFeatures(yost_bcc),
    npcs = 20,
    nfeatures.print = 15
  )

  # Define clonotypes and Get clone_size
  tmp <- yost_bcc@meta.data %>% filter(!is.na(TRA_cdr3s_aa) | !is.na(TRB_cdr3s_aa)) %>% group_by(TRA_cdr3s_aa, TRB_cdr3s_aa) %>% summarize(clone.size.tag=n()) %>% drop_na() %>% arrange(desc(clone.size.tag))
  tmp <- cbind(clonotype.tag = paste0("clonotype",1:nrow(tmp)),tmp)
  yost_bcc@meta.data <- merge(yost_bcc@meta.data, tmp, by = c("TRA_cdr3s_aa","TRB_cdr3s_aa"), all.x=T)
  rownames(yost_bcc@meta.data) <- yost_bcc@meta.data$cellname

  # Define the expanded cells
  yost_bcc@meta.data$expDegree <- yost_bcc@meta.data$clone.size.tag
  yost_bcc@meta.data$expDegree[yost_bcc@meta.data$expDegree > 1] <- "Expanded"
  yost_bcc@meta.data$expDegree[yost_bcc@meta.data$expDegree == 1] <- "Non_expanded"

  saveRDS(yost_bcc, "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/integration/seurat/yost_GSE123813/yost_bcc_seurat_object.rds")
  saveRDS(yost_bcc@meta.data, "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/integration/seurat/yost_GSE123813/yost_bcc_mdata.rds")


  ### Crater plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater.R") # source("/home/ciro/scripts/handy_functions/devel/plots_crater.R")
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("craters", showWarnings = FALSE)

  # - PBT CD8
  f1 <- read.csv(here::here("results", "dgea/CD8/CD8_Expansion/ExpandedvsNon_expanded/mastlog2cpm_results.csv"))
  f1 <- f1 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f1, here::here("results", "dgea/CD8/CD8_Expansion/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"), row.names = F)

  # - BCC CD8 responders (expanded vs non-expanded)
  f2 <- read.csv(here::here("results", "dgea/integration/bcc_CD8_PreTreatment/CD8_Expansion_Responders/ExpandedvsNon_expanded/mastlog2cpm_results.csv"))
  f2 <- f2 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f2, here::here("results", "dgea/integration/bcc_CD8_PreTreatment/CD8_Expansion_Responders/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"), row.names = F)

  fconfigs = list(
    # BCC Responders
    list(
      fnames = c(
        "PBT_Exp_vs_NonExp" = here::here("results", "dgea/CD8/CD8_Expansion/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"),
        "BCC_Exp_vs_NonExp_Responders" = here::here("results", "dgea/integration/bcc_CD8_PreTreatment/CD8_Expansion_Responders/ExpandedvsNon_expanded/.mastlog2cpm_results.csv")
      ), result_id = "craters/BCC_PBT_CD8_Responders_", plot_squared = TRUE,
        highlight_genes = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CCL3", "CCL4", "CCR1", "CCR3", "CCR5", "TNF", "IFNG", "ITGAE", "ZNF683", "CXCR6", "PDCD1", "LAG3", "HOPX", "FASLG",
        "GNLY", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "PSMB9", "PSME1", "PSME2", "CALR", "IL7R", "TCF7", "SELL", "LEF1", "CD55", "CD27", "KLF2", "MAL", "LTB", "CST7", "CD7"),
      selectss = list(c('orig.study', 'Yost_2019'), c('expDegree', 'Expanded'), c('orig.response', 'Responder'))
    )
  )
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))

  # -----------> Yost BCC
  fconfig = fconfigs[[1]]
  fc = "0.35"
  cat(" * FC:", fc, "\n");
  void <- crater_plot(
    tests_list = fconfig$fnames,
    edataf = edata,
    annotf = mdata,
    sample_filter = fconfig$selectss,
    gene_filter = degfilt,
    feature_subset = mygenes,
    topgenes = c("top10", fconfig$highlight_genes),
    lfcthresh = fc,
    # column4stats = fconfig$columns,
    outputname = fconfig$result_id,
    plot_interactive = TRUE,
    plot_squared = fconfig$plot_squared,
    limits_col = fconfig$limits_col,
    limits_size = fconfig$limits_size,
    verbose = TRUE
  )

  fconfig = fconfigs[[3]]
  fc = "0.35"
  cat(" * FC:", fc, "\n");
  void <- crater_plot(
    tests_list = fconfig$fnames,
    edataf = edata,
    annotf = mdata,
    sample_filter = fconfig$selectss,
    gene_filter = degfilt,
    feature_subset = mygenes,
    topgenes = c("top10", fconfig$highlight_genes),
    lfcthresh = fc,
    # column4stats = fconfig$columns,
    outputname = fconfig$result_id,
    plot_interactive = TRUE,
    plot_squared = fconfig$plot_squared,
    limits_col = fconfig$limits_col,
    limits_size = fconfig$limits_size,
    verbose = TRUE
  )

  # ============== Final
  source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater_v2.R")
  fconfigs <- lapply(fconfigs, function(x){
    x$result_id <- paste0(x$result_id, "final_")
    return(x)
  })

  # -----------> Yost BCC
  edata =  expm1(yost_bcc@assays$RNA@data)
  mygenes = grep("XIST|RPS4Y1|^RP11|^RP", rownames(edata), value = TRUE, invert = TRUE) # Filter unwanted genes
  mdata = yost_bcc@meta.data

  fconfig = fconfigs[[1]]
  fc = "0.35"
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))
  plot <- crater_plot(
    tests_list = fconfig$fnames,
    edataf = edata,
    annotf = mdata,
    sample_filter = fconfig$selectss,
    gene_filter = degfilt,
    feature_subset = mygenes,
    topgenes = c("top10", fconfig$highlight_genes),
    lfcthresh = fc,
    # column4stats = fconfig$columns,
    outputname = gsub("v4", "v5", fconfig$result_id),
    plot_interactive = TRUE,
    plot_squared = fconfig$plot_squared,
    limits_col = fconfig$limits_col,
    limits_size = fconfig$limits_size,
    verbose = TRUE,
    theme_extra = theme(axis.text.x = element_blank(),axis.text.y = element_blank()) #+ labs(x = NULL, y = NULL)
  )

}
