#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2022-01-15
# ---


### ==== CD4/CD8 clustering ==== ####

# source("/home/kmlanderos/pbtumor-all/scripts/object_lock_CD4_CD8.R")

## Generate all the Seurat objects with all the combinations for CD4 and for CD8
# Saved in /mnt/beegfs/kmlanderos/; e.g. /mnt/beegfs/kmlanderos/CD4_object_lock_mean0.01_pct15_pc15.rds
redu = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2"))
setwd('/home/kmlanderos/ad_hoc/pbtumor-all/results/figures/')

packages_funcs = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/filters.R",
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  "/home/ciro/scripts/handy_functions/devel/plots.R",
  "/home/ciro/scripts/figease/source.R", # fig_global_objects
  "ggplot2", "cowplot", "patchwork", "Seurat", "stringr"
)
loaded <- lapply(X = packages_funcs, FUN = function(x){
  cat("*", x, "\n")
  if(!file.exists(x)){
    suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
  }else{ source(x) }
}); theme_set(theme_cowplot())

# Read files
objects <- list()
for (i in list.files("/mnt/beegfs/kmlanderos/", pattern = "CD\\d_object_lock") ) { # , full.names = T
  filename <- gsub(".rds","",gsub("_object_lock_mean0.01","",i))
  cat("\n","reading file:", filename,"\n")
  objects[[filename]] <- readRDS(paste0("/mnt/beegfs/kmlanderos/", i))
  objects[[filename]]@meta.data$UMAP_1 <- NULL
  objects[[filename]]@meta.data$UMAP_2 <- NULL
  objects[[filename]]@meta.data$cluster <- NULL
  objects[[filename]]@meta.data$cluster <- objects[[filename]]@meta.data$RNA_snn_res.0.8
  objects[[filename]]@meta.data <- joindf(objects[[filename]]@meta.data,
    as.data.frame(objects[[filename]]@reductions$umap@cell.embeddings))
}

# This part was done before setting a locked object
{ cat(redb("### Markers: dim. red. & Dot/curtain (Scatter) - CD4/CD8 clustering ### %%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dim_reduction_markers");

  fconfigs <- list()
  for (i in names(objects)) {
    data <- list()
    if (grepl("CD4",i)) {
      data <- list(result_id = paste0("dim_reduction_features/sc_cd3p_CD4/",i,"/"),
        edata = paste0("objects[['",i,"']]@assays$RNA@data"), metadata = paste0("objects[['",i,"']]@meta.data"),
        axis_x = redu[[1]][1], axis_y = redu[[1]][2],
        features = c("CD3D", "CD4", " CD40LG", "GZMA", "GZMB", "GZMK", "CCL4", "CCL5", "CXCL13", "IL21", "CXCR3", "TBX21",
        "IFNG", "PRF1", "TNF", "HOPX", "TNFAIP3", "NFKB1", "MX1", "IRF7", "RORA", "RORC", "CCR6", "KLRB1", "FOXP3", "IL2RA",
        "ENTPD1", "TIGIT", "PDCD1", "HLA-DRB1", "ITGAE", "TOX", "LAG3", "ICOS", "BATF", "PRDM1", "XCL1", "XCL2")
      )
    }else if(grepl("CD8",i)){
      data <- list(result_id = paste0("dim_reduction_features/sc_cd3p_CD8/",i,"/"),
        edata = paste0("objects[['",i,"']]@assays$RNA@data"), metadata = paste0("objects[['",i,"']]@meta.data"),
        axis_x = redu[[1]][1], axis_y = redu[[1]][2],
        features = c("CD3D", "CD8A", "CD8B", "GZMK", "GZMA", "IRF7", "ZNF683", "ITGAE", "KLRC2", "KLRC3", "FGFBP2", "FCGR3A",
        "KLRB1", "SLC4A10", "IL18RAP", "RORC", "GZMB", "PRF1", "PDCD1", "TNF", "IFNG", "HLA-DRB1", "LAG3", "TOX")
      )
    }
    fconfigs <- append(fconfigs, list(data))
  }

  pp_markers = fig_plot_base( ## -----------------------------------------------
    fconfigs, return_plot = FALSE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(size = 0.8) + labs(x = "Dim 1", y = "Dim 2") +
      scale_color_gradientn(colours = c("lightgrey", "blue"), breaks = scales::pretty_breaks())
    }
  )

}

# This part was done before setting a locked object
{ cat(redb("### Signature Analyses - CD4/CD8 clustering ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # cp ~/Documents/liai/cancer/pbtumor/signatures_cd3p*csv /Volumes/ciro/pbtumor/info/
  # Rscript /home/ciro/scripts/functions/csvCorrect.R /home/ciro/pbtumor/info/signatures_cd3p.csv

  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R") # stats_summary_table
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
  source("/home/ciro/scripts/handy_functions/devel/utilities.R")  # vlist2df
  # source("/home/ciro/scripts/handy_functions/devel/code.R") #vlist2df
  source("/home/ciro/scripts/figease/figease.R")
  source("/home/ciro/scripts/handy_functions/devel/filters.R")

  #setwd(setwd("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures")
  #dir.create("gsea")

  ## ---- GSEA ---- ###

  dir.create("gsea/sc_cd3p_CD8subclustering/")
  dir.create("gsea/sc_cd3p_CD4subclustering/")
  ## Table with the genes for the signatures
  fname <- "gsea/signatures_cd3p_human.csv"
  slists <- readfile(fname, stringsAsFactors = FALSE)

  # i = "CD4_pct15_pc15"
  for (i in names(objects)){
    cat("\n","Processing sample:", i,"\n")
    fconfigs = list(
      list( result_id = ifelse(grepl("CD4",i),paste0("gsea/sc_cd3p_CD4subclustering/", i, "/"),paste0("gsea/sc_cd3p_CD8subclustering/", i, "/")),
        edata = paste0("expm1(objects[['",i,"']]@assays$RNA@data)"),
        metadata = paste0("objects[['",i,"']]@meta.data"), # levels(objects[['CD4_pct15_pc15']]@meta.data$cluster)
        lists = lapply(X = slists[,c("TH17_Seumois2020", "CD4_cytotoxicity_Patil2018", "TFH_Locci2013", "th1_signature1_arlehamn", "CD4CTLvsTCM_Patil2018", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "tcell_activation_GO.0042110", "tcr_signaling_RHSA202403")], FUN = function(x) x[-c(1:2)] ),
        comparisons = c("cluster")
        #comparisons = c("cluster", unique(objects[[i]]@meta.data$cluster), "REST")
        #comparisons = c("cluster", "0", "3", "4", "5.8", "7", "9.8", "10.8", "11.8", "13", "14", "15", "16", "REST")
      )
    )
    pp_gsea = fig_gsea(fconfigs, features = rownames(objects[[i]]), return_plot = TRUE)
  }

}

### ================== Figures ================== ###

# Done after locking the clustering objectss (CD4 and CD8)
# CD8 = pct20_pc20_res0.2
# CD4 = pct25_pc15_0.4

# source("/home/kmlanderos/pbtumor-all/scripts/object_lock_CD4_CD8_subclustering.R")
source("/home/kmlanderos/pbtumor-all/scripts/figures/global_CD4_CD8.R")

setwd('/home/kmlanderos/ad_hoc/pbtumor-all/results/figures/CD4_CD8_subclustering')
system("ls -loh")

{ cat(redb("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # As you progress they appear in many tasks and so become more or less global

  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
    "/home/ciro/scripts/figease/source.R", # fig_global_objects
    "ggplot2", "cowplot", "patchwork", "Seurat", "stringr", "tidyverse", "data.table", "pheatmap"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())

  # Merge the current metadata with the TCRmetadata (obtained through Vicente's pipeline)
  # NOTE: Not necsesary bc we obtained this info from the previous object (CD4+CD8). Modify in case we want to use it/
  # sc_cd3p_tcr = data.frame(data.table::rbindlist(lapply(
  #   X = c("/home/kmlanderos/large/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd4/metadata.rds"),
  #   FUN = function(x){
  #     y <- readRDS(x); y$cellname = rownames(y); y
  # })), stringsAsFactors = FALSE)
  # rownames(sc_cd3p_tcr) <- sc_cd3p_tcr$cellname; sc_cd3p_tcr$cellname <- NULL
  # sc_cd3p_tcr <- joindf(sc_cd3p_tcr, sc_cd3p@meta.data)
  #
  # sc_cd3p_tcr = data.frame(data.table::rbindlist(lapply(
  #   X = c("/home/kmlanderos/large/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd8/metadata.rds"),
  #   FUN = function(x){
  #     y <- readRDS(x); y$cellname = rownames(y); y
  # })), stringsAsFactors = FALSE)
  # rownames(sc_cd3p_tcr) <- sc_cd3p_tcr$cellname; sc_cd3p_tcr$cellname <- NULL
  # sc_cd3p_tcr <- joindf(sc_cd3p_tcr, sc_cd3p@meta.data)

  # Define LowGrade and HighGrade patients
  sc_cd3p_cd8@meta.data$grade <- sc_cd3p_cd8@meta.data$orig.donor
  sc_cd3p_cd8@meta.data$grade[sc_cd3p_cd8@meta.data$grade %in% c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5", "BT24", "BT25")] <- "LG"
  sc_cd3p_cd8@meta.data$grade[sc_cd3p_cd8@meta.data$grade %in% c("BT10", "BT15", "BT22", "BT11", "BT12", "BT18", "BT21", "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23")] <- "HG"

  sc_cd3p_cd4@meta.data$grade <- sc_cd3p_cd4@meta.data$orig.donor
  sc_cd3p_cd4@meta.data$grade[sc_cd3p_cd4@meta.data$grade %in% c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5", "BT24", "BT25")] <- "LG"
  sc_cd3p_cd4@meta.data$grade[sc_cd3p_cd4@meta.data$grade %in% c("BT10", "BT15", "BT22", "BT11", "BT12", "BT18", "BT21", "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23")] <- "HG"

  # Define the Low- and Highly expanded cells for LowGrade and HighGrade patients
  sc_cd3p_cd8@meta.data$expDegree <- sc_cd3p_cd8@meta.data$clon.size.tag
  sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree > 4 & sc_cd3p_cd8@meta.data$grade == "LG"] <- "HighExp_LG"
  sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree > 4 & sc_cd3p_cd8@meta.data$grade == "HG"] <- "HighExp_HG"
  sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree <= 4 & sc_cd3p_cd8@meta.data$grade == "LG"] <- "LowExp_LG"
  sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree <= 4 & sc_cd3p_cd8@meta.data$grade == "HG"] <- "LowExp_HG"

  sc_cd3p_cd4@meta.data$expDegree <- sc_cd3p_cd4@meta.data$clon.size.tag
  sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree > 4 & sc_cd3p_cd4@meta.data$grade == "LG"] <- "HighExp_LG"
  sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree > 4 & sc_cd3p_cd4@meta.data$grade == "HG"] <- "HighExp_HG"
  sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree <= 4 & sc_cd3p_cd4@meta.data$grade == "LG"] <- "LowExp_LG"
  sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree <= 4 & sc_cd3p_cd4@meta.data$grade == "HG"] <- "LowExp_HG"

  # Errase previous UMAP info (from past clusterng)
  sc_cd3p_cd8@meta.data$UMAP_1 <- NULL; sc_cd3p_cd8@meta.data$UMAP_2 <- NULL
  sc_cd3p_cd4@meta.data$UMAP_1 <- NULL; sc_cd3p_cd4@meta.data$UMAP_2 <- NULL
  # Add the UMAP cell embeddings to the metadata
  sc_cd3p_cd8@meta.data <- joindf(sc_cd3p_cd8@meta.data,
    as.data.frame(sc_cd3p_cd8@reductions$umap@cell.embeddings))
  sc_cd3p_cd4@meta.data <- joindf(sc_cd3p_cd4@meta.data,
    as.data.frame(sc_cd3p_cd4@reductions$umap@cell.embeddings))

  # Group celltypes
  sc_cd3p_cd8@meta.data$cell_classification <- as.character(sc_cd3p_cd8@meta.data$RNA_snn_res.0.2)
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("0","1","7")] <- "GZMK"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("6","8")] <- "TRM"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("2")] <- "INNATE-LIKE_NKR"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("3")] <- "INNATE-LIKE_XCL1-2"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("4")] <- "TCM"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("5")] <- "MAIT"
  sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("9")] <- "Cell_Cycle"
  # sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("0")] <- "GZMK1"
  # sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("1")] <- "GZMK2"
  # sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("7")] <- "GZMK3_BT25"
  # sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("6")] <- "TRM1"
  # sc_cd3p_cd8@meta.data$cell_classification[sc_cd3p_cd8@meta.data$cell_classification %in% c("8")] <- "TRM2"

  sc_cd3p_cd4@meta.data$cell_classification <- as.character(sc_cd3p_cd4@meta.data$RNA_snn_res.0.4)
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("3", "6", "7", "8")] <- "TCM" #
  # sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("0", "2", "4", "9")] <- "CTL" #
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("1")] <- "TNF" #
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("5")] <- "TREG" #
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("10")] <- "TFH"
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("11")] <- "Cell_Cycle"
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("0")] <- "CTL1" #
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("2")] <- "CTL_TH1" #
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("4")] <- "CTL_TH17"
  sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("9")] <- "CTL_IFN"
  # sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("3")] <- "TCM1" #
  # sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("6")] <- "TCM2" #
  # sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("7")] <- "TCM3" #
  # sc_cd3p_cd4@meta.data$cell_classification[sc_cd3p_cd4@meta.data$cell_classification %in% c("8")] <- "TCM4" #

  # # Create cluster.celltype columns
  # sc_cd3p_cd8@meta.data$cluster.celltype <- NULL; sc_cd3p_cd8@meta.data$cluster.celltype <- factor(paste(sc_cd3p_cd8@meta.data$cluster,": ",sc_cd3p_cd8@meta.data$cell_classification))
  # sc_cd3p_cd4@meta.data$cluster.celltype <- NULL; sc_cd3p_cd4@meta.data$cluster.celltype <- factor(paste(sc_cd3p_cd4@meta.data$cluster,": ",sc_cd3p_cd4@meta.data$cell_classification))

  # Update the Diognosis from Neuroblastoma to High_grade_glioma
  sc_cd3p_cd8@meta.data$orig.Diagnosis[sc_cd3p_cd8@meta.data$orig.Diagnosis == "Neuroblastoma"] <- "High grade glioma"
  sc_cd3p_cd4@meta.data$orig.Diagnosis[sc_cd3p_cd4@meta.data$orig.Diagnosis == "Neuroblastoma"] <- "High grade glioma"
  sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass[sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass == "Neuroblastoma"] <- "High grade glioma"
  sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass[sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass == "Neuroblastoma"] <- "High grade glioma"

  # Let us recalculate the clone size of the clonotypes shared between CD4 and CD4, and only consider cells with Gex
  sc_cd3p_cd8@meta.data$clon.size.tag.old <- sc_cd3p_cd8@meta.data$clon.size.tag # Save old values

  # Get the problematic clonotypes
  cd8_table <- table(sc_cd3p_cd8@meta.data$clonotype.tag) %>% reshape2::melt() %>% arrange(desc(value))
  cd4_table <- table(sc_cd3p_cd4@meta.data$clonotype.tag) %>% reshape2::melt() %>% arrange(desc(value))
  merge_table <- merge(cd8_table, cd4_table, by = "Var1");
  problematic_clonotypes <- as.character(merge_table$Var1)

  # CD8
  clonotype.cs.cd8 <- as.data.frame(table(sc_cd3p_cd8@meta.data %>% pull(clonotype.tag)) ) %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd8) <- c("clonotype", "clone_size")
  tmp <- unlist(apply(sc_cd3p_cd8@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"]) & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd8[clonotype.cs.cd8$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  sc_cd3p_cd8@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]
  # sc_cd3p_cd8@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))

  # CD4
  sc_cd3p_cd4@meta.data$clon.size.tag.old <- sc_cd3p_cd4@meta.data$clon.size.tag # Save old values

  clonotype.cs.cd4 <- as.data.frame(table(sc_cd3p_cd4@meta.data %>% pull(clonotype.tag)) )  %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd4) <- c("clonotype", "clone_size")
  tmp <- unlist(apply(sc_cd3p_cd4@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"]) & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd4[clonotype.cs.cd4$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  sc_cd3p_cd4@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]
  # sc_cd3p_cd4@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))

  # --- Add specificity tag

  # CD4
  GLIPH2Input_CD4 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-02-07/GLIPH2Input_TCR-Data_All-TRB_CD4_wAlphaChain.tsv", sep = "\t")
  GLIPH2Input_CD4 <- GLIPH2Input_CD4 %>% mutate(id = paste(GLIPH2Input_CD4$cdr3b.aa.se, GLIPH2Input_CD4$trb.v, GLIPH2Input_CD4$trb.j, sep = "_"))

  GLIPH2Output_CD4 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/2022-02-07/GLIPH2Output_TCR-Data_CD4.csv")
  GLIPH2Output_CD4 <- GLIPH2Output_CD4 %>% mutate(id = paste(GLIPH2Output_CD4$TcRb, GLIPH2Output_CD4$V, GLIPH2Output_CD4$J, sep = "_"))
  GLIPH2Output_CD4 <- GLIPH2Output_CD4[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  merge <- merge(data.frame(GLIPH2Output_CD4), GLIPH2Input_CD4[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups <- unlist(str_split(merge$clonotype.tag, ";"))
  sc_cd3p_cd4@meta.data <- sc_cd3p_cd4@meta.data %>% mutate(sGroup_tag = case_when(
    clonotype.tag %in% cts_sGroups ~ TRUE,
    TRUE ~ FALSE
  ))
  # table(sc_cd3p_cd4@meta.data$sGroup_tag)

  # CD8
  GLIPH2Input_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-02-07/GLIPH2Input_TCR-Data_All-TRB_CD8_wAlphaChain.tsv", sep = "\t")
  GLIPH2Input_CD8 <- GLIPH2Input_CD8 %>% mutate(id = paste(GLIPH2Input_CD8$cdr3b.aa.se, GLIPH2Input_CD8$trb.v, GLIPH2Input_CD8$trb.j, sep = "_"))

  GLIPH2Output_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/2022-02-07/GLIPH2Output_TCR-Data_CD8.csv")
  GLIPH2Output_CD8 <- GLIPH2Output_CD8 %>% mutate(id = paste(GLIPH2Output_CD8$TcRb, GLIPH2Output_CD8$V, GLIPH2Output_CD8$J, sep = "_"))
  GLIPH2Output_CD8 <- GLIPH2Output_CD8[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  merge <- merge(data.frame(GLIPH2Output_CD8), GLIPH2Input_CD8[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups <- unlist(str_split(merge$clonotype.tag, ";"))
  sc_cd3p_cd8@meta.data <- sc_cd3p_cd8@meta.data %>% mutate(sGroup_tag = case_when(
    clonotype.tag %in% cts_sGroups ~ TRUE,
    TRUE ~ FALSE
  ))
  # table(sc_cd3p_cd8@meta.data$sGroup_tag)

  # Save the metadata files and objects.
  # Just run one time
  # dir.create("data")
  # saveRDS(sc_cd3p_cd8@meta.data, file = "data/sc_cd3p_cd8_mdata.rds")
  # saveRDS(sc_cd3p_cd4@meta.data, file = "data/sc_cd3p_cd4_mdata.rds")
  # saveRDS(sc_cd3p_cd8, file = "data/sc_cd3p_cd8_seurat_object.rds")
  # saveRDS(sc_cd3p_cd4, file = "data/sc_cd3p_cd4_seurat_object.rds")

  # Save an object with CD4 and CD8 cells for IA.
  # sc_cd3p_cd4$orig.celltype <- "CD4"
  # sc_cd3p_cd8$orig.celltype <- "CD8"
  # sc_cd3p <- merge(sc_cd3p_cd4, y = sc_cd3p_cd8, project = "sc_cd3p", merge.data = TRUE) # , add.cell.ids = c("CD4", "Cd8")
  # # Edata
  # edata_old_IA <- readRDS("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/raw/.upadhye_GSENNNNNN_umi.rds_2022_03_22_21_33_10")
  # saveRDS(edata_old_IA, file = "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/raw/upadhye_GSENNNNNN_umi.rds")
  # # Metadata
  # rownames(sc_cd3p@meta.data) <- paste0(gsub("[0-9]$", "", Cells(sc_cd3p)), gsub("_Gex$", "", gsub("^[0-9]{3}_", "", sc_cd3p@meta.data$origlib)))
  # saveRDS(sc_cd3p@meta.data[colnames(edata_old_IA),], file = "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/raw/upadhye_GSENNNNNN_annotation.rds")

  ## NOTE: Do not run below line. Only if you will run the integration_select_cells and integration_select_features (IA).
  ## colnames(sc_cd3p@assays$RNA) <- paste0(gsub("[0-9]$", "", Cells(sc_cd3p)), gsub("_Gex$", "", gsub("^[0-9]{3}_", "", sc_cd3p@meta.data$origlib)))
  ## saveRDS(sc_cd3p@assays$RNA, file = "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/raw/upadhye_GSENNNNNN_umi.rds")

}

{ cat(redb("### Signature Analyses ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # cp ~/Documents/liai/cancer/pbtumor/signatures_cd3p*csv /Volumes/ciro/pbtumor/info/
  # Rscript /home/ciro/scripts/functions/csvCorrect.R /home/ciro/pbtumor/info/signatures_cd3p.csv

  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R") # stats_summary_table
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
  source("/home/ciro/scripts/handy_functions/devel/utilities.R")  # vlist2df
  # source("/home/ciro/scripts/handy_functions/devel/code.R") #vlist2df
  source("/home/ciro/scripts/figease/figease.R") # fig_module_score
  source("/home/ciro/scripts/handy_functions/devel/filters.R")
  library(dplyr)

  #setwd(setwd("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures")
  dir.create("gsea")

  ## ---- GSEA ---- ###

  fname <- "../gsea/signatures_cd3p_human_v2.csv"
  slists <- readfile(fname, stringsAsFactors = FALSE)
  # Added CD4-CTL signature from Veena's paper (Patil3_COVID") to this file

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  # Modify metadata
  # CD8
  mdata_sc_cd3p_cd8$cell_classification <- as.character(mdata_sc_cd3p_cd8$RNA_snn_res.0.2)
  mdata_sc_cd3p_cd8$cell_classification[mdata_sc_cd3p_cd8$cell_classification %in% c("0","1","7")] <- "GZMK"
  mdata_sc_cd3p_cd8$cell_classification[mdata_sc_cd3p_cd8$cell_classification %in% c("6","8")] <- "TRM"
  # Group C3 and C4
  mdata_sc_cd3p_cd8$cell_classification_2 <- as.character(mdata_sc_cd3p_cd8$RNA_snn_res.0.2)
  mdata_sc_cd3p_cd8$cell_classification_2[mdata_sc_cd3p_cd8$cell_classification_2 %in% c("3","4")] <- "3_4"

  # CD4 : Setting 1
  mdata_sc_cd3p_cd4$cell_classification <- as.character(mdata_sc_cd3p_cd4$RNA_snn_res.0.4)
  mdata_sc_cd3p_cd4$cell_classification[mdata_sc_cd3p_cd4$cell_classification %in% c("3", "6", "7", "8")] <- "TCM"
  # CD4 : Setting 2
  mdata_sc_cd3p_cd4$cell_classification_2 <- as.character(mdata_sc_cd3p_cd4$RNA_snn_res.0.4)
  mdata_sc_cd3p_cd4$cell_classification_2[mdata_sc_cd3p_cd4$cell_classification_2 %in% c("3", "6", "7", "8")] <- "TCM"
  mdata_sc_cd3p_cd4$cell_classification_2[mdata_sc_cd3p_cd4$cell_classification_2 %in% c("0", "2", "4", "9")] <- "CTL"

  # Expansion Degree
  # CD4
  mdata_sc_cd3p_cd4$expDegree <- mdata_sc_cd3p_cd4$clon.size.tag
  mdata_sc_cd3p_cd4$expDegree[mdata_sc_cd3p_cd4$expDegree > 1] <- "HighExp"
  mdata_sc_cd3p_cd4$expDegree[mdata_sc_cd3p_cd4$expDegree == 1] <- "LowExp"
  # CD8
  mdata_sc_cd3p_cd8$expDegree <- mdata_sc_cd3p_cd8$clon.size.tag
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree > 1] <- "HighExp"
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree == 1] <- "LowExp"

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))


  fconfigs = list(
    # CD8
    list(result_id = "gsea/CD8/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(1:22,40)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cell_classification")
    ),
    # CD8 PDCD1p vs PDCD1n
    list(result_id = "gsea/CD8_PDCD1/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CellCycle_Best2013", "CD4CTLvsTCM_Patil2018")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("PDCD1_tag", "PDCD1n", "PDCD1p"),
      comparisons = c("PDCD1_tag")
    ),
    # CD4_1
    list(result_id = "gsea/CD4_1/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(23:26,28:32,37:39)], FUN = function(x) x[-c(1:2)] ), # CD4CTLvsTCM_Patil2018 (26)
      sample_filter = c("cell_classification", "TCM", "0", "2", "9", "4", "10", "1", "5", "11"),
      comparisons = c("cell_classification")
    ),
    # CD4_2
    list(result_id = "gsea/CD4_2/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(23:26,28:32,37:39)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification_2", "TCM", "CTL", "10", "1", "5", "11"),
      comparisons = c("cell_classification_2")
    ),
    # ---> Neo TCR
    # Neo TCR - CD4 all; clusters comparison
    list(result_id = "gsea/CD4_neoantigenTCR_all_cluster/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster")
    ),
    # Neo TCR - CD4 LG; expansion comparison
    list(result_id = "gsea/CD4_neoantigenTCR_LG_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("grade", "LG"),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD4 HG; expansion comparison
    list(result_id = "gsea/CD4_neoantigenTCR_HG_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("grade", "HG"),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD4 all; expansion comparison
    list(result_id = "gsea/CD4_neoantigenTCR_all_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD4 all; specificityGroups comparison (9)
    list(result_id = "gsea/CD4_neoantigenTCR_all_specificityGroups/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("sGroup_tag")
    ),
    # Neo TCR - CD8 all; clusters comparison
    list(result_id = "gsea/CD8_neoantigenTCR_all_cluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster")
    ),
    # Neo TCR - CD8 LG; expansion comparison
    list(result_id = "gsea/CD8_neoantigenTCR_LG_expansion/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("grade", "LG"),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD8 HG; expansion comparison
    list(result_id = "gsea/CD8_neoantigenTCR_HG_expansion/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("grade", "HG"),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD8 all; PDCD1 comparison
    list(result_id = "gsea/CD8_neoantigenTCR_all_PDCD1/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("PDCD1_tag")
    ),
    # Neo TCR - CD8 all; specificityGroups comparison (14)
    list(result_id = "gsea/CD8_neoantigenTCR_all_specificityGroups/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("sGroup_tag")
    ),
    # ---> TSCM
    # TSCM - CD4 C3,C6,C7,C8; cluster comparison
    list(result_id = "gsea/CD4_TSCM_all_cluster/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(45)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cluster", "3", "6", "7", "8"),
      comparisons = c("cluster", "3", "6", "7", "8", "REST")
    ),
    # TSCM - CD8 all; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster", "3", "4", "REST")
    ),
    # TSCM - CD8 all; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/ClustervsCluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cluster", "3", "4"),
      comparisons = c("cluster", "3", "4")
    ),
    # TSCM - CD8 all-C4; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/ClustervsREST_minus_C/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cluster", "-4"),
      comparisons = c("cluster", "3", "REST")
    ),
    # TSCM - CD8 all-C3; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/ClustervsREST_minus_C/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cluster", "-3"),
      comparisons = c("cluster", "4", "REST")
    ),
    # TSCM - CD8 all; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/C3C4vsREST/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cell_classification_2", "3_4", "REST")
    )
  )
  pp_gsea = fig_gsea(fconfigs[2], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[3:4], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  # ---> Neo TCR
  pp_gsea = fig_gsea(fconfigs[5:9], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[10:14], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # ---> TSCM
  pp_gsea = fig_gsea(fconfigs[15], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[16:20], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)

  pp_gsea = fig_gsea(fconfigs[c(9,14)], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)


  # --> Customized Radar Plots
  library(dplyr)
  library(fmsb)
  table <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/CD4_1/summary_108tests_NES.txt", sep = "\t")
  # Select our signature of interest
  table <- table %>% select(comparison, pathway, padj, NES) %>% filter(pathway == "CD4CTLvsTCM_Patil2018") #  & padj < 0.05 & NES > 0
  table_pos <- table %>% mutate(NES_m = case_when(
    padj >= 0.05 | NES < 0 ~ 0,
    TRUE ~ NES
  ))
  table_neg <- table %>% mutate(NES_m = case_when(
    padj >= 0.05 | NES > 0 ~ 0,
    TRUE ~ abs(NES)
  ))

  data <- rbind(table_pos$NES_m, table_neg$NES_m); rownames(data) <- c("positive_score", "negative_score")
  data <- rbind(rep(4,nrow(table)), rep(0,nrow(table)), data.frame(data))
  colnames(data) <- table$comparison
  # Color vector
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/CD4_1/summary_NES1_padj0.05_radar_custom.pdf")
  radarchart(data, axistype=1 ,
    pcol=colors_border[1:2] , pfcol=colors_in[1:2] , plwd=4 , plty=1, #custom polygon
    cglty=1, caxislabels=seq(0,4,1), cglwd=0.8, #custom the grid
    vlcex=0.8, title = "CD4CTLvsTCM_Patil2018" #custom labels
    )
  legend(x=0.5, y=1.35, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in, cex=1.2, pt.cex=3) # Add a legend
  dev.off()
  # blank version
  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/CD4_1/summary_NES1_padj0.05_radar_custom_blank.pdf")
  radarchart(data, axistype = 0,
    pcol=colors_border[1:2] , pfcol=colors_in[1:2] , plwd=4 , plty=1, #custom polygon
    cglty=1, cglwd=0.8, #custom the grid
    vlabels = rep("", ncol(data))  #custom labels
    )
  dev.off()

  ## ---- Module scoring analysis ---- ###

  source("/home/ciro/scripts/handy_functions/R/gsea_signature.R") # clean_feature_list, signature_scoring
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  source("/home/kmlanderos/scripts/handy_functions/R/moduleScore_functions.R") # signature_scoring (updated)
  dir.create("modulescores");
  library(ggplot2); library(dplyr)
  # sc_cd3p_lists <- readfile("/home/ciro/pbtumor/large/results/figures/gsea/sc_cd3p_mod_signatures.csv", stringsAsFactors = FALSE)
  sc_cd3p_lists <- readfile("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/gsea/signatures_cd3p_human_v2.csv", stringsAsFactors = FALSE)

  # Modify metadata
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  mdata_sc_cd3p_cd4$cell_classification <- as.character(mdata_sc_cd3p_cd4$RNA_snn_res.0.4)
  mdata_sc_cd3p_cd4$cell_classification[mdata_sc_cd3p_cd4$cell_classification %in% c("3", "6", "7", "8")] <- "TCM"
  mdata_sc_cd3p_cd4$cell_classification[mdata_sc_cd3p_cd4$cell_classification %in% c("0", "2", "4", "9")] <- "CTL"

  fconfigs = list(
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4", object = "sc_cd3p_cd4",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2),c("CD4CTLvsTCM_Patil2018", "TFH_Locci2013", "Treg_Schmiedel2018", "TH17_Seumois2020", "tcell_activation_GO.0042110.1", "tcr_signaling_RHSA202403.1", "CellCycle_Best2013.1", "TypeIandII_IFNsignaling_Seumois2020.1", "th1_signature1_arlehamn", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")]),
      axis_x = list(col = "cell_classification" , order=c("TCM", "CTL", "1", "10", "11", "5")),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    ),
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "sc_cd3p_cd8@meta.data", object = "sc_cd3p_cd8",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2), c("CD8_GZMK_Guo2018", "TRM_192_UP_HobitScience", "Gutierrez_Innateness2019", "CYTOTOXICITY.SIGNATURE_CD161.paper_matthewson..Cell")]),
      axis_x = list(col = "cell_classification"),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    )
  )
  pp_ms = fig_module_score(fconfigs[1], reductions = redu[1], violins_color = 'mean', couls = c("steelblue1", "yellow", "red",  "#670000")) # calls signature_scoring
  pp_ms = fig_module_score(fconfigs[2], reductions = redu[1], violins_color = 'mean', couls = c("steelblue1", "yellow", "red",  "#670000")) # calls signature_scoring
  pp_ms = fig_module_score(fconfigs, reductions = redu[1], violins_color = 'pct')

  #   couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000"),

  # Adjust scale of the CD4-CTL signature (0-0.2)
  signature_scores <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/modulescores/sc_cd3p_cd4_test1/signatures.csv", row.names=1)
  if( all(rownames(mdata_sc_cd3p_cd4) == rownames(signature_scores)) ) {
    mdata_sc_cd3p_cd4$CD4CTL_signature <- signature_scores[,"CD4CTLvsTCM_Patil2018"]
  } # range(mdata_sc_cd3p_cd4$CD4CTL_signature)

  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(CD4CTL_signature_adjusted = case_when(
      CD4CTL_signature > 0.2 ~ 0.2, # 0.20, Good 0.15
      TRUE ~ CD4CTL_signature)
    )
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(CD4CTL_signature_adjusted = case_when(
      CD4CTL_signature > 0.15 ~ 0.15, # 0.20, Good 0.15
      TRUE ~ CD4CTL_signature)
    )

  dir <- "modulescores/sc_cd3p_cd4_test3/"
  dir.create(dir)
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.2_umap.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c("steelblue1", "yellow", "red")) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()
  # Version 2
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.2_umap_v2.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c("#fffffa", "#f9ec8f", "#ffe080", "orange", "#FF5E0E", "red")) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()
  #
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.15_umap_v2.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c("#fffffa", "#f9ec8f", "#ffe080", "orange", "#FF5E0E", "red")) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()



}

{ cat(redb("### Markers: Dot/curtain ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dotplots");

  fconfigs = list(
   list(result_id = "dotplots/CD4/sc_cd3p_cd4_",
       edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd4", axis_x = list(col = "RNA_snn_res.0.4",
         order = c("0", "2", "9", "4", "10", "1", "5", "3", "6", "7", "8", "11")),
       features = rev(c("CD3D", "CD4", "CD8A", "CD8B", "CD40LG", "GZMA", "GZMK", "PRF1", "TBX21", "IFNG", "CXCR3", "ISG15", "XAF1",
       "IFI6", "RORA", "RORC", "CCR6", "KLRB1", "CXCL13", "BTLA", "SH2D1A", "PDCD1", "TNF", "TNFAIP3", "NFKB1", "FOXP3", "CTLA4",
       "IL2RA", "ENTPD1", "S1PR1", "SELL", "CCR7", "TCF7", "STMN1", "TOP2A", "TUBA1B"))
     ),
   list(result_id = "dotplots/CD4/sc_cd3p_cd4_ctl",
       edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
       size = c(7,5),
       object = "sc_cd3p_cd4", axis_x = list(col = "RNA_snn_res.0.4",
         order = c("0", "2", "9", "4", "10", "1", "5", "3", "6", "7", "8", "11")),
       features = rev(c("GZMB", "EOMES", "SLAMF7", "CRTAM"))
     ),
    list(result_id = "dotplots/CD8/sc_cd3p_cd8_",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "3", "2", "5", "4", "9")),
       features = rev(c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "GZMA", "ARAP2", "ISG15", "XAF1", "IFI6", "ZNF683", "ITGAE", "ITGA1",
       "IKZF2", "KLRC2", "KLRC3", "XCL1", "XCL2", "FCGR3A", "FGFBP2", "KLRD1", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "SELL", "CCR7",
       "TCF7", "STMN1", "TOP2A", "TUBA1B"))
     ),
    list(result_id = "dotplots/CD8/sc_cd3p_cd8_v2",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")),
       features = c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "GZMA", "PRF1", "ISG15", "XAF1", "IFI6", "ZNF683", "ITGAE", "ITGA1",
       "FCGR3A", "FGFBP2", "KLRD1", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "KLRG1", "CX3CR1", "PDCD1", "CXCR5", "TCF7", "SELL", "CCR7",
       "STMN1", "TOP2A", "TUBA1B")
     ),
    list(result_id = "dotplots/CD8/sc_cd3p_cd8_v3",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("7", "1", "0", "8", "6", "9", "2", "5", "3", "4")),
       features = c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "GZMA", "PRF1", "ISG15", "XAF1", "IFI6", "ZNF683", "ITGAE", "ITGA1", "STMN1",
       "TOP2A", "TUBA1B", "FCGR3A", "FGFBP2", "KLRD1", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "KLRG1", "CX3CR1", "PDCD1", "TOX", "CXCR5",
       "TCF7", "BCL2", "SLAMF6", "SELL", "CCR7")
     ),
    list(result_id = "dotplots/CD8/sc_cd3p_cd8_v4",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")),
       features = c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "GZMA", "PRF1", "ISG15", "XAF1", "IFI6", "ZNF683", "ITGAE", "ITGA1", "FCGR3A",
       "FGFBP2", "KLRD1", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "KLRC3", "KLRG1", "CX3CR1", "TCF7", "SELL", "CCR7", "STMN1", "TOP2A", "TUBA1B")
     ),
    list(result_id = "dotplots/CD8/sc_cd3p_cd8_v5",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")),
       features = c("CD3D", "CD8A", "CD8B", "CD4", "GZMK", "GZMA", "PRF1", "ISG15", "XAF1", "IFI6", "ZNF683", "ITGAE", "ITGA1", "FCGR3A", "FGFBP2",
       "KLRD1", "KLRG1", "CX3CR1", "KLRC3", "KLRB1", "SLC4A10", "IL18RAP", "RORC", "TCF7", "SELL", "CCR7", "STMN1", "TOP2A", "TUBA1B")
     ),
    list(result_id = "dotplots/CD8/hilde_paper",
       edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
       size = c(7,10),
       object = "sc_cd3p_cd8", axis_x = list(col = "RNA_snn_res.0.2",
         order = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")),
       features = c("FCER1G", "IL2RB", "KLRB1", "GZMB")
     )
   )
  pp_curtains = fig_plot_curtain(fconfigs[8], dot.scale = 12, col.min = -1.5, col.max = 1.5)
}

{ cat(redb("### Markers: dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dim_reduction_features", showWarnings = FALSE)
  fconfigs = list(
    list(result_id = "dim_reduction_features/sc_cd3p_cd4_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "PRF1", "GZMA", "IFNG", "TBX21", "EOMES", "GZMB", "SLAMF7", "CRTAM", "TNF", "GZMH", "GNLY", "CCL4", "CCL5",
      "CTSW", "ZNF683", "FCRL6", "SPON2", "CX3CR1", "S1PR5", "NKG7", "CD244", "CXCR3", "S1PR1", "SELL", "TCF7", "FOXP3", "RORA", "RORC",
      "CCR6", "KLRB1", "IRF7", "CXCL13", "BTLA", "CD40LG"),
      col = c("lightgrey", "blue")
    ),
    list(result_id = "dim_reduction_features/sc_cd3p_cd8_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd8@assays$RNA@data",
      metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "GZMK", "GZMA", "PRF1", "LAG3", "IFNG", "TNF", "CCL4", "CCL3", "FASLG", "HLA-DRB1", "TOX"),
      col = c("lightgrey", "blue")
    )
  )

  pp_markers = fig_plot_scatters(fconfigs[2])
}

{ cat(redb("### Markers: violins ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  require(tidyverse)
  dir.create("violins")

  # NOTE: Already in "Secondary global variables"; Here we just choose another threshold to establish clonally expanded cells. (Only for CD8)
  # CLONAL EXPANSION (cs == 1 low and cs > 1 high)
  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  mdata_sc_cd3p_cd8$expDegree <- mdata_sc_cd3p_cd8$clon.size.tag
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree > 1 & mdata_sc_cd3p_cd8$grade == "LG"] <- "HighExp_LG"
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree > 1 & mdata_sc_cd3p_cd8$grade == "HG"] <- "HighExp_HG"
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree == 1 & mdata_sc_cd3p_cd8$grade == "LG"] <- "LowExp_LG"
  mdata_sc_cd3p_cd8$expDegree[mdata_sc_cd3p_cd8$expDegree == 1 & mdata_sc_cd3p_cd8$grade == "HG"] <- "LowExp_HG"

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  mdata_sc_cd3p_cd4$expDegree <- mdata_sc_cd3p_cd4$clon.size.tag
  mdata_sc_cd3p_cd4$expDegree[mdata_sc_cd3p_cd4$expDegree > 1] <- "HighExp"
  mdata_sc_cd3p_cd4$expDegree[mdata_sc_cd3p_cd4$expDegree == 1] <- "LowExp"

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  fconfigs = list(
    # CD8 thold = 5
    list(result_id = "violins/sc_cd3p_cd8_thold5_TumorGrade_ExpansionGrade/violin_",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data[!is.na(sc_cd3p_cd8@meta.data$clon.size.tag) & !is.na(sc_cd3p_cd8@meta.data$orig.donor),]",
      axis_x = list(col = "expDegree", order = c("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")),
      # sample_filter = c("celltype", "CD8"),
      features = c("CD2", "COTL1", "HLA-DRB1", "PDCD1", "TOX", "GZMK", "GZMA", "CCL3", "CCL4", "FASLG", "CTSW", "LAG3", "TNF",
      "SLAMF7", "IFNG", "PRF1", "BATF", "XCL2", "CXCR6", "CXCR3", "CCL5", "FGFBP2", "GNLY", "TBX21", "ZEB2", "ZNF683", "GZMB",
      "ARAP2", "HOPX", "NFATC2", "PRDM1", "CX3CR1", "MYH9", "NFATC2", "IL2RB", "IL2RG")
    ),
    # CD8 thold = 2
    list(result_id = "violins/sc_cd3p_cd8_thold2_TumorGrade_ExpansionGrade/violin_",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor),]",
      axis_x = list(col = "expDegree", order = c("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")),
      features = c("CD2", "COTL1", "HLA-DRB1", "PDCD1", "TOX", "GZMK", "GZMA", "CCL3", "CCL4", "FASLG", "CTSW", "LAG3", "TNF",
      "SLAMF7", "IFNG", "PRF1", "BATF", "XCL2", "CXCR6", "CXCR3", "CCL5", "FGFBP2", "GNLY", "TBX21", "ZEB2", "ZNF683", "GZMB",
      "ARAP2", "HOPX", "NFATC2", "PRDM1", "CX3CR1", "MYH9", "NFATC2", "IL2RB", "IL2RG")
    ),
    # CD4 thold = 5; Exclude TREGs (cluster 5)
    list(result_id = "violins/sc_cd3p_cd4_thold5_TumorGrade_ExpansionGrade/violin_",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data[!is.na(sc_cd3p_cd4@meta.data$clon.size.tag) & !is.na(sc_cd3p_cd4@meta.data$orig.donor),]",
      axis_x = list(col = "expDegree", order = c("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")),
      sample_filter = c("cluster", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")),
      features = c("PRF1", "IFNG", "CD40LG", "GZMK", "GZMA", "TBX21", "HLA-DRB1", "CCL5", "CCL4", "TOX", "LAG3", "PDCD1", "NKG7",
      "ITGA1", "ITGAE", "EOMES", "SPON2", "ZEB2", "RBPJ", "BATF", "XCL2", "CXCL13", "XCL1", "BTLA", "STAT1", "ICOS", "TIGIT", "IL21",
      "SH2D1A", "HOPX")
    ),
    # CD4 thold = 2; Exclude TREGs (cluster 5)
    list(result_id = "violins/sc_cd3p_cd4_thold2_ExpansionGrade/violin_",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag),]",
      axis_x = list(col = "expDegree", order = c("LowExp", "HighExp")),
      sample_filter = c("cluster", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")),
      features = c("PRF1", "IFNG", "CD40LG", "GZMK", "GZMA", "TBX21", "HLA-DRB1", "CCL5", "CCL4", "TOX", "LAG3", "PDCD1", "NKG7",
      "ITGA1", "ITGAE", "EOMES", "SPON2", "ZEB2", "RBPJ", "BATF", "XCL2", "CXCL13", "XCL1", "BTLA", "STAT1", "ICOS", "TIGIT", "IL21",
      "SH2D1A", "HOPX")
    ),
    # CD4; PDCD1p vs PDCD1n
    list(result_id = "violins/sc_cd3p_cd8_PDCD1p_PDCD1n/violin_",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "mdata_sc_cd3p_cd8",
      axis_x = list(col = "PDCD1_tag", order = c("PDCD1n", "PDCD1p")),
      features = c("CCL4", "GZMK", "GZMA", "HLA-DRB1", "HLA-DPA1", "HLA-DRA", "CCL5", "COTL1", "TOX", "ITM2A", "IFNG", "CD2", "TIGIT", "PRDM1", "ICOS",
      "LAG3", "CCL3", "EOMES", "XCL2", "TNF", "ITGAE", "CCR5", "FASLG", "ARAP2", "CD70", "IL32", "CD44", "TANK", "PRF1", "SLAMF7", "CRTAM", "TNFRSF9", "TNFSF9",
      "CCR2", "CCR4", "CTLA4", "MAF")
    )
  )

  pp_violins = fig_plot_violins(
    fconfigs[4],
    theme_extra = labs(x = NULL, y = "Seurat Normalized"),
    colour_by = "pct", couls = couls_opt$red_gradient$white
  )

  # -----> PDCD1 expression by donors

  # Add expresiosn to the metadata
  genes <- c("PDCD1", "IFNG", "GZMA", "CD4", "CD8A", "CD8B")
  if(all(colnames(sc_cd3p_cd4@assays$RNA@data) == rownames(mdata_sc_cd3p_cd4))){
    mdata_sc_cd3p_cd4 <- cbind(mdata_sc_cd3p_cd4, t(as.matrix(sc_cd3p_cd4@assays$RNA@data[genes,])))
  }
  if(all(colnames(sc_cd3p_cd8@assays$RNA@data) == rownames(mdata_sc_cd3p_cd8))){
    mdata_sc_cd3p_cd8 <- cbind(mdata_sc_cd3p_cd8, t(as.matrix(sc_cd3p_cd8@assays$RNA@data[genes,])))
  }
  # Set donor order
  mdata_sc_cd3p_cd4$orig.donor <- factor(mdata_sc_cd3p_cd4$orig.donor, levels = c("BT1","BT3","BT4","BT7","BT8","BT9","BT19","BT27","BT5","BT24","BT25",
    "BT10","BT15","BT22","BT11","BT12","BT18","BT21","BT26","BT2","BT13","BT17_brain","BT20","BT23"))
  mdata_sc_cd3p_cd8$orig.donor <- factor(mdata_sc_cd3p_cd8$orig.donor, levels = c("BT1","BT3","BT4","BT7","BT8","BT9","BT19","BT27","BT5","BT24","BT25",
    "BT10","BT15","BT22","BT11","BT12","BT18","BT21","BT26","BT2","BT13","BT17_brain","BT20","BT23"))

  # PDCD1p and PDCD1n cells
  pdf("violins/sc_cd3p_cd4_expression_by_donor/violin_PDCD1.pdf", 10, 7)
   mdata_sc_cd3p_cd4 %>% filter(!is.na(orig.donor)) %>% ggplot(aes(x = orig.donor, y = PDCD1))+
    geom_violin(fill='#A4A4A4', color="darkred") + geom_jitter(size = 0.7, width = .3) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) + labs(x = "Donor", y = "PDCD1 expression") +
    stat_summary(fun.y=median, geom="point", size=0.7, color="red")
  dev.off()

  pdf("violins/sc_cd3p_cd8_expression_by_donor/violin_PDCD1.pdf")
   mdata_sc_cd3p_cd8 %>% filter(!is.na(orig.donor)) %>% ggplot(aes(x = orig.donor, y = PDCD1))+
    geom_violin(fill='#A4A4A4', color="darkred") + geom_jitter(size = 0.7, width = .3) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) + labs(x = "Donor", y = "PDCD1 expression") +
    stat_summary(fun.y=median, geom="point", size=0.7, color="red")
  dev.off()

  # PDCD1p cells
  pdf("violins/sc_cd3p_cd4_expression_by_donor/violin_PDCD1_pos.pdf", 10, 7)
   mdata_sc_cd3p_cd4 %>% filter(!is.na(orig.donor) & PDCD1>0) %>% ggplot(aes(x = orig.donor, y = PDCD1))+
    geom_violin(fill='#A4A4A4', color="darkred") + geom_jitter(size = 0.7, width = .3) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) + labs(x = "Donor", y = "PDCD1 expression") +
    stat_summary(fun.y=median, geom="point", size=0.7, color="red")
  dev.off()

  pdf("violins/sc_cd3p_cd8_expression_by_donor/violin_PDCD1_pos.pdf")
   mdata_sc_cd3p_cd8 %>% filter(!is.na(orig.donor) & PDCD1>0) %>% ggplot(aes(x = orig.donor, y = PDCD1))+
    geom_violin(fill='#A4A4A4', color="darkred") + geom_jitter(size = 0.7, width = .3) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) + labs(x = "Donor", y = "PDCD1 expression") +
    stat_summary(fun.y=median, geom="point", size=0.7, color="red")
  dev.off()


  # -----> Gen expression for BT23 clones (CD8)

  library(ggplot2)
  genes <- c("IFNG", "PDCD1", "GZMA", "GZMK", "PRF1", "TNF", "CCL4", "HLA-DRB1", "LAG3", "TOX")
  donor <- "BT23"
  genes_expr <- sc_cd3p_cd8@assays$RNA@data[genes,]
  # Get the cells for each category in the violin plot.
  BT23_cells <- list(
    All = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clon.size.tag),]),
    Expanded = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clon.size.tag) & mdata_sc_cd3p_cd8$clon.size.tag > 1,]),
    Non_expanded = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clon.size.tag) & mdata_sc_cd3p_cd8$clon.size.tag == 1,]),
    Chosen = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clon.size.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype22088", "clonotype22421", "clonotype30069", "clonotype30163", "clonotype30266", "clonotype22574", "clonotype23396", "clonotype29218", "clonotype30408", "clonotype25268", "clonotype30240", "clonotype26322", "clonotype30105", "clonotype30074", "clonotype30190") ,]),
    TCR1 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype22088", "clonotype22421") ,]),
    TCR2 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype30069", "clonotype30163", "clonotype30266") ,]),
    TCR3 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype22574", "clonotype23396") ,]),
    TCR4 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype29218", "clonotype30408") ,]),
    TCR5 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype25268", "clonotype30240") ,]),
    TCR6 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype26322", "clonotype30105") ,]),
    TCR7 = rownames(mdata_sc_cd3p_cd8[mdata_sc_cd3p_cd8$orig.donor == donor & !is.na(mdata_sc_cd3p_cd8$orig.donor) & !is.na(mdata_sc_cd3p_cd8$clonotype.tag) & mdata_sc_cd3p_cd8$clonotype.tag %in% c("clonotype30074", "clonotype30190") ,])
  )

  # Expression values
  # expr_TCRtag_list <- lapply(BT23_cells, FUN = function(cells) { as.matrix(genes_expr[, cells])
  # } )
  # df <- reshape2::melt(as.data.frame(expr_TCRtag_list)); df$genes <- genes;

  expr_TCRtag_list_test <- lapply(BT23_cells, FUN = function(cells) { as.vector(as.matrix(genes_expr[, cells]))
  } );
  df <- data.frame(value = unlist(expr_TCRtag_list_test))
  # Add genes
  df$genes <- rep(genes, times = sum(sapply(BT23_cells, length)))

  # Add TCRtag
  # df$TCRtag <- factor(rep(names(BT23_cells), times = sapply(BT23_cells, length)), levels = c("All", "Expanded", "Non_expanded", "Chosen", "TCR1", "TCR2", "TCR3", "TCR4", "TCR5", "TCR6", "TCR7"))
  df$TCRtag <- factor(rep(names(BT23_cells), times = 10*sapply(BT23_cells, length)), levels = c("All", "Expanded", "Non_expanded", "Chosen", "TCR1", "TCR2", "TCR3", "TCR4", "TCR5", "TCR6", "TCR7"))

  # Add pct of post
  pct_TCRtag_list <- lapply(BT23_cells, FUN = function(cells) {
    apply(as.matrix(genes_expr[, cells]), 1, function(x) rep(100*sum(x>0)/length(x), times = length(cells)) )
  } );
  pct_TCRtag_list <- lapply(pct_TCRtag_list_tmp, function(x) as.vector(t(x)))
  df$pct <- unlist(pct_TCRtag_list_tmp)

  # Pct of Pos
  # lapply(BT23_cells, FUN = function(cells) {
  #   apply(as.matrix(genes_expr[, cells]), 1, function(x) 100*sum(x>0)/length(x) )
  # } );


  # ---- Plot
  dir.create("tcr")
  pdf("tcr/TCRtag_GeneExpression.pdf", 10, 7)
  # gene = "GZMK"
  for (gene in genes){
    p <- ggplot(df %>% filter(genes == gene), aes(x=TCRtag, y=value, fill = pct)) +
      scale_fill_gradientn(colours = c("yellow", "orange", "red", "red4"), limits = c(0,100)) +
      geom_violin(trim=FALSE) + labs(title = gene, x = "TCR group", y = "Normalized Seurat Expression") +
      geom_boxplot(width=0.1) + # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) + # geom_jitter()
      theme(axis.text.x = element_text(angle = 30, hjust = 1))#rotatedAxisElementText(30,'x'))
    print(p)
  }
  dev.off()

}

{ cat(redb("### Blanks volcano ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R") # volplot
  source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
  source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
  source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
  source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes
  source("/home/ciro/scripts/handy_functions/devel/volcano.R")
  require(ggplot2)

  # NOTE: Run for the 2 files.
  # file <- "/home/kmlanderos/large/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/HighGrade/highvslow/mastlog2cpm_results.csv"
  # file <- "/home/kmlanderos/large/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/LowGrade/highvslow/mastlog2cpm_results.csv"
  file <- "/home/kmlanderos/large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD8_PDCD1/All/PDCD1pvsPDCD1n/mastlog2cpm_results.csv"
  # file <- "/home/kmlanderos/large/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD4_subclustering/LowGrade_and_HighGrade_thold2/highvslow/mastlog2cpm_results.csv"
  results <- read.csv(file)

  group1 = "PDCD1p" #"high" #
  group2 = "PDCD1n" #"low" #
  padjthr = 0.05
  fcthr = 0.25
  cols_filt = "minExp"
  output = "./"
  verbose = TRUE
  return_report = FALSE
  ##############
  pseuc = 1
  dtype = "CST"
  # volcano parameters
  showgenes = NULL
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
    genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
  if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
  resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
  if("pct_diff" %in% colnames(resdf))
    resdf$pct_diff <- abs(resdf$pct_diff)
    resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
  if(is.null(showgenes)){
    showgenes <- unique(c(
      bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
      head(resdf[genes2plot, "gene"], 10)))
  }

   cat("---------------------- Volcano -------------------------\n")

  # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
  source("/home/kmlanderos/scripts/handy_functions/devel/volplot.R")

  # pdf("CD8_covsex_volcano_blankq.pdf") #Change name for CD8
  # pdf(paste0(dirname(file),"/mastlog2cpm_volcano_blank.pdf"), 10, 10)
  pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots_blank.pdf"), 10, 10)
    volplot(
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
    )
    # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
    #   title = paste(group2, "(-) vs ", group1, "(+)"))
  dev.off()

  # With names (big dots)
  pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots.pdf"), 10, 10)
    volplot_wnames(
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
    )
    # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
    #   title = paste(group2, "(-) vs ", group1, "(+)"))
  dev.off()

{ cat(redb("### Main dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  dir.create("dim_reduction")
  source("/home/ciro/scripts/clustering/R/plotting.R")

  sc_cd3p_cd4@meta.data$orig.Diagnosis <- gsub(" ", "_", sc_cd3p_cd4@meta.data$orig.Diagnosis)
  sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass <- gsub(" ", "_", sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass)
  sc_cd3p_cd8@meta.data$orig.Diagnosis <- gsub(" ", "_", sc_cd3p_cd8@meta.data$orig.Diagnosis)
  sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass <- gsub(" ", "_", sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass)

  fconfigs = list(
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = c("grade"),
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = c("grade"), facet = "grade",
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = c("orig.Diagnosis"),
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = c("cell_classification"), facet = "orig.Diagnosis",
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "grade",
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "grade", facet = "grade",
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "orig.Diagnosis",
      redu = redu[[1]]
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cell_classification", facet = "orig.Diagnosis",
      redu = redu[[1]]
    ),
    # ----- FINAL FIGURES ----- #
    # CD4
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = c("cluster"),
      colors = c("0" = "#c5b8ff", "1" = "#5fcfab", "7" = "#e57f08", "6" = "#ffaf4d", "8" = "#b35806", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#369af7",
       "9" = "#5d1d88", "11" = "#9c0006")
      # V1
      # colors = c("0" = "#bca0dc", "2" = "#b491c8", "4" = "#7c5295", "9" = "#663a82", "5" = "#fdfd96", "3" = "#c7ea46", "6" = "#98fb98", "7" = "#abc32f", "8" = "#00a572",
      #  "1" = "#d6deff", "10" = "lightsalmon1", "11" = "gray70")
      # colors = c("0" = "#7373F3", "2" = "#54D8FD", "4" = "#02AEBD", "9" = "#015DBB", "5" = "#FF6F00", "3" = "#FFB94F", "6" = "#D655A7", "7" = "#44A306", "8" = "#C994C7",
      #  "1" = "#969696", "10" = "#969696", "11" = "#525252")
    ),
    # CD8
    list(result_id = "dim_reduction/sc_cd3p_cd8",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster",
      colors = c("0" = "#72e5d2", "1" = "#3bb8ae", "7" = "#007e74", "6" = "#ffa7d7", "8" = "#ff6ede", "2" = "#fcee62", "5" = "#9e89e8", "4" = "#4da7fd", "3" = "#ff971f",
       "9" = "#9c0006")
       # V1
       # colors = c("0" = "cadetblue1", "1" = "deepskyblue", "7" = "dodgerblue2", "6" = "magenta", "8" = "pink", "2" = "khaki1", "5" = "lightgreen", "4" = "mediumpurple1", "3" = "tan1",
       #  "9" = "gray75")
       # V2
       # colors = c("0" = "#80CDC1", "1" = "#35978F", "7" = "#01665E", "6" = "#F1B6DA", "8" = "#DE77AE", "2" = "#FED976", "5" = "#5AAE61", "4" = "#C2A5CF", "3" = "#FC8D59",
       #  "9" = "#878787")
       # V3
       # colors = c("0" = "#80CDC1", "1" = "#35978F", "7" = "#007E74", "6" = "#E4A4CB", "8" = "#D37EAB", "2" = "#F8DA7F", "5" = "#A18DE0", "4" = "#72B0EA", "3" = "#A8D792",
       #  "9" = "#878787")
    )
  )
  pp_main_dr = fig_plot_base(
    fconfigs[-c(4,8)],
    theme_extra = function(x){
      p <- x + geom_point(size = 0.4, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2", color = NULL) +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left")) +
      xlim(-10, 10) + ylim(-10, 10)
      # Seurat::LabelClusters(plot = p, id = "cluster", color = "black") # , fontface = "bold"  , point.padding = 1, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc"))
    }
  )
  pp_main_dr = fig_plot_base(
    fconfigs[c(4,8)],
    theme_extra = function(x){
      p <- x + geom_point(size = 1.2, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2", color = NULL) +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left")) +
      xlim(-10, 10) + ylim(-10, 10)
      # Seurat::LabelClusters(plot = p, id = "cluster", color = "black") # , fontface = "bold"  , point.padding = 1, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc"))
    }
  )
  pp_main_dr = fig_plot_base(
    fconfigs[c(9)],
    theme_extra = function(x){
      p <- x + geom_point(size = 0.4, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2", color = NULL) +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left")) #+
      # xlim(-10, 10) + ylim(-10, 10)
      # Seurat::LabelClusters(plot = p, id = "cluster", color = "black") # , fontface = "bold"  , point.padding = 1, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc"))
    }
  )

}

{ cat(redb("### NK T-cells receptors ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  library(data.table)
  # Read GLIPH2Input table that has VDJ genes' info
  path = "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-02-13/"
  GLIPH_input_cd8  <- data.table(read.csv(paste0(path, ".intermidiate_file_TCRsimilarity_CD8.tsv"), sep = "\t"))
  GLIPH_input_cd4  <- data.table(read.csv(paste0(path, ".intermidiate_file_TCRsimilarity_CD4.tsv"), sep = "\t"))

  # CD8
  beta.info_cd8 <- GLIPH_input_cd8[
    ,
    .(
      clonotype.tag=paste0(
        unique(
          sort(paste(clonotype.tag, sep='|'))
        ), collapse=';'
      ),
      orig.donor=paste0(
        unique(
          sort(paste(orig.donor, sep='|'))
        ), collapse=';'
      ),
      size=sum(size)
     ),
    by=.(cdr3b.aa.seq, trb.v, trb.j)
  ]
  # CD4
  beta.info_cd4 <- GLIPH_input_cd4[
    ,
    .(
      clonotype.tag=paste0(
        unique(
          sort(paste(clonotype.tag, sep='|'))
        ), collapse=';'
      ),
      orig.donor=paste0(
        unique(
          sort(paste(orig.donor, sep='|'))
        ), collapse=';'
      ),
      size=sum(size)
     ),
    by=.(cdr3b.aa.seq, trb.v, trb.j)
  ]

  beta.info_cd4[beta.info_cd4$tra.v == "TRAV10" & beta.info_cd4$tra.j == "TRAJ18"  & beta.info_cd4$trb.v == "TRBV25",]
  beta.info_cd8[beta.info_cd8$tra.v == "TRAV10" & beta.info_cd8$tra.j == "TRAJ18"  & beta.info_cd8$trb.v == "TRBV25",]

}

{ cat(redb("### Markers: contour ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # /home/ciro/asthma_airways/scripts/final_figures_cd4/init.R
  dir.create("scatter_contour")
  require(dplyr)

  # ----> CD4

  # Adjust the metadata
  sc_cd3p_cd4_mdata <- sc_cd3p_cd4@meta.data
  sc_cd3p_cd4_mdata <- sc_cd3p_cd4_mdata[!is.na(sc_cd3p_cd4_mdata$clon.size.tag),]
  # sc_cd3p_cd4_mdata$ExpDegree <- ifelse(sc_cd3p_cd4_mdata$clon.size.tag > 4, "high", "low")
  sc_cd3p_cd4_mdata$ExpDegree <- ifelse(sc_cd3p_cd4_mdata$clon.size.tag > 1, "high", "low")

  # Gene comparisons
  tmp <- list(
    list(x = "PRF1", y = c("GZMA", "GZMB", "TNF", "IFNG")),
    list(x = "GZMA", y = c("TNF", "IFNG")),
    list(x = "TNF", y = c("IFNG"))
  )

  # sc_cd4_HighlyExpanded. Excluding TREGs (cluster 5).
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd4_HighlyExpanded/",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster %in% c('0', '1', '2', '3', '4', '6', '7', '8', '9', '10', '11'),]", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
      features = x,
      # sample_filter = c("cluster", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")),
      sample_filter = list(c("ExpDegree", "high"))
      # sample_filter = setNames(object = list("high", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")), c("ExpDegree", "cluster"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )
  # sc_cd4_LowlyExpanded. Excluding TREGs (cluster 5).
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd4_LowlyExpanded/",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster %in% c('0', '1', '2', '3', '4', '6', '7', '8', '9', '10', '11'),]", # [sc_cd3p_cd4_mdata$ExpDegree=='low',]
      features = x,
      # sample_filter = c("cluster", "0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11"),
      sample_filter = list(c("ExpDegree", "low"))
      # sample_filter = setNames(object = list("low", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")), c("ExpDegree", "cluster"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )

  # sc_cd4_Expanded. Excluding TREGs (cluster 5).
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd4_Expanded/",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster %in% c('0', '1', '2', '3', '4', '6', '7', '8', '9', '10', '11'),]", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
      features = x,
      # sample_filter = c("cluster", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")),
      sample_filter = list(c("ExpDegree", "high"))
      # sample_filter = setNames(object = list("high", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")), c("ExpDegree", "cluster"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )
  # sc_cd4_NonExpanded. Excluding TREGs (cluster 5).
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd4_NonExpanded/",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster %in% c('0', '1', '2', '3', '4', '6', '7', '8', '9', '10', '11'),]", # [sc_cd3p_cd4_mdata$ExpDegree=='low',]
      features = x,
      # sample_filter = c("cluster", "0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11"),
      sample_filter = list(c("ExpDegree", "low"))
      # sample_filter = setNames(object = list("low", c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11")), c("ExpDegree", "cluster"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )

  # Other way of doing the same
    # fconfigs = list(
    #   list(result_id = "scatter_contour/sc_cd4_HighlyExpanded/",
    #     edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata", filters = list(c("ExpDegree", "high")),
    #     axis_x = "cluster", features = list(x = c("PRF1"), y = c("GZMA", "GZMB", "TNF", "IFNG"))
    #   ),
    #   list(result_id = "scatter_contour/sc_cd4_LowlyExpanded/", filters = list(c("celltype", "CD8"), c("ExpDegree", "high")),
    #     edata = "sc_cd3p@assays$RNA@data", metadata = "sc_cd3p_tcr",
    #     axis_x = "cluster", features = list(x = c("IFNG"), y = c("TNF", "PRF1"))
    #   )
    # )

  # ----> CD8

  # Adjust the metadata
  sc_cd3p_cd8_mdata <- sc_cd3p_cd8@meta.data
  sc_cd3p_cd8_mdata <- sc_cd3p_cd8_mdata[!is.na(sc_cd3p_cd8_mdata$clon.size.tag),]
  sc_cd3p_cd8_mdata$ExpDegree <- ifelse(sc_cd3p_cd8_mdata$clon.size.tag > 4, "high", "low")
  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  sc_cd3p_cd8_mdata[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  sc_cd3p_cd8_mdata <- sc_cd3p_cd8_mdata %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  # Gene comparisons
  tmp <- list(
    list(x = "PDCD1", y = c("IFNG", "TNF", "PRF1", "TOP2A"))
  )
  # sc_cd8_HighlyExpanded.
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd8_HighlyExpanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
      features = x,
      sample_filter = list(c("ExpDegree", "high"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )
  # sc_cd8_LowlyExpanded.
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd8_LowlyExpanded/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
      features = x,
      sample_filter = list(c("ExpDegree", "low"))
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )
  # sc_cd8_All.
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd8_All/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
      features = x
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma"), type = "percent")
  )
  tmp <- list(
    list(x = "IFNG", y = c("TNF", "PRF1")),
    list(x = "PRF1", y = c("GZMA"))
  )
  # PDCD1p
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd8_PDCD1p/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata[sc_cd3p_cd8_mdata$PDCD1_tag=='PDCD1p',]",
      features = x
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)), type = "percent")
  )
  # PDCD1n
  fconfigs = lapply(tmp, function(x){
    list(result_id = "scatter_contour/sc_cd8_PDCD1n/",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8_mdata[sc_cd3p_cd8_mdata$PDCD1_tag=='PDCD1n',]",
      features = x
    )
  })
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,1)), type = "percent")
  )

  # -----> Coexpression UMAPs
  # CD8 (test)
  sc_cd3p_cd8@meta.data$PDCD1_tag <- sc_cd3p_cd8_mdata$PDCD1_tag
  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/scatter_contour/coexpression_UMAP.pdf", 12,10)
  FeaturePlot(sc_cd3p_cd8, features = c("PRF1", "GZMA"), blend = TRUE, , split.by = "PDCD1_tag", order = T)
  dev.off()

  # CD4
  sc_cd3p_cd4@meta.data$Exp_tag <- ifelse(sc_cd3p_cd4@meta.data$clon.size.tag > 1, "exp", "non-exp")
  sc_cd3p_cd4@meta.data$PDCD1_tag <- sc_cd3p_cd4_mdata$PDCD1_tag
  sc_cd3p_cd4_subset <- subset(x = sc_cd3p_cd4, subset = Exp_tag %in% c("exp", "non-exp"))
  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/scatter_contour/coexpression_UMAPs_CD4.pdf", 12,10)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("PRF1", "GZMA"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("PRF1", "GZMB"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("PRF1", "TNF"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("PRF1", "IFNG"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("GZMA", "TNF"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("GZMA", "IFNG"), blend = TRUE, , split.by = "Exp_tag", order = T)
  FeaturePlot(sc_cd3p_cd4_subset, features = c("TNF", "IFNG"), blend = TRUE, , split.by = "Exp_tag", order = T)
  dev.off()

}

{ cat(redb("### TCR plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source('/home/ciro/scripts/handy_functions/devel/overlap.R') # fig_plot_overlaps
  #source("/home/kmlanderos/scripts/figease/figease.R")
  suppressPackageStartupMessages(library(dplyr)); library(stringr); library(data.table)
  dir.create("tcr")

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

  # Modify the metadata
  # Put all TH1 together, so they have the same color.
  # mdata_sc_cd3p_cd4$cell_classification[mdata_sc_cd3p_cd4$cell_classification %in% c("TH1_CTL_IFN", "TH1_TFH")] <- "TH1" #

  # sc_mdata_tcr = data.frame(data.table::rbindlist(lapply(
  #   X = c("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd4_mdata.rds",
  #         "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd8_mdata.rds"),
  #   FUN = function(x){
  #     y <- readRDS(x); y$cellname = rownames(y); y
  # })), stringsAsFactors = FALSE)
  # rownames(sc_mdata_tcr) <- sc_mdata_tcr$cellname; sc_mdata_tcr$cellname <- NULL
  # #sc_mdata_tcr <- joindf(sc_mdata_tcr, sc_cd3p@meta.data)

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
      axis_x = "cluster", axis_y = "clonotype.tag", # filters = c("celltype", "CD8"),
      col = c("GZMK" = "#54D8FD", "TRM" = "#FFB94F", "INNATE-LIKE_NKR" = "#D655A7", "INNATE-LIKE_XCL1-2" = "#44A306", "TCM" = "#969696", "MAIT" = "#C994C7", "Cell_Cycle" = "#464E2E"),
      vars2col = c("cell_classification"),
      element_type = c("clones", "cells")),
    list(result_id = "tcr/sc_cd3p_cd4/", edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag), ]",
      axis_x = "cluster", axis_y = "clonotype.tag", # filters = c("celltype", "CD4"),
      col = c("TCM" = "#969696", "CTL1" = "#66C2A4", "CTL_TH17" = "#B4A7D6", "CTL_TH1" = "#B25785", "CTL_IFN" = "#703352", "TNF" = "#EA9999", "TREG" = "#FFE599", "TFH" = "#DFB9CD", "Cell_Cycle" = "#464E2E"),
      vars2col = c("cell_classification"),
      element_type = c("clones", "cells"))
  )
  ## 6) Upsetplot: Clones shared between clusters ## ---------------------------
  pp_overlaps = fig_plot_overlaps(fconfigs)
  # pp_overlaps = fig_plot_overlaps(fconfigs, nintersects = 70, empty.intersections = "on") Complete plot
  pp_overlaps = fig_plot_overlaps(lapply(fconfigs, function(x) {x$vars2col <- c("cluster"); return(x)})) # Color by cluster
  # file.remove(list.files("tcr", pattern = "_features", full.names = TRUE, recursive = TRUE))
  # ## 1) Boxplot for clone size ; X: clusters ## --------------------------------
  # for (i in 1:length(fconfigs)){ fconfigs[[i]]$sufix <- "boxplot_"; fconfigs[[i]]$axis_y = "clon.size.tag" }
  # pp_clones = fig_plot_base(
  #   fconfigs, return_plot = TRUE, verbose = TRUE,
  #   theme_extra = function(x){
  #     dots <- 100; set.seed(27);
  #     d2show <- x$data[unname(sample_even(x$data, rlang::as_name(x$mapping$x), -dots)), ]
  #     x <- x + geom_boxplot(outlier.shape = NA) +
  #       geom_jitter(data = d2show, size = 0.3, position = position_jitter(0.2), color = 'black') +
  #       labs(x = NULL, y = "Clone size") + theme(legend.position = "none")
  #   }
  # )
  # ## 4) Piechart: Clonally expanded and non-expanded per cluster ## ------------
  # for (i in 1:length(fconfigs)){ fconfigs[[i]]$sufix <- "piechart_"; fconfigs[[i]]$axis_y = "TCR.tag" }
  # pp_clones = fig_plot_base(
  #   fconfigs, return_plot = TRUE, verbose = TRUE,
  #   theme_extra = pie_fun,
  #   plot_blank_fun = function(x){
  #     plot_rm_layer(x) + theme(legend.position = "none", strip.text = element_blank())
  #   }
  # )
  #
  # ## 7) Barplots: Clonally expanded and non-expanded per cluster ## ------------
  # # eval(expr = parse(text = grep("coord_polar.*|facet_wrap.*", deparse(pie_fun), value = TRUE, invert = TRUE)))
  # for (i in 1:length(fconfigs)) fconfigs[[i]]$sufix <- "barplot_cells_"
  # pp_clones = fig_plot_base(fconfigs, theme_extra = bar_fun)
  # for (i in 1:length(fconfigs)) fconfigs[[i]]$sufix <- "barplot_clones_"
  # pp_clones = fig_plot_base(fconfigs, theme_extra = function(x){ x$data <- x$data[!duplicated(x$data$clonotype.tag), ]; bar_fun(x) })
  #
  # ## 3) Boxplot for clonal size per donor; X: donors ## ------------------------
  # for (i in 1:length(fconfigs)){
  #   fconfigs[[i]]$sufix <- "boxplot_";
  #   fconfigs[[i]]$axis_x = list(col = "orig.donor",
  #     order = c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5",
  #     "BT24", "BT25", "BT10", "BT15", "BT22", "BT11", "BT12", "BT18", "BT21",
  #     "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23"))
  #   fconfigs[[i]]$vars2col = "orig.donor"; fconfigs[[i]]$axis_y = "clon.size.tag"
  # }
  # pp_clones = fig_plot_base(
  #   fconfigs, return_plot = TRUE, verbose = TRUE,
  #   theme_extra = function(x){
  #     x$data <- x$data[!is.na(x$data$Identity), ]; dots <- 100; set.seed(27);
  #     print(rlang::as_name(x$mapping$x))
  #     d2show <- x$data[unname(sample_even(x$data, rlang::as_name(x$mapping$x), -dots)), ]
  #     x <- x + geom_boxplot(outlier.shape = NA) +
  #     geom_jitter(data = d2show, size = 0.3, position = position_jitter(0.2), color = 'black') +
  #     labs(x = NULL, y = "Clone size") +
  #     theme(legend.position = "none", axis.text.x = element_text(angle = 90))
  #   }
  # )
  # ## 5) Piechart: Clonally expanded and non-expanded per donor ## --------------
  # for (i in 1:length(fconfigs)){ fconfigs[[i]]$sufix <- "piechart_"; fconfigs[[i]]$axis_y = "TCR.tag" }
  # suppressPackageStartupMessages(library(dplyr))
  # pp_clones = fig_plot_base(
  #   fconfigs, return_plot = TRUE, verbose = TRUE,
  #   theme_extra = pie_fun,
  #   plot_blank_fun = function(x){
  #     plot_rm_layer(x) + theme(legend.position = "none", strip.text = element_blank())
  #   }
  # )
  ## 2) Scatter plot (umap) for clone size ## ----------------------------------
  # Change names of  just so they are not too large.
  mdata_sc_cd3p_cd8$cell_classification[mdata_sc_cd3p_cd8$cell_classification == "INNATE-LIKE_NKR"] <- "IL-NKR"
  mdata_sc_cd3p_cd8$cell_classification[mdata_sc_cd3p_cd8$cell_classification == "INNATE-LIKE_XCL1-2"] <- "IL-XCL1/2"

  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$sufix <- "clone_size_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification"
  }
  pp_clones = fig_plot_base(
    fconfigs, return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = clon.size.tag)) +
        geom_point(aes(size = clon.size.tag), shape = 1, color = "gray50", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6)) +
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )

  ## Scatter plot (umap) for cluster clone size ## ----------------------------------
  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$sufix <- "cluster_clone_size_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification"
  }

  # Create the "cluster_clones_size" column in the metadata of the 2 objects.
  clonotype_cd4 <- sc_cd3p_cd4@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd4) <- clonotype_cd4$clonotype.tag; clonotype_cd4$clonotype.tag <- NULL

  clonotype_cd8 <- sc_cd3p_cd8@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
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
    cluster_clone_size >= 40 ~ 40,
    TRUE ~ cluster_clone_size
  ))

  # Update CD8 mdata
  cluster_clone_size <- c()
  for (i in 1:dim(mdata_sc_cd3p_cd8)[1]){
    clonotype <- mdata_sc_cd3p_cd8[i,"clonotype.tag"]
    c <- as.character(mdata_sc_cd3p_cd8[i,"cluster"])
    cluster_clone_size <- c(cluster_clone_size, clonotype_cd8[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  mdata_sc_cd3p_cd8$cluster_clone_size <- as.numeric(cluster_clone_size) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(cluster_clone_size = case_when(
    cluster_clone_size >= 40 ~ 40,
    TRUE ~ cluster_clone_size
  ))

  pp_clones = fig_plot_base(
    fconfigs[1], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray45", alpha = 0.1, stroke = 0.4) + # gray47
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  pp_clones = fig_plot_base(
    fconfigs[2], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray20", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )

  ## Cluster Clone Size (facet LG - HG)
  fconfigs2 <- fconfigs
  for (i in 1:length(fconfigs2)) {
    fconfigs2[[i]]$sufix <- "cluster_clone_size_HG"; fconfigs2[[i]]$axis_x = "UMAP_1"; fconfigs2[[i]]$axis_y = "UMAP_2"
    fconfigs2[[i]]$vars2col = "cell_classification"
  }
  fconfigs2[[1]]$metadata <- "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag) & !is.na(mdata_sc_cd3p_cd8$grade), ]"
  fconfigs2[[2]]$metadata <- "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag) & !is.na(mdata_sc_cd3p_cd4$grade), ]"

  pp_clones = fig_plot_base(
    fconfigs2, return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray50", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") +
        ggforce::facet_wrap_paginate(~ grade, ncol = 1, nrow = 1, page = 1)
        # facet_grid(cols = vars(grade)) # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )

  for (i in 1:length(fconfigs2)) {
    fconfigs2[[i]]$sufix <- "cluster_clone_size_LG"; fconfigs2[[i]]$axis_x = "UMAP_1"; fconfigs2[[i]]$axis_y = "UMAP_2"
    fconfigs2[[i]]$vars2col = "cell_classification"
  }

  pp_clones = fig_plot_base(
    fconfigs2, return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "gray50", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") +
        ggforce::facet_wrap_paginate(~ grade, ncol = 1, nrow = 1, page = 2)
        # facet_grid(cols = vars(grade)) # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )

  ## Scatter plot (umap) for cluster clone size; Colored by gene expression ## ----------------------------------
  dir.create("tcr/sc_cd3p_cd4/markers/LG-HG", recursive = T); dir.create("tcr/sc_cd3p_cd8/markers/LG-HG", recursive = T);

  genes <- c("ZNF683", "FCGR3A", "FGFBP2", "IFNG", "PDCD1", "GZMA", "GZMK", "PRF1", "TNF", "CCL4", "HLA-DRB1", "LAG3", "TOX")
  if(all(colnames(sc_cd3p_cd4@assays$RNA@data) == rownames(mdata_sc_cd3p_cd4))){
    mdata_sc_cd3p_cd4 <- cbind(mdata_sc_cd3p_cd4, t(as.matrix(sc_cd3p_cd4@assays$RNA@data[genes,])))
    colnames(mdata_sc_cd3p_cd4) <- gsub("HLA-DRB1", "HLA_DRB1", colnames(mdata_sc_cd3p_cd4))
  }
  if(all(colnames(sc_cd3p_cd8@assays$RNA@data) == rownames(mdata_sc_cd3p_cd8))){
    mdata_sc_cd3p_cd8 <- cbind(mdata_sc_cd3p_cd8, t(as.matrix(sc_cd3p_cd8@assays$RNA@data[genes,])))
    colnames(mdata_sc_cd3p_cd8) <- gsub("HLA-DRB1", "HLA_DRB1", colnames(mdata_sc_cd3p_cd8))
  }
  fconfigs_ <- rep(fconfigs, each = length(genes))
  for (i in 1:length(genes)){
    fconfigs_[[i]]$vars2col <-  gsub("-", "_", genes[i]) #NULL #
    fconfigs_[[length(genes) + i]]$vars2col <-  gsub("-", "_", genes[i]) #NULL #
  }

  for (fconfig in fconfigs_){ # fconfig <- fconfigs_[[2]]
    cat(fconfig$vars2col, "\n")
    df <- eval(parse(text=fconfig$metadata)) %>% arrange(eval(parse(text=fconfig$vars2col)))
    pdf(paste0(fconfig$result_id, "markers/", fconfig$sufix, fconfig$vars2col, ".pdf"))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col))), alpha = 0.8) +
      geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
      scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
      # scale_colour_gradient(colorRampPalette(c("red", "yellow", "white"))) +
      scale_colour_gradient2(low = "steelblue1", mid = "yellow", high = "red", midpoint = 2.2) +
      labs(x = "Dim 1", y = "Dim 2", size = "Clone\nSize", color = fconfig$vars2col)
    print(p)
    dev.off()
    pdf(paste0(fconfig$result_id, "markers/", fconfig$sufix, fconfig$vars2col, "_blank.pdf"))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col)) ), alpha = 0.8, show.legend = FALSE) +
      geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4, show.legend = FALSE) +
      scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
      # scale_colour_gradient(colorRampPalette(c("red", "yellow", "white"))) +
      scale_colour_gradient2(low = "steelblue1", mid = "yellow", high = "red", midpoint = 2.2) +
      labs(x = "", y = "") + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
    print(p)
    dev.off()
  }
  # Devided by LG-HG
  for (fconfig in fconfigs_){
    cat(fconfig$vars2col, "\n")
    df <- eval(parse(text=fconfig$metadata)) %>% arrange(eval(parse(text=fconfig$vars2col)))
    pdf(paste0(fconfig$result_id, "markers/LG-HG/", fconfig$sufix, fconfig$vars2col, ".pdf"))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col))), alpha = 0.8) +
      geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
      scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
      scale_colour_gradient(low = "yellow", high = "red") +
      labs(x = "Dim 1", y = "Dim 2", size = "Clone\nSize", color = fconfig$vars2col) +
      facet_wrap( ~ grade, ncol=2)
    print(p)
    dev.off()
    pdf(paste0(fconfig$result_id, "markers/LG-HG/", fconfig$sufix, fconfig$vars2col, "_blank.pdf"))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col)) ), alpha = 0.8, show.legend = FALSE) +
      geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4, show.legend = FALSE) +
      scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
      scale_colour_gradient(low = "yellow", high = "red") +
      labs(x = "", y = "") + theme(axis.text.x=element_blank(), axis.text.y=element_blank() ) +
      facet_wrap( ~ grade, ncol=2) # facet_grid(. ~ grade)
    print(p)
    dev.off()
  }

  # ----> Stacked Barplot - clonotypes per donor

{ cat(redb("### GG-paired plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  library(ggpubr)
  dir.create("heatmaps")

  # ---> PDCD1n and PDCD1p Gene expression (CD8)

  tmp <- as.matrix(sc_cd3p_cd8@assays$RNA["PDCD1",])
  cells <- colnames(sc_cd3p_cd8@assays$RNA["PDCD1",])
  PDCD1p_cells <- cells[tmp>0]
  PDCD1n_cells <- cells[tmp==0]

  genes <- c("CCL4", "GZMK", "GZMA", "HLA-DRB1", "HLA-DPA1", "HLA-DRA", "CCL5", "COTL1", "TOX", "ITM2A", "IFNG", "CD2", "TIGIT",
  "PRDM1", "ICOS", "LAG3", "CCL3", "EOMES", "XCL2", "TNF", "ITGAE", "CCR5", "FASLG", "ARAP2", "CD70", "IL32", "CD44", "TANK", "PRF1",
  "SLAMF7", "CRTAM", "TNFRSF9", "TNFSF9", "CCR2", "CCR4", "CTLA4", "MAF")
  donors <- c("BT5", "BT1", "BT7", "BT4", "BT3", "BT8", "BT2", "BT18", "BT12", "BT11", "BT15", "BT17_brain", "BT13", "BT9", "BT10", "BT21", "BT25", "BT26", "BT23", "BT20", "BT19", "BT24", "BT27", "BT22")

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
  pdf("heatmaps/paired_dots_percentage.pdf", 15, 15)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", facet.by = "genes", palette = "jco", ylab = "Percentage")
  dev.off()

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

  library(tidyr)
  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor)) %>% group_by(orig.donor, PDCD1_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, PDCD1_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = PDCD1_tag, values_from = pct_exp_cells)

  write.csv(df, file = "heatmaps/paired_dots_PDCD1_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_PDCD1_expansion_pct.pdf", 15, 15)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # sc_cd3p_cd8_donor_mdata <- sc_cd3p_cd8@meta.data[intersect(PDCD1p_cells,donor_cells),]

  # df <- c()
  # # donor <- "BT5"
  # for (donor in donors){
  #   donor_cells <- rownames(sc_cd3p_cd8@meta.data[sc_cd3p_cd8@meta.data$orig.donor == donor & !is.na(sc_cd3p_cd8@meta.data$orig.donor),])
  #   # PDCD1p cells
  #   sc_cd3p_cd8_donor_mdata <- sc_cd3p_cd8@meta.data[intersect(PDCD1p_cells,donor_cells),]
  #
  #   # PDCD1p_exp <- apply(PDCD1p_exp, MARGIN=1, FUN=function(x){
  #   #   pCells <- 100*sum(x>0)/length(x)
  #   # })
  #   # PDCD1n cells
  #   PDCD1n_exp <- sc_cd3p_cd8@assays$RNA[genes,intersect(PDCD1n_cells,donor_cells)]
  #   PDCD1n_exp <- apply(PDCD1n_exp, MARGIN=1, FUN=function(x){
  #     pCells <- 100*sum(x>0)/length(x)
  #   })
  #   tmp <- data.frame(donor,PDCD1p_exp,PDCD1n_exp, genes)
  #   df <- rbind(df,tmp)
  # }
  # colnames(df) <- c("donor","PDCD1p","PDCD1n","genes")

  ########################
  # Clonotyper per donor #
  ########################

  # Get PDCD1 %pct of the clonotypes per donor
  sc_cd3p_cd4@meta.data$tag_PDCD1 <- add_gene_tag(c("PDCD1"), sc_cd3p_cd4@meta.data, sc_cd3p_cd4@assays$RNA@data, thresh = 0, tag = c('tag', 'p', 'n'))
  tmp <- sc_cd3p_cd4@meta.data %>% filter(!is.na(orig.donor), !is.na(clonotype.tag)) %>% group_by(orig.donor, clonotype.tag) %>% summarize(PDCD1_pct_byDonor = 100*sum(tag_PDCD1 == "PDCD1p")/n() ) # , .groups = c("orig.donor", "clonotype.tag")
  sc_cd3p_cd4@meta.data <- merge(sc_cd3p_cd4@meta.data, tmp, by = c("orig.donor", "clonotype.tag"), all.x = TRUE)

  sc_cd3p_cd8@meta.data$tag_PDCD1 <- add_gene_tag(c("PDCD1"), sc_cd3p_cd8@meta.data, sc_cd3p_cd8@assays$RNA@data, thresh = 0, tag = c('tag', 'p', 'n'))
  tmp <- sc_cd3p_cd8@meta.data %>% filter(!is.na(orig.donor), !is.na(clonotype.tag)) %>% group_by(orig.donor, clonotype.tag) %>% summarize(PDCD1_pct_byDonor = 100*sum(tag_PDCD1 == "PDCD1p")/n() ) # , .groups = c("orig.donor", "clonotype.tag")
  sc_cd3p_cd8@meta.data <- merge(sc_cd3p_cd8@meta.data, tmp, by = c("orig.donor", "clonotype.tag"), all.x = TRUE)

  output.dir <- "tcr"
  srt.objs.list <- list(sc_cd3p_cd4, sc_cd3p_cd8); names(srt.objs.list) <- c("CD4", "CD8")
  top.clones <- 50

  # --------> Proportions

  # obj <- sc_cd3p_cd4
  top.clonotypes <- lapply(X = srt.objs.list, FUN = function(obj){
    meta.data <- as.data.table(obj@meta.data)

    # Identify x top clonally expanded clonotypes.
    tmp.data.1 <- meta.data[!is.na(orig.donor) & !is.na(clonotype.tag), .(cell.count=.N, PDCD1_pct_byDonor=unique(PDCD1_pct_byDonor)), by=.(donor=orig.donor, clonotype=clonotype.tag)]
    setorderv(x=tmp.data.1, cols='cell.count', order=-1) # Maybe want to order by PDCD1_pct_byDonor as a 2nd variable
    uniq.donors <- tmp.data.1[, unique(donor)]
    tmp.top.clonotypes <- lapply(X=uniq.donors, FUN=function(tmp.donor){
      tmp.data <- tmp.data.1[donor==tmp.donor]
      tmp.data <- tmp.data[1:top.clones, .(donor, clonotype, cell.count, PDCD1_pct_byDonor)]
      tmp.data <- na.omit(tmp.data)
      # if(tmp.data[])
      return(tmp.data)
    })
    tmp.top.clonotypes <- rbindlist(l=tmp.top.clonotypes, use.names=TRUE, idcol=NULL)
    tmp.top.clonotypes[,donor:=factor(x = donor,
      levels = c("BT5", "BT1", "BT7", "BT4", "BT3", "BT8", "BT2", "BT18", "BT12", "BT11", "BT15", "BT17_brain", "BT13", "BT9", "BT10", "BT21", "BT25", "BT26", "BT23", "BT20", "BT19",
      "BT24", "BT27", "BT22"))]
    # @ Items for x top clonotypes.
    # Get proportions of final classes per donor.
    tmp.top.clonotypes[,proportion:=(cell.count*100)/sum(cell.count),by = "donor"]
  })

  for (i in 1:length(top.clonotypes)) {
    celltype <- str_extract(string = names(top.clonotypes[i]), pattern = "CD(4|8)")
    pdf(paste0("tcr/proportions_top50_mostExpandedClones_PBT_",celltype,".pdf"))
    tmp.y.axis <- "Proportion of T cell clones (%)"
    tmp.x.axis <- "Patients"
    tmp.title <- paste("Top 50 most-expanded",celltype,"clones in brain tumors")
    tmp.clones <- top.clonotypes[[i]]
    print(ggplot(tmp.clones, aes(x = donor, y = proportion, fill = PDCD1_pct_byDonor, group = -proportion)) + #, colour = PDCD1_pct_byDonor
      geom_bar(stat = "identity", color = "white") +     # scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
      labs(x = tmp.x.axis, y = tmp.y.axis, fill = "PDCD1(%)", title = tmp.title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()))
    dev.off()
    pdf(paste0("tcr/proportions_top50_mostExpandedClones_PBT_",celltype,"_noColor.pdf"))
    print(ggplot(tmp.clones, aes(x = donor, y = proportion, group = -proportion)) + #, colour = PDCD1_pct_byDonor
      geom_bar(stat = "identity", color = "white") +     # scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      labs(x = tmp.x.axis, y = tmp.y.axis, fill = "PDCD1(%)", title = tmp.title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()))
    dev.off()
  }

  # --------> Clone Size


  # obj <- sc_cd3p_cd4
  all.clonotypes <- lapply(X = srt.objs.list, FUN = function(obj){
    meta.data <- as.data.table(obj@meta.data)

    # Identify x top clonally expanded clonotypes.
    tmp.data.1 <- meta.data[!is.na(orig.donor) & !is.na(clonotype.tag), .(cell.count=.N, PDCD1_pct_byDonor=unique(PDCD1_pct_byDonor)), by=.(donor=orig.donor, clonotype=clonotype.tag)]
    setorderv(x=tmp.data.1, cols='cell.count', order=-1) # Maybe want to order by PDCD1_pct_byDonor as a 2nd variable
    uniq.donors <- tmp.data.1[, unique(donor)]
    # This function might not be necessary.
    tmp.all.clonotypes <- lapply(X=uniq.donors, FUN=function(tmp.donor){
      tmp.data <- tmp.data.1[donor==tmp.donor]
      tmp.data <- tmp.data[, .(donor, clonotype, cell.count, PDCD1_pct_byDonor)]
      tmp.data <- na.omit(tmp.data)
      # if(tmp.data[])
      return(tmp.data)
    })
    tmp.all.clonotypes <- rbindlist(l=tmp.all.clonotypes, use.names=TRUE, idcol=NULL)
    tmp.all.clonotypes[,donor:=factor(x = donor,
      levels = c("BT5", "BT1", "BT7", "BT4", "BT3", "BT8", "BT2", "BT18", "BT12", "BT11", "BT15", "BT17_brain", "BT13", "BT9", "BT10", "BT21", "BT25", "BT26", "BT23", "BT20", "BT19",
      "BT24", "BT27", "BT22"))]
  })

  # tmp.clones <- all.clonotypes[[1]]
  for (i in 1:length(all.clonotypes)) {
    celltype <- str_extract(string = names(all.clonotypes[i]), pattern = "CD(4|8)")
    pdf(paste0("tcr/cellCount_allClones_PBT_",celltype,".pdf"))
    tmp.y.axis <- "Number of cells"
    tmp.x.axis <- "Patients"
    tmp.title <- paste("",celltype,"clones in brain tumors")
    tmp.clones <- all.clonotypes[[i]]
    print(ggplot(tmp.clones, aes(x = donor, y = cell.count, fill = PDCD1_pct_byDonor, group = -cell.count)) + #, colour = PDCD1_pct_byDonor
      geom_bar(stat = "identity") +     # scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
      labs(x = tmp.x.axis, y = tmp.y.axis, fill = "PDCD1(%)", title = tmp.title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()))
    dev.off()
    pdf(paste0("tcr/cellCount_allClones_PBT_",celltype,"_noColor.pdf"))
    print(ggplot(tmp.clones, aes(x = donor, y = cell.count, color = "gray35", group = -cell.count)) + #, colour = PDCD1_pct_byDonor
      geom_bar(stat = "identity", color = "gray60") +     # scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      labs(x = tmp.x.axis, y = tmp.y.axis, title = tmp.title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()))
    dev.off()
  }


}

{ cat(redb("### Heatmaps ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("heatmaps")
  require(tidyverse)
  require(pheatmap)

  # CD8 Heatmap - Highly Expanded DEGs
  # ---> 1. Get Genes

  # Obtain DEGs from Highly Expanded cells (LG and HG separately).
  # Filter significant genes by p-value < 0.05 and log2FC > 0.25
  # Clonal expansion (cs>=5)
  # thold = ""
  # genes_LG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/LowGrade/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  # genes_HG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/HighGrade/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  # Clonal expansion (cs>=2)
  thold = "_thold2"
  genes_LG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/LowGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  genes_HG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD8_subclustering/HighGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)

  genes_LG_HG_shared <- intersect(genes_LG_all, genes_HG_all)
  genes_LG_unique <- setdiff(genes_LG_all, genes_LG_HG_shared)
  genes_HG_unique <- setdiff(genes_HG_all, genes_LG_HG_shared)
  genes <- c(genes_LG_HG_shared, genes_LG_unique, genes_HG_unique)

  # Iterate over the 4 combination of grade and expansion ("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")
  for (i in tmp){
    # ---> 2. Get Cells
    cells <- rownames(sc_cd3p_cd8@meta.data[!is.na(sc_cd3p_cd8@meta.data$expDegree) & sc_cd3p_cd8@meta.data$expDegree == i,])
    edata <- sc_cd3p_cd8@assays$RNA@data[genes, cells]
    # ---> 3. Get Mean/Pct
    # Mean
    mean <- rowMeans(as.matrix(edata), na.rm = T)
    mean_matrix <- cbind(mean_matrix, mean)
    # Pct
    pct <- apply(edata, MARGIN = 1, FUN = function(x){
      100*sum(x>0)/length(x)
    })
    pct_matrix <- cbind(pct_matrix, pct)
  }
  mean_matrix <- as.data.frame(round(mean_matrix, 2)); colnames(mean_matrix) <- tmp
  pct_matrix <- as.data.frame(round(pct_matrix, 2)); colnames(pct_matrix) <- tmp

  # ---> 3. Heatmap

  row_tags = data.frame("Group" = c(rep(c("shared-DEGs", "LG-DEGs", "HG-DEGs"), times=c(length(genes_LG_HG_shared), length(genes_LG_unique), length(genes_HG_unique))) )); rownames(row_tags) <- genes

  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_meanMatrix', thold,'.pdf'))
  pheatmap(mean_matrix,main = "Mean Expression HighExp DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         scale = "row",
         fontsize_number = 7)
  dev.off()

  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix,main = "% Expression HighExp DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        annotation_row = row_tags,
        scale = "row",
        fontsize_number = 7)
  dev.off()

  # CD4 (-TREGs) Heatmap - Highly Expanded DEGs
  # ---> 1. Get Genes

  # Obtain DEGs from Highly Expanded cells (LG and HG separately).
  # Filter significant genes by p-value < 0.05 and log2FC > 0.25.
  # Remove TREGs (cluster 5).
  # thold = ""
  # genes_LG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD4_subclustering/LowGrade/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  # genes_HG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD4_subclustering/HighGrade/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  thold = "_thold2"
  genes_LG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD4_subclustering/LowGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  genes_HG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size_CD4_subclustering/HighGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)

  genes_LG_HG_shared <- intersect(genes_LG_all, genes_HG_all)
  genes_LG_unique <- setdiff(genes_LG_all, genes_LG_HG_shared)
  genes_HG_unique <- setdiff(genes_HG_all, genes_LG_HG_shared)
  genes <- c(genes_LG_HG_shared, genes_LG_unique, genes_HG_unique)

  # Iterate over the 4 combination of grade and expansion ("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("LowExp_LG", "LowExp_HG", "HighExp_LG", "HighExp_HG")
  for (i in tmp){
    # ---> 2. Get Cells
    cells <- rownames(sc_cd3p_cd4@meta.data[!is.na(sc_cd3p_cd4@meta.data$expDegree) & sc_cd3p_cd4@meta.data$expDegree == i & sc_cd3p_cd4@meta.data$cluster %in% c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11"),])
    edata <- sc_cd3p_cd4@assays$RNA@data[genes, cells]
    # ---> 3. Get Mean/Pct
    # Mean
    mean <- rowMeans(as.matrix(edata))
    mean_matrix <- cbind(mean_matrix, mean)
    # Pct
    pct <- apply(edata, MARGIN = 1, FUN = function(x){
      100*sum(x>0)/length(x)
    })
    pct_matrix <- cbind(pct_matrix, pct)
  }
  mean_matrix <- as.data.frame(round(mean_matrix, 2)); colnames(mean_matrix) <- tmp
  pct_matrix <- as.data.frame(round(pct_matrix, 2)); colnames(pct_matrix) <- tmp

  # ---> 3. Heatmap

  row_tags = data.frame("Group" = c(rep(c("shared-DEGs", "LG-DEGs", "HG-DEGs"), times=c(length(genes_LG_HG_shared), length(genes_LG_unique), length(genes_HG_unique))) )); rownames(row_tags) <- genes

  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_meanMatrix', thold,'.pdf'))
  pheatmap(mean_matrix,main = "Mean Expression HighExp DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         scale = "row",
         fontsize_number = 7)
  dev.off()

  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix,main = "% Expression HighExp DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        annotation_row = row_tags,
        scale = "row",
        fontsize_number = 7)
  dev.off()


}
