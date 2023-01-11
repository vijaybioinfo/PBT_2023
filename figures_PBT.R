#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2023-01-10
# ---

### ==== CD4/CD8 clustering ==== ####

### ================== Figures ================== ###

# Done after locking the clustering objectss (CD4 and CD8)
# CD8 = pct20_pc20_res0.2
# CD4 = pct25_pc15_0.4

# Save Object
# source("/home/kmlanderos/pbtumor-all/scripts/object_lock_CD4_CD8_subclustering.R")
# source("/home/kmlanderos/pbtumor-all/scripts/figures/global_CD4_CD8.R")

# Read CD4 and CD8 Seurat Objects
# setwd('/home/kmlanderos/ad_hoc/pbtumor-all/results/figures/CD4_CD8_subclustering')
setwd('/home/kmlanderos/tmp_large/pbtumor-all/results/figures/CD4_CD8_subclustering')
sc_cd3p_cd8 <- readRDS(file = "data/sc_cd3p_cd8_seurat_object.rds")
sc_cd3p_cd4 <- readRDS(file = "data/sc_cd3p_cd4_seurat_object.rds")
system("ls -loh")

{ cat(redb("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # As you progress they appear in many tasks and so become more or less global

  redu = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2"))
  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
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

  sc_cd3p_cd4@meta.data <- sc_cd3p_cd4@meta.data %>% mutate(cluster.new = case_when(
    RNA_snn_res.0.4 == "1" & RNA_snn_res.0.8 %in% c("3","0","2","9","10") ~ "1.1", # CTL
    RNA_snn_res.0.4 == "1" & RNA_snn_res.0.8 %in% c("11","1","5","6","7","8") ~ "1.2", # TCM
    RNA_snn_res.0.4 == "1" & RNA_snn_res.0.8 == "4" ~ "5", # Assing to Treg (1 cell)
    RNA_snn_res.0.4 == "1" & RNA_snn_res.0.8 == "14" ~ "11", # Assign to Cell Clycle (1 cell)
    TRUE ~ as.character(RNA_snn_res.0.4)
  )) %>% mutate(cell_classification.new = case_when(
    cluster.new  %in% c("0","2","4","9","1.1") ~ "CTL",
    cluster.new  %in% c("3","6","7","8","1.2") ~ "TCM_TN",
    cluster.new  %in% c("11") ~ "Cell_Cycle",
    cluster.new  %in% c("5") ~ "TREG",
    cluster.new  %in% c("10") ~ "TFH"
    )
  )
  sc_cd3p_cd4@meta.data$cluster.new <- factor(sc_cd3p_cd4@meta.data$cluster.new)
  # table(sc_cd3p_cd4@meta.data$cell_classification.new, sc_cd3p_cd4@meta.data$cluster.new)
  # table(sc_cd3p_cd4@meta.data$RNA_snn_res.0.4, sc_cd3p_cd4@meta.data$RNA_snn_res.0.8)
  # table(sc_cd3p_cd4@meta.data$cell_classification, sc_cd3p_cd4@meta.data$cell_classification.new)
  # table(sc_cd3p_cd4@meta.data$cluster, sc_cd3p_cd4@meta.data$cluster.new)

  sc_cd3p_cd4@meta.data$cellname <- rownames(sc_cd3p_cd4@meta.data); cd4_cell_order <- rownames(sc_cd3p_cd4@meta.data)
  sc_cd3p_cd8@meta.data$cellname <- rownames(sc_cd3p_cd8@meta.data); cd8_cell_order <- rownames(sc_cd3p_cd8@meta.data)

  # Let us recalculate the clone size only consider cells with Gex.
  # CD8
  sc_cd3p_cd8@meta.data$clon.size.tag.old <- sc_cd3p_cd8@meta.data$clon.size.tag # Save old values
  clonotype.cs.cd8 <- as.data.frame(table(sc_cd3p_cd8@meta.data %>% pull(clonotype.tag)) ) %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd8) <- c("clonotype", "clone_size")
  sc_cd3p_cd8@meta.data$clon.size.tag <- unlist(apply(sc_cd3p_cd8@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"])) clonotype.cs.cd8[clonotype.cs.cd8$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  # sc_cd3p_cd8@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))

  # CD4
  sc_cd3p_cd4@meta.data$clon.size.tag.old <- sc_cd3p_cd4@meta.data$clon.size.tag # Save old values
  clonotype.cs.cd4 <- as.data.frame(table(sc_cd3p_cd4@meta.data %>% pull(clonotype.tag)) )  %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd4) <- c("clonotype", "clone_size")
  sc_cd3p_cd4@meta.data$clon.size.tag <- unlist(apply(sc_cd3p_cd4@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"])) clonotype.cs.cd4[clonotype.cs.cd4$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  # sc_cd3p_cd4@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))

  # Define the Low- and Highly expanded cells for LowGrade and HighGrade patients
   sc_cd3p_cd4@meta.data$expDegree <- sc_cd3p_cd4@meta.data$clon.size.tag
   sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree > 1] <- "Expanded"
   sc_cd3p_cd4@meta.data$expDegree[sc_cd3p_cd4@meta.data$expDegree == 1] <- "Non_expanded"
   sc_cd3p_cd4@meta.data$ExpDegree <- sc_cd3p_cd4@meta.data$clon.size.tag
   sc_cd3p_cd4@meta.data$ExpDegree[sc_cd3p_cd4@meta.data$ExpDegree > 1 & sc_cd3p_cd4@meta.data$grade == "LG"] <- "Exp_LG"
   sc_cd3p_cd4@meta.data$ExpDegree[sc_cd3p_cd4@meta.data$ExpDegree > 1 & sc_cd3p_cd4@meta.data$grade == "HG"] <- "Exp_HG"
   sc_cd3p_cd4@meta.data$ExpDegree[sc_cd3p_cd4@meta.data$ExpDegree == 1 & sc_cd3p_cd4@meta.data$grade == "LG"] <- "NonExp_LG"
   sc_cd3p_cd4@meta.data$ExpDegree[sc_cd3p_cd4@meta.data$ExpDegree == 1 & sc_cd3p_cd4@meta.data$grade == "HG"] <- "NonExp_HG"

   # Define the Low- and Highly expanded cells for LowGrade and HighGrade patients
   sc_cd3p_cd8@meta.data$expDegree <- sc_cd3p_cd8@meta.data$clon.size.tag
   sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree > 1] <- "Expanded"
   sc_cd3p_cd8@meta.data$expDegree[sc_cd3p_cd8@meta.data$expDegree == 1] <- "Non_expanded"
   sc_cd3p_cd8@meta.data$ExpDegree <- sc_cd3p_cd8@meta.data$clon.size.tag
   sc_cd3p_cd8@meta.data$ExpDegree[sc_cd3p_cd8@meta.data$ExpDegree > 1 & sc_cd3p_cd8@meta.data$grade == "LG"] <- "Exp_LG"
   sc_cd3p_cd8@meta.data$ExpDegree[sc_cd3p_cd8@meta.data$ExpDegree > 1 & sc_cd3p_cd8@meta.data$grade == "HG"] <- "Exp_HG"
   sc_cd3p_cd8@meta.data$ExpDegree[sc_cd3p_cd8@meta.data$ExpDegree == 1 & sc_cd3p_cd8@meta.data$grade == "LG"] <- "NonExp_LG"
   sc_cd3p_cd8@meta.data$ExpDegree[sc_cd3p_cd8@meta.data$ExpDegree == 1 & sc_cd3p_cd8@meta.data$grade == "HG"] <- "NonExp_HG"

  # --- Add specificity tag

  # CD4
  GLIPH2Input_CD4 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-07-15/GLIPH2Input_TCR-Data_All-TRB_CD4_hto_wAlphaChain.tsv", sep = "\t") # 2022-02-07
  GLIPH2Input_CD4 <- GLIPH2Input_CD4 %>% mutate(id = paste(GLIPH2Input_CD4$cdr3b.aa.se, GLIPH2Input_CD4$trb.v, GLIPH2Input_CD4$trb.j, sep = "_"))

  GLIPH2Output_CD4 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/2022-07-15/GLIPH2Output_TCR-Data_CD4_hto.csv") # 2022-02-07
  GLIPH2Output_CD4 <- GLIPH2Output_CD4 %>% mutate(id = paste(GLIPH2Output_CD4$TcRb, GLIPH2Output_CD4$V, GLIPH2Output_CD4$J, sep = "_"))
  GLIPH2Output_CD4 <- GLIPH2Output_CD4[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  merge <- merge(data.frame(GLIPH2Output_CD4), GLIPH2Input_CD4[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups <- unlist(str_split(merge$clonotype.tag, ";"))
  sc_cd3p_cd4@meta.data <- sc_cd3p_cd4@meta.data %>% mutate(sGroup_tag = case_when(
    clonotype.tag %in% cts_sGroups ~ TRUE,
    TRUE ~ FALSE
  ))
  # table(sc_cd3p_cd4@meta.data$sGroup_tag)

  # Add GLIPH2 pattern and n_sG
  tmp <- merge %>% mutate(clonotype.tag = as.character(clonotype.tag), n_sG = str_count(merge$pattern, "_")+1) %>% select(clonotype.tag, pattern, n_sG)
  tmp.2 <- c()
  for (i in 1:dim(tmp)[1]) {
    x <- tmp[i,]
    clonotypes <- unlist(str_split(x["clonotype.tag"], ";"))
    n_clonotypes <- length(clonotypes)
    for(clonotype in clonotypes){ tmp.2 <- rbind(tmp.2, c(clonotype, as.character(x[-1]))) }
  }; tmp.2 <- as.data.frame(tmp.2); colnames(tmp.2) <- colnames(tmp)
  tmp.meta.data <- merge(sc_cd3p_cd4@meta.data, tmp.2, by = "clonotype.tag", all.x = TRUE); rownames(tmp.meta.data) <- tmp.meta.data$cellname
  sc_cd3p_cd4@meta.data <- tmp.meta.data[cd4_cell_order, ]

  # CD8
  GLIPH2Input_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-07-15/GLIPH2Input_TCR-Data_All-TRB_CD8_hto_wAlphaChain.tsv", sep = "\t") # 2022-02-07
  GLIPH2Input_CD8 <- GLIPH2Input_CD8 %>% mutate(id = paste(GLIPH2Input_CD8$cdr3b.aa.se, GLIPH2Input_CD8$trb.v, GLIPH2Input_CD8$trb.j, sep = "_"))

  GLIPH2Output_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/2022-07-15/GLIPH2Output_TCR-Data_CD8_hto.csv") # 2022-02-07
  GLIPH2Output_CD8 <- GLIPH2Output_CD8 %>% mutate(id = paste(GLIPH2Output_CD8$TcRb, GLIPH2Output_CD8$V, GLIPH2Output_CD8$J, sep = "_"))
  GLIPH2Output_CD8 <- GLIPH2Output_CD8[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  # Add sGroup_tag
  merge <- merge(data.frame(GLIPH2Output_CD8), GLIPH2Input_CD8[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups <- unlist(str_split(merge$clonotype.tag, ";"))
  sc_cd3p_cd8@meta.data <- sc_cd3p_cd8@meta.data %>% mutate(sGroup_tag = case_when(
    clonotype.tag %in% cts_sGroups ~ TRUE,
    TRUE ~ FALSE
  ))# table(sc_cd3p_cd8@meta.data$sGroup_tag)

  # Add GLIPH2 pattern and n_sG
  tmp <- merge %>% mutate(clonotype.tag = as.character(clonotype.tag), n_sG = str_count(merge$pattern, "_")+1) %>% select(clonotype.tag, pattern, n_sG)
  tmp.2 <- c()
  for (i in 1:dim(tmp)[1]) {
    x <- tmp[i,]
    clonotypes <- unlist(str_split(x["clonotype.tag"], ";"))
    n_clonotypes <- length(clonotypes)
    for(clonotype in clonotypes){ tmp.2 <- rbind(tmp.2, c(clonotype, as.character(x[-1]))) }
  }; tmp.2 <- as.data.frame(tmp.2); colnames(tmp.2) <- colnames(tmp)
  tmp.meta.data <- merge(sc_cd3p_cd8@meta.data, tmp.2, by = "clonotype.tag", all.x = TRUE); rownames(tmp.meta.data) <- tmp.meta.data$cellname
  sc_cd3p_cd8@meta.data <- tmp.meta.data[cd8_cell_order, ]

  # Adding TRV TRJ to the metadata

  # CD8
  cells.clons.info <- read.csv(file='/home/kmlanderos/large/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd8/filtered_contig_annotations_aggr.csv', stringsAsFactors=FALSE)
  cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
  feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
  cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]

  cells.clons.info <- as.data.table(cells.clons.info)
  vdj.gene.info <- cells.clons.info[
    ,
    .(
      v.gene=.SD[,
                 .(count=.N),
                 by=v_gene
      ][count==max(count), paste0(unique(v_gene), collapse='|')],
      j.gene=.SD[,
                 .(count=.N),
                 by=j_gene
      ][count==max(count), paste0(unique(j_gene), collapse='|')]
    ),
    by=.(
      raw_clonotype_id,
      chain
    )
  ]
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=tibble(vdj.gene.info), col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=tibble(vdj.gene.info), col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)

  # Merge VDJ data with object metadata
  sc_cd3p_cd8@meta.data <- merge(sc_cd3p_cd8@meta.data, vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id")
  rownames(sc_cd3p_cd8@meta.data) <- sc_cd3p_cd8@meta.data$cellname

  # CD4
  cells.clons.info <- read.csv(file='/home/kmlanderos/large/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd4/filtered_contig_annotations_aggr.csv', stringsAsFactors=FALSE)
  cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
  feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
  cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]

  cells.clons.info <- as.data.table(cells.clons.info)
  vdj.gene.info <- cells.clons.info[
    ,
    .(
      v.gene=.SD[,
                 .(count=.N),
                 by=v_gene
      ][count==max(count), paste0(unique(v_gene), collapse='|')],
      j.gene=.SD[,
                 .(count=.N),
                 by=j_gene
      ][count==max(count), paste0(unique(j_gene), collapse='|')]
    ),
    by=.(
      raw_clonotype_id,
      chain
    )
  ]
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=tibble(vdj.gene.info), col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=tibble(vdj.gene.info), col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)

  # Merge VDJ data with object metadata
  sc_cd3p_cd4@meta.data <- merge(sc_cd3p_cd4@meta.data, vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id")
  rownames(sc_cd3p_cd4@meta.data) <- sc_cd3p_cd4@meta.data$cellname

  # Save the metadata files and objects.
  # Just run one time
  # dir.create("data")
  # saveRDS(sc_cd3p_cd8@meta.data, file = "data/sc_cd3p_cd8_mdata.rds")
  # saveRDS(sc_cd3p_cd4@meta.data, file = "data/sc_cd3p_cd4_mdata.rds")
  # saveRDS(sc_cd3p_cd8, file = "data/sc_cd3p_cd8_seurat_object.rds")
  # saveRDS(sc_cd3p_cd4, file = "data/sc_cd3p_cd4_seurat_object.rds")

}

{ cat(redb("### Signature Analyses ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # cp /Users/kmlanderos/Documents/PBT/Data/signatures_cd3p*csv /Volumes/kmlanderos/pbtumor-all/info/
  # Rscript /home/ciro/scripts/functions/csvCorrect.R /home/kmlanderos/pbtumor-all/info/signatures_cd3p.csv

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

  # colnames(slists)[grep("ctl",colnames(slists), ignore.case = T)]

  fconfigs = list(
    ## ------ CD8
    list(result_id = "gsea/CD8/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(1:22,40)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("cell_classification", "GZMK", "TRM", "2", "3", "4", "5", "9"),
      comparisons = c("cell_classification")
    ),
    ## CD8 PDCD1p vs PDCD1n
    list(result_id = "gsea/CD8_PDCD1/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c("CellCycle_Best2013", "CD4CTLvsTCM_Patil2018")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("PDCD1_tag", "PDCD1n", "PDCD1p"),
      comparisons = c("PDCD1_tag")
    ),
    ## Neo TCR - CD8 all; expansion comparison
    list(result_id = "gsea/CD8_neoantigenTCR_all_expansion/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
    ),
    ## TSCM - CD8 all; cluster comparison
    list(result_id = "gsea/CD8_TSCM_all_cluster/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(43,44)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster", "3", "4", "REST")
    ),
    ## -------- CD4
    list(result_id = "gsea/CD4_Newclustering_TCM/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(23:26,28:32,37:39)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cell_classification.new")
    ),
    ## Neo TCR - CD4 all; expansion comparison
    list(result_id = "gsea/CD4_neoantigenTCR_all_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
    ),
    ## CD4_nonTREG_Expansion
    list(result_id = "gsea/CD4_nonTREG_Expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "th1_signature1_arlehamn", "th2_signature0_arlehamn", "TH17_Seumois2020")], FUN = function(x) x[-c(1:2)] ), # 46 is the one that matters
      sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("expDegree")
    ),
    ## Identify CTL and TCM/TN
    list(result_id = "gsea/CD4_clust1_subclustering/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Hu.CD4.Naive.CM")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("RNA_snn_res.deconv", "1.0", "1.1", "1.2"),
      comparisons = c("RNA_snn_res.deconv")
    ),
    ## CTL
    list(result_id = "gsea/CD4_Newclustering_CTL/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Hu.CD4.Naive.CM")], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cell_classification.new", "CTL", "REST")
    ),
    ## TCM
    list(result_id = "gsea/CD4_Newclustering_TCM/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Hu.CD4.Naive.CM")], FUN = function(x) x[-c(1:2)] ),
      # sample_filter = c("RNA_snn_res.deconv", "1.0", "1.1", "1.2"),
      comparisons = c("cell_classification.new", "TCM_TN", "REST")
    )
  )

  # CD8
  pp_gsea = fig_gsea(fconfigs[c(1:4)], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # CD4
  pp_gsea = fig_gsea(fconfigs[c(5:10)], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)

  ## ---- Module scoring analysis ---- ###

  source("/home/ciro/scripts/handy_functions/R/gsea_signature.R") # clean_feature_list, signature_scoring
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  source("/home/kmlanderos/scripts/handy_functions/R/moduleScore_functions.R") # signature_scoring (updated)
  dir.create("modulescores");

  # sc_cd3p_lists <- readfile("/home/ciro/pbtumor/large/results/figures/gsea/sc_cd3p_mod_signatures.csv", stringsAsFactors = FALSE)
  sc_cd3p_lists <- readfile("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/gsea/signatures_cd3p_human_v2.csv", stringsAsFactors = FALSE)

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

  fconfigs = list(
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4", object = "sc_cd3p_cd4",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2),c("CD4CTLvsTCM_Patil2018", "TFH_Locci2013", "Treg_Schmiedel2018", "TH17_Seumois2020", "tcell_activation_GO.0042110.1", "tcr_signaling_RHSA202403.1", "CellCycle_Best2013.1", "TypeIandII_IFNsignaling_Seumois2020.1", "th1_signature1_arlehamn", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")]),
      axis_x = list(col = "cell_classification_2" , order=c("TCM", "CTL", "1", "10", "11", "5")),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    ),
    list(result_id = "modulescores/", sufix = "_test/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "sc_cd3p_cd8@meta.data", object = "sc_cd3p_cd8",
      lists = gsea_process_list(sc_cd3p_lists[-c(1:2), c("CD8_GZMK_Guo2018", "TRM_192_UP_HobitScience", "Gutierrez_Innateness2019", "CYTOTOXICITY.SIGNATURE_CD161.paper_matthewson..Cell")]),
      axis_x = list(col = "cell_classification_2"),
      # sample_filter = c("cluster", "0", "3", "9.8", "15", "14", "16", "4"),
      variables = "Identity"
    )
  )
  pp_ms = fig_module_score(fconfigs[1], reductions = redu[1], violins_color = 'mean', couls = c("steelblue1", "yellow", "red",  "#670000")) # calls signature_scoring
  pp_ms = fig_module_score(fconfigs[2], reductions = redu[1], violins_color = 'mean', couls = c("steelblue1", "yellow", "red",  "#670000")) # calls signature_scoring

  # Adjust scale of the CD4-CTL signature (0-0.2)
  signature_scores <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/modulescores/sc_cd3p_cd4_test1/signatures.csv", row.names=1)
  if( all(rownames(mdata_sc_cd3p_cd4) == rownames(signature_scores)) ) {
    mdata_sc_cd3p_cd4$CD4CTL_signature <- signature_scores[,"CD4CTLvsTCM_Patil2018"]
  }

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
  # Vicente's color scale 0.15
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.15_umap_v3.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()


}

{ cat(redb("### V and J genes usage ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  # ==================== CD8

  # ---> Beta Chain
  # trb.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.v)) %>% select(cell_classification, trb.v) %>% group_by(cell_classification, trb.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(trb.j)) %>% select(cell_classification, trb.j) %>% group_by(cell_classification, trb.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.v)) %>% select(cell_classification, tra.v) %>% group_by(cell_classification, tra.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRV_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd8@meta.data %>% filter(!is.na(tra.j)) %>% select(cell_classification, tra.j) %>% group_by(cell_classification, tra.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd8/TRJ_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)

  # ==================== CD4

  # ---> Beta Chain
  # trb.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.v)) %>% select(cell_classification, trb.v) %>% group_by(cell_classification, trb.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # trb.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(trb.j)) %>% select(cell_classification, trb.j) %>% group_by(cell_classification, trb.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRJ_beta_usage_perCluster.tsv", sep = "\t", row.names = F)
  # ---> Alpha Chain
  # tra.v
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.v)) %>% select(cell_classification, tra.v) %>% group_by(cell_classification, tra.v) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
  write.table(df, "tcr/sc_cd3p_cd4/TRV_alpha_usage_perCluster.tsv", sep = "\t", row.names = F)
  # tra.j
  df <- sc_cd3p_cd4@meta.data %>% filter(!is.na(tra.j)) %>% select(cell_classification, tra.j) %>% group_by(cell_classification, tra.j) %>% summarize(n_cells = n()) %>% arrange(cell_classification, desc(n_cells))
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

{ cat(redb("### Markers: dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dim_reduction_features", showWarnings = FALSE)

  fconfigs = list(
    list(result_id = "dim_reduction_features/sc_cd3p_cd4_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "PRF1", "GZMA", "IFNG", "TBX21", "EOMES", "GZMB", "SLAMF7", "CRTAM", "TNF", "GZMH", "GNLY", "CCL4", "CCL5",
      "CTSW", "ZNF683", "FCRL6", "SPON2", "CX3CR1", "S1PR5", "NKG7", "CD244", "CXCR3", "S1PR1", "SELL", "TCF7", "FOXP3", "RORA", "RORC",
      "CCR6", "KLRB1", "IRF7", "CXCL13", "BTLA", "CD40LG", "ITGAE", "CD69"),
      col = c("lightgrey", "blue")
    ),
    list(result_id = "dim_reduction_features/sc_cd3p_cd8_", sufix = "exm_markers/",
      edata = "sc_cd3p_cd8@assays$RNA@data",
      metadata = "sc_cd3p_cd8@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("PDCD1", "GZMK", "GZMA", "PRF1", "LAG3", "IFNG", "TNF", "CCL4", "CCL3", "FASLG", "HLA-DRB1", "TOX"),
      col = c("lightgrey", "blue")
    )

  pp_markers = fig_plot_scatters(fconfigs)

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

  # ------------------- Inset plot
  # ---> CD4
  # PDCD1
  pdf("dim_reduction_features/sc_cd3p_cd4_exm_markers/PDCD1_inset_plot.pdf")
  mdata <- sc_cd3p_cd4@meta.data; edata <- sc_cd3p_cd4@assays$RNA
  tmp <- as.vector(as.matrix(edata["PDCD1",rownames(mdata)]))
  mdata <- cbind(mdata, PDCD1_expr = tmp)
  mdata$tag_PDCD1 <- add_gene_tag(c("PDCD1"), mdata, edata@data, thresh = 0, tag = c('tag', 'p', 'n'))
  a <- mdata %>% group_by(cell_classification.new) %>% summarize(mean_expr = mean(PDCD1_expr),
    pct.expr = 100*sum(tag_PDCD1 == "PDCD1p")/n())
  a$category <- "Expansion"; a$cell_classification.new <- factor(a$cell_classification.new, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM_TN"))
  p <- a %>% ggplot(aes(x = cell_classification.new, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30, 45)) +
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
  a <- mdata %>% group_by(cell_classification.new) %>% summarize(mean_expr = mean(LAG3_expr),
    pct.expr = 100*sum(tag_LAG3 == "LAG3p")/n())
  a$category <- "Expansion"; a$cell_classification.new <- factor(a$cell_classification.new, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM_TN"))
  p <- a %>% ggplot(aes(x = cell_classification.new, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 5, 15, 25)) +
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
  a <- mdata %>% group_by(cell_classification) %>% summarize(mean_expr = mean(PDCD1_expr),
    pct.expr = 100*sum(tag_PDCD1 == "PDCD1p")/n())
  a$category <- "Expansion"; a$cell_classification <- factor(a$cell_classification, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification, y = category, color = mean_expr)) +
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
  a <- mdata %>% group_by(cell_classification) %>% summarize(mean_expr = mean(LAG3_expr),
    pct.expr = 100*sum(tag_LAG3 == "LAG3p")/n())
  a$category <- "Expansion"; a$cell_classification <- factor(a$cell_classification, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
  p <- a %>% ggplot(aes(x = cell_classification, y = category, color = mean_expr)) +
    geom_point(aes(size = pct.expr), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(1, 10, 20, 30, 40, 50)) +
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
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD8/CD8_Expansion/ExpandedvsNon_expanded/mastlog2cpm_results.csv",
      group1 = "Expanded", group2 = "Non_expanded"
    ),
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD8/CD8_PDCD1_tag/PDCD1pvsPDCD1n/mastlog2cpm_results.csv",
      group1 = "PDCD1p", group2 = "PDCD1n",
      showgenes = c("GZMA", "GZMK", "GZMH", "CCL4", "CCL5", "CD74", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "RGS1", "CST7", "NKG7", "IFNG", "CXCR6", "TIGIT", "PRDM1", "ICOS", "LAG3", "CCL3", "CXCR3",
        "TNFSF9", "EOMES", "ITGAE", "XCL2", "TNF", "ITGA1", "CCR5", "FASLG", "NFATC2", "RBPJ", "CD44", "PRF1", "BATF", "SLAMF7", "FAS", "CXCR4", "IL7R", "KLF2", "TCF7", "CD55", "CD7", "TMEM123", "FLT3LG", "LEF1",
        "SELL", "FAM65B", "CCR7")
    ),
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD4/CD4_nonTREG_Expansion/ExpandedvsNon_expanded/mastlog2cpm_results.csv",
      group1 = "Expanded", group2 = "Non_expanded",
      showgenes = c("CCL4","CCL5","GZMA","GZMH","GZMK","IFNG","TBX21","PDCD1","TOX","LAG3","RGS1","CXCR6","RBPJ","ITGAE","PRDM1","RUNX2","HOPX","BATF","RUNX3","HLA-DRB1","NKG7","THEMIS","TCF7","LEF1","IL7R","SELL","KLF2","CD7","CD55")
      # showgenes = c("GZMA","GZMH","GZMK","PRDM1","HLA-DRB1","IFNG","CCL4","CCL5","IFNG","CXCR6","PDCD1","RGS1","CST7","CD2","CXCR3","CD74","KLRB1","IL2RB","RBPJ","CTLA4","CD40LG","CXCL13","BATF","KLF2","SELL","FAM65B","CCR7","TCF7","MAL","LEF1","CD55","CD7","SOCS3","IL7R","LTB","FOS","JUNB","IER2")
    ),
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD4/CD4_nonTREG_PDCD1_tag/PDCD1pvsPDCD1n/mastlog2cpm_results.csv",
      group1 = "PDCD1p", group2 = "PDCD1n",
      showgenes = c("GZMA", "CCL4", "CCL5", "GZMH", "HLA-DRB1", "TOX", "CXCR6", "PRDM1", "IFNG", "ITGAE", "RBPJ", "PRF1", "ZEB2", "ICOS")
    ),
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/dgea/CD4_CD8_subclustering/CD8/CD8_LAG3_tag/LAG3pvsLAG3n/mastlog2cpm_results.csv",
      group1 = "LAG3p", group2 = "LAG3n",
      showgenes = c("GZMA", "CCL4", "CCL5", "GZMH", "HLA-DRB1", "TOX", "CXCR6", "PRDM1", "IFNG", "ITGAE", "RBPJ", "PRF1", "ZEB2", "ICOS", "GZMB")
    )
  )

  # fconfig=fconfigs[[5]]
  for (fconfig in fconfigs[5]){
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
    pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots_blank.pdf"), 10, 10)
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
      )+ xlim(-3,3) )
      # + xlim(-5,5)
      # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
      #   title = paste(group2, "(-) vs ", group1, "(+)"))
    dev.off()

    # With names (big dots)
    pdf(paste0(dirname(file),"/mastlog2cpm_volcano_biggerDots.pdf"), 10, 10)
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
      ))
       # + xlim(-5,5)
      # + labs(size = "Delta %", color = paste0("Mean (", dtype, ")"),
      #   title = paste(group2, "(-) vs ", group1, "(+)"))
    dev.off()

  }

}

{ cat(redb("### Main dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  dir.create("dim_reduction")
  source("/home/ciro/scripts/clustering/R/plotting.R")

  sc_cd3p_cd4@meta.data$orig.Diagnosis <- gsub(" ", "_", sc_cd3p_cd4@meta.data$orig.Diagnosis)
  sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass <- gsub(" ", "_", sc_cd3p_cd4@meta.data$orig.Diagnosis.subclass)
  sc_cd3p_cd8@meta.data$orig.Diagnosis <- gsub(" ", "_", sc_cd3p_cd8@meta.data$orig.Diagnosis)
  sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass <- gsub(" ", "_", sc_cd3p_cd8@meta.data$orig.Diagnosis.subclass)

  fconfigs = list(
    # CD4
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster",
      colors = c("0" = "#c5b8ff", "1" = "#5fcfab", "7" = "#e57f08", "6" = "#ffaf4d", "8" = "#b35806", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#fcd093",
       "9" = "#5d1d88", "10" = "#ff717c", "11" = "#9c0006"),
      redu = redu[[1]]
    ),
    # CD8
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster.new",
      colors = c("0" = "#c5b8ff", "1.1" = "#d26dec", "1.2" = "#52edf4", "7" = "#379ee4", "6" = "#97ccf5", "8" = "#1069b9", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#c0f1f5",
       "9" = "#5d1d88", "10" = "#ff717c", "11" = "#9c0006"),
      redu = redu[[1]]
    )
  )

  pp_main_dr = fig_plot_base(
    fconfigs,
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

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data
  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

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
      # col = c("GZMK" = "#54D8FD", "TRM" = "#FFB94F", "INNATE-LIKE_NKR" = "#D655A7", "INNATE-LIKE_XCL1-2" = "#44A306", "TCM" = "#969696", "MAIT" = "#C994C7", "Cell_Cycle" = "#464E2E"),
      col = c("Cell_Cycle" = "#9c0006", "GZMK_HI" = "#72e5d2", "TRM" = "#ffa7d7", "CD16p_Effector" = "#fcee62", "Effector" = "#ff971f", "MAIT" = "#9e89e8", "TCF7_HI" = "#4da7fd"),
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
  pp_overlaps = fig_plot_overlaps(fconfigs[1])

  if(!exists("gcolours")) gcolours = NULL
  fconfigs = list(
    list(result_id = "tcr/sc_cd3p_cd8/", edata = "sc_cd3p_cd8@assays$RNA@data",
      metadata = "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag), ]",
      axis_x = "cluster", axis_y = "clonotype.tag", # filters = c("celltype", "CD8"),
      col = c("Cell_Cycle" = "#9c0006", "GZMK_HI" = "#3bb8ae", "TRM" = "#ff6ede", "CD16p_Effector" = "#fcee62", "Effector" = "#ff971f", "MAIT" = "#9e89e8", "TCF7_HI" = "#4da7fd"),
      vars2plot = c("cell_classification"),
      element_type = c("clones", "cells")),
    list(result_id = "tcr/sc_cd3p_cd4/", edata = "sc_cd3p_cd4@assays$RNA@data",
      metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag), ]",
      axis_x = "cluster", axis_y = "clonotype.tag", # filters = c("celltype", "CD4"),
      col = c("TCM" = "#969696", "CTL1" = "#66C2A4", "CTL_TH17" = "#B4A7D6", "CTL_TH1" = "#B25785", "CTL_IFN" = "#703352", "TNF" = "#EA9999", "TREG" = "#FFE599", "TFH" = "#DFB9CD", "Cell_Cycle" = "#464E2E"),
      vars2plot = c("cell_classification"),
      element_type = c("clones", "cells"))
  )
  dfplot <- mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag), ]
  verbose=T
  for(i in 1:length(fconfigs[1])){ # fconfig=fconfigs[[1]]; i=1
      # Setting data
      fconfig = fig_config_check(fconfigs[[i]])
      fig_data = fig_set_data(fconfig,
        return_pdata = TRUE, verbose = verbose)
      ddf = fig_data$pdata
      axes = c(fconfig$axis_x[[1]][[1]], fconfig$axis_y[[1]][[1]])
      # feats <- fconfig$features[fconfig$features %in% fig_data$features]
    for(i in c(fconfig$vars2plot)){ # vars2plot_i=fconfig$vars2plot[1]; i="cell_classification"
      if(verbose) cat("-", i, "\n")
      fname0 <- paste0(fig_data$name, gsub("orig\\.", "", i))
      if(!is.numeric(ddf[, i])) ddf[, i] <- droplevels(factor(ddf[, i]))
      elist = base::split(x = ddf[[axes[2]]], f = ddf[[i]])
      cols_i = v2cols(names(elist), fconfig$col)
      if(is.null(fconfig$element_type)) fconfig$element_type = c("collapsed", "complete")

      for(type in fconfig$element_type){ # type = fconfig$element_type[2]
        if(verbose) cat(" *", type, "\n")
        elist_i = if(grepl("cells|complete", type)){ # clones size as the overlap
          lapply(X = setNames(nm = names(elist)), FUN = function(x){
            cells_with_clones <- ddf[[axes[2]]] %in% elist[[x]]
            rownames(ddf)[which(cells_with_clones)]
          })
        }else{ elist }

        fname <- paste0(fconfig$result_id, gsub("orig|\\.", "", i), "_")
        cat("* Upset\n")
        uddf <- as.data.frame(ComplexHeatmap::list_to_matrix(elist_i))
        pdf(paste0(fname, "upset.pdf"), onefile = FALSE, width = 10)
        print(UpSetR::upset(data = uddf, sets.bar.color = rev(fconfig$col), sets = rev(c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI")), keep.order = TRUE))
        graphics.off()
      }
    }
  }


  # pp_overlaps = fig_plot_overlaps(fconfigs, nintersects = 70, empty.intersections = "on") Complete plot
  pp_overlaps = fig_plot_overlaps(fconfigs)
  pp_overlaps = fig_plot_overlaps(lapply(fconfigs, function(x) {x$vars2col <- c("cluster"); return(x)})) # Color by cluster

  ## Scatter plot (umap) for cluster clone size ## ----------------------------------
  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$sufix <- "cluster_clone_size_v2_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification"; fconfigs[[i]]$vars2col = "cluster"
  }
  fconfigs[[2]]$vars2col = "cluster.new"

  # CD8 colors
  fconfigs[[1]]$col <- c("0" = "#72e5d2", "1" = "#72e5d2", "7" = "#72e5d2", "6" = "#ffa7d7", "8" = "#ffa7d7", "2" = "#fcee62", "5" = "#9e89e8", "4" = "#4da7fd", "3" = "#ff971f",
   "9" = "#9c0006")
  # CD4 colors
  fconfigs[[2]]$col <- c("0" = "#a084ff", "1.1" = "#a084ff", "1.2" = "#97ccf5", "7" = "#97ccf5", "6" = "#97ccf5", "8" = "#97ccf5", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#a084ff", "3" = "#97ccf5",
   "9" = "#a084ff", "11" = "#9c0006", "10" = "#ff717c")

  # Create the "cluster_clones_size" column in the metadata of the 2 objects.
  clonotype_cd4 <- sc_cd3p_cd4@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster.new) %>% group_by(clonotype.tag, cluster.new) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster.new) %>%
    pivot_wider(names_from = cluster.new, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd4) <- clonotype_cd4$clonotype.tag; clonotype_cd4$clonotype.tag <- NULL

  clonotype_cd8 <- sc_cd3p_cd8@meta.data %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd8) <- clonotype_cd8$clonotype.tag; clonotype_cd8$clonotype.tag <- NULL

  # Update CD4 mdata
  cluster_clone_size <- c()
  for (i in 1:dim(mdata_sc_cd3p_cd4)[1]){
    clonotype <- mdata_sc_cd3p_cd4[i,"clonotype.tag"]
    c <- as.character(mdata_sc_cd3p_cd4[i,"cluster.new"])
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
    geom_point(aes(size = pct.exp), alpha = 0.8, show.legend = TRUE) +
    scale_radius(breaks = c(10, 30, 50, 70)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
    scale_colour_manual(values = c("Cell_Cycle" = "#9c0006", "GZMK_HI" = "#3bb8ae", "TRM" = "#ff6ede", "CD16p_Effector" = "#fcee62", "Effector" = "#ff971f", "MAIT" = "#9e89e8", "TCF7_HI" = "#4da7fd"))
  print(p)
  dev.off()
  # --- CD4
  pdf(paste0(fconfigs[[2]]$result_id, "inset_plot_clone_size.pdf"))
  a <- sc_cd3p_cd4@meta.data %>% select(clonotype.tag, cell_classification.new, expDegree) %>% filter(!is.na(clonotype.tag)) %>% group_by(cell_classification.new) %>%
    summarize(pct.exp = 100*sum(expDegree == "Expanded")/length(expDegree)) %>% data.frame()
  a$category <- "Expansion"; a$cell_classification.new <- factor(a$cell_classification.new, levels = c("Cell_Cycle", "CTL", "TFH", "TREG", "TCM_TN"))
  p <- a %>% ggplot(aes(x = cell_classification.new, y = category, color = cell_classification.new)) +
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

  # --------> Clone Size

  all.clonotypes <- lapply(X = srt.objs.list, FUN = function(obj){
    meta.data <- as.data.table(obj@meta.data)

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
      levels = c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5", "BT24", "BT25", "BT10", "BT15", "BT22", "BT11", "BT12", "BT18", "BT21", "BT26", "BT2", "BT13",
       "BT17_brain", "BT20", "BT23"))]
  })

  # Save CD4 and CD8 XL tables
  write.csv(all.clonotypes[["CD4"]], file = "tcr/cellCount_allClones_PBT_CD4.csv", row.names = FALSE)
  write.csv(all.clonotypes[["CD8"]], file = "tcr/cellCount_allClones_PBT_CD8.csv", row.names = FALSE)

  ################ - Collapsing clonotypes with cs == 1

  # df <- all.clonotypes[["CD4"]]
  all.clonotypes.1 <- lapply(X = all.clonotypes, FUN = function(df){
    df.new <- df %>% filter(cell.count != 1)
    df.tmp <- df %>% filter(cell.count == 1) %>% group_by(donor) %>% summarize(
      cell.count = sum(cell.count), PDCD1_pct_byDonor = sum(PDCD1_pct_byDonor)/n()
    )
    df.tmp$clonotype <- "clonotypes_nonExp"; df.tmp <- df.tmp[,c("donor", "clonotype", "cell.count", "PDCD1_pct_byDonor")]
    df.final <- rbind(df.new, df.tmp)
    df.final
  })

  # Plot

  donors <- c("BT5", "BT25", "BT24", "BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT26", "BT21", "BT11", "BT18", "BT12", "BT22", "BT15", "BT10", "BT23", "BT13", "BT20", "BT17_brain", "BT2")
  cd4.matrix$donor <- factor(cd4.matrix$donor, levels = rev(donors))
  cd8.matrix$donor <- factor(cd8.matrix$donor, levels = rev(donors))

  pdf(paste0("tcr/cellCount_allClones_PBT_CD4_CD8_nonExpCollapsed_v5.pdf"))
  tmp.y.axis <- "CD4 cells"
  tmp.x.axis <- "Patient"
  n_non_Exp <- sum(grepl("clonotypes_nonExp", cd4.matrix$donor_clonotype))
  p1 <- ggplot(cd4.matrix, aes(x = donor, fill = donor_clonotype)) + geom_bar(position = "stack", color = "cadetblue4", size = 0.0001) + # gray60
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_reverse(limits = c(2000,0)) + coord_flip() +
    scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"), legend.position="none")

  tmp.y.axis <- "CD8 cells"
  p2 <- ggplot(cd8.matrix, aes(x = donor, fill = donor_clonotype)) + geom_bar(position = "stack", color = "royalblue4", size = 0.0001) + # gray60
    labs(x = tmp.x.axis, y = tmp.y.axis, fill = "PDCD1 %") +
    scale_y_continuous(limits = c(0,2000)) + coord_flip() +
    scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd8.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd8.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y = element_blank(), axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"), legend.position="none")

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
  # print(p1)
  dev.off()

}

{ cat(redb("### GLIPH2 plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("tcr/gliph")

  donor_order <- c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5", "BT24", "BT25", "BT10", "BT15", "BT22", "BT11", "BT12", "BT18",
  "BT21", "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23")

  # Version 5 - Back2back plot. Stacked barplot of EXPANDED CELLS. Each box is a different specificity group and the size is the cell within it. color Cd8 and Cd4 as in Fig1A

  donor_order <- rev(c("BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT22", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT23", "BT13", "BT20", "BT17_brain", "BT2"))

  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8_Expanded.pdf'))
  # CD4
  tmp <- sc_cd3p_cd4@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  p1 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "cadetblue4", fill = "white", size = 0.15) + # gray60
      labs(x = "Patient", y = "CD4 Cells") +
      scale_y_reverse(limits = c(100,0)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))
  # CD8
  tmp <- sc_cd3p_cd8@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  p2 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "royalblue4", fill = "white", size = 0.15) + # gray60
      labs(x = "Patient", y = "CD8 Cells") +
      scale_y_continuous(limits = c(0,400)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
            axis.line.y = element_blank(), axis.ticks.y=element_blank(),
            plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()


  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8.pdf'))
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
  p1 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "cadetblue4", fill = "white", size = 0.15) + # gray60
      labs(x = "Patient", y = "CD4 Cells") +
      scale_y_reverse(limits = c(400,0)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))
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
  p2 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "royalblue4", fill = "white", size = 0.15) + # gray60
      labs(x = "Patient", y = "CD8 Cells") +
      scale_y_continuous(limits = c(0,600)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
            axis.line.y = element_blank(), axis.ticks.y=element_blank(),
            plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()


  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8_v2.pdf'))
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
  p1 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "cadetblue4", fill = "white", size = 0.05) + # gray60
      labs(x = "Patient", y = "CD4 Cells") +
      scale_y_reverse(limits = c(400,0)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))
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
  p2 <- df %>% group_by(orig.donor, pattern) %>% summarize(cells = n()) %>% arrange(cells) %>% unite(col = donor_pattern, orig.donor, pattern, sep = "_", remove = FALSE) %>%
    mutate(donor_pattern = factor(donor_pattern, levels = donor_pattern)) %>%
    ggplot(aes(x = orig.donor, y = cells, fill = donor_pattern)) + geom_bar(stat="identity", position = "stack", color = "royalblue4", fill = "white", size = 0.05) + # gray60
      labs(x = "Patient", y = "CD8 Cells") +
      scale_y_continuous(limits = c(0,600)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
            axis.line.y = element_blank(), axis.ticks.y=element_blank(),
            plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()

  ############################################################
  # NUMBER of Specificity groups & clones per donor (tables)
  ############################################################

  # ---> Expanded

  # CD4
  tmp <- sc_cd3p_cd4@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  df %>% group_by(orig.donor) %>% summarize(n_SG = length(unique(pattern))) %>% as.data.frame()

  # CD8
  tmp <- sc_cd3p_cd8@meta.data %>% filter(sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  df %>% group_by(orig.donor) %>% summarize(n_SG = length(unique(pattern))) %>% as.data.frame()

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

  # Plot

  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8_simplified.pdf'))
  p1 <- df_cd4 %>%
    ggplot(aes(x = orig.donor, y = n_SG)) + geom_bar(stat="identity", fill = "cadetblue4", size = 0.15) + # gray60 , color = "cadetblue4"
      labs(x = "Patient", y = "CD4 Specificity Groups") +
      scale_y_reverse(limits = c(250,0)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))

  p2 <- df_cd8 %>%
    ggplot(aes(x = orig.donor, y = n_SG)) + geom_bar(stat="identity", fill = "royalblue4", size = 0.15) + # gray60 , color = "royalblue4"
      labs(x = "Patient", y = "CD8 Specificity Groups") +
      scale_y_continuous(limits = c(0,200)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
            axis.line.y = element_blank(), axis.ticks.y=element_blank(),
            plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
  dev.off()


  # ------------------------- Number of Expanded specificity groups -------------------------
  # NOTE: Definded by the presence of at least 1 expanded clonotype in that group

  date <- "2022-07-15" # date <- "2022-02-07" #
  HTO <- TRUE
  suffix <- c("CD8", "CD4"); if(HTO) suffix <- paste0(suffix,"_hto")

  # Read GLIPH Output
  results.path = paste0('/home/kmlanderos/tmp_large/pbtumor-all/results/tcr/gliph/',date)
  tmp.file.name <- paste0(results.path, '/GLIPH2Output_TCR-Data_', suffix, '.csv')
  gliph_outs <- lapply(tmp.file.name, function(x) {
    gliph.res <- fread(file=x, blank.lines.skip=TRUE)
    gliph.res <- gliph.res[pattern!='single']
    gliph.res$Sample <- gsub(":Exp1", "", gliph.res$Sample)
    return(gliph.res)
  })
  names(gliph_outs) <- gsub("_hto","",suffix)

  # Get number of Expanded specificity groups (presence of at least 1 expanded clonotype in that group)
  # gliph_outs[[1]] %>% select(index, Freq) %>% data.frame() %>% head(100)
  gliph_outs[[1]] %>% group_by(index) %>% summarize(expansion_tag = any(Freq > 1)) %>% pull(expansion_tag) %>% sum()
  gliph_outs[[2]] %>% group_by(index) %>% summarize(expansion_tag = any(Freq > 1)) %>% pull(expansion_tag) %>% sum()


}

{ cat(redb("### GG-paired plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("heatmaps")

  # ---> PDCD1n and PDCD1p Gene expression (CD8) - MANY GENES

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

  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd8$orig.donor)) %>% group_by(orig.donor, PDCD1_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, PDCD1_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = PDCD1_tag, values_from = pct_exp_cells)

  write.csv(df, file = "heatmaps/paired_dots_CD8_PDCD1_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD8_PDCD1_expansion_pct.pdf", 7, 7)
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

  write.csv(df, file = "heatmaps/paired_dots_CD8_expansion_PDCD1_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD8_expansion_PDCD1_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "Non_expanded", cond2 = "Expanded", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of cells expressing PDCD1")
  dev.off()

  # CD4

  # ---> PDCD1n and PDCD1p clonal expansion percentage (%) (CD4)

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd4@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd4@assays$RNA@data)
  mdata_sc_cd3p_cd4[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))
  mdata_sc_cd3p_cd4$expansion_tag <- ifelse(mdata_sc_cd3p_cd4$clon.size.tag > 1, TRUE, FALSE)

  df <- mdata_sc_cd3p_cd4 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd4$orig.donor)) %>% group_by(orig.donor, PDCD1_tag, expansion_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, PDCD1_tag) %>% summarize(pct_exp_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = PDCD1_tag, values_from = pct_exp_cells)

  write.csv(df, file = "heatmaps/paired_dots_CD4_PDCD1_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD4_PDCD1_expansion_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # ---> Expanded and Non-Expanded cells PDCD1 expression (CD4)

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

  df <- mdata_sc_cd3p_cd4 %>% filter(!is.na(expansion_tag) & !is.na(mdata_sc_cd3p_cd4$orig.donor)) %>% group_by(orig.donor, expansion_tag, PDCD1_tag) %>% summarize(n=n()) %>% ungroup() %>%
    group_by(orig.donor, expansion_tag) %>% summarize(pct_PDCD1_cells = 100*n[2]/sum(n))  %>%
    pivot_wider(names_from = expansion_tag, values_from = pct_PDCD1_cells)

  write.csv(df, file = "heatmaps/paired_dots_CD4_expansion_PDCD1_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD4_expansion_PDCD1_pct.pdf", 7, 7)
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
  donors <- c("BT5", "BT1", "BT7", "BT4", "BT3", "BT8", "BT2", "BT18", "BT12", "BT11", "BT15", "BT17_brain", "BT13", "BT9", "BT10", "BT21", "BT25", "BT26", "BT23", "BT20", "BT19", "BT24", "BT27", "BT22")

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

  # save table
  write.csv(df, file = "heatmaps/paired_dots_CD8_LAG3_percentageExpr.csv", row.names = FALSE)

  # ---> PDCD1n and PDCD1p clonal expansion percentage (%) (CD8)

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

  write.csv(df, file = "heatmaps/paired_dots_CD8_LAG3_expansion_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD8_LAG3_expansion_pct.pdf", 7, 7)
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

  write.csv(df, file = "heatmaps/paired_dots_CD8_expansion_LAG3_pct.csv", row.names = FALSE)

  # Plot
  pdf("heatmaps/paired_dots_CD8_expansion_LAG3_pct.pdf", 7, 7)
  ggpaired(df, cond1 = "Non_expanded", cond2 = "Expanded", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of cells expressing LAG3")
  dev.off()


}

{ cat(redb("### Cluster Markers ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("cluster_markers")

  fnames <- paste0(
    "/home/kmlanderos/tmp_large/pbtumor-all/results/clustering/CD45pCD3p_clean2/",
    c("CD4/seurat_mean0.01_pct25_pc15_res0.4/dgea_MAST_fc0.25_padj0.05_summary_stats.csv",
      "CD8/seurat_mean0.01_pct20_pc20_res0.2/dgea_MAST_fc0.25_padj0.05_summary_stats.csv")
  );names(fnames) <- c("CD4", "CD8")

  ## ---> CD8

  celltype <- "CD8"
  n_features <- 50

  df <- read.csv(fnames[celltype], stringsAsFactors = FALSE, row.names = 1, check.names = FALSE) %>%
    select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  head(df); dim(df)

  top10 <- df %>%
    group_by(cluster) %>%
    top_n(n = n_features, wt = avg_logFC)
  write.table(top10, paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features.csv"), sep = ",", row.names = F, quote = F)

  pdf(paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"))
  sc_cd3p_cd8$cluster <- factor(sc_cd3p_cd8$cluster, levels = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9"))
  top10 <- top10 %>% arrange(factor(cluster, levels = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")))
  DoHeatmap(sc_cd3p_cd8, features = top10$gene, slot = "scale.data", group.by = "cluster", draw.lines = F, group.colors = c("#72e5d2", "#3bb8ae", "#007e74", "#ffa7d7", "#ff6ede", "#fcee62", "#9e89e8", "#ff971f", "#4da7fd", "#9c0006")) +
  scale_fill_gradientn(colors = c("blue", "black", "yellow")) #+ NoLegend()
  dev.off()

  ## ---> CD4

  celltype <- "CD4"
  n_features <- 50

  # df <- read.csv(fnames[celltype], stringsAsFactors = FALSE, row.names = 1, check.names = FALSE) %>%
  #   select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  # head(df); dim(df)

  Idents(sc_cd3p_cd4) <- "cluster.new"
  df <- FindAllMarkers(
          object = sc_cd3p_cd4,
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

  # Idents(sc_cd3p_cd4) <- "cluster"
  # df <- FindAllMarkers(sc_cd3p_cd4, test.use = "MAST")

  df %>% select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  top10 <- df %>%
    group_by(cluster) %>%
    top_n(n = n_features, wt = avg_logFC)
  write.table(top10, paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features.csv"), sep = ",", row.names = F, quote = F)

  pdf(paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"))
  sc_cd3p_cd4$cluster.new <- factor(sc_cd3p_cd4$cluster.new, levels = c("0", "1.1", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11"))
  top10 <- top10 %>% arrange(factor(cluster, levels = c("0", "1.1", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")))
  DoHeatmap(sc_cd3p_cd4, features = top10$gene, slot = "scale.data", group.by = "cluster.new", draw.lines = F, group.colors = c("#c5b8ff", "#d26dec", "#a084ff", "#8c50cf", "#5d1d88", "#52edf4", "#c0f1f5", "#97ccf5", "#379ee4", "#1069b9", "#ffa7d7", "#ff717c", "#9c0006")) +
  scale_fill_gradientn(colors = c("blue", "black", "yellow")) #+ NoLegend()
  dev.off()


}

{ cat(redb("### Supplementary QC Plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("QC")

  # --- CD8
  df <- sc_cd3p_cd8@meta.data %>%
    mutate( cluster=factor(cluster,levels=c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")) )

  p1 <- ggplot(df, aes(x=cluster, y=nCount_RNA, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#72e5d2", "#3bb8ae", "#007e74", "#ffa7d7", "#ff6ede", "#fcee62", "#9e89e8", "#ff971f", "#4da7fd", "#9c0006")) + # c("#72e5d2", "#3bb8ae", "#fcee62", "#ff971f", "#4da7fd", "#9e89e8", "#ffa7d7", "#007e74", "#ff6ede", "#9c0006")
    theme_classic() +
    theme(legend.position="none")

  p2 <- ggplot(df, aes(x=cluster, y=nFeature_RNA, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#72e5d2", "#3bb8ae", "#007e74", "#ffa7d7", "#ff6ede", "#fcee62", "#9e89e8", "#ff971f", "#4da7fd", "#9c0006")) + # c("#72e5d2", "#3bb8ae", "#fcee62", "#ff971f", "#4da7fd", "#9e89e8", "#ffa7d7", "#007e74", "#ff6ede", "#9c0006")
    theme_classic() +
    theme(legend.position="none")

  p3 <- ggplot(df, aes(x=cluster, y=percent.mt, fill = cluster)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#72e5d2", "#3bb8ae", "#007e74", "#ffa7d7", "#ff6ede", "#fcee62", "#9e89e8", "#ff971f", "#4da7fd", "#9c0006")) + # c("#72e5d2", "#3bb8ae", "#fcee62", "#ff971f", "#4da7fd", "#9e89e8", "#ffa7d7", "#007e74", "#ff6ede", "#9c0006")
    theme_classic() +
    theme(legend.position="none")

  df <- sc_cd3p_cd8@meta.data %>% select(cluster, origlib) %>% group_by(cluster, origlib) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")) )
  p4 <- ggplot(df, aes(fill=origlib, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Libraries", fill = "Library")

  df <- sc_cd3p_cd8@meta.data %>% select(cluster, orig.Run) %>% group_by(cluster, orig.Run) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")) )
  p5 <- ggplot(df, aes(fill=orig.Run, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Runs", fill = "Run")

  df <- sc_cd3p_cd8@meta.data %>% select(cluster, orig.donor) %>% group_by(cluster, orig.donor) %>% summarize(n = n()) %>%
    mutate( cluster=factor(cluster,levels=c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")) )
  p6 <- ggplot(df, aes(fill=orig.donor, y=n, x=cluster)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Donors", fill = "Donor")

  # Change violin plot colors by groups
  pdf("QC/qc_cd8.pdf", 10, 7)
  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
  dev.off()

  # --- CD4

  df <- sc_cd3p_cd4@meta.data %>%
    mutate( cluster.new=factor(cluster.new,levels=c("0", "1.1", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")) )

  p1 <- ggplot(df, aes(x=cluster.new, y=nCount_RNA, fill = cluster.new)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#c5b8ff", "#d26dec", "#a084ff", "#8c50cf", "#5d1d88", "#52edf4", "#c0f1f5", "#97ccf5", "#379ee4", "#1069b9", "#ffa7d7", "#ff717c", "#9c0006")) +
    theme_classic() +
    theme(legend.position="none")

  p2 <- ggplot(df, aes(x=cluster.new, y=nFeature_RNA, fill = cluster.new)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#c5b8ff", "#d26dec", "#a084ff", "#8c50cf", "#5d1d88", "#52edf4", "#c0f1f5", "#97ccf5", "#379ee4", "#1069b9", "#ffa7d7", "#ff717c", "#9c0006")) +
    theme_classic() +
    theme(legend.position="none")

  p3 <- ggplot(df, aes(x=cluster.new, y=percent.mt, fill = cluster.new)) +
    geom_jitter(width = 0.10) +
    geom_violin(trim=FALSE, color="black", alpha = 0.8) +
    geom_boxplot(width=0.1) +
    scale_fill_manual(values = c("#c5b8ff", "#d26dec", "#a084ff", "#8c50cf", "#5d1d88", "#52edf4", "#c0f1f5", "#97ccf5", "#379ee4", "#1069b9", "#ffa7d7", "#ff717c", "#9c0006")) +
    theme_classic() +
    theme(legend.position="none")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster.new, origlib) %>% group_by(cluster.new, origlib) %>% summarize(n = n()) %>%
    mutate( cluster.new=factor(cluster.new,levels=c("1.1", "0", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")) )
  p4 <- ggplot(df, aes(fill=origlib, y=n, x=cluster.new)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Libraries", fill = "Library")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster.new, orig.Run) %>% group_by(cluster.new, orig.Run) %>% summarize(n = n()) %>%
    mutate( cluster.new=factor(cluster.new,levels=c("1.1", "0", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")) )
  p5 <- ggplot(df, aes(fill=orig.Run, y=n, x=cluster.new)) +
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Cluster", y = "Runs", fill = "Run")

  df <- sc_cd3p_cd4@meta.data %>% select(cluster.new, orig.donor) %>% group_by(cluster.new, orig.donor) %>% summarize(n = n()) %>%
    mutate( cluster.new=factor(cluster.new,levels=c("1.1", "0", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")) )
  p6 <- ggplot(df, aes(fill=orig.donor, y=n, x=cluster.new)) +
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
