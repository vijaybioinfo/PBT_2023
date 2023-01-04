#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2022-01-15
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


  { # Change clone size of problematic clonotypes
    # # Let us recalculate the clone size of the clonotypes shared between CD4 and CD4, and only consider cells with Gex
    # sc_cd3p_cd8@meta.data$clon.size.tag.old <- sc_cd3p_cd8@meta.data$clon.size.tag # Save old values
    #
    # # Get the problematic clonotypes
    # cd8_table <- table(sc_cd3p_cd8@meta.data$clonotype.tag) %>% reshape2::melt() %>% arrange(desc(value))
    # cd4_table <- table(sc_cd3p_cd4@meta.data$clonotype.tag) %>% reshape2::melt() %>% arrange(desc(value))
    # merge_table <- merge(cd8_table, cd4_table, by = "Var1");
    # problematic_clonotypes <- as.character(merge_table$Var1)
    #
    # # CD8
    # clonotype.cs.cd8 <- as.data.frame(table(sc_cd3p_cd8@meta.data %>% pull(clonotype.tag)) ) %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd8) <- c("clonotype", "clone_size")
    # tmp <- unlist(apply(sc_cd3p_cd8@meta.data, MARGIN = 1, function(x){
    #   if(!is.na(x["clonotype.tag"]) & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd8[clonotype.cs.cd8$clonotype == x["clonotype.tag"],"clone_size"] else NA
    # }))
    # sc_cd3p_cd8@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]
    # # sc_cd3p_cd8@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))
    #
    # # CD4
    # sc_cd3p_cd4@meta.data$clon.size.tag.old <- sc_cd3p_cd4@meta.data$clon.size.tag # Save old values
    #
    # clonotype.cs.cd4 <- as.data.frame(table(sc_cd3p_cd4@meta.data %>% pull(clonotype.tag)) )  %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd4) <- c("clonotype", "clone_size")
    # tmp <- unlist(apply(sc_cd3p_cd4@meta.data, MARGIN = 1, function(x){
    #   if(!is.na(x["clonotype.tag"]) & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd4[clonotype.cs.cd4$clonotype == x["clonotype.tag"],"clone_size"] else NA
    # }))
    # sc_cd3p_cd4@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]
    # # sc_cd3p_cd4@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))
  }

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
  saveRDS(sc_cd3p_cd8@meta.data, file = "data/sc_cd3p_cd8_mdata.rds")
  saveRDS(sc_cd3p_cd4@meta.data, file = "data/sc_cd3p_cd4_mdata.rds")
  saveRDS(sc_cd3p_cd8, file = "data/sc_cd3p_cd8_seurat_object.rds")
  saveRDS(sc_cd3p_cd4, file = "data/sc_cd3p_cd4_seurat_object.rds")

}

{ cat(redb("### Signature Analyses ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # cp /Users/kmlanderos/Documents/PBT/Data/signatures_cd3p*csv /Volumes/kmlanderos/pbtumor-all/info/
  # Rscript /home/ciro/scripts/functions/csvCorrect.R /home/kmlanderos/pbtumor-all/info/signatures_cd3p.csv

  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R") # stats_summary_table
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
  source("/home/ciro/scripts/handy_functions/devel/utilities.R")  # vlist2df
  # source("/home/ciro/scripts/handy_functions/devel/code.R") #vlist2df
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
    # Neo TCR - CD4 all; expansion comparison
    list(result_id = "gsea/CD4_neoantigenTCR_all_expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
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
    # Neo TCR - CD8 all; expansion comparison
    list(result_id = "gsea/CD8_neoantigenTCR_all_expansion/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("expDegree")
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
    ),
    # Supplement
    list(result_id = "gsea/Supplement/CD8_CX3CR1_C2vsREST/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(45,46)], FUN = function(x) x[-c(1:2)] ), # 46 is the one that matters
      comparisons = c("cell_classification", "INNATE-LIKE_NKR", "REST")
    ),
    # CD4_nonTREG_Expansion
    list(result_id = "gsea/CD4_nonTREG_Expansion/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "TRM_CLARKE", "Trm_Savas2018", "TRM_192_UP_HobitScience", "th1_signature1_arlehamn", "th2_signature0_arlehamn", "TH17_Seumois2020")], FUN = function(x) x[-c(1:2)] ), # 46 is the one that matters
      sample_filter = c("cell_classification", "-TREG"),
      comparisons = c("expDegree")
    ),
    # Supplement
    list(result_id = "gsea/Supplement/CD8_TSCM_C3vsREST/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(47,48)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster", "3", "REST")
    ),
    # Supplement
    list(result_id = "gsea/Supplement/CD8_TSCM_C4vsREST/",
      edata = "expm1(sc_cd3p_cd8@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd8",
      lists = lapply(X = slists[,c(47,48)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster", "4", "REST")
    ),
    # Supplement
    list(result_id = "gsea/Supplement/CD4_TSCM_TCMvsREST/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(49,50)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cell_classification", "TCM", "REST")
    ),
    # Supplement
    list(result_id = "gsea/Supplement/CD4_TSCM_allvsall/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c(49,50)], FUN = function(x) x[-c(1:2)] ),
      comparisons = c("cluster")
    ),
    # ---> CD4 INTERMEDIATE sub-clustering
    # Identify CTL and TCM/TN
    list(result_id = "gsea/CD4_clust1_subclustering/",
      edata = "expm1(sc_cd3p_cd4@assays$RNA@data)",
      metadata = "mdata_sc_cd3p_cd4",
      lists = lapply(X = slists[,c("CD4CTLvsTCM_Patil2018", "Hu.CD4.Naive.CM")], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("RNA_snn_res.deconv", "1.0", "1.1", "1.2"),
      comparisons = c("RNA_snn_res.deconv")
    )
  )

  pp_gsea = fig_gsea(fconfigs[2], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[3:4], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  # ---> Neo TCR
  pp_gsea = fig_gsea(fconfigs[5:9], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[10:15], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # ---> TSCM
  pp_gsea = fig_gsea(fconfigs[16], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[17:21], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)

  pp_gsea = fig_gsea(fconfigs[c(6,11)], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # Supplement
  pp_gsea = fig_gsea(fconfigs[22], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # CD4_nonTREG_Expansion
  pp_gsea = fig_gsea(fconfigs[23], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # Supplement
  pp_gsea = fig_gsea(fconfigs[c(24,25)], features = rownames(sc_cd3p_cd8), return_plot = TRUE, verbose = T)
  # Supplement
  pp_gsea = fig_gsea(fconfigs[c(26,27)], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)
  # CD4 INTERMEDIATE sub-clustering
  pp_gsea = fig_gsea(fconfigs[c(28)], features = rownames(sc_cd3p_cd4), return_plot = TRUE, verbose = T)

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
  pp_ms = fig_module_score(fconfigs, reductions = redu[1], violins_color = 'pct')

  #   couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000"),

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
  # Vicente's color scale 0.15
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.15_umap_v3.pdf"))
  mdata_sc_cd3p_cd4 %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = CD4CTL_signature_adjusted)) + geom_point(size = 0.1) + scale_color_gradientn(colours = c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')) + #, "#670000"
    labs(colour = NULL, x = NULL, y = NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
  graphics.off()
  # Vicente's color scale 0.2
  pdf(paste0(dir, "cd4ctlvstcm_patil2018_0.2_umap_v3.pdf"))
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
     ),
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
     ),
    list(result_id = "dotplots/CD4/greg_ThIFNRs_",
       edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
       size = c(7,12),
       object = "sc_cd3p_cd4", axis_x = list(col = "RNA_snn_res.0.4",
         order = c("0", "2", "9", "4", "1", "10", "5", "3", "6", "7", "8", "11")),
       features = c("CD3D", "CD4", "CD8A", "CD8B", "GZMA", "GZMK", "PRF1", "MX1", "IFI44L", "IFI6",  "ISG15", "LY6E", "TNFSF10", "CXCL10", "TNF", "NFKB1", "CXCL13", "BTLA", "SH2D1A", "PDCD1", "FOXP3", "CTLA4", "IL2RA", "ENTPD1", "S1PR1", "SELL", "CCR7", "TCF7", "STMN1", "TOP2A", "TUBA1B")
     )
   )
  pp_curtains = fig_plot_curtain(fconfigs[11], dot.scale = 12, col.min = -1.5, col.max = 1.5)
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
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/LAG3_new_colorScale.pdf")
  FeaturePlot(sc_cd3p_cd8, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) )
  dev.off()
  # Blank version
  pdf("dim_reduction_features/sc_cd3p_cd8_exm_markers/LAG3_new_colorScale_blank.pdf")
  plot_blank_fun(FeaturePlot(sc_cd3p_cd8, features = c("LAG3"), order = T) +
    scale_colour_gradientn(colours = colorRampPalette(c("gray90", "#5800FF"))(100) ))
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

{ cat(redb("### Markers: violins ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("violins")

  # --- CD8

  mdata_sc_cd3p_cd8 <- sc_cd3p_cd8@meta.data

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd8@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd8@assays$RNA@data)
  mdata_sc_cd3p_cd8[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))


  fconfigs = list(
    # CD8 Expanded vs Non_Expanded
    list(result_id = "violins/sc_cd3p_cd8_expansion/violin_",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "mdata_sc_cd3p_cd8[!is.na(mdata_sc_cd3p_cd8$clon.size.tag),]",
      axis_x = list(col = "expDegree", order = c("Non_expanded", "Expanded")),
      # sample_filter = c("celltype", "CD8"),
      features = c("PRF1", "GZMA", "GZMK", "GZMH", "IFNG", "TNF", "CCL4", "CCL5", "HLA-DRB1", "LAG3", "EOMES", "PRDM1", "TOX", "ITGAE", "CCR5", "FASLG", "CCL3", "ZNF683", "CXCR6", "PDCD1")
    ),  # CD8 PDCD1p vs PDCD1n
    list(result_id = "violins/sc_cd3p_cd8_PDCD1p_PDCD1n/violin_",
      edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "mdata_sc_cd3p_cd8",
      axis_x = list(col = "PDCD1_tag", order = c("PDCD1n", "PDCD1p")),
      # sample_filter = c("celltype", "CD8"),
      features = c("PRF1", "GZMA", "GZMK", "GZMH", "IFNG", "TNF", "CCL4", "CCL5", "HLA-DRB1", "LAG3", "EOMES", "PRDM1", "TOX", "ITGAE", "CCR5", "FASLG", "GZMB")
    )
  )

  pp_violins = fig_plot_violins(
    fconfigs,
    theme_extra = labs(x = NULL, y = "Seurat Normalized"),
    colour_by = "pct", couls = couls_opt$red_gradient$white
  )

  # Create table with info per donor to make paired plots
  genes <- c("PRF1", "GZMA", "GZMK", "GZMH", "IFNG", "TNF", "CCL4", "CCL5", "HLA-DRB1", "LAG3", "EOMES", "PRDM1", "TOX", "ITGAE", "CCR5", "FASLG", "GZMB")
  gene_tmp_tags <- as.matrix(sc_cd3p_cd8@assays$RNA@data[genes,]) > 0
  mdata_sc_cd3p_cd8 <- cbind(mdata_sc_cd3p_cd8, t(gene_tmp_tags)[rownames(mdata_sc_cd3p_cd8),])
  df <- mdata_sc_cd3p_cd8 %>% filter(!is.na(orig.donor)) %>% group_by(orig.donor, PDCD1_tag) %>% summarize( across(genes, ~ 100*sum(.x, na.rm = TRUE)/length(.x)) ) %>%
    pivot_longer(genes, names_to = "gene", values_to = "pct_expr") %>% data.frame()

  write.table(df, file = "violins/sc_cd3p_cd8_PDCD1p_PDCD1n/table_donorwise_effector_molecules_pct_expr.csv", sep = ",", row.names = F, quote = F)

  # --- CD4

  mdata_sc_cd3p_cd4 <- sc_cd3p_cd4@meta.data

  # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  PDCD1_tag_tmp <- as.matrix(sc_cd3p_cd4@assays$RNA@data["PDCD1",]) > 0
  PDCD1_cells <- colnames(sc_cd3p_cd4@assays$RNA@data)
  mdata_sc_cd3p_cd4[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  mdata_sc_cd3p_cd4 <- mdata_sc_cd3p_cd4 %>% mutate(PDCD1_tag = case_when(
    PDCD1_tag_tmp > 0 ~ "PDCD1p",
    PDCD1_tag_tmp == 0 ~ "PDCD1n"
  ))

  fconfigs = list(
    # CD4 Expanded vs Non_Expanded (non-TREG)
    list(result_id = "violins/sc_cd3p_cd4_nonTREG_expansion/violin_",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "mdata_sc_cd3p_cd4[!is.na(mdata_sc_cd3p_cd4$clon.size.tag),]",
      axis_x = list(col = "expDegree", order = c("Non_expanded", "Expanded")),
      sample_filter = c("cell_classification", "-TREG"),
      features = c("CCL4","CCL5","GZMA","GZMH","GZMK","IFNG","TBX21","PDCD1","TOX","LAG3","RGS1","CXCR6","RBPJ","ITGAE","PRDM1","RUNX2","HOPX","BATF",
        "RUNX3","HLA-DRB1","NKG7","THEMIS","PRF1")
      # features = c("PDCD1", "TBX21", "PRF1", "GZMA", "HOPX", "PRDM1", "GZMK", "GZMH", "TNF", "TOX", "IFNG", "CCL5", "CXCR6", "CXCR3", "TIGIT", "ICOS",
      #   "TNFRSF4", "CTLA4", "BATF", "RBPJ", "HLA-DRB1", "IL2RB", "CD40LG")
    ),
    # CD4 PDCD1p vs PDCD1n (non-TREG)
    list(result_id = "violins/sc_cd3p_cd4_nonTREG_PDCD1p_PDCD1n/violin_",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "mdata_sc_cd3p_cd4",
      axis_x = list(col = "PDCD1_tag", order = c("PDCD1n", "PDCD1p")),
      sample_filter = c("cell_classification", "-TREG"),
      features = c("CCL4", "CCL5", "GZMA", "GZMH", "GZMK", "IFNG", "TBX21", "PDCD1", "TOX", "LAG3", "RGS1", "CXCR6", "RBPJ", "ITGAE", "PRDM1", "RUNX2",
        "HOPX", "BATF", "RUNX3", "HLA-DRB1", "NKG7", "THEMIS", "GZMB", "CD69", "SELL", "CCR7", "S1PR1", "PRF1", "TNF")
    )
  )

  pp_violins = fig_plot_violins(
    fconfigs[2],
    theme_extra = labs(x = NULL, y = "Seurat Normalized"),
    colour_by = "pct", couls = couls_opt$red_gradient$white
  )

  # Create table with info per donor to make paired plots
  genes <- c("CCL4", "CCL5", "GZMA", "GZMH", "GZMK", "IFNG", "TBX21", "PDCD1", "TOX", "LAG3", "RGS1", "CXCR6", "RBPJ", "ITGAE", "PRDM1", "RUNX2", "HOPX", "BATF", "RUNX3", "HLA-DRB1", "NKG7", "THEMIS", "GZMB", "CD69", "SELL", "CCR7", "S1PR1", "PRF1", "TNF")
  gene_tmp_tags <- as.matrix(sc_cd3p_cd4@assays$RNA@data[genes,]) > 0
  mdata_sc_cd3p_cd4 <- cbind(mdata_sc_cd3p_cd4, t(gene_tmp_tags)[rownames(mdata_sc_cd3p_cd4),])
  df <- mdata_sc_cd3p_cd4 %>% filter(!is.na(orig.donor)) %>% group_by(orig.donor, PDCD1_tag) %>% summarize( across(genes, ~ 100*sum(.x, na.rm = TRUE)/length(.x)) ) %>%
    pivot_longer(genes, names_to = "gene", values_to = "pct_expr") %>% data.frame()

  write.table(df, file = "violins/sc_cd3p_cd4_nonTREG_PDCD1p_PDCD1n/table_donorwise_effector_molecules_pct_expr.csv", sep = ",", row.names = F, quote = F)



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


}

{ cat(redb("### Markers: contour ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("scatter_contour")

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

  # Gene comparisons
  # tmp <- list(
  #   list(x = "GZMA", y = c("PRF1"))
  # )

  # Expanded
  # fconfigs = lapply(tmp, function(x){
  #   list(result_id = "scatter_contour/sc_cd4_Expanded/",
  #     edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4_mdata[sc_cd3p_cd4_mdata$cluster!='5',]", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
  #     features = list(x = "GZMA", y = c("PRF1")),
  #     sample_filter = list(c("expDegree", "Expanded"))
  #   )
  # })
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

  # ------------------------ LAG3 CO-EXPRESSION (CD8)

  fconfigs = list(
      list(result_id = "scatter_contour/sc_cd8_LAG3/",
        edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
        features = list(x = c("TNF", "IFNG", "GZMA", "GZMB", "PDCD1", "GZMK", "CCL4", "CCL5", "HLA-DRB1", "ITGAE", "PRF1"), y = c("LAG3"))#,
        # sample_filter = list(c("expDegree", "Expanded"))
      )
  )
  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent") #  + scale_x_continuous(limits = c(0, )) + scale_y_continuous(limits = c(0, 5)) # , type = "percent"
  )

  # Donor wise plots
  dir.create("scatter_contour/sc_cd8_LAG3/donor_wise")

  donors <- unique(sc_cd3p_cd8$orig.donor); donors <- donors[!is.na(donors)]
  fconfigs <- rep(list(
      list(result_id = "scatter_contour/sc_cd8_LAG3/donor_wise/",
        edata = "sc_cd3p_cd8@assays$RNA@data", metadata = "sc_cd3p_cd8@meta.data", # [sc_cd3p_cd4_mdata$ExpDegree=='high',]
        features = list(x = c("TNF", "IFNG", "GZMA", "GZMB", "PDCD1", "GZMK", "CCL4", "CCL5", "HLA-DRB1", "ITGAE", "PRF1"), y = c("LAG3")),
        sample_filter = list(c("orig.donor", "donor"))
      )
  ), length(donors))

  fconfigs <- lapply(1:length(donors), function(i){
    fconfigs[[i]]$result_id <- paste0(fconfigs[[i]]$result_id, donors[i],  "/")
    fconfigs[[i]]$sample_filter <- list(c("orig.donor", donors[i]))
    return(fconfigs[[i]])
  })

  pp_contour = fig_plot_contour(fconfigs,
    theme_extra = function(x) plot_add_quadrants(x + viridis::scale_color_viridis(option = "magma", limits = c(0,0.8)) + scale_y_continuous(limits = c(0, 5)) + geom_point(size=2), limits = list(0.5, 0.5), type = "percent") #  + scale_x_continuous(limits = c(0, 5)) + scale_y_continuous(limits = c(0, 5))
  )
  # for i in $(ls | grep pdf); do mv $i "${i%.pdf}_percent.pdf"; done
  # for i in $(ls); do for j in $(ls $i | grep pdf); do mv $i/$j "$i/${j%.pdf}_percent.pdf"; done; done


  scatter_contour <- function(
    data, axis_x, axis_y, axis_z = NULL, dp = TRUE) {
    if(is.null(data$Density) && is.null(axis_z)){
      data$Density <- 0
      dps = if(isTRUE(dp)){
        rowSums(data[, c(axis_x, axis_y)] > 0) > 1
      }else{ rownames(data) }
      x.sp = setdiff(rownames(data)[data[, c(axis_x)] > 0], rownames(data)[dps==T])
      y.sp = setdiff(rownames(data)[data[, c(axis_y)] > 0], rownames(data)[dps==T])
      dns <- rownames(data)[rowSums(data[, c(axis_x, axis_y)] > 0) == 0]
      tmp = try(MASS_kde2d( # Add density to mdata
        x = data[dps, axis_x],
        y = data[dps, axis_y]
      ), silent = TRUE)
      if(class(tmp) != "try-error"){
        data[dps, ]$Density <- tmp; axis_z = "Density" # Add density to mdata
      }else{ data$Density = NULL }
    }
    aesy = if(!is.null(axis_z)){
      aes_string(x = axis_x, y = axis_y, color = axis_z)
      aes(x = !!rlang::sym(axis_x), y = !!rlang::sym(axis_y), color = !!rlang::sym(axis_z))
    }else{ aes(x = !!rlang::sym(axis_x), y = !!rlang::sym(axis_y)) }
    p <- ggplot() + geom_point(data = data[dps, ], mapping = aesy, size = 2) + geom_point(data = data[y.sp, ], mapping = aesy, size = 1.5, position = position_jitterdodge(jitter.width = 0.2)) + # + geom_jitter(width = 0.20)
      geom_point(data = data[x.sp, ], mapping = aesy, size = 1.5, position = position_jitterdodge(jitter.height = 0.2)) +
      geom_point(data = data[dns, ], mapping = aesy, size = 2) # , position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05)
    return(p)
  }

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
    )
  )

  # fconfig=fconfigs[[4]]
  for (fconfig in fconfigs[4]){
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
    ),
    list(result_id = "dim_reduction/sc_cd3p_cd4",
      edata = "sc_cd3p_cd4@assays$RNA@data", metadata = "sc_cd3p_cd4@meta.data",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      vars2col = "cluster",
      colors = c("0" = "#c5b8ff", "1" = "#5fcfab", "7" = "#e57f08", "6" = "#ffaf4d", "8" = "#b35806", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#fcd093",
       "9" = "#5d1d88", "10" = "#ff717c", "11" = "#9c0006"),
      redu = redu[[1]]
    ),
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
  pp_main_dr = fig_plot_base(
    fconfigs[11],
    theme_extra = function(x){
      p <- x + geom_point(size = 0.4, na.rm = TRUE) + labs(x = "Dim 1", y = "Dim 2", color = NULL) +
      guides(col = guide_legend(ncol = 1, override.aes = list(size = 6), label.position = "left")) #+
      # xlim(-10, 10) + ylim(-10, 10)
      # Seurat::LabelClusters(plot = p, id = "cluster", color = "black") # , fontface = "bold"  , point.padding = 1, min.segment.length = 0, arrow = arrow(length = unit(0.015, "npc"))
    }
  )
  pp_main_dr = fig_plot_base(
    fconfigs[12],
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
  # c("0" = "#72e5d2", "1" = "#72e5d2", "7" = "#72e5d2", "6" = "#ffa7d7", "8" = "#ffa7d7", "2" = "#fcee62", "5" = "#9e89e8", "4" = "#4da7fd", "3" = "#ff971f",
  #  "9" = "#9c0006")
  #  CD8 GZMK-HI, CD8 TRM, CD8 16+ Effector (C2), CD8 Effector (C3), CD8 TCF7-HI (C4), CD8 MAIT (C5), CD8 cell cycle (C9).
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
    fconfigs[[i]]$sufix <- "cluster_clone_size_v2_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification"; fconfigs[[i]]$vars2col = "cluster"
  }
  fconfigs[[2]]$vars2col = "cluster.new"

  # CD8 colors
  # fconfigs[[1]]$col <- c("0" = "#c5f0e9", "1" = "#72e5d2", "7" = "#3bb8ae", "6" = "#ffa7d7", "8" = "#ff6ede", "2" = "#fcee62", "5" = "#9e89e8", "4" = "#4da7fd", "3" = "#ff971f",
  #  "9" = "#9c0006")
  fconfigs[[1]]$col <- c("0" = "#72e5d2", "1" = "#72e5d2", "7" = "#72e5d2", "6" = "#ffa7d7", "8" = "#ffa7d7", "2" = "#fcee62", "5" = "#9e89e8", "4" = "#4da7fd", "3" = "#ff971f",
   "9" = "#9c0006")
  # CD4 colors
  fconfigs[[2]]$col <- c("0" = "#a084ff", "1.1" = "#a084ff", "1.2" = "#97ccf5", "7" = "#97ccf5", "6" = "#97ccf5", "8" = "#97ccf5", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#a084ff", "3" = "#97ccf5",
   "9" = "#a084ff", "11" = "#9c0006", "10" = "#ff717c")
  # fconfigs[[2]]$col <- c("0" = "#c5b8ff", "1" = "#5fcfab", "1" = "#5fcfab", "7" = "#e57f08", "6" = "#ffaf4d", "8" = "#b35806", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#fcd093",
  #  "9" = "#5d1d88", "11" = "#9c0006", "10" = "#ff717c")

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
  # c("0" = "#a084ff", "1.1" = "#a084ff", "1.2" = "#97ccf5", "7" = "#97ccf5", "6" = "#97ccf5", "8" = "#97ccf5", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#a084ff", "3" = "#97ccf5",
  #  "9" = "#a084ff", "11" = "#9c0006", "10" = "#ff717c")
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

  ## Scatter plot (umap) for GEX + TCR + HTO clone size ## ----------------------------------

  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$sufix <- "cluster_clone_size_hto_"; fconfigs[[i]]$axis_x = "UMAP_1"; fconfigs[[i]]$axis_y = "UMAP_2"
    fconfigs[[i]]$vars2col = "cell_classification"; fconfigs[[i]]$metadata = paste0(basename(fconfigs[[i]]$result_id), "_mdata_hto")
  }

  sc_cd3p_cd4_mdata_hto <- sc_cd3p_cd4@meta.data %>% filter(!is.na(orig.donor))
  sc_cd3p_cd8_mdata_hto <- sc_cd3p_cd8@meta.data %>% filter(!is.na(orig.donor))

  sc_cd3p_cd8_mdata_hto$cell_classification[sc_cd3p_cd8_mdata_hto$cell_classification == "INNATE-LIKE_NKR"] <- "IL-NKR"
  sc_cd3p_cd8_mdata_hto$cell_classification[sc_cd3p_cd8_mdata_hto$cell_classification == "INNATE-LIKE_XCL1-2"] <- "IL-XCL1/2"

  # Create the "cluster_clones_size" column in the metadata of the 2 objects.
  clonotype_cd4 <- sc_cd3p_cd4_mdata_hto %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd4) <- clonotype_cd4$clonotype.tag; clonotype_cd4$clonotype.tag <- NULL

  clonotype_cd8 <- sc_cd3p_cd8_mdata_hto %>% filter(!is.na(clonotype.tag)) %>% select(clonotype.tag, cluster) %>% group_by(clonotype.tag, cluster) %>% summarize(clone_size = n()) %>% arrange(clonotype.tag, cluster) %>%
    pivot_wider(names_from = cluster, values_from = clone_size, values_fill = 0) %>% as.data.frame() #%>% head()
  rownames(clonotype_cd8) <- clonotype_cd8$clonotype.tag; clonotype_cd8$clonotype.tag <- NULL

  # Update CD4 mdata
  cluster_clone_size_hto <- c()
  for (i in 1:dim(sc_cd3p_cd4_mdata_hto)[1]){
    clonotype <- sc_cd3p_cd4_mdata_hto[i,"clonotype.tag"]
    c <- as.character(sc_cd3p_cd4_mdata_hto[i,"cluster"])
    cluster_clone_size_hto <- c(cluster_clone_size_hto, clonotype_cd4[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  sc_cd3p_cd4_mdata_hto$cluster_clone_size_hto <- as.numeric(cluster_clone_size_hto) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  sc_cd3p_cd4_mdata_hto <- sc_cd3p_cd4_mdata_hto %>% mutate(cluster_clone_size_hto = case_when(
    cluster_clone_size_hto >= 40 ~ 40,
    TRUE ~ cluster_clone_size_hto
  ))

  # Update CD8 mdata
  cluster_clone_size_hto <- c()
  for (i in 1:dim(sc_cd3p_cd8_mdata_hto)[1]){
    clonotype <- sc_cd3p_cd8_mdata_hto[i,"clonotype.tag"]
    c <- as.character(sc_cd3p_cd8_mdata_hto[i,"cluster"])
    cluster_clone_size_hto <- c(cluster_clone_size_hto, clonotype_cd8[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  sc_cd3p_cd8_mdata_hto$cluster_clone_size_hto <- as.numeric(cluster_clone_size_hto) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  sc_cd3p_cd8_mdata_hto <- sc_cd3p_cd8_mdata_hto %>% mutate(cluster_clone_size_hto = case_when(
    cluster_clone_size_hto >= 40 ~ 40,
    TRUE ~ cluster_clone_size_hto
  ))

  pp_clones = fig_plot_base(
    fconfigs[1], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size_hto)) +
        geom_point(aes(size = cluster_clone_size_hto), shape = 1, color = "gray45", alpha = 0.1, stroke = 0.4) + # gray47
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )
  pp_clones = fig_plot_base(
    fconfigs[2], return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size_hto)) +
        geom_point(aes(size = cluster_clone_size_hto), shape = 1, color = "gray20", alpha = 0.1, stroke = 0.4) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6) # c(1, 10, 20, 40, 50) # breaks = c(1, 5, 10, 20, 40), range = c(0.5, 6)
        guides(colour = guide_legend(override.aes = list(size = 6))) +
        labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
        # xlim(-10, 10) + ylim(-10, 10)
    }
  )

  ## Scatter plot (umap) for cluster clone size; Colored by gene expression ## ----------------------------------
  # NOTE: Add cluster_clone_size column by running code above
  dir.create("tcr/sc_cd3p_cd4/markers/LG-HG", recursive = T); dir.create("tcr/sc_cd3p_cd8/markers/LG-HG", recursive = T);

  # genes <- c("ZNF683", "FCGR3A", "FGFBP2", "IFNG", "PDCD1", "GZMA", "GZMK", "PRF1", "TNF", "CCL4", "HLA-DRB1", "LAG3", "TOX")
  genes <- c("ZNF683", "PDCD1")
  mdata_sc_cd3p_cd4 <- cbind(mdata_sc_cd3p_cd4[colnames(sc_cd3p_cd4@assays$RNA@data),], t(as.matrix(sc_cd3p_cd4@assays$RNA@data[genes,])))
  colnames(mdata_sc_cd3p_cd4) <- gsub("HLA-DRB1", "HLA_DRB1", colnames(mdata_sc_cd3p_cd4))
  mdata_sc_cd3p_cd8 <- cbind(mdata_sc_cd3p_cd8[colnames(sc_cd3p_cd8@assays$RNA@data),], t(as.matrix(sc_cd3p_cd8@assays$RNA@data[genes,])))
  colnames(mdata_sc_cd3p_cd8) <- gsub("HLA-DRB1", "HLA_DRB1", colnames(mdata_sc_cd3p_cd8))
  fconfigs_ <- rep(fconfigs, each = length(genes))
  for (i in 1:length(genes)){
    fconfigs_[[i]]$vars2col <-  gsub("-", "_", genes[i]) #NULL #
    fconfigs_[[length(genes) + i]]$vars2col <-  gsub("-", "_", genes[i]) #NULL #
  }

  # Pink = "#FFACC7" "#FF8DC7" "#FF597B"; "#FF597B" "#FA2FB5"
  # Turqouise = "#83c5be" "#00a896" "#006d77"
  # Red = "#d81159" "#8f2d56"
  # Blue = "#0899ba" "#1c558e" "#1d4e89"
  # Purple = "#9368b7" "#aa3e98"
  # Orange = "#fe621d" "#fd5200"
  palettes <- list("pink" = c("#FF8DC7", "#FF597B"), "turquoise" = c("#00a896", "#006d77"), "red" = c("#d81159", "#8f2d56"), "blue" = c("#0899ba", "#1c558e"), "purple" = c("#9368b7", "#aa3e98"), "orange" = c("#fe621d", "#fd5200"),
    "bright_blue" = c("gray", "blue", "blue"))

  for (fconfig in fconfigs_){ # fconfig <- fconfigs_[[2]]
    for (color in names(palettes)){ # color <- names(palettes)[1]; color <- "bright_blue"
      palette <- palettes[[color]]
      cat(fconfig$vars2col, "\n")
      df <- eval(parse(text=fconfig$metadata)) %>% arrange(eval(parse(text=fconfig$vars2col)))
      pdf(paste0(fconfig$result_id, "markers/", fconfig$sufix, fconfig$vars2col, "_", color, ".pdf"))
      p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
        geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col))), alpha = 0.8) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.45) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
        # scale_colour_gradient(colorRampPalette(c("red", "yellow", "white"))) +
        scale_colour_gradient2(low = "gray", mid = palette[1], high = palette[2], midpoint = 2.2) + # low = "steelblue1", mid = "yellow", high = "red"
        labs(x = "Dim 1", y = "Dim 2", size = "Clone\nSize", color = fconfig$vars2col)
      print(p)
      dev.off()
      pdf(paste0(fconfig$result_id, "markers/", fconfig$sufix, fconfig$vars2col, "_", color, "_blank.pdf"))
      p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
        geom_point(aes(size = cluster_clone_size, color = eval(parse(text=fconfig$vars2col)) ), alpha = 0.8, show.legend = FALSE) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.45, show.legend = FALSE) +
        scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) +
        # scale_colour_gradient(colorRampPalette(c("red", "yellow", "white"))) +
        scale_colour_gradient2(low = "gray", mid = palette[1], high = palette[2], midpoint = 2.2) +
        labs(x = "", y = "") + theme(axis.text.x=element_blank(), axis.text.y=element_blank() )
      print(p)
      dev.off()
      # Inset plot
      # NOTE: Only for cd8 for the moment
      # pdf(paste0(fconfig$result_id, "markers/PDCD1_", color, "_inset_plot.pdf"))
      # a <- DotPlot(object = sc_cd3p_cd8, features = c("PDCD1"), group.by = "cell_classification")
      # a <- a$data %>% select(features.plot, avg.exp.scaled, pct.exp, id); a$id <- factor(a$id, levels = c("Cell_Cycle", "GZMK_HI", "TRM", "CD16p_Effector", "Effector", "MAIT", "TCF7_HI"))
      # p <- a %>% ggplot(aes(x = id, y = features.plot, color = avg.exp.scaled )) +
      #   geom_point(aes(size = pct.exp), alpha = 0.8, show.legend = TRUE) +
      #   scale_radius(breaks = c(18, 19, 20, 21), range = c(2, 6)) +
      #   theme_classic() +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10), panel.background = element_blank()) +
      #   scale_colour_gradient(low = "gray", high = palette[2])
      # print(p)
      # dev.off()
    }
  }
  # Pink = "#FFACC7" "#FF8DC7" "#FF597B"; "#FF597B" "#FA2FB5"
  # Turqouise = "#83c5be" "#00a896" "#006d77"
  # Red = "#d81159" "#8f2d56"
  # Blue = "#0899ba" "#1c558e" "#1d4e89"
  # Purple = "#9368b7" "#aa3e98"
  # Orange = "#fe621d" "#fd5200"


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

  # Version 5   -- Change donor order

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


  # ------------------------------------------------------------------

  ################ - Adding cells with HTO but no Gex (Fig 1A)

  lib_pattern <- "3N|scripts|cd3x$"
  invert_grep <- TRUE
  dirs <- list.dirs("/home/ciro/large/pbtumor/results/ab_demux/cd3x")
  # dirs <- list.dirs("/home/kmlanderos/kmlanderos/pbtumor-all/results/ab_demux/cd3x")

  files <- paste0(dirs[grep(lib_pattern, dirs, invert = invert_grep)], "/step_6_object.rdata")

  HTO_matrixes_mdata <- lapply(files, function(file) {load(file); return(ht_object@meta.data)} )
  names(HTO_matrixes_mdata) <- basename(gsub("/step_6_object.rdata", "", files))

  # Get important HTO info
  tmp <- lapply(HTO_matrixes_mdata, function(x){
    tmp.var <- x %>% filter(in_gex == FALSE & !is.na(HT_ID) & !HT_ID %in% c("Doublet", "Negative")) %>%
      select(HT_ID, in_gex)
    tmp.var$cellname <- rownames(tmp.var)
    tmp.var$donor <- gsub("-TSC[0-9]+-C0[0-9]+$", "", tmp.var$HT_ID)
    # tmp.var$HT_ID <- NULL
    tmp.var
  })

  table_HTO <- Reduce(rbind, tmp)
  # Remove unwanted donors
  table_HTO <- table_HTO %>% filter(!donor %in% c("BT17_spinal", "BN8"))

  # Get TCR data
  # celltype <- "cd8"
  # for (celltype in c("cd4", "cd8")){
  table_summary_list <- lapply( c("cd4", "cd8"), function(celltype) {
    fname <- paste0("/home/kmlanderos/ad_hoc/pbtumor-all/results/tcr/CD45pCD3p_clean2_", celltype, "/filtered_contig_annotations_aggr.csv")
    contig_annotations_aggr <- read.csv(fname)
    cells.to.keep <- contig_annotations_aggr$is_cell=='True' & contig_annotations_aggr$high_confidence=='True' & contig_annotations_aggr$productive=='True'
    contig_annotations_aggr <- contig_annotations_aggr[cells.to.keep, ]

    table_TCR <- contig_annotations_aggr %>% select(barcode, raw_clonotype_id) %>% distinct()

    # ---> Merge HTO & TCR
    table_merge <- merge(table_HTO, table_TCR, by.x = "cellname", by.y = "barcode")
    # Remove shared clonotypes between CD4 and CD8
    table_merge <- table_merge %>% filter(!raw_clonotype_id %in% problematic_clonotypes) # problematic_clonotypes was obtained in secondary variables section

    # Only keep clonotypes exclusive from the celltype being analyzed
    text.tmp <- paste0("sc_cd3p_", celltype, "@meta.data$clonotype.tag")
    table_merge <- table_merge %>% filter(raw_clonotype_id %in% unique(eval(parse(text = text.tmp))))

    # Summarize cells by donor and clonotype
    table_summary <- table_merge %>% group_by(donor, raw_clonotype_id) %>% summarize(cells = n()) %>% arrange(desc(cells))
    table_summary
  }); names(table_summary_list) <- c("CD4", "CD8")

  # lapply(table_summary_list, head)
  #
  # table_summary <- table_summary_list[["CD4"]]
  # sum(table_summary$cells)
  # dim(table_summary)
  # length(unique(table_summary$raw_clonotype_id))
  # sum(table_summary$cells > 1)
  # unique(table_summary$donor)
  #
  # table_summary <- table_summary_list[["CD8"]]
  # sum(table_summary$cells)
  # dim(table_summary)
  # length(unique(table_summary$raw_clonotype_id))
  # sum(table_summary$cells > 1)
  # unique(table_summary$donor)

  #  Plot

  # Eliminate PDCD1_pct_byDonor column
  all.clonotypes.tmp <- lapply(all.clonotypes, function(x){x$PDCD1_pct_byDonor <- NULL; x})

  all.clonotypes.tmp <- lapply(c("CD4", "CD8"), function(celltype){
    merge.tmp <- merge(all.clonotypes.tmp[[celltype]], table_summary_list[[celltype]], by.x=c("donor", "clonotype"), by.y = c("donor", "raw_clonotype_id"), all.x = TRUE)
    merge.tmp$cell.count <- rowSums(cbind(merge.tmp$cell.count,merge.tmp$cells), na.rm=TRUE)
    all.clonotypes.tmp <- merge.tmp %>% select(donor, clonotype, cell.count)
    all.clonotypes.tmp
  }); names(all.clonotypes.tmp) <- c("CD4", "CD8")

  all.clonotypes.tmp <- lapply(all.clonotypes.tmp, function(x){
    x$donor <- factor(x$donor, levels = c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT5", "BT24", "BT25", "BT10", "BT15", "BT22", "BT11", "BT12", "BT18",
    "BT21", "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23"))
    x
  })

  # Back-to-back barplot
  pdf(paste0("tcr/cellCount_allClones_PBT_CD4_CD8_noColor_extended.pdf"))
  tmp.y.axis <- "CD4 cells"
  tmp.x.axis <- "Patient"
  p1 <- ggplot(all.clonotypes.tmp[["CD4"]], aes(x = donor, y = cell.count, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_reverse(limits = c(2000,0)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))

  tmp.y.axis <- "CD8 cells"
  p2 <- ggplot(all.clonotypes.tmp[["CD8"]], aes(x = donor, y = cell.count, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_continuous(limits = c(0,2000)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y = element_blank(), axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()

  # Save CD4 and CD8 XL tables
  write.csv(all.clonotypes.tmp[["CD4"]], file = "tcr/cellCount_allClones_PBT_CD4_extended.csv", row.names = FALSE)
  write.csv(all.clonotypes.tmp[["CD8"]], file = "tcr/cellCount_allClones_PBT_CD8_extended.csv", row.names = FALSE)

  # ----> Change scale
  # log
  pdf(paste0("tcr/cellCount_allClones_PBT_CD4_CD8_noColor_log_extended.pdf"))
  tmp.y.axis <- "CD4 cells"
  tmp.x.axis <- "Patient"
  p1 <- ggplot(all.clonotypes.tmp[["CD4"]] %>% mutate(cell.count.log2 = log2(cell.count)), aes(x = donor, y = cell.count.log2, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_reverse(limits = c(500,0)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))

  tmp.y.axis <- "CD8 cells"
  p2 <- ggplot(all.clonotypes.tmp[["CD8"]] %>% mutate(cell.count.log2 = log2(cell.count)), aes(x = donor, y = cell.count.log2, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_continuous(limits = c(0,500)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y = element_blank(), axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()

  # exp
  pdf(paste0("tcr/cellCount_allClones_PBT_CD4_CD8_noColor_exp_extended.pdf"))
  tmp.y.axis <- "CD4 cells"
  tmp.x.axis <- "Patient"
  p1 <- ggplot(all.clonotypes.tmp[["CD4"]] %>% mutate(cell.count.exp = exp(cell.count)), aes(x = donor, y = cell.count.exp, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_reverse(limits = c(500,0)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))

  tmp.y.axis <- "CD8 cells"
  p2 <- ggplot(all.clonotypes.tmp[["CD8"]] %>% mutate(cell.count.exp = exp(cell.count)), aes(x = donor, y = cell.count.exp, color = "gray35", group = -cell.count)) + geom_bar(stat = "identity", color = "gray60") +
    labs(x = tmp.x.axis, y = tmp.y.axis) +
    scale_y_continuous(limits = c(0,500)) + coord_flip() +
    theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y = element_blank(), axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 15.5, 5.5, -10), "pt"))

  grid.newpage()
  grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

  dev.off()

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


}

{ cat(redb("### GLIPH2 Heatmap ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

  date <- "2022-07-15" # date <- "2022-02-07" #
  HTO <- TRUE
  suffix <- c("CD8"); if(HTO) suffix <- paste0(suffix,"_hto")

  # Read donor HLA info. # NOTE: No HLA data for BT2 and BT9
  HLA_info <- fread(paste0('/home/kmlanderos/tmp_large/pbtumor-all/results/tcr/gliph/',date, "/HLA_info.csv"), header=F)
  colnames(HLA_info)<- c("orig.donor", "A1", "A2", "B1", "B2", "C1", "C2")
  # Get Number of patients with each HLA
  hla_perDonor <- HLA_info %>% pivot_longer(!orig.donor, names_to = "HLA_gene", values_to = "allele") %>% select(orig.donor, allele) %>% distinct() %>%
    group_by(allele) %>% summarize(count = n()) %>% arrange(desc(count)) %>% as.data.frame()
  rownames(hla_perDonor) <- hla_perDonor$allele

  # Read GLIPH Input
  inputs.path = paste0("/home/kmlanderos/large/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_", date)
  tmp.file.name <- paste0(inputs.path, '/GLIPH2Input_TCR-Data_All-TRB_', suffix, '.tsv')
  gliph_inputs <- lapply(tmp.file.name, function(x) {
    gliph.input <- fread(file=x, blank.lines.skip=TRUE)
    gliph.input <- merge(gliph.input, HLA_info, by.x = "donor.id", by.y="orig.donor", all.x=T)
  })


  # Read GLIPH Output
  results.path = paste0('/home/kmlanderos/tmp_large/pbtumor-all/results/tcr/gliph/',date)
  tmp.file.name <- paste0(results.path, '/GLIPH2Output_TCR-Data_', suffix, '.csv')
  gliph_outs <- lapply(tmp.file.name, function(x) {
    gliph.res <- fread(file=x, blank.lines.skip=TRUE)
    gliph.res <- gliph.res[pattern!='single']
    gliph.res$Sample <- gsub(":Exp1", "", gliph.res$Sample)
    gliph.res[expansion_score < 0.05,.(index, pattern,TcRb,Sample,Freq)] # Filter only Expanded specificity groups and keep variables of interest
  })
  names(gliph_outs) <- gsub("_hto","",suffix)


  # i = 1
  for(i in 1:length(gliph_outs)){
    x <- gliph_outs[[i]]
    ctype <- names(gliph_outs)[i]
    cat("\n\nCelltype: ", ctype, "\n")

    cat("Number of specificity groups before filtering: ", length(unique(x$index)), "\n")

    index_donor <- x[!Sample %in% c("BT2","BT9"),.(SG_cs = sum(Freq), nClonotypes = .N, orig.donor=unlist(str_split(Sample, ";"))), by=.(index,pattern)]
    # index_donor <- x[!grepl("BT2;|BT2$|BT9;|BT9$", Sample),.(SG_cs = sum(Freq), nClonotypes = .N, orig.donor=unlist(str_split(Sample, ";"))), by=.(index,pattern)] # Whenever a clonotype is present in >1 donor, we generate a row per each of the donor. # Remove TCRS from donors BT2 and BT9 due to missing HLA data
    index_donor <- index_donor[!orig.donor %in% c("BT2","BT9"),]
    tmp <- merge(index_donor, HLA_info, by="orig.donor", all.x=T) %>% # Add HLA data
            group_by(index, pattern) %>% filter(SG_cs > nClonotypes) %>% mutate(nDonor = sum(length(unique(orig.donor)))) %>% ungroup() %>% # Keep only clonally expanded SGs. # NOTE: Should not be necessary, but keeping it should affect the results
            unite("HLA", c("A1","A2","B1","B2","C1","C2"), remove = TRUE, sep = "_") # Then, substitute the donor for the HLA data of that donor

    # tmp.2 <- as.data.table(tmp)[,.(HLA=unlist(str_split(HLA, "_"))), by=.(index,pattern)]
    tmp.2 <- as.data.table(tmp)[,.(HLA=unlist(str_split(HLA, "_"))), by=.(index,pattern, orig.donor)] %>% distinct() %>% as.data.table() # We want to get the amount fo donors that have a certain HLA in the specificty group (it does not matter if it is present in both alleles)
    tmp.2 <- tmp.2 %>% group_by(index) %>% mutate(sgDonors = length(unique(orig.donor))) %>% as.data.table() # Get number of donors per SG
    tmp.2 <- tmp.2[, .(sgDonors, hlaDonors = .N), by = .(index,pattern,HLA)] %>% distinct() %>% as.data.table() # Get number of donors per HLA, per SG

    # Get number of sepecity groups after filtering for expansion
    # cat("Number of specificity groups after expansion filtering: ", length(unique(tmp.2$index)), "\n")

    # Hypergeometric tests
    hyper.test.res <- tmp.2 %>% arrange(index) %>%
            mutate(p_val = phyper(hlaDonors, hla_perDonor[HLA, "count"], length(HLA_info$orig.donor) - hla_perDonor[HLA, "count"], sgDonors, lower.tail = FALSE, log.p = FALSE)) %>%  # Apply Hypergeometric test for each HLA gene for each sG
            mutate(p_val = ifelse(p_val == 0, 3.1416*10**-30, p_val)) %>% # P-values of 0 are assigned a really low p-value instead. Should be eliminated . This happens when we keep SG coming from a single donor
            mutate(neg.log.pval = -log10(p_val)) %>% ungroup() %>% arrange(index,desc(neg.log.pval))
    hyper.test.res %>% arrange(desc(neg.log.pval)) %>% as.data.frame() %>% head(n=30)

    # Get the most enriched HLA gene per specificity group.
    # If >1 HLA genes are equally enriched, choose one randomly
    sg_order <- as.data.table(hyper.test.res)[,.(top_SG = HLA[which(p_val < 0.05 & hlaDonors == max(hlaDonors))]), by=index] %>% #%>% as.data.frame() %>% head(n=30)
     distinct(index, .keep_all= TRUE) %>% arrange(top_SG)

    # Creating matrix to plot.
    # Capping the p-value to 7. Also setting all non significant p-values to 0
    matrix <- hyper.test.res %>%
                mutate(neg.log.pval = ifelse(neg.log.pval > 7, 7, neg.log.pval)) %>%
                mutate(neg.log.pval = ifelse(neg.log.pval < 1.302, 0, neg.log.pval)) %>%
                select(index,HLA,neg.log.pval) %>% pivot_wider(names_from = HLA, values_from = neg.log.pval, values_fill = 0)
    indexes <- matrix$index; matrix$index <- NULL
    matrix <- as.matrix(matrix); rownames(matrix) <- indexes
    matrix <- matrix[as.character(sg_order$index),unique(as.character(sg_order$top_SG))]
    # matrix[1:10,1:10]

    # Keep the sepecity groups that show enrichment for at least one HLA allele.
    tmp.sig <- apply(matrix, MARGIN = 1, function(x){
      ifelse(any(x >= 1.302), TRUE, FALSE)
    })
    significant_matrix <- matrix[tmp.sig,]
    cat("Number of specificity groups after expansion filtering: ", nrow(significant_matrix), "\n")

    # Proportions of specificity groups that are enriched for a single HLA allele.
    tmp.mult.sig <- apply(significant_matrix, MARGIN = 1, function(x){
      ifelse(sum(x >= 1.302) == 1, TRUE, FALSE)
    })
    cat("Number of specificity groups enriched for a single HLA allele: ", sum(tmp.mult.sig), "\n")

    # Amount of SGs that are significant for each allele
    cat("Amount of SGs that are significant for each allele: \n")
    tmp.var <- colSums(significant_matrix >= 1.302); sort(tmp.var)
    cat(names(tmp.var), "\n")
    cat(tmp.var, "\n")

    # ---> Plot

    # Heatmap
    palette <- colorRampPalette(c("white", "red")) #, "darkred"
    palette_alleles <- pal_ucscgb("default")(ncol(significant_matrix)+1)[-1]; names(palette_alleles) <- colnames(significant_matrix)
    # table(as.matrix(HLA_info[,-1]))
    # table(as.matrix(HLA_info[,-1]))[colnames(significant_matrix)]
    donors_per_allele <- as.vector(table(as.matrix(HLA_info[,-1]))[colnames(significant_matrix)])

    pdf(paste0("/home/kmlanderos/tmp_large/pbtumor-all/results/figures/CD4_CD8_subclustering/tcr/gliph/HLA_heatmap_", ctype, ".pdf"))
    print(Heatmap(significant_matrix, name = "-log10(p-value)", #border = grid::gpar(col = "black", lty = 2),
      col = palette(20),
      column_title = "HLA enrichment in GLIPH2 Expanded Specificity Groups",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      top_annotation = HeatmapAnnotation(donors = anno_barplot(donors_per_allele, gp = gpar(fill = palette_alleles))),
      left_annotation = rowAnnotation(allele = sg_order$top_SG[tmp.sig], col = list(allele = palette_alleles)),
      width = unit(8, "cm")
    ))
    dev.off()
  }

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

}

{ cat(redb("### Heatmaps ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("heatmaps")

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

  # Iterate over the 4 combination of grade and expansion ("NonExp_LG", "NonExp_HG", "Exp_LG", "Exp_HG")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("NonExp_LG", "NonExp_HG", "Exp_LG", "Exp_HG")
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

  pdf(paste0("heatmaps/",'CD8_Expanded_DEGs_meanMatrix', thold,'.pdf'))
  pheatmap(mean_matrix,main = "Mean Expression Expanded DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         scale = "row",
         fontsize_number = 7)
  dev.off()

  pdf(paste0("heatmaps/",'CD8_Expanded_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix,main = "% Expression Expanded DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
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

  # Iterate over the 4 combination of grade and expansion ("NonExp_LG", "NonExp_HG", "Exp_LG", "Exp_HG")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("NonExp_LG", "NonExp_HG", "Exp_LG", "Exp_HG")
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

  pdf(paste0("heatmaps/",'CD4_Expanded_DEGs_meanMatrix', thold,'.pdf'))
  pheatmap(mean_matrix,main = "Mean Expression Expanded DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         scale = "row",
         fontsize_number = 7)
  dev.off()

  pdf(paste0("heatmaps/",'CD4_Expanded_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix,main = "% Expression Expanded DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        annotation_row = row_tags,
        scale = "row",
        fontsize_number = 7)
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

  # Expression Matrix
  # cells_df <- sc_cd3p_cd8@meta.data %>% select(cellname, cluster) %>% arrange(factor(cluster, levels = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9")))#; cells_df$cluster <- factor(cells_df$cluster, levels = c("0", "1", "7", "6", "8", "2", "5", "3", "4", "9"))
  # mat <- t(scale(t(as.matrix(sc_cd3p_cd8@assays$RNA@data[top10$gene, cells_df$cellname]))))
  # tmp <- ScaleData(sc_cd3p_cd8, features = top10$gene, vars.to.regress = c("nCount_RNA", "percent.mt"), block.size = 2000, scale.max = 6)
  # mat <- as.matrix(tmp@assays$RNA@scale.data[top10$gene, cells_df$cellname])

  # Color`
  # palettes <- c(pal1 = colorRampPalette(c("yellow", "red")), pal2 = colorRampPalette(c("purple", "yellow")), pal3 = colorRampPalette(c("blue", "yellow")), pal4 = colorRampPalette(c("blue", "red")), pal5 = colorRampPalette(c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')),
  #   pal6 = colorRampPalette(c("blue", "yellow", "red")), pal7 = colorRampPalette(c("blue", "white", "red"))
  # )
  # palette <- palettes[["pal7"]]
  # pheatmap(mat, fontsize_col = 7, main = "Cluster markers", annotation_col = cells_df %>% select(cluster), cluster_rows=F, cluster_cols=F, filename = paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"), color = palette(100), show_rownames = F, show_colnames = F)



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

  Idents(sc_cd3p_cd4) <- "cluster"
  df <- FindAllMarkers(sc_cd3p_cd4, test.use = "MAST")
  df %>% select(p_val, avg_logFC, pct.1, pct.2, p_val_adj, cluster, gene) %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
  top10 <- df %>%
    group_by(cluster) %>%
    top_n(n = n_features, wt = avg_logFC)
  write.table(top10, paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features.csv"), sep = ",", row.names = F, quote = F)

  pdf(paste0("/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/cluster_markers/cluster_markers_", celltype, "_top", as.character(n_features), "features_heatmap.pdf"))
  sc_cd3p_cd4$cluster <- factor(sc_cd3p_cd4$cluster, levels = c("0", "1.1", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11"))
  top10 <- top10 %>% arrange(factor(cluster, levels = c("1.1", "0", "2", "4", "9", "1.2", "3", "6", "7", "8", "5", "10", "11")))
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
  # c("0" = "#c5b8ff", "1.1" = "#d26dec", "1.2" = "#52edf4", "7" = "#379ee4", "6" = "#97ccf5", "8" = "#1069b9", "2" = "#a084ff", "5" = "#ffa7d7", "4" = "#8c50cf", "3" = "#c0f1f5",
  #  "9" = "#5d1d88", "10" = "#ff717c", "11" = "#9c0006")

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

{ cat(redb("### Clonality dot plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")=

  donor = "BT1" # donor = "BT5"
  donors <- unique(sc_cd3p_cd8@meta.data$orig.donor); donors <- donors[!is.na(donors)]

  for (donor in donors){
    cat("\nDonor: ", donor)

    # -----------> CD8
    df <- sc_cd3p_cd8@meta.data %>% filter(orig.donor == donor & !is.na(orig.donor) & !is.na(clonotype.tag)) %>% select(clonotype.tag, clon.size.tag) %>%
      distinct() %>% arrange(clon.size.tag)
    if(max(df$clon.size.tag) != 1){
      pdf(paste0("tcr/sc_cd3p_cd8/", donor, "_clones_dotplot.pdf"))
      net <- make_empty_graph(nrow(df))
      V(net)$size <- ifelse(df$clon.size.tag == 1, df$clon.size.tag, df$clon.size.tag+1.5) # df$clon.size.tag # Set node size based on clone size
      min.val=min(df$clon.size.tag); max.val=max(df$clon.size.tag); min_dotsize=1.5; max_dotsize=25
      V(net)$size <- ( (V(net)$size - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # Scale so that it goes from 2 to 10.
      plot(net, vertex.label=NA, vertex.color=c("royalblue4", "gray80")[1+(V(net)$size== 1.5)], vertex.frame.color	= c("white", "white")[1+(V(net)$size== 1.5)]) # c("none", "gray")[1+(V(net)$size== 1.5)] # "gray"
      dev.off()
    }
    # -----------> CD4
    df <- sc_cd3p_cd4@meta.data %>% filter(orig.donor == donor & !is.na(orig.donor) & !is.na(clonotype.tag)) %>% select(clonotype.tag, clon.size.tag) %>%
      distinct() %>% arrange(clon.size.tag)
    if(max(df$clon.size.tag) != 1){
      pdf(paste0("tcr/sc_cd3p_cd4/", donor, "_clones_dotplot.pdf"))
      net <- make_empty_graph(nrow(df))
      V(net)$size <- ifelse(df$clon.size.tag == 1, df$clon.size.tag, df$clon.size.tag+1.5) # df$clon.size.tag  # Set node size based on clone size
      min.val=min(df$clon.size.tag); max.val=max(df$clon.size.tag); min_dotsize=1.5; max_dotsize=25
      V(net)$size <- ( (V(net)$size - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # Scale so that it goes from 2 to 10.
      plot(net, vertex.label=NA, vertex.color=c("cadetblue4", "gray80")[1+(V(net)$size== 1.5)], vertex.frame.color	= c("white", "white")[1+(V(net)$size== 1.5)])
      dev.off()
    }
  }

}

{ cat(redb("### TCR Diversity Tables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  # ---  Generate tables V1 - (picking randomly 1 beta chain, when 2)

  # PBT

  donors <- c("BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT22", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT23", "BT13", "BT20", "BT17_brain",
  "BT2")

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd8_pbt_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/PBT_", donor, "_cd8_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd8.tsv"), sep = "\t", quote = F, row.names = F)


  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
   tmp <- sc_cd3p_cd4_pbt_mdata %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4") %>%
     select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
     filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
     ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
     rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
     mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
   colnames(tmp) <- gsub("count", "#count",  colnames(tmp))
   tmp$d <- "."
   tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
   # Write donor table
   fname = paste0("tcr/vdjtools/PBT_", donor, "_cd4_TCR_table.tsv")
   write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
   # Add info to metadata
   metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd4.tsv"), sep = "\t", quote = F, row.names = F)

 # LUNG

  donors <- sc_cd3p@meta.data %>% filter(orig.site == "lung") %>% pull(orig.donor) %>% unique()
  donors <- donors[!is.na(donors)]

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
   tmp <- sc_cd3p_cd8_lung_mdata %>%
     select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
     filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
     ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
     rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
     mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
   colnames(tmp) <- gsub("count", "#count", colnames(tmp))
   tmp$d <- "."
   tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
   # Write donor table
   fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd8_TCR_table.tsv")
   write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
   # Add info to metadata
   metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd8.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd4_lung_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd4_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd4.tsv"), sep = "\t", quote = F, row.names = F)


  # ---  Generate tables V2 - (picking randomly 1 beta chain, when 2)

  # PBT

  donors <- c("BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT22", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT23", "BT13", "BT20", "BT17_brain", "BT2")

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd8_pbt_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      filter(str_detect(cdr3aa, ";", negate = TRUE))
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/PBT_", donor, "_cd8_TCR_table_V2.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd8_V2.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd4_pbt_mdata %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4") %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      filter(str_detect(cdr3aa, ";", negate = TRUE))
    colnames(tmp) <- gsub("count", "#count",  colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/PBT_", donor, "_cd4_TCR_table_V2.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd4_V2.tsv"), sep = "\t", quote = F, row.names = F)

  # LUNG

  donors <- sc_cd3p@meta.data %>% filter(orig.site == "lung") %>% pull(orig.donor) %>% unique()
  donors <- donors[!is.na(donors)]

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd8_lung_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      filter(str_detect(cdr3aa, ";", negate = TRUE))
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd8_TCR_table_V2.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd8_V2.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd4_lung_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = TRB.nt.chains.tag, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
      filter(str_detect(cdr3aa, ";", negate = TRUE))
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd4_TCR_table_V2.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd4_V2.tsv"), sep = "\t", quote = F, row.names = F)
}

{ cat(redb("### TCR Diversity Plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/", recursive = T)

  # CD8
  # PBT_cd8 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/PBT_cd8.diversity.strict.exact.txt", sep = "\t")
  # LUNG_cd8 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/LUNG_cd8.diversity.strict.exact.txt", sep = "\t")
  PBT_cd8 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_cd8.diversity.strict.exact.txt", sep = "\t")
  LUNG_cd8 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/LUNG_cd8.diversity.strict.exact.txt", sep = "\t")
  tmp.1 <- PBT_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "brain")
  tmp.2 <- LUNG_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "lung")
  PBT_LUNG_cd8 <- rbind(tmp.1, tmp.2)

  # CD4
  # PBT_cd4 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/PBT_cd4.diversity.strict.exact.txt", sep = "\t")
  # LUNG_cd4 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/LUNG_cd4.diversity.strict.exact.txt", sep = "\t")
  PBT_cd4 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_cd4.diversity.strict.exact.txt", sep = "\t")
  LUNG_cd4 <- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/LUNG_cd4.diversity.strict.exact.txt", sep = "\t")
  tmp.1 <- PBT_cd4 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD4", tissue = "brain")
  tmp.2 <- LUNG_cd4 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD4", tissue = "lung")
  PBT_LUNG_cd4 <- rbind(tmp.1, tmp.2)

  df <- rbind(PBT_LUNG_cd8, PBT_LUNG_cd4)
  # ---> Write info
  # Shannon Index
  df_brain_cd4 <- df %>% filter(celltype == "CD4" & tissue == "brain" & statistic == "shannonWienerIndex_mean")
  df_brain_cd8 <- df %>% filter(celltype == "CD8" & tissue == "brain" & statistic == "shannonWienerIndex_mean")
  df_lung_cd4 <- df %>% filter(celltype == "CD4" & tissue == "lung" & statistic == "shannonWienerIndex_mean")
  df_lung_cd8 <- df %>% filter(celltype == "CD8" & tissue == "lung" & statistic == "shannonWienerIndex_mean")
  write.table(df_brain_cd4, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_CD4_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_brain_cd8, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_CD8_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd4, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/NSCLC_CD4_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd8, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/NSCLC_CD8_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)

  # inverseSimpson Index
  df_brain_cd4 <- df %>% filter(celltype == "CD4" & tissue == "brain" & statistic == "inverseSimpsonIndex_mean")
  df_brain_cd8 <- df %>% filter(celltype == "CD8" & tissue == "brain" & statistic == "inverseSimpsonIndex_mean")
  df_lung_cd4 <- df %>% filter(celltype == "CD4" & tissue == "lung" & statistic == "inverseSimpsonIndex_mean")
  df_lung_cd8 <- df %>% filter(celltype == "CD8" & tissue == "lung" & statistic == "inverseSimpsonIndex_mean")
  write.table(df_brain_cd4, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_CD4_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_brain_cd8, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_CD8_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd4, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/NSCLC_CD4_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd8, file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/NSCLC_CD8_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)

  # pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/shannon_diversity_mean.pdf")
  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/shannon_diversity_mean_v2.pdf")
  df1 <- df %>% filter(statistic == "shannonWienerIndex_mean") %>% group_by(tissue, statistic, celltype) %>% summarize(mean = mean(value), error = sd(value)) %>% ungroup()
  df1 %>% ggplot(aes(x = celltype, y = mean, fill = tissue)) + geom_col(position = "dodge") + geom_errorbar(aes(ymin = mean, ymax = mean+error, col= tissue),
   position = position_dodge(0.9), width = .3) + labs(x = "Celltype", y = "Shannon-Wiener Index", fill = "Tissue") #+
   # ggpubr::theme_pubclean() + scale_y_continuous(breaks = seq(50,250,50))
  dev.off()

  # pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/inverse_Simpson_diversity_mean.pdf")
  pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/inverse_Simpson_diversity_mean_v2.pdf")
  df2 <- df %>% filter(statistic == "inverseSimpsonIndex_mean") %>% group_by(tissue, statistic, celltype) %>% summarize(mean = mean(value), error = sd(value)) %>% ungroup() %>%
  df2 %>% ggplot(aes(x = celltype, y = mean, fill = tissue)) + geom_col(position = "dodge") + geom_errorbar(aes(ymin = mean, ymax = mean+error, col= tissue),
   position = position_dodge(0.9), width = .3) + labs(x = "Celltype", y = "Inverse Simpson Index", fill = "Tissue") #+
   # ggpubr::theme_pubclean() + scale_y_continuous(breaks = seq(50,250,50))
  dev.off()

}
