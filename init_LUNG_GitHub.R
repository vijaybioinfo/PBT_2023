#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2023-01-10
# ---

### ================== Figures ================== ###

# source("/home/kmlanderos/pbtumor-all/scripts/DICE_lung/figures/object_lock.R")
# source("/home/kmlanderos/pbtumor-all/scripts/DICE_lung/figures/global.R")
# sc_cd3p <- readRDS("/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_filtered_PBT_LUNG_seurat_object.rds") # file: f6.r we create this object # Good one

# NOTE: Instead of loading global, load this object:
# sc_cd3p <- readRDS("/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_filtered_PBT_LUNG_seurat_object.rds")

# Load final object
sc_cd3p <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_seurat_object.rds")

setwd('/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/')
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
    "ggplot2", "cowplot", "patchwork", "Seurat", "stringr", "tidyverse", "data.table", "pheatmap", "fmsb", "grid", "ggpubr"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())

  # Merge the current metadata with the TCRmetadata (obtained through Vicente's pipeline)
  sc_cd3p_tcr = data.frame(data.table::rbindlist(lapply(
    X = list(c("/home/kmlanderos/large/pbtumor-all/results/DICE_lung/tcr/PBT/metadata.rds", "y$orig.site == 'brain'"),
      c("/home/kmlanderos/large/pbtumor-all/results/DICE_lung/tcr/LUNG_CD4/metadata.rds", "y$orig.site == 'lung' & y$orig.Cell.Type == 'CD4'"),
      c("/home/kmlanderos/large/pbtumor-all/results/DICE_lung/tcr/LUNG_CD8/metadata.rds", "y$orig.site == 'lung' & y$orig.Cell.Type == 'CD8'")),
    FUN = function(x){
      y <- readRDS(x[1]); y$cellname = rownames(y); y
      y <- y[eval(parse(text=x[2])),]
  })), stringsAsFactors = FALSE)
  rownames(sc_cd3p_tcr) <- sc_cd3p_tcr$cellname; sc_cd3p_tcr$cellname <- NULL
  sc_cd3p@meta.data <- joindf(sc_cd3p@meta.data, sc_cd3p_tcr)

  sc_cd3p@meta.data <- sc_cd3p@meta.data[!is.na(sc_cd3p@meta.data$orig.site),] # Sanity check
  sc_cd3p@meta.data$cellname <- rownames(sc_cd3p@meta.data)

  # Errase previous UMAP info (from past clusterng)
  sc_cd3p@meta.data$UMAP_1 <- NULL; sc_cd3p@meta.data$UMAP_2 <- NULL
  # Add the UMAP cell embeddings to the metadata
  sc_cd3p@meta.data <- joindf(sc_cd3p@meta.data,
    as.data.frame(sc_cd3p@reductions$umap@cell.embeddings))

  # ------------->
  # Add celltypes to PBT CD4 and CD8 cells.
  change_barcodes <- function(
    old_barcodes
  ){
    new_barcodes <- old_barcodes %>% str_replace_all(., "1$", "7") %>%
      str_replace_all(., "2$", "8") %>% str_replace_all(., "3$", "9") %>%
      str_replace_all(., "4$", "10") %>% str_replace_all(., "5$", "11") %>%
      str_replace_all(., "6$", "12")
    return(new_barcodes)
  }
  # PBT
  pbt_cd8_cell_celltype <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd8_mdata.rds") %>% select(cell_classification)
  pbt_cd4_cell_celltype <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd4_mdata.rds") %>% select(cell_classification)
  sc_cd3p@meta.data$cell_classification <- NA
  sc_cd3p@meta.data[change_barcodes(rownames(pbt_cd8_cell_celltype)), "cell_classification"] <- pbt_cd8_cell_celltype$cell_classification
  sc_cd3p@meta.data[change_barcodes(rownames(pbt_cd4_cell_celltype)), "cell_classification"] <- pbt_cd4_cell_celltype$cell_classification
  table(sc_cd3p@meta.data$orig.site, sc_cd3p@meta.data$cell_classification)

  # LUNG
  lung_mdata <- readRDS("/home/kmlanderos/large/pbtumor-all/results/DICE_lung/clustering/DICElung/.object_meta.data_seurat_mean0.01_pct25_pc20.rds")
  lung_treg_cells <- rownames(lung_mdata[lung_mdata$RNA_snn_res.0.4 == "4" & lung_mdata$orig.Cell.Type == "CD4",])
  sc_cd3p@meta.data[lung_treg_cells[lung_treg_cells %in% rownames(sc_cd3p@meta.data)], "cell_classification"] <- "TREG"
  table(sc_cd3p@meta.data$orig.site, sc_cd3p@meta.data$cell_classification)
  table(sc_cd3p@meta.data$cell_classification)

  # Add TREG_tag
  sc_cd3p@meta.data$TREG_tag <- FALSE
  sc_cd3p@meta.data[!is.na(sc_cd3p@meta.data$cell_classification) & sc_cd3p@meta.data$cell_classification == "TREG", "TREG_tag"] <- TRUE
  table(sc_cd3p@meta.data$TREG_tag)

  # # Assing orig.site metadata variable to the lung samples
  # sc_cd3p@meta.data[sc_cd3p@meta.data$orig.Project == "R24G", "orig.site"] <- "lung"
  # sc_cd3p@meta.data[sc_cd3p@meta.data$orig.Project %in% c("AdUp01","AdUp02","AdUp04"), "orig.site"] <- "brain"

  # ------------->
  # Let us recalculate the clone size of of cells with Gex
  sc_cd3p@meta.data$clon.size.tag.old <- sc_cd3p@meta.data$clon.size.tag # Save old values

  # Get the problematic clonotypes
  cd8_table <- table(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD8") %>% pull(clonotype.tag) ) %>% reshape2::melt() %>% arrange(desc(value))
  cd4_table <- table(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4") %>% pull(clonotype.tag) ) %>% reshape2::melt() %>% arrange(desc(value))


  # Add sGroup_tag
  merge_table <- merge(cd8_table, cd4_table, by = "Var1");
  problematic_clonotypes <- as.character(merge_table$Var1)

  # PBT - CD8
  clonotype.cs.cd8.pbt <- as.data.frame(table(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD8") %>% pull(clonotype.tag)) ) %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd8.pbt) <- c("clonotype", "clone_size")
  tmp <- unlist(apply(sc_cd3p@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"]) & x["orig.Cell.Type"] == "CD8" & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd8.pbt[clonotype.cs.cd8.pbt$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  sc_cd3p@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]
  # sc_cd3p@meta.data[,c("clonotype.tag", "clon.size.tag", "clon.size.tag")] %>% arrange(desc(clon.size.tag))

  # PBT - CD4
  clonotype.cs.cd4.pbt <- as.data.frame(table(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4") %>% pull(clonotype.tag)) ) %>% arrange(desc(Freq),Var1); colnames(clonotype.cs.cd4.pbt) <- c("clonotype", "clone_size")
  tmp <- unlist(apply(sc_cd3p@meta.data, MARGIN = 1, function(x){
    if(!is.na(x["clonotype.tag"]) & x["orig.Cell.Type"] == "CD4" & x["clonotype.tag"] %in% problematic_clonotypes) clonotype.cs.cd4.pbt[clonotype.cs.cd4.pbt$clonotype == x["clonotype.tag"],"clone_size"] else NA
  }))
  sc_cd3p@meta.data$clon.size.tag[!is.na(tmp)] <- tmp[!is.na(tmp)]

  # ------------->
  # Expansion Degree
  sc_cd3p@meta.data$expDegree <- sc_cd3p@meta.data$clon.size.tag
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree > 1] <- "Expanded"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree == 1] <- "Non_expanded"

  # --- Add specificity tag (LUNG ONLY)

  # CD4
  GLIPH2Input_CD4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/tcr/preproces_4_Gliph2/ag_specificity_assessment/reports/ag_specificity_assessment_2022-09-12/GLIPH2Input_TCR-Data_All-TRB_LUNG_CD4_hto_wAlphaChain.tsv", sep = "\t") # 2022-02-07
  GLIPH2Input_CD4 <- GLIPH2Input_CD4 %>% mutate(id = paste(GLIPH2Input_CD4$cdr3b.aa.se, GLIPH2Input_CD4$trb.v, GLIPH2Input_CD4$trb.j, sep = "_"))

  GLIPH2Output_CD4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/tcr/gliph/2022-09-12/GLIPH2Output_TCR-Data_LUNG_CD4_hto.csv") # 2022-02-07
  GLIPH2Output_CD4 <- GLIPH2Output_CD4 %>% mutate(id = paste(GLIPH2Output_CD4$TcRb, GLIPH2Output_CD4$V, GLIPH2Output_CD4$J, sep = "_"))
  GLIPH2Output_CD4 <- GLIPH2Output_CD4[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  # Add sGroup_tag
  merge <- merge(data.frame(GLIPH2Output_CD4), GLIPH2Input_CD4[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups_cd4 <- unlist(str_split(merge$clonotype.tag, ";"))
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
  tmp.meta.data.lung.cd4 <- merge(sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD4"), tmp.2, by = "clonotype.tag", all.x = TRUE); rownames(tmp.meta.data.lung.cd4) <- tmp.meta.data.lung.cd4$cellname
  tmp.meta.data.lung.cd4$sGroup_tag <- NULL

  # CD8
  GLIPH2Input_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_2022-07-15/GLIPH2Input_TCR-Data_All-TRB_CD8_hto_wAlphaChain.tsv", sep = "\t") # 2022-02-07
  GLIPH2Input_CD8 <- GLIPH2Input_CD8 %>% mutate(id = paste(GLIPH2Input_CD8$cdr3b.aa.se, GLIPH2Input_CD8$trb.v, GLIPH2Input_CD8$trb.j, sep = "_"))

  GLIPH2Output_CD8 <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/2022-07-15/GLIPH2Output_TCR-Data_CD8_hto.csv") # 2022-02-07
  GLIPH2Output_CD8 <- GLIPH2Output_CD8 %>% mutate(id = paste(GLIPH2Output_CD8$TcRb, GLIPH2Output_CD8$V, GLIPH2Output_CD8$J, sep = "_"))
  GLIPH2Output_CD8 <- GLIPH2Output_CD8[,c("id", "pattern")] %>% group_by(id) %>% summarize(pattern = paste(pattern, collapse = "_")) %>% mutate(sGroup = TRUE)

  # Add sGroup_tag
  merge <- merge(data.frame(GLIPH2Output_CD8), GLIPH2Input_CD8[,c("id", "clonotype.tag")], by = "id", all.x = TRUE)
  cts_sGroups_cd8 <- unlist(str_split(merge$clonotype.tag, ";"))

  # Add GLIPH2 pattern and n_sG
  tmp <- merge %>% mutate(clonotype.tag = as.character(clonotype.tag), n_sG = str_count(merge$pattern, "_")+1) %>% select(clonotype.tag, pattern, n_sG)
  tmp.2 <- c()
  for (i in 1:dim(tmp)[1]) {
    x <- tmp[i,]
    clonotypes <- unlist(str_split(x["clonotype.tag"], ";"))
    n_clonotypes <- length(clonotypes)
    for(clonotype in clonotypes){ tmp.2 <- rbind(tmp.2, c(clonotype, as.character(x[-1]))) }
  }; tmp.2 <- as.data.frame(tmp.2); colnames(tmp.2) <- colnames(tmp)

  tmp.meta.data.lung.cd8 <- merge(sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD8"), tmp.2, by = "clonotype.tag", all.x = TRUE); rownames(tmp.meta.data.lung.cd8) <- tmp.meta.data.lung.cd8$cellname

  # PUT TOGETHER EVERYTHING
  # Put toguether cd4, cd8 metadata and PBT metadata
  tmp.meta.data.brain <- sc_cd3p@meta.data %>% filter(orig.site == "brain"); tmp.meta.data.brain$pattern <- NA; tmp.meta.data.brain$n_sG <- NA
  sc_cd3p@meta.data <- rbind(tmp.meta.data.lung.cd4, tmp.meta.data.lung.cd8, tmp.meta.data.brain)[rownames(sc_cd3p@meta.data),]
  # colnames(tmp.meta.data.lung.cd4)[!colnames(tmp.meta.data.lung.cd4) %in% colnames(tmp.meta.data.lung.cd8)]

  # Add sGroup_tag
  sc_cd3p@meta.data <- sc_cd3p@meta.data %>% mutate(sGroup_tag = case_when(
    clonotype.tag %in% cts_sGroups_cd4 & orig.site == "lung" & orig.Cell.Type == "CD4" ~ TRUE,
    clonotype.tag %in% cts_sGroups_cd8 & orig.site == "lung" & orig.Cell.Type == "CD8" ~ TRUE,
    TRUE ~ FALSE
  ))# table(sc_cd3p_cd8@meta.data$sGroup_tag)

  # --- Add V and J genes information
  # Take the cells from the same clonotype and define V and J genes based on the UMI support

  # PBT
  cells.clons.info <- read.csv(file='/home/ciro/large/pbtumor/results/tcr/CD45pCD3p_clean2_cd8/filtered_contig_annotations_aggr.csv', stringsAsFactors=FALSE)
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
      ][count==max(count), paste0(unique(j_gene), collapse='|')],
      cdr3.aa=.SD[,
        .(count=.N),
        by=cdr3
      ][count==max(count), paste0(unique(cdr3), collapse='|')],
      cdr3.nt=.SD[,
        .(count=.N),
        by=cdr3_nt
      ][count==max(count), paste0(unique(cdr3_nt), collapse='|')]
    ),
    by=.(
      raw_clonotype_id,
      chain
    )
  ]
  # Amount of them that have multiple equally supported v and j gene
  cat("\n",dim(vdj.gene.info[str_detect(string=v.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported v gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=j.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported j gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'),])[1], "clonotypes have multiple equally supported CDR3 (alpha or beta).\n")
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'), cdr3.aa:=str_replace(string=cdr3.aa, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.nt, pattern='\\|'), cdr3.nt:=str_replace(string=cdr3.nt, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, cdr3.aa, cdr3.nt, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j', 'cdr3a.aa', 'cdr3a.nt'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j', 'cdr3b.aa', 'cdr3b.nt'), sep=',', convert=FALSE)
  # Merge VDJ data with object metadata
  sc_cd3p_cd8_pbt_mdata <- merge(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD8"), vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id"); rownames(sc_cd3p_cd8_pbt_mdata) <- sc_cd3p_cd8_pbt_mdata$cellname
  sc_cd3p_cd4_pbt_mdata <- merge(sc_cd3p@meta.data %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4"), vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id"); rownames(sc_cd3p_cd4_pbt_mdata) <- sc_cd3p_cd4_pbt_mdata$cellname

  # LUNG
  # CD8
  cells.clons.info <- read.csv(file='/home/kmlanderos/large/pbtumor-all/results/DICE_lung/tcr/LUNG_CD8/filtered_contig_annotations_aggr.csv', stringsAsFactors=FALSE)
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
      ][count==max(count), paste0(unique(j_gene), collapse='|')],
      cdr3.aa=.SD[,
        .(count=.N),
        by=cdr3
      ][count==max(count), paste0(unique(cdr3), collapse='|')],
      cdr3.nt=.SD[,
        .(count=.N),
        by=cdr3_nt
      ][count==max(count), paste0(unique(cdr3_nt), collapse='|')]
    ),
    by=.(
      raw_clonotype_id,
      chain
    )
  ]
  # Amount of them that have multiple equally supported v and j gene
  cat("\n",dim(vdj.gene.info[str_detect(string=v.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported v gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=j.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported j gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'),])[1], "clonotypes have multiple equally supported CDR3 (alpha or beta).\n")
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'), cdr3.aa:=str_replace(string=cdr3.aa, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.nt, pattern='\\|'), cdr3.nt:=str_replace(string=cdr3.nt, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, cdr3.aa, cdr3.nt, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j', 'cdr3a.aa', 'cdr3a.nt'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j', 'cdr3b.aa', 'cdr3b.nt'), sep=',', convert=FALSE)
  # Merge VDJ data with object metadata
  sc_cd3p_cd8_lung_mdata <- merge(sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD8"), vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id"); rownames(sc_cd3p_cd8_lung_mdata) <- sc_cd3p_cd8_lung_mdata$cellname

  # CD4
  cells.clons.info <- read.csv(file='/home/kmlanderos/large/pbtumor-all/results/DICE_lung/tcr/LUNG_CD4/filtered_contig_annotations_aggr.csv', stringsAsFactors=FALSE)
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
      ][count==max(count), paste0(unique(j_gene), collapse='|')],
      cdr3.aa=.SD[,
        .(count=.N),
        by=cdr3
      ][count==max(count), paste0(unique(cdr3), collapse='|')],
      cdr3.nt=.SD[,
        .(count=.N),
        by=cdr3_nt
      ][count==max(count), paste0(unique(cdr3_nt), collapse='|')]
    ),
    by=.(
      raw_clonotype_id,
      chain
    )
  ]
  # Amount of them that have multiple equally supported v and j gene
  cat("\n",dim(vdj.gene.info[str_detect(string=v.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported v gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=j.gene, pattern='\\|'),])[1], "clonotypes have multiple equally supported j gene (alpha or beta).\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'),])[1], "clonotypes have multiple equally supported CDR3 (alpha or beta).\n")
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.aa, pattern='\\|'), cdr3.aa:=str_replace(string=cdr3.aa, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=cdr3.nt, pattern='\\|'), cdr3.nt:=str_replace(string=cdr3.nt, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, cdr3.aa, cdr3.nt, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j', 'cdr3a.aa', 'cdr3a.nt'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j', 'cdr3b.aa', 'cdr3b.nt'), sep=',', convert=FALSE)
  # Merge VDJ data with object metadata
  sc_cd3p_cd4_lung_mdata <- merge(sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD4"), vdj.gene.info, all.x = TRUE, by.x = "clonotype.tag", by.y = "raw_clonotype_id"); rownames(sc_cd3p_cd4_lung_mdata) <- sc_cd3p_cd4_lung_mdata$cellname

  # PUT TOGETHER EVERYTHING
  sc_cd3p@meta.data <- rbind(sc_cd3p_cd8_pbt_mdata, sc_cd3p_cd4_pbt_mdata, sc_cd3p_cd8_lung_mdata, sc_cd3p_cd4_lung_mdata)[rownames(sc_cd3p@meta.data),]

  # ---- Save Files.

  # Just run one time
  # dir.create("data")
  # This comes from file: sc_cd3p_filtered_PBT_LUNG_seurat_object.rds; but the TCR data is present and the clone size for PBT has been done.
  # saveRDS(sc_cd3p@meta.data, file = "data/sc_cd3p_mdata.rds")
  # saveRDS(sc_cd3p, file = "data/sc_cd3p_seurat_object.rds")

}

{ cat(redb("### TCR Diversity Tables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # ---> Save tables for TCR Diversity Estimation (vdjtools)
  dir.create("tcr/vdjtools/")

  # ---  Generate tables V1 - (picking randomly 1 beta chain, when 2)

  # PBT

  donors <- c("BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT22", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT23", "BT13", "BT20", "BT17_brain", "BT2")

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p@meta.data %>% filter(orig.Cell.Type == "CD8" & orig.site == "brain") %>%
      select(clonotype.tag, orig.donor, cdr3b.nt, cdr3a.aa, cdr3b.aa, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = cdr3b.aa, v = trb.v, j = trb.j) %>% arrange(desc(count)) #%>%
      # mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/PBT_", donor, "_cd8_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd8.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd4_pbt_mdata %>% filter(orig.site == "brain" & orig.Cell.Type == "CD4") %>%
      select(clonotype.tag, orig.donor, cdr3b.nt, cdr3a.aa, cdr3b.aa, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = cdr3b.aa, v = trb.v, j = trb.j) %>% arrange(desc(count)) #%>%
      # mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count",  colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/PBT_", donor, "_cd4_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
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
      select(clonotype.tag, orig.donor, cdr3b.nt, cdr3a.aa, cdr3b.aa, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = cdr3b.aa, v = trb.v, j = trb.j) %>% arrange(desc(count)) #%>%
      # mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd8_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd8.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd4_lung_mdata %>%
      select(clonotype.tag, orig.donor, cdr3b.nt, cdr3a.aa, cdr3b.aa, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, cdr3b.aa, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = cdr3b.aa, v = trb.v, j = trb.j) %>% arrange(desc(count)) #%>%
      # mutate(cdr3nt = unlist(str_split(cdr3nt, ";"))[1], cdr3aa = unlist(str_split(cdr3nt, ";"))[1])
    colnames(tmp) <- gsub("count", "#count", colnames(tmp))
    tmp$d <- "."
    tmp <- tmp[, c("#count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j")]
    # Write donor table
    fname = paste0("tcr/vdjtools/LUNG_", donor, "_cd4_TCR_table.tsv")
    write.table(tmp, fname, sep = "\t", quote = F, row.names = F)
    # Add info to metadata
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd4.tsv"), sep = "\t", quote = F, row.names = F)


  # ---  Generate tables V2 - (Eliminate Clonotypes with more than 1 beta chain)

  # PBT

  donors <- c("BT1", "BT7", "BT19", "BT3", "BT27", "BT9", "BT4", "BT8", "BT5", "BT25", "BT24", "BT22", "BT15", "BT10", "BT26", "BT21", "BT11", "BT18", "BT12", "BT23", "BT13", "BT20", "BT17_brain", "BT2")

  ##### CD8 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p_cd8_pbt_mdata %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, cdr3b.nt, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
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
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, cdr3b.nt, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor) %>% group_by(clonotype.tag, cdr3b.nt, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
      ungroup() %>% mutate(freq = count/sum(count)) %>% select(count, freq, cdr3b.nt, TRB.aa.chains.tag, trb.v, trb.j) %>%
      rename(cdr3nt = cdr3b.nt, cdr3aa = TRB.aa.chains.tag, v = trb.v, j = trb.j) %>% arrange(desc(count)) %>%
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

{ cat(redb("### GLIPH2 plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("tcr/gliph")

  # Version 5 - Back2back plot. Stacked barplot of EXPANDED CELLS. Each box is a different specificity group and the size is the cell within it. color Cd8 and Cd4 as in Fig1A

  donor_order <- rev(paste0("Hashtag", c("05","06","09","07","02","01","03","10","04","08")))
  # donor_order <- sc_cd3p@meta.data %>% filter(orig.site == "lung" & !is.na(orig.donor)) %>% pull(orig.donor) %>% unique()

  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8_Expanded.pdf'))
  # CD4
  tmp <- sc_cd3p@meta.data %>% filter(orig.Cell.Type == "CD4" & orig.site == "lung" & sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
      scale_y_reverse(limits = c(600,0)) + coord_flip() +
      scale_x_discrete(drop=FALSE) +
      # scale_fill_manual(values = c(rep("white", sum(!grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype)))), rep("gray90", sum(grepl("clonotypes_nonExp",levels(cd4.matrix$donor_clonotype))))) ) + # scale_fill_gradient(low = "yellow", high = "red") +
      theme(panel.spacing.x = unit(0, "mm")) + theme_classic() +
      theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))
  # CD8
  tmp <- sc_cd3p@meta.data %>% filter(orig.Cell.Type == "CD8" & orig.site == "lung" & sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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

  ############################################################
  # NUMBER of Specificity groups & clones per donor (tables)
  ############################################################

  # CD4
  tmp <- sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD4" & sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  tmp <- sc_cd3p@meta.data %>% filter(orig.site == "lung" & orig.Cell.Type == "CD8" & sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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

}

{ cat(redb("### Crater plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater.R") # /home/ciro/scripts/handy_functions/devel/plots_crater.R
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("craters", showWarnings = FALSE)

  edata =  expm1(sc_cd3p@assays$RNA@data)
  mygenes = grep("XIST|RPS4Y1|^RP11|^RP", rownames(edata), value = TRUE, invert = TRUE) # Filter unwanted genes
  mdata = sc_cd3p@meta.data

  # CD8
  f1 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpandedvsNon_expanded/mastlog2cpm_results.csv")
  f1 <- f1 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f1, "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv", row.names = F)

  f2 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/ExpandedvsNon_expanded/mastlog2cpm_results.csv")
  f2 <- f2 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f2, "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv", row.names = F)

  # NOTE: Change Version if necessary
  fconfigs = list(list(
      fnames = c(
        "PBT_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv",
        "LUNG_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"
      ), result_id = "craters/PBT_LUNG_CD8_v4_", plot_squared = TRUE, #, columns = c("sex_disease")
      highlight_genes = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CCL3", "CCL4", "CCR1", "CCR3", "CCR5", "TNF", "IFNG", "ITGAE", "ZNF683", "CXCR6", "PDCD1", "LAG3", "HOPX", "FASLG",
        "GNLY", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "PSMB9", "PSME1", "PSME2", "CALR", "IL7R", "TCF7", "SELL", "LEF1", "CD55", "CD27", "KLF2", "MAL", "LTB", "CST7", "CD7 "),
      selectss = list(c('orig.site', 'brain'), c('orig.Cell.Type', 'CD8')) # brain, lung; # CD4, CD8
  )#,
  # list(fnames = c(
  #   "PBT_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv",
  #   "LUNG_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"
  # ), result_id = "craters/PBT_LUNG_CD8_v4_", plot_squared = TRUE, #, columns = c("sex_disease")
  # highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL"),
  # selectss = list(c('orig.site', 'lung'), c('orig.Cell.Type', 'CD8')) # brain, lung; # CD4, CD8
  # )
  )
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))

  for(fconfig in fconfigs){
    cat(c("- ", fconfig$result_id, "\n"), sep = "")
    for(fc in as.character(c(0.35))){
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
    }
  }

  # Final
  source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater_v2.R")
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
    outputname = fconfig$result_id,
    plot_interactive = TRUE,
    plot_squared = fconfig$plot_squared,
    limits_col = fconfig$limits_col,
    limits_size = fconfig$limits_size,
    verbose = TRUE,
    theme_extra = theme(axis.text.x = element_blank(),axis.text.y = element_blank()) #+ labs(x = NULL, y = NULL)
  )

  # CD4

  f1 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD4_expansion_thold2/ExpvsNonExp/mastlog2cpm_results.csv")
  f1 <- f1 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f1, "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv", row.names = F)

  f2 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD4_expansion_thold2/ExpvsNonExp/mastlog2cpm_results.csv")
  f2 <- f2 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f2, "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv", row.names = F)

  # NOTE: Change Version if necessary
  fconfigs = list(list(
      fnames = c(
        "PBT_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv",
        "LUNG_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv"
      ), result_id = "craters/PBT_LUNG_CD4_v3", plot_squared = TRUE, #, columns = c("sex_disease")
      highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL"),
      selectss = list(c('orig.site', 'brain'), c('orig.Cell.Type', 'CD4')) # brain, lung; # CD4, CD8
  ),
  list(
      fnames = c(
        "PBT_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv",
        "LUNG_Exp_vs_NonExp" = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD4_expansion_thold2/ExpvsNonExp/.mastlog2cpm_results.csv"
      ), result_id = "craters/PBT_LUNG_CD4_v4", plot_squared = TRUE, #, columns = c("sex_disease")
      highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL"),
      selectss = list(c('orig.site', 'lung'), c('orig.Cell.Type', 'CD4')) # brain, lung; # CD4, CD8
  )
  )
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))

  for(fconfig in fconfigs){
    cat(c("- ", fconfig$result_id, "\n"), sep = "")
    for(fc in as.character(c(0.35))){
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
    }
  }



}
