#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2022-05-22
# ---

# ls -holt /home/ciro/*/scripts/*tables*

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

  # function
  meanData <- function(data, column){
    mean <- data %>%
      group_by(eval(parse(text=column))) %>%
      summarise(mean = mean(value), error = sd(value)) %>% as.data.frame()
    colnames(mean) <- c(column, "mean", "error")
    return(mean)
  }

  # CD8
  PBT_cd8 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/PBT_cd8.diversity.strict.exact.txt", sep = "\t")
  LUNG_cd8 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/LUNG_cd8.diversity.strict.exact.txt", sep = "\t")
  tmp.1 <- PBT_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "brain")
  tmp.2 <- LUNG_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "lung")
  PBT_LUNG_cd8 <- rbind(tmp.1, tmp.2)

  # CD4
  PBT_cd4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/PBT_cd4.diversity.strict.exact.txt", sep = "\t")
  LUNG_cd4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/LUNG_cd4.diversity.strict.exact.txt", sep = "\t")
  tmp.1 <- PBT_cd4 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD4", tissue = "brain")
  tmp.2 <- LUNG_cd4 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD4", tissue = "lung")
  PBT_LUNG_cd4 <- rbind(tmp.1, tmp.2)

  df <- rbind(PBT_LUNG_cd8, PBT_LUNG_cd4)

  #  Plot

  pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/shannon_diversity_mean_v2.pdf")
  tmp.df <- df %>% filter(statistic == "shannonWienerIndex_mean") %>% mutate(celltype_tissue = paste(celltype, tissue, sep="_"))
  plot <- tmp.df %>% ggplot(aes(x = celltype_tissue, y = value, fill = celltype_tissue)) +
    geom_bar(stat="identity", position=position_dodge(.5), color = "black", data = meanData(tmp.df, "celltype_tissue"), aes(x=celltype_tissue, y=mean))  +
    geom_errorbar(data = meanData(tmp.df, "celltype_tissue"), aes(y = mean, ymin = mean, ymax = mean+error), color= "black", position = position_dodge(0.9), width = .3) +
    scale_fill_manual(values = c("cadetblue4", "cadetblue4", "royalblue4", "royalblue4"))+
    geom_jitter(width = 0.30, height = 0.5, color = 'black')+
    labs(title="Clonal diversity", x="Celltype", y = "Shannon-Wiener Index")+
    scale_y_continuous(breaks = seq(0,2000,500)) + # ggpubr::theme_pubclean() +
    theme(legend.position="none")

  print(plot)
  dev.off()

  # pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/inverse_Simpson_diversity_mean.pdf")
  pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/inverse_Simpson_diversity_mean_v2.pdf")
  tmp.df <- df %>% filter(statistic == "inverseSimpsonIndex_mean") %>% mutate(celltype_tissue = paste(celltype, tissue, sep="_"))
  plot <- tmp.df %>% ggplot(aes(x = celltype_tissue, y = value, fill = celltype_tissue)) +
    geom_bar(stat="identity", position=position_dodge(.5), color = "black", data = meanData(tmp.df, "celltype_tissue"), aes(x=celltype_tissue, y=mean))  +
    geom_errorbar(data = meanData(tmp.df, "celltype_tissue"), aes(y = mean, ymin = mean, ymax = mean+error), color= "black", position = position_dodge(0.9), width = .3) +
    scale_fill_manual(values = c("cadetblue4", "cadetblue4", "royalblue4", "royalblue4"))+
    geom_jitter(width = 0.30, height = 0.5, color = 'black')+
    labs(title="Clonal diversity", x="Celltype", y = "Inverse Simpson Index")+
    scale_y_continuous(breaks = seq(0,1000,250)) + # ggpubr::theme_pubclean() +
    theme(legend.position="none")

  print(plot)
  dev.off()

  # GRANT (CD8)
  pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/shannon_diversity_mean_CD8_v2.pdf")
  tmp.df <- df %>% filter(statistic == "shannonWienerIndex_mean" & celltype == "CD8")
  plot <- tmp.df %>% ggplot(aes(x = tissue, y = value, fill = tissue)) +
    geom_bar(stat="identity", position=position_dodge(.5), color = "black", data = meanData(tmp.df, "tissue"), aes(x=tissue, y=mean))  +
    geom_errorbar(data = meanData(tmp.df, "tissue"), aes(y = mean, ymin = mean, ymax = mean+error), color= "black", position = position_dodge(0.9), width = .3) +
    scale_fill_manual(values = c("chocolate1", "cadetblue1"))+
    geom_jitter(width = 0.30, height = 0.5, color = 'black')+
    labs(title="Clonal diversity", x="Celltype", y = "Shannon-Wiener Index")+
    scale_y_continuous(breaks = seq(0,500,100)) + # ggpubr::theme_pubclean() +
    theme(legend.position="none")

  print(plot)
  dev.off()

  # pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_randomBchain/inverse_Simpson_diversity_mean.pdf")
  pdf("/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/figures/tcr/vdjtools/outs/noJgene_single_ab/inverse_Simpson_diversity_mean_CD8_v2.pdf")
  tmp.df <- df %>% filter(statistic == "inverseSimpsonIndex_mean" & celltype == "CD8")
  plot <- tmp.df %>% ggplot(aes(x = tissue, y = value, fill = tissue)) +
    geom_bar(stat="identity", position=position_dodge(.5), color = "black", data = meanData(tmp.df, "tissue"), aes(x=tissue, y=mean))  +
    geom_errorbar(data = meanData(tmp.df, "tissue"), aes(y = mean, ymin = mean, ymax = mean+error), color= "black", position = position_dodge(0.9), width = .3) +
    scale_fill_manual(values = c("chocolate1", "cadetblue1"))+
    geom_jitter(width = 0.30, height = 0.5, color = 'black')+
    labs(title="Clonal diversity", x="Celltype", y = "Inverse Simpson Index")+
    scale_y_continuous(breaks = seq(0,250,50)) + # ggpubr::theme_pubclean() +
    theme(legend.position="none")

  print(plot)
  dev.off()



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

  # setwd(setwd("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures")
  dir.create("gsea")

  ## ---- GSEA ---- ###

  fname <- "/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/gsea/signatures_cd3p_human_v2.csv"
  slists <- readfile(fname, stringsAsFactors = FALSE)
  # Added CD4-CTL signature from Veena's paper (Patil3_COVID") to this file

  mdata_sc_cd3p <- sc_cd3p@meta.data

  # # Create PDCD1_tag to compare PDCD1p vs PDCD1n
  # PDCD1_tag_tmp <- as.matrix(sc_cd3p@assays$RNA@data["PDCD1",]) > 0
  # PDCD1_cells <- colnames(sc_cd3p@assays$RNA@data)
  # mdata_sc_cd3p[PDCD1_cells, "PDCD1_tag_tmp"] <- PDCD1_tag_tmp
  # mdata_sc_cd3p <- mdata_sc_cd3p %>% mutate(PDCD1_tag = case_when(
  #   PDCD1_tag_tmp > 0 ~ "PDCD1p",
  #   PDCD1_tag_tmp == 0 ~ "PDCD1n"
  # ))

  fconfigs = list(
    # ---> Neo TCR
    # Neo TCR - CD4 all; expansion comparison
    list(result_id = "gsea/LUNG_CD4_all_neoantigenTCR_expansion/",
      edata = "expm1(sc_cd3p@assays$RNA@data)",
      metadata = "mdata_sc_cd3p %>% filter(orig.site == 'lung')",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = list(c("orig.Cell.Type", "CD4"), c("TREG_tag", "FALSE")),
      comparisons = c("expDegree")
    ),
    # Neo TCR - CD8 all; expansion comparison
    list(result_id = "gsea/LUNG_CD8_all_neoantigenTCR_expansion/",
      edata = "expm1(sc_cd3p@assays$RNA@data)",
      metadata = "mdata_sc_cd3p %>% filter(orig.site == 'lung')",
      lists = lapply(X = slists[,c(41,42)], FUN = function(x) x[-c(1:2)] ),
      sample_filter = c("orig.Cell.Type", "CD8"),
      comparisons = c("expDegree")
    )
  )

  # ---> Neo TCR
  pp_gsea = fig_gsea(fconfigs[1], features = rownames(sc_cd3p), return_plot = TRUE, verbose = T)
  pp_gsea = fig_gsea(fconfigs[2], features = rownames(sc_cd3p), return_plot = TRUE, verbose = T)


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
     )
   )
  pp_curtains = fig_plot_curtain(fconfigs[7], dot.scale = 12, col.min = -1.5, col.max = 1.5)
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

  fconfigs = list(
    # FOR GRANT
    list(file = "/home/kmlanderos/tmp_large/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpandedvsNon_expanded//mastlog2cpm_results.csv",
      group1 = "Expanded", group2 = "Non_expanded",
      showgenes = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CCL3", "CCL4", "CCR1", "CCR3", "CCR5", "TNF", "IFNG", "ITGAE", "ZNF683", "CXCR6", "PDCD1", "LAG3", "HOPX", "FASLG", "GNLY", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "PSMB9", "PSME1", "PSME2", "CALR", "IL7R", "TCF7", "SELL", "LEF1", "CD55", "CD27", "KLF2", "MAL", "LTB", "CST7", "CD7")
    )
  )

  # fconfig=fconfigs[[1]]
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
      ))
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

{ cat(redb("### NK T-cells receptors ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

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
        geom_point(aes(size = clon.size.tag), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
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
  clonotype_f_cd4 = read.csv("/home/kmlanderos/large/pbtumor-all/results/figures/CD4_CD8_subclustering/supplementary_tables/single-cell_tcr_cd4_clonotypes.csv", header = TRUE)
  clonotype_f_cd4 = clonotype_f_cd4[,c(1, 8:19)]; rownames(clonotype_f_cd4) <- clonotype_f_cd4$Clonotype_ID
  clonotype_f_cd8 = read.csv("/home/kmlanderos/large/pbtumor-all/results/figures/CD4_CD8_subclustering/supplementary_tables/single-cell_tcr_cd8_clonotypes.csv", header = TRUE)
  clonotype_f_cd8 = clonotype_f_cd8[,c(1, 8:17)]; rownames(clonotype_f_cd8) <- clonotype_f_cd8$Clonotype_ID
  rownames(clonotype_f_cd4) <- clonotype_f_cd4$Clonotype.ID; clonotype_f_cd4$Clonotype.ID <- NULL
  rownames(clonotype_f_cd8) <- clonotype_f_cd8$Clonotype.ID; clonotype_f_cd8$Clonotype.ID <- NULL
  colnames(clonotype_f_cd4) <- gsub('X', '', colnames(clonotype_f_cd4)); colnames(clonotype_f_cd8) <- gsub('X', '', colnames(clonotype_f_cd8))

  # Update CD4 mdata
  cluster_clone_size <- c()
  for (i in 1:dim(mdata_sc_cd3p_cd4)[1]){
    clonotype <- mdata_sc_cd3p_cd4[i,"clonotype.tag"]
    c <- as.character(mdata_sc_cd3p_cd4[i,"cluster"])
    cluster_clone_size <- c(cluster_clone_size, clonotype_f_cd4[clonotype, c])
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
    cluster_clone_size <- c(cluster_clone_size, clonotype_f_cd8[clonotype, c])
  } # NOTE: Some cells do not have HT so, there will be NAs
  mdata_sc_cd3p_cd8$cluster_clone_size <- as.numeric(cluster_clone_size) # Add column to mdata
  # All the cells with cluster clone size greater than 40 will have the same dot size_feature
  mdata_sc_cd3p_cd8 <- mdata_sc_cd3p_cd8 %>% mutate(cluster_clone_size = case_when(
    cluster_clone_size >= 40 ~ 40,
    TRUE ~ cluster_clone_size
  ))

  pp_clones = fig_plot_base(
    fconfigs, return_plot = TRUE, verbose = TRUE,
    theme_extra = function(x){
      x + geom_point(aes(size = cluster_clone_size)) +
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
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
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
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
        geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
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

  # pp_clones = fig_plot_base(
  #   fconfigs, return_plot = TRUE, verbose = TRUE,
  #   theme_extra = function(x){
  #     x + geom_point(aes(size = cluster_clone_size)) +
  #       geom_point(aes(size = cluster_clone_size), shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
  #       scale_radius(breaks = c(1, 5, 10, 20, 40), range = c(0, 5)) + # breaks = scales::pretty_breaks(n=7,min.n=7), range = c(0, 6)
  #       guides(colour = guide_legend(override.aes = list(size = 6))) +
  #       labs(x = "Dim 1", y = "Dim 2", color = NULL, size = "Clone\nSize") # +
  #       # xlim(-10, 10) + ylim(-10, 10)
  #   }
  # )

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

{ cat(redb("### GG-paired plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
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


  # Plot
  pdf("heatmaps/paired_dots_PDCD1_expansion_pct.pdf", 15, 15)
  ggpaired(df, cond1 = "PDCD1n", cond2 = "PDCD1p", line.color = "gray", #color = "condition",
    fill = "condition", palette = "jco", ylab = "Percentage of Expansion")
  dev.off()

  # sc_cd3p_cd8_donor_mdata <- sc_cd3p_cd8@meta.data[intersect(PDCD1p_cells,donor_cells),]
  #
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

}

{ cat(redb("### Heatmaps DGEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dir.create("heatmaps")
  require(tidyverse)
  require(pheatmap)
  library(RColorBrewer)

  ###############################
  # PBT + LUNG  (CD8)
  ###############################

  # CD8 Heatmap - Highly Expanded DEGs
  # ---> 1. Get Genes

  # Change the metadata
  sc_cd3p@meta.data$expDegree <- sc_cd3p@meta.data$clon.size.tag
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree > 1 & sc_cd3p@meta.data$orig.site == "brain"]  <- "HighExp_brain"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree == 1 & sc_cd3p@meta.data$orig.site == "brain"] <- "LowExp_brain"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree > 1 & sc_cd3p@meta.data$orig.site == "lung"]  <- "HighExp_lung"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree == 1 & sc_cd3p@meta.data$orig.site == "lung"] <- "LowExp_lung"


  thold = "_thold2"
  genes_group1_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/ExpvsNonExp/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  genes_group2_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/ExpvsNonExp/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)

  genes_group1_group2_shared <- intersect(genes_group1_all, genes_group2_all)
  genes_group1_unique <- setdiff(genes_group1_all, genes_group1_group2_shared)
  genes_group2_unique <- setdiff(genes_group2_all, genes_group1_group2_shared)
  genes <- c(genes_group1_group2_shared, genes_group1_unique, genes_group2_unique)

  # Iterate over the 4 combination of grade and expansion ("LowExp_group1", "LowExp_group2", "HighExp_group1", "HighExp_group2")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("LowExp_brain", "HighExp_brain", "LowExp_lung", "HighExp_lung")
  for (i in tmp){
    # ---> 2. Get Cells
    cells <- rownames(sc_cd3p@meta.data[!is.na(sc_cd3p@meta.data$expDegree) & sc_cd3p@meta.data$expDegree == i,])
    edata <- sc_cd3p@assays$RNA@data[genes, cells]
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
  mean_matrix <- as.data.frame(round(mean_matrix, 2)); colnames(mean_matrix) <- c("NoExp_brain", "Exp_brain", "NoExp_lung", "Exp_lung")
  pct_matrix <- as.data.frame(round(pct_matrix, 2)); colnames(pct_matrix) <- c("NoExp_brain", "Exp_brain", "NoExp_lung", "Exp_lung")

  scale <- TRUE
  if (scale){
    mean_matrix_scaled <- as.data.frame(t(scale(t(mean_matrix))))
    pct_matrix_scaled <- as.data.frame(t(scale(t(pct_matrix))))
  }

  # ---> 3. Heatmap

  groups =c(rep(c("shared-DEGs", "brain-DEGs", "lung-DEGs"), times=c(length(genes_group1_group2_shared), length(genes_group1_unique), length(genes_group2_unique))) )

  pct_matrix_scaled$groups <- groups; # pct_matrix_scaled$gene <- rownames(pct_matrix_scaled)
  tmp1 <- pct_matrix_scaled[pct_matrix_scaled$groups == "shared-DEGs",]; tmp1 <- tmp1 %>% arrange(desc(Exp_brain))
  tmp2 <- pct_matrix_scaled[pct_matrix_scaled$groups == "brain-DEGs",]; tmp2 <- tmp2 %>% arrange(desc(Exp_brain))
  tmp3 <- pct_matrix_scaled[pct_matrix_scaled$groups == "lung-DEGs",]; tmp3 <- tmp3 %>% arrange(desc(Exp_lung))
  pct_matrix_scaled <- rbind(tmp1, tmp2, tmp3); pct_matrix_scaled$groups <- NULL

  mean_matrix_scaled$groups <- groups;
  tmp1 <- mean_matrix_scaled[mean_matrix_scaled$groups == "shared-DEGs",]; tmp1 <- tmp1 %>% arrange(desc(Exp_brain))
  tmp2 <- mean_matrix_scaled[mean_matrix_scaled$groups == "brain-DEGs",]; tmp2 <- tmp2 %>% arrange(desc(Exp_brain))
  tmp3 <- mean_matrix_scaled[mean_matrix_scaled$groups == "lung-DEGs",]; tmp3 <- tmp3 %>% arrange(desc(Exp_lung))
  mean_matrix_scaled <- rbind(tmp1, tmp2, tmp3); mean_matrix_scaled$groups <- NULL

  row_tags = data.frame("Group" = c(rep(c("shared-DEGs", "brain-DEGs", "lung-DEGs"), times=c(length(genes_group1_group2_shared), length(genes_group1_unique), length(genes_group2_unique))) )); rownames(row_tags) <- genes
  # colors <- colorRampPalette(c("purple", "black", "yellow"))(100)
  # colors <- colorRampPalette(c("blue", "red"))(100)
  colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)[16:85] # default: RdYlBu
  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_meanMatrix', thold,'.pdf'), 5)
  pheatmap(as.matrix(mean_matrix_scaled) ,main = "Mean Expression HighExp DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
         cellwidth = 30,
         color = colors,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         fontsize_number = 7 #, scale = "row",
       )
  dev.off()
  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_meanMatrix', thold,'_blank.pdf'),5)
  pheatmap(as.matrix(mean_matrix_scaled), cluster_cols = F, cluster_rows = F,
         cellwidth = 30,
         color = colors,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         legend = FALSE,
         show_colnames = FALSE,
         fontsize_number = 7 #, scale = "row",
       )
  dev.off()

  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix_scaled,main = "% Expression HighExp DEGs (CD8+)", cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        annotation_row = row_tags,
        fontsize_number = 7 #, scale = "row",
      )
  dev.off()
  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_pctMatrix', thold,'_blank.pdf'))
  pheatmap(pct_matrix_scaled, cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        legend = FALSE,
        show_colnames = FALSE,
        fontsize_number = 7 #, scale = "row",
      )
  dev.off()

  # ----> loxXlog co-expression

  thold = "_thold2"
  # Up-regulated genes un PBT and LUNG
  genes_group1_up <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  genes_group2_up <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  # Down-regulated genes un PBT and LUNG
  genes_group1_down <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange < -0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  genes_group2_down <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange < -0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)

  # Shared up- and down-regulated genes in both conditions
  genes_up <- intersect(genes_group1_up %>% pull(gene), genes_group2_up %>% pull(gene))
  genes_down <- intersect(genes_group1_down %>% pull(gene), genes_group2_down %>% pull(gene))

  table_genes_up <- merge(genes_group1_up[genes_group1_up$gene %in% genes_up,], genes_group2_up[genes_group2_up$gene %in% genes_up,], by = "gene", suffixes = c("_PBT", "_LUNG"))
  table_genes_down <- merge(genes_group1_down[genes_group1_down$gene %in% genes_down,], genes_group2_down[genes_group2_down$gene %in% genes_down,], by = "gene", suffixes = c("_PBT", "_LUNG"))
  table_genes <- rbind(table_genes_up, table_genes_down)

  pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_logXlog', thold,'.pdf'))
  p <- ggplot(table_genes, aes(x=log2FoldChange_PBT, y=log2FoldChange_LUNG)) +
  geom_point() +
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
    labs(x = "log2FC PBT", y = "log2FC LUNG") +
    geom_hline(yintercept = 0.25, linetype="dashed", color = "gray") +
    geom_hline(yintercept = -0.25, linetype="dashed", color = "gray") +
    geom_vline(xintercept = 0.25, linetype="dashed", color = "gray") +
    geom_vline(xintercept = -0.25, linetype="dashed", color = "gray") +
    # theme_classic()
    theme_light()
  print(p)
  dev.off()

  ###############################
  # PBT + LUNG  (CD4)
  ###############################

  # CD4 Heatmap - Highly Expanded DEGs
  # ---> 1. Get Genes

  # Change the metadata
  sc_cd3p@meta.data$expDegree <- sc_cd3p@meta.data$clon.size.tag
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree > 1 & sc_cd3p@meta.data$orig.site == "brain"]  <- "HighExp_brain"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree == 1 & sc_cd3p@meta.data$orig.site == "brain"] <- "LowExp_brain"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree > 1 & sc_cd3p@meta.data$orig.site == "lung"]  <- "HighExp_lung"
  sc_cd3p@meta.data$expDegree[sc_cd3p@meta.data$expDegree == 1 & sc_cd3p@meta.data$orig.site == "lung"] <- "LowExp_lung"


  thold = "_thold2"
  genes_group1_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD4_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  genes_group2_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD4_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% arrange(desc(log2FoldChange)) %>% pull(gene)

  genes_group1_group2_shared <- intersect(genes_group1_all, genes_group2_all)
  genes_group1_unique <- setdiff(genes_group1_all, genes_group1_group2_shared)
  genes_group2_unique <- setdiff(genes_group2_all, genes_group1_group2_shared)
  genes <- c(genes_group1_group2_shared, genes_group1_unique, genes_group2_unique)

  # Iterate over the 4 combination of grade and expansion ("LowExp_group1", "LowExp_group2", "HighExp_group1", "HighExp_group2")
  mean_matrix <- c()
  pct_matrix <- c()

  tmp <- c("LowExp_brain", "HighExp_brain", "LowExp_lung", "HighExp_lung")
  for (i in tmp){
    # ---> 2. Get Cells
    cells <- rownames(sc_cd3p@meta.data[!is.na(sc_cd3p@meta.data$expDegree) & sc_cd3p@meta.data$expDegree == i,])
    edata <- sc_cd3p@assays$RNA@data[genes, cells]
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
  mean_matrix <- as.data.frame(round(mean_matrix, 2)); colnames(mean_matrix) <- c("NoExp_brain", "Exp_brain", "NoExp_lung", "Exp_lung")
  pct_matrix <- as.data.frame(round(pct_matrix, 2)); colnames(pct_matrix) <- c("NoExp_brain", "Exp_brain", "NoExp_lung", "Exp_lung")

  scale <- TRUE
  if (scale){
    mean_matrix_scaled <- as.data.frame(t(scale(t(mean_matrix))))
    pct_matrix_scaled <- as.data.frame(t(scale(t(pct_matrix))))
  }

  # ---> 3. Heatmap

  groups =c(rep(c("shared-DEGs", "brain-DEGs", "lung-DEGs"), times=c(length(genes_group1_group2_shared), length(genes_group1_unique), length(genes_group2_unique))) )

  pct_matrix_scaled$groups <- groups; # pct_matrix_scaled$gene <- rownames(pct_matrix_scaled)
  tmp1 <- pct_matrix_scaled[pct_matrix_scaled$groups == "shared-DEGs",]; tmp1 <- tmp1 %>% arrange(desc(Exp_brain))
  tmp2 <- pct_matrix_scaled[pct_matrix_scaled$groups == "brain-DEGs",]; tmp2 <- tmp2 %>% arrange(desc(Exp_brain))
  tmp3 <- pct_matrix_scaled[pct_matrix_scaled$groups == "lung-DEGs",]; tmp3 <- tmp3 %>% arrange(desc(Exp_lung))
  pct_matrix_scaled <- rbind(tmp1, tmp2, tmp3); pct_matrix_scaled$groups <- NULL

  mean_matrix_scaled$groups <- groups;
  tmp1 <- mean_matrix_scaled[mean_matrix_scaled$groups == "shared-DEGs",]; tmp1 <- tmp1 %>% arrange(desc(Exp_brain))
  tmp2 <- mean_matrix_scaled[mean_matrix_scaled$groups == "brain-DEGs",]; tmp2 <- tmp2 %>% arrange(desc(Exp_brain))
  tmp3 <- mean_matrix_scaled[mean_matrix_scaled$groups == "lung-DEGs",]; tmp3 <- tmp3 %>% arrange(desc(Exp_lung))
  mean_matrix_scaled <- rbind(tmp1, tmp2, tmp3); mean_matrix_scaled$groups <- NULL

  row_tags = data.frame("Group" = c(rep(c("shared-DEGs", "brain-DEGs", "lung-DEGs"), times=c(length(genes_group1_group2_shared), length(genes_group1_unique), length(genes_group2_unique))) )); rownames(row_tags) <- genes
  colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)[16:85]
  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_meanMatrix', thold,'.pdf'), 5)
  pheatmap(as.matrix(mean_matrix_scaled) ,main = "Mean Expression HighExp DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
         cellwidth = 30,
         color = colors,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         annotation_row = row_tags,
         fontsize_number = 7 #, scale = "row",
       )
  dev.off()
  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_meanMatrix', thold,'_blank.pdf'),5)
  pheatmap(as.matrix(mean_matrix_scaled), cluster_cols = F, cluster_rows = F,
         cellwidth = 30,
         color = colors,
         display_numbers = FALSE,
         show_rownames = FALSE,
         number_color = "black",
         legend = FALSE,
         show_colnames = FALSE,
         fontsize_number = 7 #, scale = "row",
       )
  dev.off()

  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_pctMatrix', thold,'.pdf'))
  pheatmap(pct_matrix_scaled,main = "% Expression HighExp DEGs (CD4+)", cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        annotation_row = row_tags,
        fontsize_number = 7 #, scale = "row",
      )
  dev.off()
  pdf(paste0("heatmaps/",'CD4_HighExp_DEGs_pctMatrix', thold,'_blank.pdf'))
  pheatmap(pct_matrix_scaled, cluster_cols = F, cluster_rows = F,
        display_numbers = FALSE,
        show_rownames = FALSE,
        number_color = "black",
        legend = FALSE,
        show_colnames = FALSE,
        fontsize_number = 7 #, scale = "row",
      )
  dev.off()

  # ----> loxXlog co-expression

  # thold = "_thold2"
  # # Up-regulated genes un PBT and LUNG
  # genes_group1_up <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  # genes_group2_up <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  # # Down-regulated genes un PBT and LUNG
  # genes_group1_down <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/PBT_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange < -0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  # genes_group2_down <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/dgea/clone_size/LUNG_CD8_expansion_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange < -0.25) %>% arrange(desc(log2FoldChange)) %>% select(gene, log2FoldChange)
  #
  # # Shared up- and down-regulated genes in both conditions
  # genes_up <- intersect(genes_group1_up %>% pull(gene), genes_group2_up %>% pull(gene))
  # genes_down <- intersect(genes_group1_down %>% pull(gene), genes_group2_down %>% pull(gene))
  #
  # table_genes_up <- merge(genes_group1_up[genes_group1_up$gene %in% genes_up,], genes_group2_up[genes_group2_up$gene %in% genes_up,], by = "gene", suffixes = c("_PBT", "_LUNG"))
  # table_genes_down <- merge(genes_group1_down[genes_group1_down$gene %in% genes_down,], genes_group2_down[genes_group2_down$gene %in% genes_down,], by = "gene", suffixes = c("_PBT", "_LUNG"))
  # table_genes <- rbind(table_genes_up, table_genes_down)
  #
  # pdf(paste0("heatmaps/",'CD8_HighExp_DEGs_logXlog', thold,'.pdf'))
  # p <- ggplot(table_genes, aes(x=log2FoldChange_PBT, y=log2FoldChange_LUNG)) +
  # geom_point() +
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  #   labs(x = "log2FC PBT", y = "log2FC LUNG") +
  #   geom_hline(yintercept = 0.25, linetype="dashed", color = "gray") +
  #   geom_hline(yintercept = -0.25, linetype="dashed", color = "gray") +
  #   geom_vline(xintercept = 0.25, linetype="dashed", color = "gray") +
  #   geom_vline(xintercept = -0.25, linetype="dashed", color = "gray") +
  #   # theme_classic()
  #   theme_light()
  # print(p)
  # dev.off()


  ###############################
  # PBT
  ###############################

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
  genes_LG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size/LowGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
  genes_HG_all <- read.csv("/home/kmlanderos/kmlanderos/pbtumor-all/results/dgea/CD4_CD8_subclustering/clone_size/HighGrade_thold2/highvslow/mastlog2cpm_results.csv") %>% filter(padj < 0.05 & log2FoldChange > 0.25) %>% arrange(desc(log2FoldChange)) %>% pull(gene)

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
