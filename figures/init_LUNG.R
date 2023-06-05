#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Kevin Meza Landeros
# Date: 2023-01-10
# ---

### ================== Figures ================== ###

setwd('/home/kmlanderos/tmp_large/pbtumor-all-Batch2'); here::here()

# Load final object
sc_cd3p <- readRDS("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_pbt_lung_Batch1_Batch2_seurat_object.rds")
# sc_cd3p_cd8_pbt_mdata <- readRDS(file = "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_cd3p_cd8_Batch1_Batch2_seurat_object_subset.rds")
# sc_cd3p_cd4_pbt_mdata <- readRDS(file = "/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/data/sc_cd3p_cd4_Batch1_Batch2_seurat_object_subset.rds")
sc_cd3p <- NormalizeData(sc_cd3p, normalization.method = "LogNormalize", scale.factor = 10000)

figures_path <- here::here("results", "figures")
if(!dir.exists(figures_path)) dir.create(figures_path)
setwd(figures_path);
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

}

{ cat(redb("### TCR Diversity Tables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # ---> Save tables for TCR Diversity Estimation (vdjtools)
  dir.create("tcr/vdjtools/")

  # ---  Generate tables V2 - (Eliminate Clonotypes with more than 1 beta chain)

  # PBT

  ##### CD8 #####
  donors <- c("BT8", "BT21", "BT20", "BT44", "BT12", "BT1", "BT9", "BT25", "BT23", "BT38", "BT3", "BT5", "BT7", "BT24", "BT19", "BT46", "BT4", "BT13", "BT35", "BT31", "BT30", "BT32", "BT37", "BT26", "BT27", "BT40", "BT17_brain", "BT15", "BT11", "BT18", "BT22", "BT2")
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p@meta.data %>% filter(celltype == "CD8" & orig.site == "brain") %>%
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
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/PBT_metadata_cd8_V2.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  donor_order <- c("BT5", "BT12", "BT26", "BT8", "BT25", "BT1", "BT21", "BT17_brain", "BT40", "BT46", "BT9", "BT15", "BT23", "BT3", "BT20", "BT24", "BT19", "BT27", "BT11", "BT7", "BT4", "BT32", "BT13", "BT18", "BT2", "BT22")
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p@meta.data %>% filter(celltype == "CD4" & orig.site == "brain") %>%
      select(clonotype.tag, orig.donor, TRA.nt.chains.tag, TRA.aa.chains.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, tra.v, tra.j, trb.v, trb.j) %>%
      filter(!is.na(clonotype.tag)) %>% filter(orig.donor == donor)
      if(nrow(tmp) == 0) next
    tmp <- tmp  %>% group_by(clonotype.tag, TRB.nt.chains.tag, TRB.aa.chains.tag, trb.v, trb.j) %>% summarize(count = n()) %>%
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
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/", fname)), ncol=2) )
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
    tmp <- tmp <- sc_cd3p@meta.data %>% filter(celltype == "CD8" & orig.site == "lung") %>%
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
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd8_V2.tsv"), sep = "\t", quote = F, row.names = F)

  ##### CD4 #####
  metadata_df <- c()
  # donor = "BT1"
  for (donor in donors){
    tmp <- sc_cd3p@meta.data %>% filter(celltype == "CD4" & orig.site == "lung") %>%
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
    metadata_df <- rbind(metadata_df, matrix(c(donor, paste0("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/", fname)), ncol=2) )
  }
  # Write metadata file
  colnames(metadata_df) <- c("sample.id", "#file.name")
  write.table(metadata_df[,c("#file.name", "sample.id")], paste0("tcr/vdjtools/LUNG_metadata_cd4_V2.tsv"), sep = "\t", quote = F, row.names = F)

}

{ cat(redb("### TCR Diversity Plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/", recursive = T)

  # CD8
  PBT_cd8 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_cd8.diversity.strict.resampled.txt", sep = "\t")
  LUNG_cd8 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/LUNG_cd8.diversity.strict.resampled.txt", sep = "\t")
  tmp.1 <- PBT_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "brain")
  tmp.2 <- LUNG_cd8 %>% select(sample_id, shannonWienerIndex_mean, inverseSimpsonIndex_mean) %>% pivot_longer(!sample_id, names_to = "statistic", values_to = "value") %>% mutate(celltype = "CD8", tissue = "lung")
  PBT_LUNG_cd8 <- rbind(tmp.1, tmp.2)

  # CD4
  PBT_cd4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/PBT_cd4.diversity.strict.resampled.txt", sep = "\t")
  LUNG_cd4 <- read.csv("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain/LUNG_cd4.diversity.strict.resampled.txt", sep = "\t")
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
  write.table(df_brain_cd4, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/PBT_CD4_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_brain_cd8, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/PBT_CD8_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd4, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/NSCLC_CD4_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd8, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/NSCLC_CD8_shannon_diversity_mean_v2.csv", sep = ",", row.names = F)

  # inverseSimpson Index
  df_brain_cd4 <- df %>% filter(celltype == "CD4" & tissue == "brain" & statistic == "inverseSimpsonIndex_mean")
  df_brain_cd8 <- df %>% filter(celltype == "CD8" & tissue == "brain" & statistic == "inverseSimpsonIndex_mean")
  df_lung_cd4 <- df %>% filter(celltype == "CD4" & tissue == "lung" & statistic == "inverseSimpsonIndex_mean")
  df_lung_cd8 <- df %>% filter(celltype == "CD8" & tissue == "lung" & statistic == "inverseSimpsonIndex_mean")
  write.table(df_brain_cd4, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/PBT_CD4_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_brain_cd8, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/PBT_CD8_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd4, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/NSCLC_CD4_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)
  write.table(df_lung_cd8, file = "tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/NSCLC_CD8_inverseSimpson_diversity_mean_v2.csv", sep = ",", row.names = F)

  # pdf("/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_randomBchain/shannon_diversity_mean.pdf")
  pdf("tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/shannon_diversity_mean_v2.pdf")
  df1 <- df %>% filter(statistic == "shannonWienerIndex_mean") %>% group_by(tissue, statistic, celltype) %>% summarize(mean = mean(value), error = sd(value)) %>% ungroup()
  df1 %>% ggplot(aes(x = celltype, y = mean, fill = tissue)) + geom_col(position = "dodge") + geom_errorbar(aes(ymin = mean, ymax = mean+error, col= tissue),
   position = position_dodge(0.9), width = .3) + labs(x = "Celltype", y = "Shannon-Wiener Index", fill = "Tissue")
  dev.off()
  pdf("tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/shannon_diversity_mean_dots_v2.pdf")
  df1 <- df %>% filter(statistic == "shannonWienerIndex_mean")
  df1$celltype_tissue <- paste(df1$celltype, df1$tissue)
  df1 <- df1 %>% filter(celltype == "CD8")
  df1 %>% ggplot(aes(x = celltype_tissue, y = value, fill = tissue)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
    stat_compare_means(paired = F)
  df1 <- df %>% filter(statistic == "shannonWienerIndex_mean")
  df1$celltype_tissue <- paste(df1$celltype, df1$tissue)
  df1 <- df1 %>% filter(celltype == "CD4")
  df1 %>% ggplot(aes(x = celltype_tissue, y = value, fill = tissue)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
    stat_compare_means(paired = F)
  dev.off()

  # pdf("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_randomBchain/inverse_Simpson_diversity_mean.pdf")
  pdf("tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/inverse_Simpson_diversity_mean_v2.pdf")
  df2 <- df %>% filter(statistic == "inverseSimpsonIndex_mean") %>% group_by(tissue, statistic, celltype) %>% summarize(mean = mean(value), error = sd(value)) %>% ungroup()
  df2 %>% ggplot(aes(x = celltype, y = mean, fill = tissue)) + geom_col(position = "dodge") + geom_errorbar(aes(ymin = mean, ymax = mean+error, col= tissue),
   position = position_dodge(0.9), width = .3) + labs(x = "Celltype", y = "Inverse Simpson Index", fill = "Tissue") #+
   # ggpubr::theme_pubclean() + scale_y_continuous(breaks = seq(50,250,50))
  dev.off()
  pdf("tcr/vdjtools/outs/noJgene_removeMultipleBchain/summary/inverse_Simpson_diversity_mean_dots_v2.pdf")
  df1 <- df %>% filter(statistic == "inverseSimpsonIndex_mean")
  df1$celltype_tissue <- paste(df1$celltype, df1$tissue)
  df1 <- df1 %>% filter(celltype == "CD8")
  df1 %>% ggplot(aes(x = celltype_tissue, y = value, fill = tissue)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
    stat_compare_means(paired = F)
  df1 <- df %>% filter(statistic == "inverseSimpsonIndex_mean")
  df1$celltype_tissue <- paste(df1$celltype, df1$tissue)
  df1 <- df1 %>% filter(celltype == "CD4")
  df1 %>% ggplot(aes(x = celltype_tissue, y = value, fill = tissue)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
    stat_compare_means(paired = F)
  dev.off()

 stat_compare_means(paired = F)

}

{ cat(redb("### GLIPH2 plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("tcr/gliph")

  # Version 5 - Back2back plot. Stacked barplot of EXPANDED CELLS. Each box is a different specificity group and the size is the cell within it. color Cd8 and Cd4 as in Fig1A

  donor_order <- rev(paste0("Hashtag", c("05","06","09","07","02","01","03","10","04","08")))
  # donor_order <- sc_cd3p@meta.data %>% filter(orig.site == "lung" & !is.na(orig.donor)) %>% pull(orig.donor) %>% unique()

  pdf(paste0("tcr/gliph/",'specificityGroups_per_Donor_CD4_CD8_Expanded.pdf'))
  # CD4
  tmp <- sc_cd3p@meta.data %>% filter(celltype == "CD4" & orig.site == "lung" & sGroup_tag == TRUE & !is.na(orig.donor) & expDegree == "Expanded") %>% select(orig.donor, pattern) %>% mutate(orig.donor = factor(orig.donor, levels = rev(donor_order)))
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
  source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater.R") # source("/home/ciro/scripts/handy_functions/devel/plots_crater.R")
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("craters", showWarnings = FALSE)

  edata =  expm1(sc_cd3p@assays$RNA@data)
  mygenes = grep("XIST|RPS4Y1|^RP11|^RP", rownames(edata), value = TRUE, invert = TRUE) # Filter unwanted genes
  mdata = sc_cd3p@meta.data

  # CD8

  f1 <- read.csv(here::here("results", "dgea/PBT_LUNG/CD8_Expansion_PBT/ExpandedvsNon_expanded/mastlog2cpm_results.csv"))
  f1 <- f1 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))

  write.csv(f1, here::here("results", "dgea/PBT_LUNG/CD8_Expansion_PBT/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"), row.names = F)

  f2 <- read.csv(here::here("results", "dgea/PBT_LUNG/CD8_Expansion_LUNG/ExpandedvsNon_expanded/mastlog2cpm_results.csv"))
  f2 <- f2 %>% mutate(padj = case_when(
    padj < 1*10**-100 ~ 1*10**-100,
    TRUE ~ padj
  ))
  write.csv(f2, here::here("results", "dgea/PBT_LUNG/CD8_Expansion_LUNG/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"), row.names = F)

  fconfigs = list(list(
      fnames = c(
        "PBT_Exp_vs_NonExp" = here::here("results", "dgea/PBT_LUNG/CD8_Expansion_PBT/ExpandedvsNon_expanded/.mastlog2cpm_results.csv"),
        "LUNG_Exp_vs_NonExp" = here::here("results", "dgea/PBT_LUNG/CD8_Expansion_LUNG/ExpandedvsNon_expanded/.mastlog2cpm_results.csv")
      ), result_id = "craters/PBT_LUNG_CD8_v4_", plot_squared = TRUE, #, columns = c("sex_disease")
      highlight_genes = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CCL3", "CCL4", "CCR1", "CCR3", "CCR5", "TNF", "IFNG", "ITGAE", "ZNF683", "CXCR6", "PDCD1", "LAG3", "HOPX", "FASLG",
        "GNLY", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "PSMB9", "PSME1", "PSME2", "CALR", "IL7R", "TCF7", "SELL", "LEF1", "CD55", "CD27", "KLF2", "MAL", "LTB", "CST7", "CD7")#,
      # selectss = list(c('orig.site', 'brain'), c('celltype', 'CD8')) # brain, lung; # CD4, CD8
  ))
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))

  for(fconfig in fconfigs){
    cat(c("- ", fconfig$result_id, "\n"), sep = "")
    for(fc in as.character(c(0.25, 0.35))){
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
    outputname = gsub("v4", "v5", fconfig$result_id),
    plot_interactive = TRUE,
    plot_squared = fconfig$plot_squared,
    limits_col = fconfig$limits_col,
    limits_size = fconfig$limits_size,
    verbose = TRUE,
    theme_extra = theme(axis.text.x = element_blank(),axis.text.y = element_blank()) #+ labs(x = NULL, y = NULL)
  )

}
