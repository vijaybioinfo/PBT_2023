#!/usr/bin/R

#####################
# Tables compendium #
#####################

# ---
# Author: Kevin Meza Landeros
# Date: 2022-01-17
# ---

# I need to find were is the metadata od the CD4 and CD8 or crate those files to read them here.

# This script will format supplementary tables.

source("/home/ciro/scripts/handy_functions/devel/supp_table.R")
source("/home/ciro/scripts/handy_functions/devel/utilities.R") # summarise_table
library(dplyr)

fig_dir = "/home/kmlanderos/ad_hoc/pbtumor-all/results/figures/CD4_CD8_subclustering"
if(!file.exists(fig_dir)) dir.create(fig_dir, showWarnings = FALSE);
setwd(fig_dir); cat("Working at (fig_dir):", fig_dir, "\n")
if(!dir.exists("supplementary_tables/")) dir.create("supplementary_tables/")
system("ls -loh supplementary_tables/*")

{ cat("### Single-cell TCR ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  configs = list(
    list(
      clonotype_f = "/home/kmlanderos/ad_hoc/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd4/clonotypes_aggr.csv",
      metadata = readRDS("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd4_mdata.rds") %>%
        mutate(
          clusters_m = dplyr::case_when(
            RNA_snn_res.0.4 %in% c(3, 6, 7, 8) ~ "TCM",
            TRUE ~ as.character(RNA_snn_res.0.4)
          )
        ),
      # cols_extra = c("orig.donor"),
      cols_extra = c("orig.donor", "clusters_m"),
      cols = "RNA_snn_res.0.4", name = "CD4 T cells.",
      # cols = "clusters_m", name = "CD4 T cells.",
      filename = "supplementary_tables/single-cell_tcr_cd4_clonotypes"
    ),
    list(
      clonotype_f = "/home/kmlanderos/ad_hoc/pbtumor-all/results/tcr/CD45pCD3p_clean2_cd8/clonotypes_aggr.csv",
      metadata = readRDS("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd8_mdata.rds") %>%
        mutate(
          clusters_m = dplyr::case_when(
            RNA_snn_res.0.2 %in% c(0, 1, 7) ~ "GZMK",
            RNA_snn_res.0.2 %in% c(6, 8) ~ "TRM",
            TRUE ~ as.character(RNA_snn_res.0.2)
          )
        ),
      # cols_extra = c("orig.donor"),
      cols_extra = c("orig.donor", "clusters_m"),
      cols = "RNA_snn_res.0.2", name = "CD8 T cells.",
      # cols = "clusters_m", name = "CD8 T cells.",
      filename = "supplementary_tables/single-cell_tcr_cd8_clonotypes"
    )
  )
  for (i in 1:2) {
    configs[[i+2]] <- configs[[i]]
    configs[[i+2]]$cols_extra = NULL
    configs[[i+2]]$cols = c("orig.donor", "clusters_m") # , "cell_classification")
    configs[[i+2]]$name = paste("Extended table", configs[[i+2]]$name)
    configs[[i+2]]$filename = paste0(configs[[i+2]]$filename, "_extended")
  }
  # config = configs[[1]]; config$cols = "clusters_m"
  for (config in configs[1:2]) {  #configs[1:2]
    cat(config$filename, "\n")
    clonotype_df <- read.csv(config$clonotype_f, stringsAsFactors = FALSE)
    tcr_df <- config$metadata
    tcr_df$clonotype_id <- tcr_df$clonotype.tag
    tcr_df <- tcr_df[!is.na(tcr_df$orig.donor), ] # removing TCR info data    # Remove info wo hashtag
    tcr_df$main_id <- ident_combine(tcr_df, config$cols, "~")
    if(!any(grepl("[A-z]", levels(tcr_df$main_id))))
      tcr_df$main_id <- paste0("C", tcr_df$main_id)
    ddf_sum <- summarise_table(
      tcr_df[, c("clonotype_id", "main_id", config$cols_extra)],
      "clonotype_id", expand = TRUE, verbose = TRUE)
    colnames(ddf_sum) <- sub("^total$", "clone_size_Gex", colnames(ddf_sum))
    tmp <- lapply(
      X = tcr_df[, c("main_id", config$cols_extra), drop = FALSE],
      FUN = function(x) gtools::mixedsort(names(table(x)))
    ); cols <- c("clonotype_id", "clone_size_Gex", unname(unlist(tmp)))
    ddf_sum <- ddf_sum[, colnames(ddf_sum) %in% cols]
    ddf_sum[is.na(ddf_sum)] <- 0
    ddf_sum$proportion_Gex <- ddf_sum$clone_size_Gex / sum(ddf_sum$clone_size_Gex)
    colnames(clonotype_df) <- sub("frequency", "clone_size_total", colnames(clonotype_df))
    colnames(clonotype_df) <- sub("proportion", "proportion_total", colnames(clonotype_df))
    ddf_final <- dplyr::left_join(
      ddf_sum, clonotype_df[, !grepl("samples", colnames(clonotype_df))], by = "clonotype_id")
    # saveRDS(ddf_final, file = paste0(gsub("supplementary_tables", "data", config$filename), ".rds"))
    headers_i <- list(
      none = c("clonotype_id$" = "TEXT"),
      `CDR3 Sequences` = c("cdr3s_" = "TEXT"),
      `Clonotype Metrics` = c("^clone_" = "NUMBER", "^prop" = "PERCENTAGE"),
      `Clone size per cluster` = c("^C[0-9]{1,}\\.?" = "NUMBER"),
      `Clone size per donor` = c("^BT[0-9]{1,}|brain$" = "NUMBER"),
      `Clone size per grouped clusters` = c("^[0-9]{1,}|TRM|TCM|GZMK\\.?" = "NUMBER"),
      `Clone size per donor-cluster` = c("~" = "NUMBER")
    )[if(grepl("Extended", config$name)) -c(4:6) else if("clusters_m" %in% config$cols_extra) -7 else -c(6:7)]
    # )[if(grepl("Extended", config$name)) -c(4:5) else -6]
    supp_table(
      mytables = setNames(list(ddf_final), config$name),
      headers = headers_i,
      title_name = "Clonotypes.",
      filename = config$filename,
      rename = list(c("pvalue", "P-value"), c("padj", "Adjusted P-value"),
        c("clonotype_id", "Clonotype ID"), c("proportion_", "Frequency "),
        c("clone_size_", "Clone size "), c("^C([0-9]{1,})$", "\\1"),
        c("cdr3s_aa", "Amino acid"), c("cdr3s_nt", "Nucleotide")),
      body_start_h = 20, round_num = 0,
      verbose = 1
    )
  }
}
#str(read_excel(paste0(config$filename, ".xlsx"), skip = 4))
