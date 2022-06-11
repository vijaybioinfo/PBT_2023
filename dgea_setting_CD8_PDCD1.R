#!/usr/bin/R

cat("====================== Setting DGEA comparisons ======================\n")

config_f = "/home/kmlanderos/pbtumor-all/scripts/dgea_CD8_PDCD1.yaml"
comp_f = sub("scripts", "info", sub("yaml", "csv", config_f))
config = yaml::read_yaml("/home/kmlanderos/scripts/dgea/config.yaml")

{ cat("========== Comparison CD8+ T cells\n") ## ---------------------------
  config$method = "mastlog2cpm"
  config$project = "CD8_PDCD1"
  config$metadata = "/home/kmlanderos/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd8_mdata.rds"
  config$expression_data = "/home/kmlanderos/kmlanderos/pbtumor-all/results/clustering/CD45pCD3p_clean2/CD8/.object_stem_seurat_mean0.01_pct20.rds"
  config$output_dir = "/home/kmlanderos/ad_hoc/pbtumor-all/results/dgea/CD4_CD8_subclustering"
  config$comparisons <- config$comparisons[1] # only the file
  config$comparisons$file = NULL
  config$job$main$ppn = 4
  config$job$main$nodes = 1

  cat("---- Comparing cells PDCD1p vs cells PDCD1n \n") ## ----------------------------
  # config$comparisons[["LowGrade"]] <- list(
  #   context = "LowGrade", test_column = "tag_PDCD1",
  #   contrast = c("PDCD1p", "PDCD1n"),
  #   filters = setNames(list(c("BT1", "BT3", "BT4", "BT7", "BT8", "BT9", "BT19", "BT27", "BT24", "BT25", "BT5")), c("orig.donor"),
  #   job = list(mem = "90gb")
  # )
  # config$comparisons[["HighGrade"]] <- list(
  #   context = "HighGrade", test_column = "tag_PDCD1",
  #   contrast = c("PDCD1p", "PDCD1n"),
  #   filters = setNames(list(c("BT10", "BT15", "BT22", "BT11", "BT12", "BT18", "BT21", "BT26", "BT2", "BT13", "BT17_brain", "BT20", "BT23")), c("orig.donor"),
  #   job = list(mem = "90gb")
  # )
  config$comparisons[["All"]] <- list(
    context = "All", test_column = "tag_PDCD1",
    # filters = setNames(list(c("CD8")), c("celltype")),
    contrast = c("PDCD1p", "PDCD1n"),
    job = list(mem = "90gb")
  )

  cat("Creating metadata\n"); mdata = readRDS(config$metadata)
  mdata <- mdata[, c("orig.donor", "RNA_snn_res.0.2", "orig.Diagnosis.subclass", "celltype")]
  # NOTE: This file will be created
  config$metadata = metadata_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/sc_cd3p_cd8_mdata_PDCD1.rds"
  saveRDS(mdata, file = metadata_f)

  # Modifying expression_data object (seurat object's metadata).
  cat("Updating seurat object\n"); seurat_obj = readRDS(config$expression_data)
  seurat_obj@meta.data <- mdata
  # NOTE: This file will be generated.
  config$expression_data = expression_data_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/figures/CD4_CD8_subclustering/data/.object_stem_seurat_CD8_PDCD1_mean0.01_pct20.rds"
  saveRDS(seurat_obj, file = expression_data_f)

  cat("Config:", config_f, "\n")
  yaml::write_yaml(config, file = config_f)
}

# RUN
# Rscript3 /home/kmlanderos/pbtumor-all/scripts/dgea_setting_CD8_PDCD1.R
# Rscript3 ~/scripts/dgea/R/dgea_jobs.R -s TRUE -y /home/kmlanderos/pbtumor-all/scripts/dgea_CD8_PDCD1.yaml
