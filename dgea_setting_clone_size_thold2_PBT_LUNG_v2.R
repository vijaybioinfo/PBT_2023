#!/usr/bin/R

cat("====================== Setting DGEA comparisons ======================\n")

setwd("/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/clustering")
config_f = "/home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_clone_size_thold2_DICElung_PBT.yaml"
comp_f = sub("scripts", "info", sub("yaml", "csv", config_f))
config = yaml::read_yaml("/home/kmlanderos/scripts/dgea/config.yaml")

{ cat("========== Clone size comparison CD4/CD8 T cells\n") ## ---------------------------
  config$method = "mastlog2cpm"
  config$project = "clone_size"
  config$metadata = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_mdata.rds"
  config$expression_data = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_seurat_object.rds"
  config$output_dir = "/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/dgea"
  config$comparisons <- config$comparisons[1] # only the file
  config$comparisons$file = NULL
  config$job$main$ppn = 4
  config$job$main$nodes = 1

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["PBT_CD8_expansion_thold2"]] <- list(
    context = "PBT_CD8_expansion_thold2", test_column = "cloneSize_dgea",
    contrast = c("high", "low"), #c("low", "high"),
    filters = setNames(list(c("brain"),c("CD8")), c("orig.site", "orig.Cell.Type")),
    job = list(mem = "140gb")
  )

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["LUNG_CD8_expansion_thold2"]] <- list(
    context = "LUNG_CD8_expansion_thold2", test_column = "cloneSize_dgea",
    contrast = c("high", "low"), #c("low", "high"),
    filters = setNames(list(c("lung"),c("CD8")), c("orig.site", "orig.Cell.Type")),
    job = list(mem = "350gb")
  )

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["PBT_CD4_expansion_thold2"]] <- list(
   context = "PBT_CD4_expansion_thold2", test_column = "cloneSize_dgea",
   contrast = c("high", "low"), #c("low", "high"),
   filters = setNames(list(c("brain"),c("CD4")), c("orig.site", "orig.Cell.Type")),
   job = list(mem = "140gb")
  )

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["LUNG_CD4_expansion_thold2"]] <- list(
   context = "LUNG_CD4_expansion_thold2", test_column = "cloneSize_dgea",
   contrast = c("high", "low"), #c("low", "high"),
   filters = setNames(list(c("lung"),c("CD4")), c("orig.site", "orig.Cell.Type")),
   job = list(mem = "350gb")
  )

  cat("Creating metadata\n"); mdata = readRDS(config$metadata)
  mdata <- mdata[!is.na(mdata$clon.size.tag), c("orig.donor", "clon.size.tag", "orig.site", "orig.Cell.Type")]
  mdata$cloneSize_dgea <- ifelse(mdata$clon.size.tag > 1, "high", "low")
  # NOTE: This fille will be generated.
  config$metadata = metadata_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/cs_thold2_dgea_PBT_LUNG_mdata.rds"
  saveRDS(mdata, file = metadata_f)

  # Modifying expression_data object (seurat object's metadata).
  cat("Updating seurat object\n"); seurat_obj = readRDS(config$expression_data)
  seurat_obj@meta.data <- mdata
  # NOTE: This fille will be generated.
  config$expression_data = expression_data_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/figures/data/.object_stem_seurat_cs_thold2.rds"
  saveRDS(seurat_obj, file = expression_data_f)

  cat("Config:", config_f, "\n")
  yaml::write_yaml(config, file = config_f)
}

# RUN
# Rscript3 /home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_setting_clone_size_thold2_PBT_LUNG.R
# Rscript3 ~/scripts/dgea/R/dgea_jobs.R -s TRUE -y /home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_clone_size_thold2_DICElung_PBT.yaml

setwd("/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/clustering")
config_f = "/home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_clone_size_thold2_DICElung.yaml"
comp_f = sub("scripts", "info", sub("yaml", "csv", config_f))
config = yaml::read_yaml("/home/kmlanderos/scripts/dgea/config.yaml")

{ cat("========== Clone size comparison CD4/CD8 T cells\n") ## ---------------------------
  config$method = "mastlog2cpm"
  config$project = "clone_size"
  config$metadata = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_mdata.rds"
  config$expression_data = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/sc_cd3p_seurat_object.rds"
  config$output_dir = "/home/kmlanderos/ad_hoc/pbtumor-all/results/DICE_lung/dgea"
  config$comparisons <- config$comparisons[1] # only the file
  config$comparisons$file = NULL
  config$job$main$ppn = 4
  config$job$main$nodes = 1

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["LUNG_CD8_expansion_thold2"]] <- list(
    context = "LUNG_CD8_expansion_thold2", test_column = "cloneSize_dgea",
    contrast = c("high", "low"), #c("low", "high"),
    filters = setNames(list(c("lung"),c("CD8")), c("orig.site", "orig.Cell.Type")),
    job = list(mem = "350gb")
  )

  cat("---- Comparing cells w/cs > 1 vs cells w/cs < 1 \n") ## ----------------------------
  config$comparisons[["LUNG_CD4_expansion_thold2"]] <- list(
   context = "LUNG_CD4_expansion_thold2", test_column = "cloneSize_dgea",
   contrast = c("high", "low"), #c("low", "high"),
   filters = setNames(list(c("lung"),c("CD4")), c("orig.site", "orig.Cell.Type")),
   job = list(mem = "350gb")
  )

  cat("Creating metadata\n"); mdata = readRDS(config$metadata)
  mdata <- mdata[!is.na(mdata$clon.size.tag), c("orig.donor", "clon.size.tag", "orig.site", "orig.Cell.Type")]
  mdata$cloneSize_dgea <- ifelse(mdata$clon.size.tag > 1, "high", "low")
  mdata <- mdata[mdata$orig.site == "lung",]
  # NOTE: This fille will be generated.
  config$metadata = metadata_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/DICE_lung/figures/data/cs_thold2_dgea_LUNG_mdata.rds"
  saveRDS(mdata, file = metadata_f)

  # Modifying expression_data object (seurat object's metadata).
  cat("Updating seurat object\n"); seurat_obj = readRDS(config$expression_data)
  seurat_obj <- subset(x = seurat_obj, subset = orig.site == "lung")
  seurat_obj@meta.data <- mdata
  # NOTE: This fille will be generated.
  config$expression_data = expression_data_f = "/home/kmlanderos/kmlanderos/pbtumor-all/results/figures/data/.object_stem_seurat_LUNG_cs_thold2.rds"
  saveRDS(seurat_obj, file = expression_data_f)


# RUN
# Rscript3 /home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_setting_clone_size_thold2_PBT_LUNG_v2.R
# Rscript3 ~/scripts/dgea/R/dgea_jobs.R -s TRUE -y /home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_clone_size_thold2_DICElung_PBT.yaml
# Rscript3 ~/scripts/dgea/R/dgea_jobs.R -s TRUE -y /home/kmlanderos/pbtumor-all/scripts/DICE_lung/dgea_clone_size_thold2_DICElung.yaml
