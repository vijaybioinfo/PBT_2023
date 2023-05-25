
library(cBioPortalData)
library(httr)
library(dplyr)
library(stringr)
library(ggplot2)
library(MultiAssayExperiment)
library(S4Vectors)
library(UpSetR)
suppressPackageStartupMessages({library(survival); library(survminer)})
library(stringr)

# CIBERSORT
apply_CIBERSORT <- FALSE
# CIBERSORT_file <- "~/Documents/CIBERSORTx_CBTN_TPM_absoluteMode_Results.csv" # "~/Documents/CIBERSORTx_CBTN_TPM_relativeMode_Results.csv", "~/Documents/CIBERSORTx_CBTN_TPM_absoluteMode_Results.csv"

# Filter by Grade
filter_gradde = T; grade <- "HG"
filter_diagnosis = T;
# diagnosis <- c("Medulloblastoma, group 4", "Medulloblastoma, SHH-activated", "Medulloblastoma, group 3", "Medulloblastoma", "Medulloblastoma, WNT-activated")
diagnosis <- c("High-grade glioma/astrocytoma, H3 G35-mutant", "Brainstem glioma- Diffuse intrinsic pontine glioma", "Diffuse midline glioma, H3 K28-mutant", "High-grade glioma/astrocytoma (WHO grade III/IV)", "High-grade glioma/astrocytoma, H3 wildtype", "Brainstem glioma- Diffuse intrinsic pontine glioma", "Diffuse midline glioma, H3 K28-mutant", "High-grade glioma/astrocytoma (WHO grade III/IV)", "High-grade glioma/astrocytoma, H3 wildtype")
# diagnosis <- c("Pilocytic astrocytoma, KIAA1549-BRAF", "Low-grade glioma/astrocytoma, KIAA1549-BRAF", "Low-grade glioma/astrocytoma, wildtype", "Pilocytic astrocytoma, wildtype", "Low-grade glioma/astrocytoma, BRAF V600E", "Ganglioglioma, BRAF V600E", "Pilocytic astrocytoma, FGFR", "Low-grade glioma/astrocytoma, RTK", "Ganglioglioma, other MAPK", "Low-grade glioma/astrocytoma, FGFR", "Ganglioglioma, RTK", "Pilocytic astrocytoma, RTK", "Low-grade glioma/astrocytoma, NF1-germline, FGFR", "Pilocytic astrocytoma, KIAA1549-BRAF, other MAPK", "Pilocytic astrocytoma", "Ganglioglioma, other MAPK, IDH", "Ganglioglioma, MYB/MYBL1", "Ganglioglioma, FGFR", "Ganglioglioma", "Ganglioglioma, BRAF V600E, CDKN2A/B", "Ganglioglioma, H3", "Ganglioglioma, NF1-germline", "Low-grade glioma/astrocytoma, NF1-somatic", "Pilocytic astrocytoma, KIAA1549-BRAF, NF1-germline", "Pilocytic astrocytoma, NF1-somatic", "Ganglioglioma, IDH", "Low-grade glioma/astrocytoma, MYB/MYBL1", "Pilocytic astrocytoma, NF1-germline", "Ganglioglioma, KIAA1549-BRAF", "Low-grade glioma/astrocytoma, BRAF V600E, CDKN2A/B", "Low-grade glioma/astrocytoma, other MAPK", "Low-grade glioma/astrocytoma, NF1-germline", "Pilocytic astrocytoma, BRAF V600E", "Low-grade glioma/astrocytoma (WHO grade I/II)", "Pilocytic astrocytoma, other MAPK", "Ganglioglioma, wildtype", "Pilocytic astrocytoma, wildtype")

# Split Cohort
percentile_down <- 0.30; percentile_up <- 0.70

# Clinical Data
mdata_cBioPortal <- read.csv("~/Documents/Clinicaldata_CBTN_noCellLines_wSurvData_noMetastasis_under20_noSpinal.tsv", sep = "\t")
mdata <- mdata_cBioPortal[,c("Patient_ID", "GENDER", "AGE", "OS_STATUS", "OS_MONTHS", "TUMOR_TYPE", "CANCER_TYPE_DETAILED", "pathology_free_text_diagnosis", "ANALYSIS")]; colnames(mdata_cBioPortal) <- c("Patient_ID", "Sex", "age_at_diagnosis", "OS_STATUS", "OS_MONTHS", "TUMOR_TYPE", "disease_type", "pathology_free_text_diagnosis") # , "Sample.ID"
mdata$VITAL_STATUS = ifelse(mdata$OS_STATUS == '0:LIVING', 0, 1)
mdata$Patient_ID <- make.names(mdata$Patient_ID, unique=TRUE)
row.names(mdata) <- mdata$Patient_ID
if(filter_gradde) mdata <- mdata %>% filter(ANALYSIS == grade)
if(filter_diagnosis) mdata <- mdata %>% filter(CANCER_TYPE_DETAILED %in% diagnosis)
table(mdata$VITAL_STATUS)

# Gex
GeneExpression <- read.csv("~/Documents/CBTN_tpm_matrix.tsv", sep = "\t");
# GeneExpression$gene_id <- sapply(str_split(GeneExpression$gene_id, "_"), function(x){x[2]})
dim(GeneExpression) # 58347   645
genes_names <- GeneExpression$gene_name; GeneExpression$gene_name <- NULL
# apply z-score of log2(TPM+1)
GeneExpression <- data.frame(t(scale(t(log2(GeneExpression+1))))); GeneExpression$gene_name <- genes_names

# subset
signature_list <- list(
  CD3 = c("CD3D", "CD3G", "CD3D"),
  Gide_PD1_response = c("BTN3A1", "TAP1", "UBD", "RP1", "SMTNL1", "LAG3", "CD8A", "DTHD1", "TBX21", "IFNG", "HERC6", "PSMB9", "RTL1", "CD8B", "HLA-F", "BATF2", "CCL4", "FRMPD3", "GBP5", "CEP128", "BTN3A3", "SLAMF7", "CXCL9", "STAT1", "DDX60", "LAP3", "UBE2L6", "NLRC3", "GABRA2", "HLA-A", "TIGIT", "CXCL13", "FCRL6", "OTOF", "IRF1", "SAMD3", "CCL5", "UBASH3A", "PDCD1", "SLFN12L", "NKG7", "CXCL10", "CD274", "IL12RB1", "BFSP2", "CD38", "SIRPG", "GBP1", "GPR98", "IL4I1", "CCDC136", "CD7", "GBP4", "IFIH1", "CD48", "CD244", "SLA2", "IDO1", "SIT1", "TRAJ42", "IKZF3", "IFI44L", "DBH-AS1", "PRKCQ", "PSMB10", "CCL3L3", "ZNFX1", "IL2RB", "PYHIN1", "GCH1", "NMI", "SP110", "MAGI1", "CA9", "TRAJ5", "ITGAL", "JAKMIP1", "SLFN14", "IPCEF1", "HERC5", "CCL3", "PLA2G2D", "IFIT5", "BTN3A2", "PSMB8", "PRF1", "OAS2", "ITGAD", "TRIM21", "WARS", "SYTL5", "PSMB8-AS1", "LCK", "PIK3CD-AS1", "GFI1", "CRTAM", "ZNF831", "PARP14", "IFIT3", "CST7", "RFX5", "TRAJ43", "ITGAE", "CCL4L1", "AKAP5", "GRIA2", "XAF1", "GPR171", "TOX", "CXCL11", "RBM38", "RLTPR", "NLRC5", "MAP4K1", "TMPRSS3", "AIM2", "STAT4", "SIDT1", "RNF212", "HLA-E", "TRGV10", "ARHGEF25", "TRAJ27", "SYNGR3", "NCR1", "TRBV7-6", "CD247", "DBH", "TMEM246", "CD6", "KRT18", "PARP11", "SP140", "TWIST1", "RBFOX2", "PPP1R16B", "DDX58", "LY9", "IGFBP2", "TRGC2", "TAP2", "TRAJ13", "CLEC4E", "P2RY8", "SLC27A2", "TRBV5-6", "PATL2", "CD2", "TRAF3IP3", "CD3E", "TRAJ18", "TRAJ38", "TNIP3", "OPRM1", "ZAP70", "SAMD9L", "TAGAP", "PIK3CD", "IL18RAP", "AADAC", "ASB2", "CMTM8", "WNT3", "TRAJ1", "GSG2", "TRAJ17", "GPR114", "TST", "PZP", "ANXA10", "GZMA", "RARRES3", "HSH2D", "DOCK8", "TAPBP", "SLFN5", "GDPD2", "IGLL5", "TRAJ44", "TRAJ16", "MX1", "OPHN1", "IFIT2", "DUSP27", "PSME1", "ARRDC5", "GPR174", "KCNA3", "SLAMF6", "TRAJ10", "TMC8", "APOL6", "ENPP5", "ICOS", "THEMIS", "TRAC", "TRAV8-3", "SH2D1A", "TNFRSF9", "NDST3", "CARD11", "TRAJ22", "LILRA4", "FAM26F", "PIM2", "C4orf50", "GZMK", "BTN2A1", "TNFAIP3", "GP5", "APOL2", "DHH", "EPSTI1", "RNF157", "ZBP1", "PTPN22", "P2RY10", "ITIH5", "BIRC3", "RASAL3", "TRAJ26", "TRAJ36", "FAM131B", "NLRP7", "FBXO6", "1-Sep", "LAX1", "SP100", "ITK", "SEMA4D", "DSPP", "AKNA", "GRIA4", "KLRD1", "TLR3", "ADAR", "TRAT1", "SPN", "KLRC4", "NUGGC", "DNAJC22", "PTPRCAP", "KIR2DL4", "ADORA2A", "TRBV28", "CNTN1", "TRAJ45", "LUZP2", "XIRP1", "FAM196B", "CLEC6A", "IL2RG", "KLRC4-KLRK1", "FAIM3", "CMPK2", "IRF8", "RNF157-AS1", "BST2", "TRBC2", "XRN1", "ACOXL", "KLHL6", "NUP43", "PTPRC", "C5orf56", "TRBV27", "ARHGAP9", "TRAJ14", "DOCK10", "SLAMF1", "ZNF101", "TRBV25-1", "CCR5", "TRAV4", "HAPLN3", "HLA-B", "OAF", "SELL", "PRSS3P2", "NDUFA4L2", "TRAJ2", "KIR3DX1", "CXCR3", "DUSP2", "PTCRA", "PSME2", "LINC01515", "CORO1A", "GRIK3", "ZNF683", "LYZ", "MYBPC1", "TRAJ12", "LAMB4", "TRIM69", "TNF", "KPNA5", "TRAJ52", "SPOCK2", "DENND1C", "TRAJ50", "CD79A", "EDIL3", "DAPK1", "OAS3", "NCF1", "ZDHHC18", "CYLD", "C1orf74", "CECR1"),
  NK_C2_Upadhye = c("FGFBP2", "GZMB", "GNLY", "FCGR3A", "NKG7", "EFHD2", "GPR56", "PRF1", "KLF2", "CX3CR1", "SPON2", "PLEK", "KLRD1", "S1PR5", "PRSS23", "FLNA", "LITAF", "GZMH", "TBX21", "C1orf21", "FGR", "S100A4", "ASCL2", "CST7", "TTC38", "ABHD17A", "KLF3", "EMP3", "MTSS1", "ARL4C", "AES", "PLAC8", "PXN", "SSBP3", "AHNAK", "ANXA1", "ITGB1", "TGFB1", "CLIC3", "S1PR1", "ITGB2", "RAP1B", "TNFRSF1B", "KLRG1", "FAM65B", "C12orf75", "SH3BP5", "CEP78", "RASGRP2", "CD247", "ADRB2", "LYAR", "SYNE2", "UCP2", "DSTN", "TPST2", "SUN2", "BHLHE40", "STK38", "SYNE1", "SELPLG", "RORA", "MYO1G", "PPP2R5C", "FCRL6", "MYO1F", "IQGAP2", "PRDM1", "IFITM2", "LPCAT1", "TSC22D3", "METRNL", "RHOB", "FBXW5", "YWHAQ", "KLRB1", "HOPX", "ADD3", "ICAM2", "MBP", "CYTH1", "SPN", "MYBL1", "SSBP4", "P2RY8", "RASA3", "CTSW", "SYTL3", "GZMM", "CDC25B", "SYTL1", "SH2D2A", "ADAM8", "RAB7L1", "TAGLN2", "SLC9A3R1", "CDC42SE1", "ITGB7", "LGALS1", "PTP4A2", "PTPN18", "GPR65", "PYHIN1", "GLUL", "CAPN2", "S100A10", "RASSF1", "MRPL10", "S1PR4", "IGF2R", "CD47", "BATF", "ITGAL", "C9orf142", "CCDC88C", "TPM4", "XBP1", "PRKCB", "CAPNS1", "FMNL1", "ATP2B4", "DOK2", "CDKN2D", "SRPK2", "RNF166", "NFKBIA"),# , "ZEB2"
  NK_CD8_CX3CR1_GUO= c("FGFBP2", "CX3CR1", "FCGR3A", "ADGRG1", "PLEK", "FCGR3B", "KLRD1", "S1PR1", "LITAF", "GZMH", "FCRL6", "GNLY", "KLRG1", "NKG7", "S1PR5", "PRF1", "PLAC8", "A2M", "FGR", "PXN", "CST7", "EFHD2", "PRSS23", "KLF2", "TGFBR3", "SYNE1", "FLNA", "BIN2", "PATL2", "LILRB1", "TBX21", "TARP", "MYO1F", "ITGB2", "UCP2", "C1orf21", "PYHIN1", "SPON2", "PZP", "ADRB2", "PXN-AS1", "LAIR2", "ITGAL", "APMAP", "CTSW", "HLA-C", "TTC38", "FCRL3", "EMP3", "ADD3", "ZEB2-AS1", "RAP1B", "CYTH1", "RASSF1", "TPST2", "MTSS1", "TSPAN32", "RAP2B", "RASA3", "ITGAM", "MYO1G", "CD47", "C12orf75", "SAMD3", "RAP1GAP2", "MYOM2", "PFN1", "SPN", "ITGA4", "SLC9A3R1", "LPCAT1", "USP28", "SH3BP5", "CD99", "SLC44A2", "SLAMF7", "PTP4A2", "CCND3", "RAB29", "PIK3R5", "CEP78", "ARHGAP25", "RAB37", "SH3BP5-AS1", "OSBPL5", "LINC00612", "CYBA", "SUN2", "PLEKHG3", "SYTL1", "TTC16", "LOC100130872", "CYTH4", "PRKCB", "VCL", "ITGB1", "CD52", "DOK2", "DIP2A", "PTPN4", "KLRF1"), #, "ZEB2"
  TRM_CLARKE = c("ITGAE","GSG2","GZMB","MYO7A","GPR25","LAYN","KRT86","STMN1","PDCD1","RBPJ","ENTPD1","ZNF683","SPRY1","KLRC2","ALOX5AP","TOX","TNS3","SRGAP3","CLNK","AFAP1L2","KLRC1","GOLIM4","CXCL13","AKAP5","HAVCR2","KIR2DL4","CD63","PHLDA1","CHRM3-AS2","ATP8B4","TNFRSF9","CSF1","RAB27A","DAPK2","FAM3C","ARHGAP11A","RRM2","ETV1","CD109","CD7","DBH-AS1","TMIGD2","SIRPG","DOCK5","CXCR6","CCL3","KIF2C","CD101","PAQR4","XCL1","CAPG","WBP4","IVNS1ABP","ADAM19","CTLA4","CLIC3","IFITM10","GEM","RAD51AP1","KIFC1","PTMS","UBASH3B","NUSAP1","TPX2","AURKA","KIF5C","VDR","SYNJ2","ATP10A","ANKRD35","KLRC3","SCCPDH","KIAA0101","CHN1","GTSE1","TTC24","ZC3H12C","SARDH","ZBED2","DPF3","ARHGEF12","CHEK1","SLC2A8","PLAGL1","LILRP2","UAP1L1","XCL2","FANCI","CDT1","TNFSF4","ABAT","AHI1","ASB2","DBN1","PTGIS","DFNB31","ABCB1","ATP10D","GCNT1","SPNS3","RUNX2","KLRB1","LINC00963","AMICA1","CAMK1","ANKS1B","TMEM200A","PGLYRP2","SUOX","KCNK5","DIXDC1","KIF14","SNAP47","SYNGR3","RBBP9","C1orf106","FASLG","RGS16","CCNA2","SLC16A6","CCRL2","GINS1","ACP5","BCL2L11","USP14","MIR155HG","TK1","BIRC5","TTYH3","CASC5","HELLS","NHS","CENPU","BRCA1","PDLIM7","C15orf53","INPP5F","APEX2","MCM2","HJURP","RDH10","FUT8","MKI67","RYBP","MYO1E","TOP2A","NDFIP2","MAD2L2","BRCA2","ITGA1","MAP3K6","KRT81","GALNT2","MPST","CEP41","ELK1","PPP1R21","DHFR","CCDC50","PDE7B","KIAA1524","ATL2","CARD6","AGPS","CDKN2C","CD200R1","SLC4A2","ACSL4","LOC101928988","CDC6","LRRN3","TET2","ZWINT","TNFRSF18","RHOC","PALB2","PMM2","EPSTI1","WIPF3","LINC00539","LIMK1","MAN1A1","SAC3D1","CKAP2L","NAIF1","AMZ1","ZDHHC18","SDHAP1","SSH1","GINS3","CASP9","ITGA2","MZB1","VCAM1","TIAM2","SLC27A2","XYLT1","ADAMTS17","IL18RAP","ACOT7","SNX9","GPA33","UHRF1","GMNN","CDCA3","RMND1","PDE4A","RIC1","CDKN3","KCTD9","EMC9","NUDT14","CD226","AZIN2","CDK1"),
  TRM_SAVAS = c("NBL1", "RP4-728D4.2", "LMO4", "VCAM1", "MLLT11", "SEMA4A", "CD244", "FASLG", "TNFSF4", "C1orf21", "RGS13", "NR5A2", "SNAP47", "GALNT2", "LYST", "AC092580.4", "FAM49A", "FAM179A", "CRIM1", "GNLY", "CD8A", "CD8B", "GPAT2", "NR4A2", "GPD2", "CERS6", "TTN", "CCDC141", "ITM2C", "UBE2F", "PDCD1", "SRGAP3", "RP11-222K16.2", "EOMES", "CMC1", "CCR1", "CCR5", "ABHD6", "DZIP3", "CD200R1", "GTPBP8", "HEG1", "GOLIM4", "GNB4", "CD38", "DTHD1", "STAP1", "TNIP3", "TMEM155", "ITGA1", "ITGA2", "GZMA", "PLPP1", "RASA1", "PDLIM4", "TIMD4", "HAVCR2", "DBN1", "GFOD1", "AIF1", "HLA-DRA", "HLA-DRB5", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "SOBP", "VNN2", "CHST12", "ETV1", "CHN2", "WIPF3", "TRGV9", "TRG-AS1", "LIMK1", "PON3", "PON2", "LRRN3", "FAM3C", "ADAM28", "TOX", "ASPH", "MSC", "FABP5", "MTSS1", "PIP5K1B", "RP11-305L7.1", "RP11-305L7.3", "SLC2A8", "RP11-492E3.2", "SLC2A6", "DBH-AS1", "CLIC3", "PNPLA7", "PRR5L", "PLA2G16", "CTSW", "TPCN2", "MYO7A", "PRSS23", "DIXDC1", "FXYD2", "CRTAM", "BARX2", "NEBL", "SPAG6", "PRF1", "ENTPD1", "PIK3AP1", "AFAP1L2", "CD9", "LPAR5", "PTMS", "LAG3", "CLECL1", "CLEC2B", "KLRD1", "RP11-277P12.20", "KLRK1", "KLRC4", "KLRC2", "KLRC1", "RP11-291B21.2", "NELL2", "KRT86", "KRT81", "ITGB7", "CD63", "IFNG", "TMTC2", "HVCN1", "KATNAL1", "SERP2", "TSC22D1", "SPRY2", "KDELC1", "GZMH", "GZMB", "AKAP5", "FUT8", "DPF3", "ATP8B4", "DAPK2", "PIF1", "SEMA7A", "MIR9-3HG", "MCTP2", "LRRC28", "GTF3C1", "NOD2", "CPNE2", "ADGRG5", "ADGRG1", "GSG2", "RASD1", "AC069363.1", "CCL3", "CCL4", "CCL3L3", "CCL4L2", "LINC00672", "RP11-357H14.17", "GNGT2", "ABI3", "SCPEP1", "PECAM1", "FBF1", "GALNT1", "BCAS4", "ZFP82", "LINC01480", "ATP1A3", "NKG7", "KIR2DL4", "DGCR6", "KIAA1671", "HMOX1", "APOBEC3H", "LINC00158", "MIR155HG"),
  TRM_UP_HOBIT = c("HSPD1", "QPCT", "CD69", "CDH1", "FOS", "LAD1", "VDAC1", "RHOB", "HPGDS", "NEURL3", "INPP4B", "ISG20", "SPSB1", "XCL1", "DUSP6", "EHD1", "P2RY10", "ITGAE", "PPP1R15A", "GPR34", "PNRC1", "P4HB", "HMGCS1", "RNF149", "SKIL", "HILPDA", "EGR1", "8430419L09RIK", "ABI3", "OSGIN1", "IFNG", "GADD45B", "LITAF", "DUSP5", "GLRX", "ZFP36", "TRAF4", "CSRNP1", "LY6G5B", "AMICA1", "GPR55", "EYA2", "SMIM3", "GPR171", "NR4A1", "FRMD4B", "ARRDC3", "DGAT1", "NR4A2", "RGS2", "FOSL2", "NFKBID", "GPR56", "JUNB", "SC4MOL", "LDLRAD4", "ATF3", "DUSP1", "PER1", "ODC1", "CTNNA1", "PYGL", "SIK1", "MAPKAPK3", "CISH", "IRF4", "HOBIT", "CD244", "B4GALNT4", "TNFAIP3", "STARD4", "PLK3", "NEDD4", "HSPA5", "JUN", "FOSB", "CXCR6", "BTG2", "KLF6", "PPP1R16B", "DHCR24", "DDX3X", "GSG2", "INSIG1"),
  GUO_GZMK_CD8 = c("GZMK","CCL4L1","ITM2C","CD74","CCL4","AOAH","CXCR4","DTHD1","CCL3L3","CCL3L1","CLDND1","CD44","SH2D1A","TRAT1","EOMES","CCL5","F2R","TC2N","FAM102A","PVRIG","PDE4DIP","CMC1","CRTAM","NEK7","EPHA1","WIPF1","MS4A1","CD84"),
  CD4CTLvsCD4TCM_Patil2018 = c("KLRD1","ASCL2","CMKLR1","CD8B","NMUR1","SPON2","GPR56","PRF1","CCL4","GZMB","NKG7","LGR6","PCDH1","FRMPD3","RP11-81H14.2","PTPN3","FASLG","FCGR3A","TRGC2","TRGV9","CST3","AC009951.1","EFHD2","CTSW","FGR","TTC38","FAM49A","TBX21","FGFBP2","OGN","CST7","TRGV2","KLRK1","METRN","EFNA5","PRSS23","GPR114","GNLY","GZMH","ERBB2","CD8A","MATK","CCL3","FCRL6","TTC16","PLEK","PTGDS","S1PR5","TYROBP","OSBPL5","FCGR3B","SLAMF7","TGFBR3","LILRB1","TRGC1","F2R","CD244","C1orf21","RASSF4","BZRAP1","SLC27A3","ENC1","SETBP1","USP28","CX3CR1","SH3RF2","RGS9","PROK2","ABI3","ZNF683","LAIR2","ARHGAP26","IFNG","PLOD1","FADS1","AKR1C3","ADAP1","CCL5","COLGALT2","C12orf75","SAMD3","C1orf177","CLIC3","SBK1","ST8SIA6","GZMA","TRGV8","SH2D1B","PRR5L","KLRG1","GAB3","PPP2R2B","PIK3AP1","SLC2A8","CHST12","C17orf66","SCD5","ADRB2","GPR141"), # ZEB2,
  TFH_Locci2013 = c("CXCR5","CXCL13","TOX2","ICA1","PDCD1","FAM43A","SCGB3A1","ASCL2","CHI3L2","CDK5R1","CEBPA","GNG4","CXXC5","GFOD1","TIGIT","POU2AF1","ATP9A","KIAA1671","SGPP2","BCL6","EGR2","LGMN","CD79A","KIAA1324","NFATC1","MYL6B","BCAT1","CTLA4","KCNK5","CTTN","LAG3","ICOS","THADA","LHFPL2","MAF","ID3","TRIB1","BTLA","SLC9A9","CAV1","SH2D1A","SH3TC1","MAGEH1","ST8SIA1","IKZF3","IGFBP4","NMB","TBC1D4","TRIM8","CCDC50","MYH10","PTPN11","FAM46C","ITM2A","FABP5","SPSB1","ODC1","MERTK","FBLN7","SERTAD2","TOX","CORO1B","SCD","SHISA2","RAB27A","FAM167A","STX11","LRRC1","PTTG1","CD200","KIAA0101","LIMS1","SLC7A5","PON2","FKBP5","NAB1","PFKFB3","DUSP6","AIM2","GADD45G","COL6A3","FAM179A","PHGDH","NUDT7","RNF19A","DUSP2","PAQR4","PPP1CC","SLC25A46","CARHSP1","HES6","H1F0","SRGN","CDCA7","SIRPG","MT1E","CD38","ASB2","KIFC1","RDH10","LRMP","DDIT4","P2RY11","NUCB2","TMEM2","ACTA2","SERPINE2","CASP9","INPP1","SGCE","MTUS1","GRAMD3","AFF3","NCALD","ANKRD55","AGMAT","TRIM32","CPA5","MAPK6","NDFIP2","UBE2E3","SEPN1","LIMS2","NETO2","GLCCI1","FAAH2","RPL39L","SSH2","PHACTR2","ERMP1","TNFRSF18","ANKS1B","SOCS1","FEN1","PYHIN1","PTPN13","QPRT","ACTN1","MYO5A","POMT1","LAT","SYT11","BCAS4","INSIG1","WSB2","RAB11FIP1","NR3C1","TRAF3IP2","GINS2","RFC5","FAM110A","RILPL2","IL2RB","ASAP1","GK","TSPAN5","ALDH5A1","TIAM1","CBLB","CLINT1","SMPDL3A","SIPA1L2","JARID2","ORMDL3","TSHR","BAZ2B","ATP6V1D","ZFPM1","TM2D3","RAB37","GPR19","F2R","SRPK2","BIK","FBXO33","MAPKAPK3","PTPN2","LBH","PTPN7","UBE2G1","FHL2","TJP2","IFNAR2","VOPP1","IRF4","CDK6","HIF1A","P2RX5","CXCR4","TMEM99","FYN","DCUN1D3","BATF","CDT1","GFI1","ARHGAP10","PRKCZ","IL6R","LTA","KPNA2","PPP2R5C","PPP1R16B","VWA5A","FAM160B1","TMEM64","MBOAT1","FAM134B","TP53INP1","RYR1","AKAP13","PCNA","HIST1H2BD","HES4","NAP1L4","PLCH2","GAPDH","HECW2","TRIB2","SAT1","SLC29A1","ZFP91","RACGAP1","MCM2","SFXN1","IL21"),
  NeoTCR8 = c("ATP10D","GZMB","ENTPD1","KIR2DL4","LAYN","HTRA1","CD70","CXCR6","HMOX1","ADGRG1","LRRN3","ACP5","CTSW","GALNT2","LINC01480","CARS","LAG3","TOX","PTPRCAP","ASB2","ITGB7","PTMS","CD8A","GPR68","NSMCE1","ABI3","SLC1A4","PLEKHF1","CD8B","LINC01871","CCL4","NKG7","CLIC3","NDFIP2","PLPP1","PCED1B","CXCL13","PDCD1","PRF1","HLA-DMA","GPR25","CD9","TIGIT","HLA-DRB5","SYTL3","SLF1","NEK1","CASP1","SMC4","TSEN54","PLSCR1","GNPTAB","HLA-DPB1","PLEKHA1","ARHGAP9","ALOX5AP","SH3BP1","NCF4","NELL2","GATA3","PPM1M","TNFRSF1A","AC022706.1","MCM5","HLA-DRB1","TNFSF10","TRIM21","HDLBP","ERN1","CALHM2","SASH3","ACTA2","MAST4","CAPG","MPST","IGFLR1","GZMA","CD27","ITGAE","SLA2","RHOC","COMMD8","MYO1G","SP140","PHPT1","CD2BP2","PLEKHO1","STAM","MRPL16","IL2RB","ID2","TESPA1","GOLGA8B","MIS18BP1","VAMP5","DAPK2","HLA-DPA1","TSG101","IL4R","CCND2","CTSC","TRAF3IP3","NLRC3","ORAI3","GNLY","MIR155HG","CARD16","CD82","ECH1","JAML","EEF1G","ETFB","DAXX","RBM4","HCST","RAB27A","YPEL2","CHST12","ARPC1B","PDIA4","PDIA6","AC243960.1","TBC1D10C","PTPN6","PYCARD","BST2","BTN3A2","MTG1","MLEC","DUSP4","GSDMD","SLAMF1","IFI6","PCID2","GIMAP1","ITGA1","CSNK2B","CDK2AP2","MYO1F","AC004687.1","PTTG1","APOBEC3C","TSPAN14","MOB3A","STXBP2","LCP2","PLA2G16","LINC00649","CST7","TADA3","SIT1","APOBEC3G","SUSD3","CD3G","CCL5","CDC25B","TNFRSF1B","HMGN3","THEMIS","ASF1A","CTNNB1","FIBP","CCDC85B","POLR3GL","GIMAP6","ARL6IP1","CALCOCO2","CCPG1","KLRB1","ACAA2","ISG15","EIF4A1","CAT","MANF","XAB2","GRINA","GLO1","LSM2","SLFN5","FKBP1A","AKNA","TAP1","LMO4","APEH","C12orf75","TMEM14A","DNPH1","C17orf49","NUDT5","MGAT1","CCDC69","EIF4EBP1","PDHB","ARL3","UCP2","IFI35","HSBP1","LYST","MRFAP1L1","ITGAL","AIP","RASAL3","CAPN1","ITGB1","RBPJ","LBH","DYNLL1","NME2","MT1F","SYNGR2","ABTB1","ZGPAT","CD63","ILK","SKA2","TMEM204","ACO2","HOPX","CRIP1","OXNAD1","CCS","GRAP2","GSTO1","HADHB","IL16","PIN4","CUEDC2","CALM3","SAMSN1","HM13","SNAP23","LPCAT4","FAAP20","EFHD2","PRDX3","CCM2","C22orf39","SDHA","ARRDC1","MAP4K1","NDUFA13","IL27RA","C14orf119"),
  NeoTCR4 = c("CXCL13","HMOX1","ETV7","ADGRG1","PDCD1","ENTPD1","CCDC50","TOX","CD4","TIGIT","TNFRSF18","NMB","MYL6B","AHI1","MAF","IFNG","LAG3","CXCR6","IGFLR1","DUSP4","ACP5","LINC01943","LIMS1","BATF","PCED1B","ITGAL","YPEL2","MAL","PPT1","ELMO1","MIS18BP1","TMEM173","ADI1","SLA","GALM","LBH","SECISBP2L","CTSB","C17orf49","CORO1B"),
  TRM_C8_Upadhye = c("ZNF683", "GNLY", "ITM2C", "HLA.DRB5", "GALNT2", "LGALS3", "ITGAE", "HLA.DRA", "APOBEC3C", "MIR155HG", "COTL1", "CD9", "HLA.DQA1", "IFNG", "APOBEC3G", "LRRN3", "GZMA", "UBXN11", "CD52", "GZMH", "CCL5", "HLA.DPB1", "ANXA5", "ITGB1", "HMGA1", "LINC00152", "CLIC3", "C12orf75", "CORO1B", "HLA.DPA1", "HLA.DQA2", "ITGA1", "ZYX", "PTTG1", "CD8A", "ICOS", "CD63", "PDE4DIP", "AC131056", "GSTM2", "ANKRD32", "TMEM14B", "MIR4435HG", "PTMS", "CD74", "SUSD3", "MT2A", "GSTP1", "LDLRAD4", "C7orf50", "CD81", "CAPG", "JUND", "CRIP1", "AMICA1", "IFI27L2", "OASL", "HAVCR2", "PSMB9", "CD70", "LAT", "FASLG", "IDH2", "CD8B", "HLA.DRB1", "PLA2G16", "LAG3", "CLEC2B", "RBPJ", "SIT1", "CXCR3", "RABAC1", "RP1124N14", "RPS4Y1", "PPDPF", "H1FX", "MT1F", "PPIB", "HMGN2", "CKLF", "MT1X", "GSTM1", "SLC2A4RG", "BLOC1S1", "TMEM14C", "SAMSN1", "CD27", "APOBEC3D", "PPP1R14B", "TYMP", "LSP1", "TSPO", "HLA.DMA", "GBP1", "MT1E", "TBCD", "PSMA7", "VAMP2", "ZFYVE28", "LSMD1", "C14orf2", "FKBP1A", "LIMD2", "RGL4", "TBC1D10C", "ANXA2", "PSMB10", "GNAI2", "PDIA6", "CALR", "RAB3GAP1", "GBP2", "PPP1CA", "ODC1", "HMGN1", "S100A11", "SUMO2", "SEPT1", "PRDX5", "RHOF", "CDK6", "PABPN1", "CD2", "CCND2", "GLIPR2", "TRAPPC1", "DAD1", "UQCR11", "CD82", "NDUFA4", "PSMB3", "FDFT1", "UBL5", "SHFM1", "IL2RB", "PLP2", "CPNE7", "DOK2", "GALNT1", "NDUFA3", "CAPZB", "MRPS6", "NDUFA13", "ETFB", "SEL1L3", "CSTB", "ARRDC3", "NDUFC1", "ABI3", "TPM4", "FUT8", "PPP1R35", "TBC1D2B", "LIME1", "HIST1H1E", "COMMD8", "VIM", "RASAL3", "TAP1", "ARHGAP9", "UBL3", "FABP5", "PRDX1", "PGAM1", "CHCHD10", "GNAS", "DNPH1", "COMMD7", "EMP3", "SRSF9", "NCOR1", "DDX24", "OST4", "PDIA3", "DBI", "RAB27A", "NME3", "CLIC1", "ACP5", "TSTD1", "PET100", "ARL6IP6", "AKAP11", "ARPC2", "ATP5E", "KMT2A", "TALDO1", "PTPN6", "RNF167", "NDUFAF3", "MYH9", "EPHA1", "CSK", "CAPNS1", "NMRK1", "EVL", "UBE2L6", "PPP1R18", "S100A6", "NDUFA11"),
  TRM_C6_Upadhye = c("TNFAIP3", "DUSP4", "RGCC", "CSRNP1", "NR4A2", "FOSB", "NFKBIA", "LMNA", "HSPA5", "DNAJA1", "UBE2S", "DNAJB1", "MYADM", "IER5L", "YPEL5", "NR4A3", "KLF6", "IER5", "HSP90AB1", "HSP90AA1", "SERTAD1", "PMAIP1", "JUNB", "HERPUD1", "DUSP5", "EIF4A3", "MCL1", "KPNA2", "TSPYL2", "HSPA1B", "FAM46C", "ZFAND5", "HOOK2", "IFRD1", "ZFP36", "BHLHE40", "ZNF683", "JUN", "EIF5", "NR4A1", "SRSF7", "LINC.PINT", "CITED2", "BRD2", "SBDS", "TUBB2A", "PER1", "PTGER4", "VIM", "PPP1R14B", "DNAJB6", "DDX3X", "ZBTB1", "RGS1", "SRSF5", "ZNF331", "ANXA1", "ZC3H12A", "H2AFX", "SRSF3", "PIM3", "DUSP1", "RGS2", "DDIT3", "ODC1", "JUND", "EZR", "NEU1", "CD83", "IDS", "CXCR4", "ITM2C", "SRSF2", "PNRC1", "TUBA1B", "KLF10", "ATF3", "GATA3", "CDKN1A", "TUBA1A", "KDM2A", "RBM38", "JMJD6", "TUBA1C", "TIPARP", "HSPE1", "TGFB1", "SLC7A5", "IER2", "BTG2", "TUBA4A", "HSPA8", "SIK1", "HSPA1A", "MIDN", "HSPH1", "HNRNPA0", "EGR1", "SNHG15", "TSC22D3", "RP1124N14", "PTP4A1", "ATP1B3", "SFPQ", "MAP3K8", "AMD1", "LDHA", "MARCKSL1", "UBB", "NDNL2", "GADD45B", "CCNL1", "IRF1", "MAP1LC3B", "SAT1", "TOB1", "RHOB", "NXT1", "STK17B", "GPR183", "RHOH", "AHR", "OSER1", "BCLAF1", "RP11.489E7", "REL", "ZFAND2A", "DUSP2", "PLK3", "RNF139", "ZFP36L1", "PLEC", "CD6", "TRA2B", "EIF4G2", "PHLDA1", "SLC3A2", "IRF2BP2", "PPP2CA", "BTG3", "VPS37B", "ZNF394", "NFKBIZ", "BUD31", "OASL", "UBALD2", "CLK1", "HAUS3", "YWHAZ", "H2AFZ", "SELK", "RSRC2", "TBCD", "ARHGAP9", "TBCC", "PPP1R15B", "MEF2D", "RELB", "HSPD1", "SIAH2", "CAPG", "CXCR3", "HEXIM1", "H1FX", "CTC50I14", "SERTAD3", "RUNX3", "CBX4", "SNHG12", "NASP", "DUSP10", "G3BP2", "PNP", "SAP18", "DEDD2", "HERPUD2", "CACYBP", "DDIT4", "MYH9", "GABARAPL1", "BCAS2", "EIF4A1", "FUS", "RASGEF1B", "FOS", "GSPT1", "PDE4D", "DNAAF2", "MAPK1IP1L", "FAM53C", "RBM39", "NFKB1", "SCAND1", "PDE4DIP", "NR1D2", "WHAMM", "ZC3HAV1", "MSL2", "CDK13", "YTHDF2", "CCDC59", "DOK2", "FOSL2", "HNRNPDL", "FAM102A", "CD44", "HNRNPU", "CDKN2AIP", "PIM2", "SLC1A5", "NEDD9", "RNF138", "CYCS", "C20orf24", "DYNLL1", "XCL1", "ARL4A", "OAT", "ZFYVE28", "RBM8A", "TOB2", "IL2RB", "B3GNT2", "CALR", "HMGA1", "KANSL2", "CEBPB", "ETF1", "ARF6", "CDK17", "SQSTM1", "C16orf80", "NFKBIB", "PLP2", "SRRT", "B4GALT1", "NRBP1", "TERF2IP", "LGALS3", "CHST11", "NUDT4", "RASAL3", "PTMS", "TBX21", "DNM2", "COTL1", "PNPLA8", "ARHGDIA", "FBRS", "SLC2A3", "NUFIP2", "ANKHD1", "DNAJB9", "AZIN1", "SPTY2D1", "EIF4A2", "TAGLN2", "CD82", "FOXO1", "CTSA", "ATXN2L", "PRNP", "UBE2A", "GRPEL1", "EHD1", "CAMLG", "SRRM2", "NCL", "BZW1", "SH3BP1", "SPTAN1", "DYNLL2", "PTS", "UBAP1", "ATP6V0C", "DNTTIP2", "POLR2A", "SURF4", "NGRN", "NR1H2", "AHNAK", "PPP1R16B", "PTGES3", "EIF1AX", "SIK3", "GNAS", "LCP2", "TNFSF9", "CPD", "RAB1A", "MRPL10", "MRFAP1", "ISCA1"),
  CD4CTL_C2_Upadhye = c("CCL4", "GZMK", "GZMH", "GZMA.1", "CCL5.1", "CST7", "ITM2C", "APOBEC3G", "CD74", "RGS1.1", "IFNG", "HLA.DPA1", "HLA.DPB1", "HLA.DRB1", "PDCD1", "CD81", "LINC00152", "RARRES3", "KIAA1551.1", "HLA.DRA", "APOBEC3C", "F2R", "DUSP4", "DUSP2.1", "MIR4435.1HG", "STK17A", "CXCR3.1", "HLA.DRB5", "COTL1", "CD2.1", "ANKRD32", "CTSD.1", "FAM46C.1", "SAMD3", "RASAL3", "CXCR6", "ARAP2", "APOBEC3H.1", "CNN2", "PIP4K2A", "DDX3Y", "CD27", "TGFBR2", "CD99.1", "DNAJC1", "CLEC2B", "ITGB2", "WNK1", "OASL", "PRDM1", "TUBA4A", "LSP1", "HCST.1", "CD6", "C12orf75", "TSEN54", "RUNX3", "FAM160B1", "NFATC2", "MXD4", "SLC9A3R1", "THEMIS", "TOX", "PARP8.1", "SRGN", "IKZF3", "LCP1.1", "SAMSN1", "PSMB9", "PTTG1", "RABGAP1L", "ARHGAP9", "CHST12", "IDH2", "EPS15", "PMAIP1.1", "PYHIN1", "TMEM2", "ABI3", "ITM2A", "PDE4DIP", "WIPF1", "NCOR1", "IFI16", "CLIC1.1", "SPON2", "RPS4Y1", "CASP4", "TBX21", "NCF1.1", "ARPC5L", "ZYX", "ZNF331", "CTSC", "GZMM", "CXCR4.1", "CCNH", "SIT1.1", "1-Sep", "ERAP2", "SUB1", "APMAP", "CAP1", "C9orf142.1", "GIMAP4", "TRAPPC10", "CKLF.1", "TC2N", "MYL12B.1", "NELFCD", "CD5", "PSTPIP1", "SASH3", "BLOC1S1", "UBL3", "ITGAL", "DENND2D", "LINC.PINT.1", "PPP1R18.1", "GLIPR2.1", "SYNE1", "CCDC167", "PLEKHF1.1", "ARHGAP18"),
  CD4CTL_C4_Upadhye = c("KLRB1.1", "CEBPD", "SLC4A10", "MYBL1", "TNFSF13B", "LTK", "RORC", "YWHAQ", "IL23R", "RORA", "ABCB1", "CCR6", "IL4I1", "AQP3.1", "CXCR6.1", "PTPN13", "CTSH", "OSTF1", "RUNX2", "FKBP11", "S100A4.1", "IFI44", "RBL2", "CD2.2", "AC092580.4", "S100A6", "DPP4.1", "YWHAH", "PRR5", "GYG1", "PLCB1", "TC2N.1", "AMICA1", "C4orf32", "RHOC.1", "MKNK1", "PDE4D.1", "CERK", "IL7R.1", "ERN1", "TNFSF14", "SYTL2", "MAN1A1", "PRF1", "HOPX.1", "IFNGR1.1", "PTPN4", "TMIGD2", "IFI44L", "CTSC.1", "OXNAD1.1", "CD96.1", "GSTM2", "SESN1", "IL18R1", "C10orf128", "SOCS2", "STAT4.1", "PERP", "CKLF.2", "TTC39C", "CMTM6", "MGAT4A", "CCR2", "CD40LG.1", "ACAA2", "REEP3", "OBFC1", "GPX1", "C9orf142.2", "SASH3.1", "ALOX5AP", "LCP1.2", "SLAMF1", "MDFIC.1", "GOLGA8B", "MIR142", "CAPG.1", "PTGER2.2", "ABRACL", "CCDC107", "PTMS", "SELT", "LAG3", "PTPN22", "TNF.1", "ZNRF1", "CISH.1", "ABI3.1", "GDI1", "SPOCK2", "BZW1", "ODF2L", "ID2", "PARP8.2", "PDCD4", "REEP5", "BCL2.1", "MACF1.1", "BHLHE40.2", "DYNLT3", "TNFSF12", "LGALS3", "CDC42EP3", "ANKRD28", "LPCAT4", "EFHD2", "ELK3", "TSEN54.1", "FAM50A", "NKG7", "IL12RB1", "TNFRSF25", "ADAM10", "GPR65", "PIM1", "DRAP1", "TP53I13", "LPIN2", "H2AFV", "TMEM167A", "PNP", "CHD9", "THEMIS.1", "CMC1", "ADA", "ARNTL", "VPS13C", "MIS18BP1", "ARL14EP", "STOM.1", "MGAT5", "DDIT4", "NEDD9", "SLK", "SEC61B", "NMRK1", "S100PBP", "MYO1F", "HSPA1B", "BAX", "HELZ", "INPP4A", "RPS6KA3", "HSPA1A.1")
)

genes <- unique(unlist(signature_list)); length(genes)
genes_in_Gex <- genes[genes %in% GeneExpression$gene_name]
GeneExpression <- GeneExpression[GeneExpression$gene_name %in% genes_in_Gex, ]; dim(GeneExpression) # 243 589
rownames(GeneExpression) <- GeneExpression$gene_name; GeneExpression$gene_name <- NULL

# Define Patients
patients_wMetadata <- mdata$Patient_ID; length(patients_wMetadata)
patients_wGex <- colnames(GeneExpression); length(patients_wGex)
patients_final <- intersect(patients_wMetadata, patients_wGex); length(patients_final) # 565
# subset
GeneExpression <- GeneExpression[,patients_final]; dim(GeneExpression)
mdata <- mdata[patients_final,]; dim(mdata)
# all(colnames(GeneExpression) == rownames(mdata))

## Apply CIBERSORT NORMALIZATION (OPTIONAL)
if(apply_CIBERSORT == TRUE){
  CIBERSORT <- read.csv(CIBERSORT_file, sep = ","); CIBERSORT$Patient_ID  <- rownames(CIBERSORT) <- gsub("[.]", "-", gsub(".01$", "", CIBERSORT$Mixture)); CIBERSORT$Mixture <- NULL
  desired_celltype <- c("T.cells.CD8")
  CIBERSORT <- CIBERSORT[patients_final, c("Patient_ID", desired_celltype)]
  normalization_fctr <- CIBERSORT[,desired_celltype]
} else{
  normalization_fctr <- rep(1, length(patients_final))
}

# Signature MeanExpresion (before cibersort normalization)
SignatureMean_noNorm_list <- lapply(signature_list, function(x){
  colMeans(GeneExpression[x,], na.rm = T)
})

# Calculate the mean expression of the genes in the signatures
# NOTE: Corrected seignature by multipliying by CIBERSORT immune cell fractions (if apply_CIBERSORT == TRUE)
SignatureMean_list <- lapply(signature_list, function(x){
  colMeans(GeneExpression[x,], na.rm = T)*normalization_fctr
})

# Get the median value for each signature or Percentile
SignatureMedian_down_list <- lapply(signature_list, function(x){
  tmp <- colMeans(GeneExpression[x,], na.rm = T)*normalization_fctr
  quantile(tmp, probs = c(percentile_down))
})

# Get the median value for each signature or Percentile
SignatureMedian_up_list <- lapply(signature_list, function(x){
  tmp <- colMeans(GeneExpression[x,], na.rm = T)*normalization_fctr
  quantile(tmp, probs = c(percentile_up))
})

# Divide Patients # x = "TRM_CLARKE"
signature_group_list <- lapply(names(SignatureMean_list), function(x){
  # if(SignatureMean_list[[x]] > SignatureMedian_list[[x]]) {"high"} else if(SignatureMean_list[[x]] <= SignatureMedian_list[[x]]) {"low"}
  out <- c()
  for(i in SignatureMean_list[[x]]){
    if(i > SignatureMedian_up_list[[x]]){ out <- c(out, "high")
    }else if(i <= SignatureMedian_down_list[[x]]){ out <- c(out, "low")
    }else{out <- c(out, "not_used")}
  }
  return(out)
}); names(signature_group_list) <- names(SignatureMean_list)



# Add signature mean Expression and groups to metadata
df_signature_value_noNorm <- data.frame(SignatureMean_noNorm_list); colnames(df_signature_value_noNorm) <- paste0(colnames(df_signature_value_noNorm), "_noNorm"); df_signature_value_noNorm$Patient_ID <- rownames(df_signature_value_noNorm)
df_signature_value <- data.frame(SignatureMean_list); df_signature_value$Patient_ID <- rownames(df_signature_value)
df_signature_group <- data.frame(signature_group_list); df_signature_group$Patient_ID <- df_signature_value$Patient_ID # df_signature_group$Patient_ID <- rownames(df_signature_group)
df <- merge(df_signature_value, df_signature_group, by = "Patient_ID", suffixes = c("_value", "_group"))
df <- merge(df, df_signature_value_noNorm, by = "Patient_ID")
df_norm_fctrs <- data.frame(normalization_fctr); df_norm_fctrs$Patient_ID <- df_signature_group$Patient_ID
df <- merge(df, df_norm_fctrs, by = "Patient_ID")
mdata <- merge(mdata, df, by = "Patient_ID")

# Save table
filename <- "~/Documents/CBTN_Survival_Analysis_table"
if(filter_gradde) filename <- paste0(filename, "_", grade)
write.table(mdata, paste0(filename, ".tsv"), row.names = F, sep = "\t")
# # Add diagnosis
# tmp.df <- merge(mdata_cBioPortal[,c("Patient_ID", "disease_type")], mdata_CAVATICA[,c("Kids.First.Participant.ID", "disease_type")], suffixes = c(".cBioPortal",".CAVATICA"), by.x = "Patient_ID", by.y = "Kids.First.Participant.ID")
# mdata2 <- merge(tmp.df, mdata, all.y = TRUE, by = "Patient_ID")
# write.table(mdata2, "~/Documents/CBTN_Survival_Analysis_table_wDiagnosis.tsv", row.names = F, sep = "\t")



## ------------------- Survival Analysis -------------------  ##

# --> NeoTCR8

mdata_surv <- mdata

# Select Desired signature
signatures <- c("NeoTCR8") # signatures <- c("CD8_C5_ZNF683", "CD8_C6_LAYN") # For Multiple Signature Classification
signatures_group <- paste0(signatures, "_group")
if(length(signatures_group) == 1) mdata_surv$SIGNATURE_GROUP <- mdata_surv[, signatures_group]
if(length(signatures_group) > 1)  mdata_surv <- unite(mdata_surv, col = "SIGNATURE_GROUP", signatures_group, sep = "_")

# Plot distribution of the signatures
for(signature in paste0(signatures, "_value")){
  p <- data.frame(sig = mdata_surv[,signature]) %>%
    ggplot(aes(x=sig)) +
    geom_density( fill="dodgerblue", alpha=0.5) +
    geom_vline(xintercept=median(mdata_surv[,signature]), size=0.5, color="black", linetype="dotted") +
    # geom_vline(xintercept=quantile(mdata_surv[,signature], probs = 0.25), size=0.5, color="black", linetype="dotted") +
    # geom_vline(xintercept=quantile(mdata_surv[,signature], probs = 0.75), size=0.5, color="black", linetype="dotted") +
    geom_vline(xintercept=quantile(mdata_surv[,signature], probs = percentile_down), size=0.5, color="red", linetype="dotted") +
    geom_vline(xintercept=quantile(mdata_surv[,signature], probs = percentile_up), size=0.5, color="red", linetype="dotted")
  plot(p)
}

# mdata_surv$SIGNATURE_GROUP <- mdata_surv$ANALYSIS

mdata_surv <- mdata_surv %>% filter(SIGNATURE_GROUP != "not_used")

mdata_surv$y = Surv(as.numeric(mdata_surv$OS_MONTHS), mdata_surv$VITAL_STATUS)
fit = survfit(y ~ SIGNATURE_GROUP,
              data = mdata_surv)

## and make a plot
ggsurvplot(fit = fit,
           data = mdata_surv,
           risk.table = F,
           pval = F,
           size = 1.2, # change line size
           palette =  c("royalblue4", "gray"), # cadetblue4 (CD4), royalblue4 (CD8)
           ggtheme = theme_classic(base_size = 11,
                                   base_family = "",
                                   base_line_size = 11/44,
                                   base_rect_size = 11/44)) +
  labs(x = 'Overall_survival [months]')

# Cox proportional hazards model
mdata_surv$SIGNATURE_GROUP <- factor(mdata_surv$SIGNATURE_GROUP, levels = c("low","high")) # Make the low group the reference
fit <- coxph(Surv(time = mdata_surv$OS_MONTHS, event = mdata_surv$VITAL_STATUS) ~ SIGNATURE_GROUP + AGE + GENDER, data = mdata_surv) # Looks better withou adding strata
ggadjustedcurves(fit, data = mdata_surv, method = "average", variable = "SIGNATURE_GROUP") # , palette = c("#E7B800", "#2E9FDF")
ggforest(fit)


# --> NeoTCR4

mdata_surv <- mdata

# Select Desired signature
signatures <- c("NeoTCR4") # signatures <- c("CD8_C5_ZNF683", "CD8_C6_LAYN") # For Multiple Signature Classification
signatures_group <- paste0(signatures, "_group")
if(length(signatures_group) == 1) mdata_surv$SIGNATURE_GROUP <- mdata_surv[, signatures_group]
if(length(signatures_group) > 1)  mdata_surv <- unite(mdata_surv, col = "SIGNATURE_GROUP", signatures_group, sep = "_")

# Plot distribution of the signatures
for(signature in paste0(signatures, "_value")){
  p <- data.frame(sig = mdata_surv[,signature]) %>%
    ggplot(aes(x=sig)) +
    geom_density( fill="dodgerblue", alpha=0.5) +
    geom_vline(xintercept=median(mdata_surv[,signature]), size=0.5, color="black", linetype="dotted") +
    # geom_vline(xintercept=quantile(mdata_surv[,signature], probs = 0.25), size=0.5, color="black", linetype="dotted") +
    # geom_vline(xintercept=quantile(mdata_surv[,signature], probs = 0.75), size=0.5, color="black", linetype="dotted") +
    geom_vline(xintercept=quantile(mdata_surv[,signature], probs = percentile_down), size=0.5, color="red", linetype="dotted") +
    geom_vline(xintercept=quantile(mdata_surv[,signature], probs = percentile_up), size=0.5, color="red", linetype="dotted")
  plot(p)
}

# mdata_surv$SIGNATURE_GROUP <- mdata_surv$ANALYSIS

mdata_surv <- mdata_surv %>% filter(SIGNATURE_GROUP != "not_used")

mdata_surv$y = Surv(as.numeric(mdata_surv$OS_MONTHS), mdata_surv$VITAL_STATUS)
fit = survfit(y ~ SIGNATURE_GROUP,
              data = mdata_surv)

## and make a plot
ggsurvplot(fit = fit,
           data = mdata_surv,
           risk.table = F,
           pval = F,
           size = 1.2, # change line size
           palette =  c("cadetblue4", "gray"), # cadetblue4 (CD4), royalblue4 (CD8)
           ggtheme = theme_classic(base_size = 11,
                                   base_family = "",
                                   base_line_size = 11/44,
                                   base_rect_size = 11/44)) +
  labs(x = 'Overall_survival [months]')

# Cox proportional hazards model
mdata_surv$SIGNATURE_GROUP <- factor(mdata_surv$SIGNATURE_GROUP, levels = c("low","high")) # Make the low group the reference
fit <- coxph(Surv(time = mdata_surv$OS_MONTHS, event = mdata_surv$VITAL_STATUS) ~ SIGNATURE_GROUP + AGE + GENDER, data = mdata_surv) # Looks better withou adding strata
ggadjustedcurves(fit, data = mdata_surv, method = "average", variable = "SIGNATURE_GROUP") # , palette = c("#E7B800", "#2E9FDF")
ggforest(fit)
