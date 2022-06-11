############    -  Assessment on antigen specificity  -    ############
############    --------   for the TCR data    --------    ############

# By Vicente Fajardo
# Modified: Kevin Meza Landeros

# Long version: 0.1
# Version: 0
# Version updates:
#   First non-stable version.
# Subversion: 7
#   Compared to v0.6, in here we:
#     * assess and select a clone size threshold.
#     * pick a single v and single j gene for each clonotype.
#     * provide a means to visualize the results provided by GLIPH.
#     * automize every step of the process.

### -------------------------- Description -------------------------- ###
# Based on both gene expression data (provided, for now, only as a seurat object) and TCR data (provided as standard 10x vdj files), this program will create:
#   The input files to run GLIPH2 over the TCR data of interest.

### Running Example
# Rscript3 ag_specificity_assessment_v1.4_cd4_and_cd8.R --celltype CD4 --date 2022-06-08 --LG FALSE --hinojosa FALSE --iGraph FALSE --expansionThold 1
# Rscript3 ag_specificity_assessment_v1.4_cd4_and_cd8.R --celltype CD8 --date 2022-06-08 --LG FALSE --hinojosa FALSE --iGraph FALSE --expansionThold 1
# Rscript3 ag_specificity_assessment_v1.4_cd4_and_cd8.R --celltype CD4_CD8 --date 2022-06-08 --LG FALSE --hinojosa FALSE --iGraph FALSE --expansionThold 1

cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(Seurat)
library(stringr)
library(tidyr)
library(data.table)
library(gtools)
library(ggplot2)
library(igraph)
library(optparse)
library(dplyr)
source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')
source('/home/vfajardo/scripts/functions/R_visualization_functions.1.5.R')
cat('Libraries imported!\n')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
  make_option(opt_str="--celltype", type="character", default="CD8", dest="celltype",
  help="Celltype that you want to obtain the TCRs from (CD4, CD8)"),
  make_option(opt_str="--date", type="character", default=NULL, dest="date", help="Date of the analysis."),
  make_option(opt_str="--LG", type="logical", default=FALSE, dest="LG", help="Boolean flag saying if
  you want to select only Low_grade samples."),
  make_option(opt_str="--hinojosa", type="logical", default=FALSE, dest="hinojosa", help="Boolean flag saying if
  you want to include the Hinojosa_NatComm_2021 TCRs."),
  make_option(opt_str="--iGraph", type="logical", default=FALSE, dest="iGraph", help="Boolean flag saying if
  you want to plot your GLIPH2 results. If it is set to false, it will give you GLIPH2 input files; otherwise you will
  obtain the plots (GLIPH2 results are needed)."),
  make_option(opt_str="--expansionThold", type="numeric", default=1, dest="expansion.thold", help="Numeric value indicating
  the minimum clone size you wish the clonotipes to have in order to get GLIPH2 input files.")
)

opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
celltype <- opt$celltype
date <- opt$date
LG <- opt$LG
hinojosa <- opt$hinojosa
iGraph <- opt$iGraph
expansionThold <- opt$expansionThold

if(interactive()){
  celltype <- "CD8"
  date <- "2022-03-30"
  LG <- FALSE
  hinojosa <- FALSE
  iGraph <- FALSE
  expansionThold <- 1
}

# ---> General definitions.
# @ On the project.
gen.prj <- 'pbtumor-all'
paper.lab <- '/results/tcr/preproces_4_Gliph2/CD3p/'
this.assess <- 'ag_specificity_assessment'

# celltype <- 'CD8' #'CD4' #
# expansion.thold <- 1 # greater or equal
# date <- '2022-02-07'
# LG <- FALSE #TRUE #

# ---> Path definitions.
reports.date <- date
gen.reports.path <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/', gen.prj, paper.lab)
gen.reports.path <- paste0(gen.reports.path, this.assess, '/reports')
data.path <- str_replace(string=gen.reports.path, pattern='reports$', replacement='data')
reports.path <- paste0(gen.reports.path, '/', this.assess, '_', if(is.null(reports.date)) Sys.Date() else reports.date)
if(!dir.exists(reports.path)) dir.create(reports.path, recursive=TRUE)

# ---> File definitions.
combined_celltypes <- FALSE
if ( celltype == "CD4") {
  folder <- "CD45pCD3p_clean2_cd4"
} else if ( celltype == "CD8") {
  folder <- "CD45pCD3p_clean2_cd8"
} else if ( celltype == "CD4_CD8") {
  combined_celltypes <- TRUE
  folder <- "CD45pCD3p_clean2_cd8"
}

# TCR data.
cells.clons.info.file <- paste0('/home/ciro/large/pbtumor/results/tcr/', folder, '/filtered_contig_annotations_aggr.csv')
clons.info.file <- paste0('/home/ciro/large/pbtumor/results/tcr/', folder, '/clonotypes_aggr.csv')

if(combined_celltypes){
  seurat.obj.file.1 <- paste0('/home/ciro/large/pbtumor/results/tcr/', "CD45pCD3p_clean2_cd4", '/metadata.rds')
  seurat.obj.file.2 <- paste0('/home/ciro/large/pbtumor/results/tcr/', "CD45pCD3p_clean2_cd8", '/metadata.rds')
} else{
  # Gene expression data. # We will only use the metadata of the Seurat object
  seurat.obj.file <- paste0('/home/ciro/large/pbtumor/results/tcr/', folder, '/metadata.rds') #'/mnt/beegfs/kmlanderos/CD45pCD3p_clean2_object_lock_mean0.01_pct15_pc15.rds'
}

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Gene expression info.
if(combined_celltypes){
  seurat.obj.1 <- readRDS(file=seurat.obj.file.1)
  seurat.obj.2 <- readRDS(file=seurat.obj.file.2)
  seurat.obj <- rbind(seurat.obj.1, seurat.obj.2)
} else{
  seurat.obj <- readRDS(file=seurat.obj.file)
}

if (LG == TRUE) {
  seurat.obj$diagnosis_grade <- ifelse(!is.na(seurat.obj$orig.Tumor.grade), ifelse(seurat.obj$orig.Tumor.grade %in% c("III","IV"),"High_grade", "Low_grade"), NA)
  seurat.obj <- seurat.obj[seurat.obj$diagnosis_grade == "Low_grade" & !is.na(seurat.obj$diagnosis_grade),]
}
# ---> TCR data files
# Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
# Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')

# ---> Rules to create new tags.
# new.tags.rules <- if(!is.null(new.tags.file)) read.csv(file=new.tags.file, stringsAsFactors=FALSE) else NULL


cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')

# ---> Cells-clonotypes relationships info (filtering).
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]


cat('\n\n')
### ---------------- Create individual GLIPH2 input ---------------- ###
cat('### ---------------- Create individual GLIPH2 input ---------------- ###\n')

# ---> Assess clonal expansion distribution across unique clonotypes for the datasets of interest in order to pick cutoffs.
tmp.file.name <- paste0(reports.path, '/CloneSizeDensities_',celltype,'.pdf')
# Get plot
pdf(file=tmp.file.name)
meta.data <- as.data.table(seurat.obj)
tmp.data <- meta.data[, .(clone.size=.N), by=clonotype.tag]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=clone.size)) +
  geom_density(fill='indianred', alpha=0.6) +
  geom_vline(xintercept=3, linetype='dashed') +
  scale_x_log10() + scale_y_continuous(expand=c(0, 0)) +
  labs(title=celltype, x='Clone size', y='Density') + theme_bw()
print(tmp.ggplot)
dev.off()

# ---> Program as a function.
get.gliph2.input.1 <- function(seurat.obj, cells.clons.info, clons.info){
  # ---> Fill cell-clonotypes relationships with tags of interest info.
  # Take the overall metadata
  meta.data <- seurat.obj
  meta.data$barcode <- rownames(meta.data)
  # Then, combine clones data with other meta data and update seurat object.
  tmp.data <- unique(cells.clons.info[, c('barcode', 'raw_clonotype_id')])
  # tmp.data$clonotype.tag <- tmp.data$raw_clonotype_id; tmp.data$raw_clonotype_id <- NULL # NOTE: this helps double check clonotype annotation. But if the seurat object has not been previously annotated, this will obviously not work.
  tmp.data <- merge(x=tmp.data, y=meta.data, by='barcode', all.x=FALSE, all.y=TRUE)
  row.names(tmp.data) <- tmp.data$barcode
  # Check previous clonotype annotations match the ones here. Else, break loop.
  tmp.check <- all(tmp.data$clonotype.tag==tmp.data$raw_clonotype_id, na.rm=TRUE)
  if(!tmp.check) stop('Unexpected base error.\n')
  # Check that up to this point merging has been successful in terms of the cells that were processed. If so, update seurat object's metadata.
  tmp.check <- all(rownames(seurat.obj) %in% row.names(tmp.data)) & all(row.names(tmp.data) %in% rownames(seurat.obj))
  if(!tmp.check) stop('Unpected error 1. Please check the code.\n')
  meta.data <- tmp.data[rownames(seurat.obj), ]
  tmp.data <- as.data.table(tmp.data)
  # Include clonotypes' sequences.
  gex.clons.info <- merge(x=tmp.data, y=clons.info[, c('clonotype_id', 'cdr3s_aa', 'cdr3s_nt')], by.x='clonotype.tag', by.y='clonotype_id', all=FALSE)
  # gex.clons.info = metadata but with clonotypes' sequences TBA, TRB in the same column.

  # ---> vdj gene usage info.
  # General vdj gene usage info.
  # We obtain unique v and j genes per clonotype. When cellranger has defined multiple v and j genes for a given clonotype, we take the one that is most supported. The most supported is the gene that has the largest amount of entries (contigs) in the table below.
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
  # Amount of them that have multiple equally supported v and j gene
  cat("\n",dim(vdj.gene.info[str_detect(string=v.gene, pattern='\\|'),])[1], "genes have multiple equally supported v gene.\n")
  cat("\n",dim(vdj.gene.info[str_detect(string=j.gene, pattern='\\|'),])[1], "genes have multiple equally supported j gene.\n")
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)
  # Sanity check. Confirm we have unique entries for the amount of unique clonotypes.
  tmp.check <- vdj.gene.info[, uniqueN(raw_clonotype_id)==.N]
  if(!tmp.check) stop('Unexpected error post vdj data retrieval.\n')

  # Temporal copy
  # vdj.gene.info.copy <- vdj.gene.info
  # vdj.gene.info <- vdj.gene.info.copy

  # ---> Info per clonotype (i.e., disregarding GEx and extended information).
  # We must keep only the clonotypes that meet certain criteria (as described below) at the time that we keep track of the proportions clonotypes that are kept after each filtering step.

  #   Filtering rule #1: Clonotypes must be private in the sense that they must be found with any size for a single donor or as clonally expanded in a single donor out of multiple donors that seem to express it.
  # Preprocess data.
  # NOTE: Eliminated
  clons.info <- gex.clons.info[
    ,
    .(
      cdr3a.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
      cdr3b.aa.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_aa), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
      cdr3a.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRA:\\w+[^;\\n]', simplify=TRUE), pattern='TRA:|;', replacement=''), collapse=';'),
      cdr3b.nt.seq=paste0(str_replace(string=str_extract_all(string=unique(cdr3s_nt), pattern='TRB:\\w+[^;\\n]', simplify=TRUE), pattern='TRB:|;', replacement=''), collapse=';'),
      size=.N#,
      # pr.tag=unique(pr.tag)
    ),
    by=.(
      clonotype.tag,
      orig.donor
    )
  ]
  clons.info[cdr3a.aa.seq=='', cdr3a.aa.seq:=NA]
  clons.info[cdr3b.aa.seq=='', cdr3b.aa.seq:=NA]
  clons.info[cdr3a.nt.seq=='', cdr3a.nt.seq:=NA]
  clons.info[cdr3b.nt.seq=='', cdr3b.nt.seq:=NA]
  # clons.info = clonotypes repeated per donor, each donor has its own clone size

  # Keep track of proportions.
  filt.track <- clons.info[, uniqueN(clonotype.tag)]
  # Apply rule.
  # tmp.data <- clons.info[, .(donor.count=uniqueN(orig.donor), expanded.count=sum(size>=expansion.thold), keep=TRUE), by=clonotype.tag]
  # tmp.data[donor.count>1 & expanded.count!=1, keep:=FALSE]; tmp.data <- tmp.data[keep==TRUE, clonotype.tag]
  # clons.info <- clons.info[clonotype.tag %chin% tmp.data]
  # If public but expanded for a single donor, keep the info only for the donor where expanded.
  # tmp.data <- clons.info[, .SD[size==max(size), .(orig.donor=unique(orig.donor))], by=clonotype.tag]
  # clons.info <- merge(x=clons.info, y=tmp.data, by=c('clonotype.tag', 'orig.donor'))

  #   Filtering rule #2: Clonotypes must have at least a single properly defined beta sequence.
  # Keep track of proportions.
  filt.track[2] <- clons.info[, uniqueN(clonotype.tag)]
  # Apply rule.
  clons.info <- clons.info[!is.na(cdr3b.aa.seq)] # 8226 to 7812 - CD8
  clons.info <- clons.info[!grepl(";",cdr3b.aa.seq)] # 7812 to 7769 - CD8

  # Keep track of proportions.
  filt.track[3] <- clons.info[, uniqueN(clonotype.tag)]

  # Further, include vdj gene usage info.
  clons.info <- merge(x=clons.info, y=vdj.gene.info, by.x='clonotype.tag', by.y='raw_clonotype_id', all.x=TRUE, all.y=FALSE)
  return(clons.info)
}
# clons.info.copy <- clons.info
# clons.info.copy == clons.info

# ---> Run for every dataset for all TCRs.
gliph.all <- get.gliph2.input.1(seurat.obj=seurat.obj, cells.clons.info=cells.clons.info, clons.info=clons.info)
# Save this file
fname=paste0("/.intermidiate_file_TCRsimilarity_",celltype,".tsv")
fwrite(file=paste0(reports.path,fname), x=gliph.all, sep='\t', na=NA, quote=FALSE)


cat('\n\n')
### ------------------- Get actual GLIPH2 inputs ------------------- ###
cat('### ------------------- Get actual GLIPH2 inputs ------------------- ###\n')

# # ---> HLA data.
# # Remove the asterisk prefixing some of the gene fields. Those cases mark the genes whose redundancy was resolved based on some software as explained by the provider. Check the raw results sheet if necessary.
# hla.genes <- setdiff(x=colnames(hla.data), y='orig.donor')
# for(tmp.gene in hla.genes){
#   hla.data[, tmp.gene] <- str_replace(string=hla.data[, tmp.gene], pattern='^\\*\\s*', replacement='')
# }
# # Separate alleles for a given type of HLA gene.
# hla.genes <- setdiff(x=colnames(hla.data), y='orig.donor')
# for(tmp.gene in hla.genes){
#   hla.data <- separate(data=hla.data, col=tmp.gene, into=paste(tmp.gene, 1:2, sep='.'), sep=' \\+ ')
# }
# # Keep only the information provided by the first four digits. This is the most important information according to the HLA gene nomenclature as explained here:
# #   http://hla.alleles.org/nomenclature/naming.html
# hla.genes <- setdiff(x=colnames(hla.data), y='orig.donor')
# for(tmp.gene in hla.genes){
#   hla.data[, tmp.gene] <- str_extract(string=hla.data[, tmp.gene], pattern='^[^:]+:[^:]+')
# }
# # Further include the HLA gene nomenclature prefix (dropping the 'HLA-').
# hla.genes <- setdiff(x=colnames(hla.data), y='orig.donor')
# for(tmp.gene in hla.genes){
#   hla.pfx <- str_replace(string=tmp.gene, pattern='^HLA\\.', replacement=''); hla.pfx <- str_replace(string=hla.pfx, pattern='\\.[12]$', replacement='')
#   hla.data[!is.na(hla.data[, tmp.gene]), tmp.gene] <- paste(hla.pfx, hla.data[!is.na(hla.data[, tmp.gene]), tmp.gene], sep='*')
# }
# # ---> GLIPH2 input @ HLA data.
# tmp.data <- clons.info[, .(orig.donor=unique(orig.donor))]#[!orig.donor%chin%hla.data$orig.donor]
# tmp.data <- merge(x=tmp.data, y=hla.data, by='orig.donor')
# tmp.file.name <- paste0(reports.path, '/GLIPH2Input_HLA-Data.tsv')
# fwrite(file=tmp.file.name, x=tmp.data, sep='\t', na='', quote=FALSE)

# ---> TCR GLIPH input, program as a function.
# @ Program as a function.
get.gliph2.input.2 <- function(gliph.all, reports.path, file.prefix=NA){
  # Get general data.
  tmp.data <- gliph.all #rbindlist(l=gliph.all, use.names=TRUE, idcol='data.set')
  # ---> Split info between chain types.
  # @ Beta TCR chains.
  beta.info <- tmp.data[!is.na(cdr3b.aa.seq)] # Not necessary, done previously, remove if necessary.
  # beta.info[, cdr3a.aa.seq:=NULL] # Can uncoment, only comment if you want to keep a GLIPH Input version with alpha chains
  # beta.info <- as.data.table(separate_rows(data=beta.info, cdr3b.aa.seq, sep=';', convert=FALSE)) # NOTE: We assumme double beta means it's a doublet.
  beta.info <- beta.info[
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
      size=sum(size),
      cdr3a.aa.seq,
      tra.v,
      tra.j
     ),
    by=.(cdr3b.aa.seq, trb.v, trb.j)
  ]
  beta.info[grepl("\\|",clonotype.tag)| grepl("\\|",orig.donor)] # Empty, okay
  # beta.info[, .N, by=.(cdr3b.aa.seq, trb.v, trb.j)][N>1] # Sanity check. Now we have unique entries of clonotypes where each clonotype is defined by a single TCR beta chain and a unique set of v and j genes.
  # Thus, check how many clonotypes were equally called as pR and non-pR (same amount of cells for each status).
  # beta.info[, .N, by=pr.tag]; beta.info[pr.tag=='non-pR;pR'] # NOTE: Sanity check. Remove if necessary.
  # ---> Clone size assessment.
  # Get plot. One for the whole dataset and one distinguishing between pR- and non-pR-expressing cells.
  tmp.ggplot.1 <- ggplot(data=beta.info, aes(x=size)) +
    geom_density(fill='indianred', alpha=0.6) +
    geom_vline(xintercept=3, linetype='dashed') +
    scale_x_log10() + scale_y_continuous(expand=c(0, 0)) +
    labs(title=celltype, x='Clonesize (log 10)', y='Density') + theme_bw()
  # Output plot.
  tmp.file.name <- paste0(reports.path, '/CloneSizeDensities_', celltype, '_Reclonotyping.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot.1)
  dev.off()
  # Pick expansion threshold and indicate it below.
  tmp.expansion.thold <- 1 # greater or equal
  # ---> Actual GLIPH2 inputs.
  # @ Beta chains.
  tmp.data <- beta.info[size>=tmp.expansion.thold, .(cdr3b.aa.seq, trb.v, trb.j, cdr3a.aa.seq=NA, donor.id=paste(orig.donor, 'Exp1', sep=':'), size)]
  # Separate multiple V and J genes.
  # tmp.data[str_detect(string=trb.v, pattern='\\|'), trb.v]
  # tmp.data <- as.data.table(separate_rows(data=tmp.data, trb.v, sep=';', convert=FALSE)) # Should remain the same, bc we kept only one.
  # tmp.data <- as.data.table(separate_rows(data=tmp.data, trb.v, sep='\\|', convert=FALSE)) # Should remain the same, bc we kept only one.
  # tmp.data <- as.data.table(separate_rows(data=tmp.data, trb.j, sep=';', convert=FALSE)) # Should remain the same, bc we kept only one.
  # tmp.data <- as.data.table(separate_rows(data=tmp.data, trb.j, sep='\\|', convert=FALSE)) # Should remain the same, bc we kept only one.
  tmp.data <- unique(tmp.data) # Sanity check. Should remain the same.
  # Output.
  if(!is.na(file.prefix)){
    tmp.file.name <- paste0(reports.path, '/', file.prefix, '.tsv')
    fwrite(file=tmp.file.name, x=tmp.data, sep='\t', na=NA, quote=FALSE)
    # fwrite(file=gsub("GLIPH2Input", ".GLIPH2Input", gsub(".tsv","_beta.info.tsv",tmp.file.name)), x=beta.info, sep='\t', na=NA, quote=FALSE) # Save this file for later iGraph visualization.
  }
  # Save another file that has the alpha chain and V- J- information.
  tmp.data.alpha <- beta.info[size>=tmp.expansion.thold, .(cdr3b.aa.seq, trb.v, trb.j, cdr3a.aa.seq, tra.v, tra.j, clonotype.tag, donor.id=paste(orig.donor, 'Exp1', sep=':'), size)] %>% mutate(id = paste(cdr3b.aa.seq, trb.v, trb.j, sep=","))
  tmp.data.alpha$cdr3a.aa.seq <- as.character(tmp.data.alpha$cdr3a.aa.seq); tmp.data.alpha[is.na(tmp.data.alpha$cdr3a.aa.seq),"cdr3a.aa.seq"] <- "NA"
  tmp.data.alpha <- tmp.data.alpha %>% group_by(id) %>% mutate(cdr3a.aa.seq.all = paste(cdr3a.aa.seq, collapse = ";")) %>% distinct(id, .keep_all= TRUE)
  if(!is.na(file.prefix)){
    tmp.file.name <- paste0(reports.path, '/', file.prefix, '_wAlphaChain.tsv')
    fwrite(file=tmp.file.name, x=tmp.data.alpha, sep='\t', na=NA, quote=FALSE)
    # fwrite(file=gsub("GLIPH2Input", ".GLIPH2Input", gsub(".tsv","_beta.info.tsv",tmp.file.name)), x=beta.info, sep='\t', na=NA, quote=FALSE) # Save this file for later iGraph visualization.
  }
  # Return results
  beta.info[, cdr3a.aa.seq:=NULL]
  return(beta.info)
}

# Run
if(LG == TRUE){
  gliph.all.beta.info <- get.gliph2.input.2(gliph.all=gliph.all, reports.path=reports.path, file.prefix=paste0('GLIPH2Input_TCR-Data_All-TRB_', celltype, '_LG') )
} else{
  gliph.all.beta.info <- get.gliph2.input.2(gliph.all=gliph.all, reports.path=reports.path, file.prefix=paste0('GLIPH2Input_TCR-Data_All-TRB_', celltype) ) #, sep.pr=TRUE)
}
cat('Data has been preprocessed!\n')

# ============= #

# Add Rivero Hinojosa (Nature Communications) pbtumor TCR Data_All
# cat TCRDetails_Hinojosa_pre.tsv | awk 'OFS="\t" {print $4,$6,$8,$3}' >  TCRDetails_Hinojosa.tsv
if(hinojosa){
  tcr_hinojosa <- read.csv("/home/kmlanderos/pbtumor-all/info/TCRDetails_Hinojosa.tsv", sep = "\t")
  tcr_hinojosa$cdr3a.aa.seq <- NA
  tcr_hinojosa[,"donor.id:condition:study"] <- paste(NA, "Exp1", "Hinojosa_NatComm_2021", sep=":")
  tcr_hinojosa$donor.id <- NULL
  tcr_hinojosa <- tcr_hinojosa[,c(1,2,3,5,6,4)]

  tcr_uphadye <- read.csv(paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_", date, "/GLIPH2Input_TCR-Data_All-TRB_", celltype, ".tsv"), sep = "\t")
  tcr_uphadye[,"donor.id:condition:study"] <- paste(tcr_uphadye$donor.id, "Uphadye_NNNNN_202N", sep=":")
  tcr_uphadye$donor.id <- NULL
  tcr_uphadye <- tcr_uphadye[,c(1,2,3,4,6,5)]

  tcr_all <- rbind(tcr_uphadye, tcr_hinojosa)
  fwrite(file=paste0("/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_", date, "/GLIPH2Input_TCR-Data_All-TRB_", celltype, "_Hinojosa.tsv"), x=tcr_all, sep='\t', na=NA, quote=FALSE)
  cat("TCR data from Hinojosa's paper has been added!\n")
}

if(iGraph){
  cat('\n\n')
  ### -------------------- Results visualization --------------------- ###
  cat('### -------------------- Results visualization --------------------- ###\n\n')

  library(tidyr)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(igraph)
  # ---> Run GLIPH2 using the web platform.
  input.path = paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_",date)
  results.path=paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/gliph/',date)
  # date <- "2022-02-07"
  # celltype <- "CD8"
  # hinojosa <- FALSE #TRUE#
  # LG <- FALSE#TRUE#
  # expansion.thold <- 1 #greater or equal

  # ---> Analaysis of the results provided by GLIPH2 through iGraph.
  # @ Program as a function.
  check.gliph2.output <- function(file.suffix, beta.info, hinojosa = FALSE){
    # @ Load data as removing the 'single' patterns.
    # Define results of interest.
    tmp.file.name <- paste0(results.path, '/GLIPH2Output_TCR-Data_', file.suffix, '.csv')
    gliph.res <- fread(file=tmp.file.name, blank.lines.skip=TRUE)
    gliph.res <- gliph.res[pattern!='single']

    # @ Further preprocessing.
    tmp.data.1 <- merge(x=gliph.res, y=beta.info[, .(TcRb=cdr3b.aa.seq, V=trb.v, J=trb.j)], by=c('TcRb', 'V', 'J'), all.x=TRUE, all.y=FALSE, sort=FALSE)
    # Define summary information for every pattern.
    #     Average clone size.
    tmp.data.2 <- tmp.data.1[, .(size.median=median(Freq)), by=pattern]

    # @ Define input for iGraph.
    # Edges' info.
    edge.info <- merge(
      x=tmp.data.1[, .(connection=paste(TcRb, V, J, sep=';')), by=pattern],
      y=tmp.data.1[, .(connection=paste(TcRb, V, J, sep=';')), by=pattern],
      by='connection', all=TRUE, allow.cartesian=TRUE
    )
    edge.info <- edge.info[pattern.x!=pattern.y]
    edge.info <- edge.info[, .(weight=.N), by=.(from=pattern.x, to=pattern.y)]
    # Nodes' info.
    if(hinojosa){
      gliph.res.tmp <- separate(gliph.res, col="Sample", into=c("donor","condition","study"), sep=":")
      gliph.res.tmp <- distinct(gliph.res.tmp[,c("pattern","study")])
      gliph.res.tmp$value <- 1
      gliph.res.tmp <- pivot_wider(gliph.res.tmp, names_from = study, values_from = value) #, values_fill = 0)
      gliph.res.tmp[is.na(gliph.res.tmp)] <- 0
      motif_study <- apply(gliph.res.tmp, 1, function(x){
        suma <- sum(as.numeric(x[2:3]))
        if(suma==2) study <- "shared"
        if(suma==1) study <- colnames(gliph.res.tmp)[which(x == "1")]
        study
      })
      tmp.data.2$study <- sapply(motif_study, FUN=function(x) str_split(x, "_")[[1]][1])
    }
    nodes.info <- tmp.data.2

    # @ iGraph.
    # General definition of an iGraph.
    igraph.obj <- graph_from_data_frame(d=edge.info, vertices=nodes.info, directed=FALSE)
    if(hinojosa){
      min.val <- min(E(igraph.obj)$weight); max.val <- max(E(igraph.obj)$weight)
      E(igraph.obj)$weight <- (E(igraph.obj)$weight*13/max.val)+2 # Scale so that it goes from 2 to 15.
      E(igraph.obj)$width <- E(igraph.obj)$weight
    }else{
      E(igraph.obj)$width <- E(igraph.obj)$weight
    }

    min.val <- min(V(igraph.obj)$size.median); max.val <- max(V(igraph.obj)$size.median)
    V(igraph.obj)$size <- (V(igraph.obj)$size.median*8/max.val)+2 # Scale so that it goes from 2 to 10.
    # Colors for the vertices according to their presence in the pbtumos datasets.
    if(hinojosa){
      colrs.v <- c(shared="gold", Uphadye="tomato", Hinojosa='green')
      V(igraph.obj)$color = colrs.v[V(igraph.obj)$study]# Get general things for plotting.
    }
    # Plot.
    tmp.file.name <- paste0(results.path, '/GLIPH2Results_iGraph_', file.suffix, '.pdf')
    pdf(file=tmp.file.name)
    plot(igraph.obj, alpha=0.6, vertex.label=NA, main= paste0("TCR motifs found by GLIPH2 in ", celltype, "+ Tcells ", "(cs=>", expansion.thold, ")" )) #, vertex.label.cex = 0.6, vertex.label.dist=2
    #legend(x='bottomleft', pch=21, pt.bg=tmp.legend$color, bty='n', title='pR degree') #, legend=tmp.legend$attr
    dev.off()
    # To include the size legend, see: https://stackoverflow.com/questions/38451431/add-legend-in-igraph-to-annotate-difference-vertices-size
    #
    # @ Processed GLIPH results table.
    # Get table.
    # tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='pattern', sort=FALSE)
    # setorderv(x=tmp.data, cols=c('index', 'Fisher_score', 'size.median'), order=c(1, 1, -1))
    # tmp.file.name <- paste0(reports.path, '/GLIPH2Results_Processed_', file.suffix, '.csv')
    # fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
  }

  # fname <- paste0("/.GLIPH2Input_TCR-Data_All-TRB_", celltype); if(hinojosa) paste0(fname,"_Hinojosa")
  fname <- paste0("/GLIPH2Input_TCR-Data_All-TRB_", celltype); if(hinojosa) fname <- paste0(fname,"_Hinojosa")
  gliph.all.beta.info <- read.csv(paste0(input.path, fname, ".tsv"), sep="\t")
  suffix <- celltype; if(hinojosa) suffix <- paste0(suffix,"_Hinojosa"); if(LG) suffix <- paste0(suffix,"_LG")
  check.gliph2.output(file.suffix=suffix, beta.info=data.table(gliph.all.beta.info), hinojosa = hinojosa)
}

if(!iGraph){
  ## --------- TCR selction to clone -------------- ##
  # Create files with the Top25 of expanded clones and all the public clones.
  cat("\nGetting PublicClones and Top Expanded Clones!")

  # date <- "2022-02-13"
  # celltype <- "CD8"
  # hinojosa <- FALSE #TRUE#
  # LG <- FALSE #TRUE#
  # expansion.thold <- 1 #greater or equal
  input.path = paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_",date)
  results.path=paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/gliph/',date)
  if(!dir.exists(results.path)) dir.create(results.path) #, recursive=TRUE)

  fname <- paste0("/GLIPH2Input_TCR-Data_All-TRB_", celltype, ".tsv")
  f <- read.csv(paste0(input.path, fname), sep="\t")
  f <- f %>% separate(col=donor.id, into=c("donor.id", "condition"), sep=":") %>%
        select(cdr3b.aa.seq,trb.v,trb.j,donor.id,size)
  f.temp <- f[order(f$size, decreasing=T),]

  # - Public clones
  f.public <- f.temp[str_detect(string=f.temp$donor.id , pattern=';'),]
  f.public <- f.public[!(str_detect(string=f.public$donor.id , pattern='NA') & str_count(f.public$donor.id, ";")==1 ),] # To conserve NA (no hashtag assigned), the clonotype needs to have at least 2 known donors, apart from the NA.
  # ==== We break down, the donor proportion for this important clones.
  f.public.donors <- data.table(read.csv(paste0(input.path, "/.intermidiate_file_TCRsimilarity_", celltype, ".tsv"), sep="\t"))
  f.public.donors <- f.public.donors[!is.na(cdr3b.aa.seq)]; f.public.donors[, cdr3a.aa.seq:=NULL]
  beta.info <- f.public.donors[
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
  ] # This object is to get the cellranger clonotypes
  # To get donor specific info
  f.public.donors <- f.public.donors[
    ,
    .(
      size=sum(size)
     ),
    by=.(cdr3b.aa.seq, trb.v, trb.j, orig.donor)
  ] %>% pivot_wider(names_from = orig.donor, values_from = size)
  f.public.donors[is.na(f.public.donors)] <- 0
  # NOTE: There are cases when we have the same cdr3b but different VJ genes, we expect to have the same J gene for those cases.

  # Keep only the public clonotypes's info
  selection <- paste(f.public.donors$cdr3b.aa.seq, f.public.donors$trb.v, f.public.donors$trb.j, sep= "_") %in% paste(f.public$cdr3b.aa.seq, f.public$trb.v, f.public$trb.j, sep= "_")
  f.public.donors <- f.public.donors[selection,] # %>% as.data.frame()
  f.public.donors$size <- apply(f.public.donors, MARGIN=1,function(x) {sum(as.integer(x[4:dim(f.public.donors)[2]])) } )
  f.public.donors <- f.public.donors[order(f.public.donors$size, decreasing=T),]

  f.public.donors.tmp <- merge(f.public.donors, beta.info[,-c("size")], all.x=T , all.y=F, by=c("cdr3b.aa.seq", "trb.v", "trb.j"))
  f.public.donors.tmp <- f.public.donors.tmp[order(f.public.donors.tmp$size, decreasing=T),]

  # - Most expanded clones
  f.expanded <- f.temp[f.temp$size > 4,] # f.temp[1:25,]
  f.expanded <- merge(f.expanded, beta.info[,-c("size")], all.x=T , all.y=F, by=c("cdr3b.aa.seq", "trb.v", "trb.j"))
  f.expanded <- f.expanded[order(f.expanded$size, decreasing=T),]

  # Write final files
  fwrite(file=paste0(results.path,"/PublicClones_", celltype, ".tsv"), x=f.public.donors.tmp, sep='\t', na=NA, quote=FALSE)
  fwrite(file=paste0(results.path,"/TopExpandedClones_", celltype, ".tsv"), x=f.expanded, sep='\t', na=NA, quote=FALSE)

  cat("\nDone.")
}

########## Cluster 17

# See if there is any clonotype of interest in the Cell Cycle cluster
if(0){

  library(stringr)
  library(data.table)
  require(dplyr)

  # source("/home/kmlanderos/pbtumor-all/scripts/object_lock.R")
  # source("/home/kmlanderos/pbtumor-all/scripts/figures/global.R")
  date <-"2022-02-13"

  # Read Cluster 17 metadata
  sc_cd3p_c17_mdata <- readRDS("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/CD45pCD3p_clean2_c17/metadata.rds")
  sc_cd3p_c17_mdata$barcode <- rownames(sc_cd3p_c17_mdata)
  # Read relevant clones (Top25 expanded and public) from CD4 and CD8
  # path <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/gliph/',date)
  # cd8.expanded <- read.csv(paste0(path, "/TopExpandedClones_CD8.tsv"), sep="\t")
  # cd8.public <- read.csv(paste0(path, "/PublicClones_CD8.tsv"), sep="\t")
  # cd4.expanded <- read.csv(paste0(path, "/TopExpandedClones_CD4.tsv"), sep="\t")
  # cd4.public <- read.csv(paste0(path, "/PublicClones_CD4.tsv"), sep="\t")

  input.path = paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_",date)
  # cd4
  f.public.donors <- data.table(read.csv(paste0(input.path, "/.intermidiate_file_TCRsimilarity_CD4.tsv"), sep="\t"))
  f.public.donors <- f.public.donors[!is.na(cdr3b.aa.seq)]; f.public.donors[, cdr3a.aa.seq:=NULL]
  beta.info.cd4 <- f.public.donors[
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
  # cd8
  f.public.donors <- data.table(read.csv(paste0(input.path, "/.intermidiate_file_TCRsimilarity_CD8.tsv"), sep="\t"))
  f.public.donors <- f.public.donors[!is.na(cdr3b.aa.seq)]; f.public.donors[, cdr3a.aa.seq:=NULL]
  beta.info.cd8 <- f.public.donors[
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

  cd4_cd8.clonotypes <- unique(unlist(
    lapply(list(beta.info.cd4$clonotype.tag, beta.info.cd8$clonotype.tag), function(x){
      unlist(str_split(x, pattern=";"))
    })
  ))# Get all the clonotypes from fom CD4 and CD8

  # ----- Get the overlap with the clonotypes in C17
  # c17 mdata filtered by overlaping clonotypes
  c17_mdata_cd4_cd8_clones <- sc_cd3p_c17_mdata[!is.na(sc_cd3p_c17_mdata$clonotype.tag) & sc_cd3p_c17_mdata$clonotype.tag %in% cd4_cd8.clonotypes, c("barcode", "clonotype.tag", "TRB.aa.chains.tag")]
  # cd4_cd8.Info filtered by overlaping clonotypes
  cd4_cd8.Info <- rbind(beta.info.cd4[,c("cdr3b.aa.seq", "trb.v", "trb.j", "clonotype.tag", "size")], beta.info.cd8[,c("cdr3b.aa.seq", "trb.v", "trb.j", "clonotype.tag", "size")] )
  selection <- which(grepl(paste(c(paste0(unique(c17_mdata_cd4_cd8_clones$clonotype.tag), "$"), paste0(unique(c17_mdata_cd4_cd8_clones$clonotype.tag), ";")),collapse="|"), cd4_cd8.Info$clonotype.tag))
  cd4_cd8.Info <- cd4_cd8.Info[selection,];
  # Obtain if they are repeated in CD4 or CD8
  cd4_cd8.Info$selection <- selection
  cells4 <- dim(beta.info.cd4)[1]; #cells8 <- dim(beta.info.cd4)[1];
  cd4_cd8.Info <- mutate(cd4_cd8.Info, selection = case_when(
      selection <= cells4 ~ "CD4",
      selection > cells4 ~ "CD8"
  ))
  # Order by size
  cd4_cd8.Info <- cd4_cd8.Info[order(cd4_cd8.Info$size, decreasing=T),]
  # If a clonotype is in CD4 and Cd8, display it.
  cd4_cd8.Info <- cd4_cd8.Info[
    ,.(
      clonotype.tag,
      size=paste0(
        unique(sort(size)), collapse=';'),
      selection=paste0(
        unique(sort(selection)), collapse=';')
    ),
    by=.(cdr3b.aa.seq, trb.v, trb.j)
  ] %>% distinct()

  c17_mdata_cd4_cd8_clones # cells
  cd4_cd8.Info # GLIPH clonotypes

  # Merge: Add info (v gene, j gene, merged clonotypes) to c17 metadata
  order_ <- vapply(c17_mdata_cd4_cd8_clones$clonotype.tag, function(x){grep(paste(x,c(";","$"), sep="", collapse= "|"), cd4_cd8.Info$clonotype.tag)[1]}, numeric(1))
  c17_mdata_cd4_cd8_clones$trb.v <- cd4_cd8.Info$trb.v[order_]
  c17_mdata_cd4_cd8_clones$trb.j <- cd4_cd8.Info$trb.j[order_]
  c17_mdata_cd4_cd8_clones$clonotype.tags <- cd4_cd8.Info$clonotype.tag[order_]
  c17_mdata_cd4_cd8_clones$size <- cd4_cd8.Info$size[order_]
  c17_mdata_cd4_cd8_clones$celltype <- cd4_cd8.Info$selection[order_]
  c17_mdata_cd4_cd8_clones

  results.path=paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/gliph/',date)
  fwrite(file=paste0(results.path,"/c17_cells_in_CD4_CD8.csv"), x=c17_mdata_cd4_cd8_clones, sep=',', na=NA, quote=FALSE)



########## GLIPH group per donor

if(0){
  require(tidyverse)
  library(ggplot2)
  library(scales)
  library(readr)

  celltype <- "CD8"#"CD4" #
  date <- "2022-04-08"#  "2022-03-30" "2022-02-07"

  ifelse(celltype == "CD8", iFile <- "GLIPH2Input_TCR-Data_All-TRB_CD8_wAlphaChain.tsv", iFile <- "GLIPH2Input_TCR-Data_All-TRB_CD4GLIPH2Input_TCR-Data_All-TRB_CD8_wAlphaChain.tsv")
  ifelse(celltype == "CD8", oFile <- "GLIPH2Output_TCR-Data_CD8.csv", oFile <- "GLIPH2Output_TCR-Data_CD4.csv")

  iPath <- paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_", date, "/")
  oPath <- paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/gliph/", date, "/")
  Gi <- read.csv(paste0(iPath, iFile), sep = "\t")
  Go <- read.csv(paste0(oPath, oFile), sep = ",")

  donors <- c("BT25", "BT23")
  for (donor in donors){
    cat(donor)
    # Create separate tables for each donor
    donor_Gi <- Gi %>% filter(grepl(donor, Gi$donor.id)) %>% mutate(id = paste(cdr3b.aa.seq, trb.v, trb.j, sep=","))
    # Find the specificity group (if any) for each clonotype in each donor.
    donor_Go <- Go %>% filter(grepl(donor, Go$Sample)) %>% mutate(id = paste(TcRb, V, J, sep=","))
    # Add specificity group information
    donor_df <- merge(donor_Gi, donor_Go, all.x=T, by="id") %>% select("cdr3b.aa.seq", "trb.v", "trb.j", "cdr3a.aa.seq", "tra.v", "tra.j", "clonotype.tag", "donor.id", "pattern", "index", "size") %>% arrange(index, desc(size))

    # Save
    write_csv(donor_df, paste0(oPath, "allClonotypes_wGliph_info_", donor, "_", celltype, ".csv"))

    # Plot: percentage of cells with motif
    data <- data.frame(
      group=c("YES", "NO"),
      value=c(sum(!is.na(donor_df$index)), sum(is.na(donor_df$index)))
    )
    # Basic piechart
    blank_theme <- theme_minimal() +
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
    )
    p <- ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=1) + coord_polar("y", start=0) +
      blank_theme + labs(title = paste0(donor, " clonotypes with GLIPH specificity group in ", celltype)) +
      geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]),
            label = percent(value/sum(value))), size=5)

    # Save Plot
    pdf(paste0(oPath, "Clonotypes_wGliph_info_proportion_", donor, "_", celltype, ".pdf"))
    print(p)
    dev.off()
  }

}
