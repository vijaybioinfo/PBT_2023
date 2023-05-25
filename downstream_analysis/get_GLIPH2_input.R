############    -  Assessment on antigen specificity  -    ############
############    --------   for the TCR data    --------    ############

# ---
# Author: Kevin Meza Landeros
# Date: 2023-01-10
# ---

### -------------------------- Description -------------------------- ###
# Based on both gene expression data (provided, for now, only as a seurat object) and TCR data (provided as standard 10x vdj files), this program will create:
#   The input files to run GLIPH2 over the TCR data of interest.

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
  the minimum clone size you wish the clonotipes to have in order to get GLIPH2 input files."),
  make_option(opt_str="--HTO", type="logical", default=FALSE, dest="HTO", help="Boolean flag saying if you want to calculate
  the clone size only considering cells with Gex and Hashtag (donor). Default = Cells with Gex.")
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
HTO <- opt$HTO

if(interactive()){
  celltype <- "CD8"
  date <- "2023-04-19"
  LG <- FALSE
  hinojosa <- FALSE
  iGraph <- FALSE
  expansionThold <- 1
  HTO <- TRUE
}

# ---> General definitions.
# @ On the project.
gen.prj <- 'pbtumor-all-Batch2'
paper.lab <- '/results/tcr/preproces_4_Gliph2/CD3p/'
this.assess <- 'ag_specificity_assessment'

# ---> Path definitions.
reports.date <- date
gen.reports.path <- paste0('/mnt/bioadhoc-temp/Groups/vd-vijay/kmlanderos/', gen.prj, paper.lab)
gen.reports.path <- paste0(gen.reports.path, this.assess, '/reports')
data.path <- str_replace(string=gen.reports.path, pattern='reports$', replacement='data')
reports.path <- paste0(gen.reports.path, '/', this.assess, '_', if(is.null(reports.date)) Sys.Date() else reports.date)
if(!dir.exists(reports.path)) dir.create(reports.path, recursive=TRUE)

# ---> File definitions.
if ( celltype == "CD4") {
  folder <- "CD4"
} else if ( celltype == "CD8") {
  folder <- "CD8"
}

# Gene expression data. # We will only use the metadata of the Seurat object
seurat.obj.file <- paste0('/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/tcr/', folder, '/metadata.rds') #'/mnt/beegfs/kmlanderos/CD45pCD3p_clean2_object_lock_mean0.01_pct15_pc15.rds'
# TCR data.
cells.clons.info.file <- paste0('/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/tcr/', folder, '/filtered_contig_annotations_aggr.csv')
clons.info.file <- paste0('/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/tcr/', folder, '/clonotypes_aggr.csv')

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Gene expression info.
seurat.obj <- readRDS(file=seurat.obj.file)
seurat.obj <- seurat.obj[,-c(66:67)] # NOTE: Remove if necesary
if (HTO == TRUE){
  seurat.obj <- seurat.obj %>% filter(!is.na(orig.donor))
  hto_suffix <- "_hto"
  cat("Calculating clone size by only considering cells with Gex and Hashtag (donor).\n")
} else {hto_suffix <- ""}
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
      size=.N
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
fname=paste0("/.intermidiate_file_TCRsimilarity_",celltype, hto_suffix,".tsv")
fwrite(file=paste0(reports.path,fname), x=gliph.all, sep='\t', na=NA, quote=FALSE)


cat('\n\n')
### ------------------- Get actual GLIPH2 inputs ------------------- ###
cat('### ------------------- Get actual GLIPH2 inputs ------------------- ###\n')

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
    tmp.file.name <- paste0(reports.path, '/', file.prefix, hto_suffix, '.tsv')
    fwrite(file=tmp.file.name, x=tmp.data, sep='\t', na=NA, quote=FALSE)
    # fwrite(file=gsub("GLIPH2Input", ".GLIPH2Input", gsub(".tsv","_beta.info.tsv",tmp.file.name)), x=beta.info, sep='\t', na=NA, quote=FALSE) # Save this file for later iGraph visualization.
  }
  # Save another file that has the alpha chain and V- J- information.
  tmp.data.alpha <- beta.info[size>=tmp.expansion.thold, .(cdr3b.aa.seq, trb.v, trb.j, cdr3a.aa.seq, tra.v, tra.j, clonotype.tag, donor.id=paste(orig.donor, 'Exp1', sep=':'), size)] %>% mutate(id = paste(cdr3b.aa.seq, trb.v, trb.j, sep=","))
  tmp.data.alpha$cdr3a.aa.seq <- as.character(tmp.data.alpha$cdr3a.aa.seq); tmp.data.alpha[is.na(tmp.data.alpha$cdr3a.aa.seq),"cdr3a.aa.seq"] <- "NA"
  tmp.data.alpha <- tmp.data.alpha %>% group_by(id) %>% mutate(cdr3a.aa.seq.all = paste(cdr3a.aa.seq, collapse = ";")) %>% distinct(id, .keep_all= TRUE)
  if(!is.na(file.prefix)){
    tmp.file.name <- paste0(reports.path, '/', file.prefix, hto_suffix, '_wAlphaChain.tsv')
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
