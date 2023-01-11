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
### -------------------- Results visualization --------------------- ###
cat('### -------------------- Results visualization --------------------- ###\n\n')

library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(igraph)
set.seed(42)

# ---> Run GLIPH2 using the web platform.
date <- "2022-07-15" # date <- "2022-02-07" #
celltype <- "CD8"
hinojosa <- FALSE #TRUE#
LG <- FALSE #TRUE#
expansion.thold <- 1 #greater or equal
HTO <- TRUE

if (HTO == TRUE){
  hto_suffix <- "_hto"
} else { hto_suffix <- "" }

# CD8
if(celltype == "CD8"){
  min_dotsize <- 3; max_dotsize <- 9
  legend_breaks <- c(2,20,40,60) # Do not include 1
  color_legend_breaks <- c(2,6,11) # Do not include 1
  cs_upper_limit <- 60 # All the values higher (clone size value) than this will be asigned this value
} else if(celltype == "CD4"){
    min_dotsize <- 3; max_dotsize <- 9
    legend_breaks <- c(2,10,20,40) # Do not include 1
    color_legend_breaks <- c(2,6,13) # Do not include 1
    cs_upper_limit <- NULL # All the values higher (clone size value) than this will be asigned this value
}


input.path = paste0("/home/kmlanderos/kmlanderos/pbtumor-all/results/tcr/preproces_4_Gliph2/CD3p/ag_specificity_assessment/reports/ag_specificity_assessment_",date)
results.path = paste0('/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/pbtumor-all/results/tcr/gliph/',date)
results.path = paste0('/home/kmlanderos/tmp_large/pbtumor-all/results/tcr/gliph/',date); dir.create(results.path, recursive=T)

# Get number of expanded specificity groups
suffix <- c("CD4","CD8"); if(HTO) suffix <- paste0(suffix,"_hto")
tmp.file.name <- paste0(results.path, '/GLIPH2Output_TCR-Data_', suffix, '.csv')
gliph_outs <- lapply(tmp.file.name, function(x) {
  gliph.res <- fread(file=x, blank.lines.skip=TRUE)
  gliph.res <- gliph.res[pattern!='single']
})

exp_sg <- lapply(gliph_outs, function(x) {
  tmp <- x %>% select(index, pattern, TcRb, Sample, Freq, expansion_score) %>% group_by(index) %>% filter(expansion_score < 0.05) %>% # %>% as.data.frame() %>% head(100)
    select(index, pattern) %>% distinct() %>% nrow()
}); names(exp_sg) <- suffix
exp_sg


# ---> Analaysis of the results provided by GLIPH2 through iGraph.
# @ Program as a function.
check.gliph2.output <- function(file.suffix, beta.info, hinojosa = FALSE, min_dotsize = 2, max_dotsize = 10, legend_breaks = NULL, cs_upper_limit = NULL, color_legend_breaks = NULL){
  # @ Load data as removing the 'single' patterns.
  # Define results of interest.
  tmp.file.name <- paste0(results.path, '/GLIPH2Output_TCR-Data_', file.suffix, '.csv')
  gliph.res <- fread(file=tmp.file.name, blank.lines.skip=TRUE)
  gliph.res <- gliph.res[pattern!='single']

  # @ Further preprocessing.
  tmp.data.1 <- merge(x=gliph.res, y=beta.info[, .(TcRb=cdr3b.aa.seq, V=trb.v, J=trb.j)], by=c('TcRb', 'V', 'J'), all.x=TRUE, all.y=FALSE, sort=FALSE)
  # Define summary information for every pattern.
  #     Average clone size.
  tmp.data.2 <- tmp.data.1[, .(size.median=median(Freq), size.sum=sum(Freq), nTCRs = .N), by=pattern]
  # Re-define clone size upper limit
  if(!is.null(cs_upper_limit)) tmp.data.2[size.sum>cs_upper_limit, size.sum:=cs_upper_limit]

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
  } else{
    gliph.res.tmp <- separate(gliph.res, col="Sample", into=c("donor","condition"), sep=":")
    gliph.res.tmp <- distinct(gliph.res.tmp[,c("pattern","donor")])
    gliph.res.tmp$value <- 1
    gliph.res.tmp <- pivot_wider(gliph.res.tmp, names_from = donor, values_from = value)
    gliph.res.tmp[is.na(gliph.res.tmp)] <- 0
    gliph.res.tmp[,"NA"] <- NULL # You are not sure that comes from another donor
    motif_donor <- apply(gliph.res.tmp, 1, function(x){
      suma <- sum(as.numeric(x[-1]))
      if(suma==1) donor <- "non_public"
      if(suma==2){
        tmp <- colnames(gliph.res.tmp)[x==1]
        if(tmp[2] == paste0(tmp[1], ";NA")){
          donor <- "non_public"
        } else{
          donor <- "public"
        }
      }
      if(suma>2) {
        donor <- "public"
      }
      # if(donor) cat(sum); cat(x)
      donor
    })
    tmp.data.2$donor <- motif_donor
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
  # Get min and max value of the node clone size
  min.val <- min(V(igraph.obj)$size.sum); max.val <- max(V(igraph.obj)$size.sum)
  if(is.null(legend_breaks)){
    legend_breaks <- pretty(n=5, x=nodes.info$size.sum); legend_breaks[1] <- 2 # range(nodes.info$size.sum) # max.val <- max(legend_breaks)
  }
  legend_breaks_scaled <- ( (legend_breaks - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # (legend_breaks*8/max.val)+2 # Scale so that it goes from 2 to 10.

  V(igraph.obj)$size <- ( (V(igraph.obj)$size.sum - min.val) * (max_dotsize - min_dotsize) / (max.val - min.val) ) + min_dotsize # Scale so that it goes from 2 to 10. #V(igraph.obj)$size <- (V(igraph.obj)$size.sum*8/max.val)+2

  # Colors for the vertices according to their presence in the pbtumos datasets.
  if(hinojosa){
    colrs.v <- c(shared="gold", Uphadye="tomato", Hinojosa='green')
    V(igraph.obj)$color = colrs.v[V(igraph.obj)$study]# Get general things for plotting.
  } else{
    # Coloring by amount of TCRs
    palette <- colorRampPalette(c('yellow','red'))
    if(!is.null(color_legend_breaks)) {
      V(igraph.obj)$color <- palette(length(color_legend_breaks))[as.numeric(cut(V(igraph.obj)$nTCRs,breaks=c(1,color_legend_breaks)))]
      color_tags_color <- palette(length(color_legend_breaks))
    } else{
      fine <- 100
      V(igraph.obj)$color <- palette(fine)[as.numeric(cut(V(igraph.obj)$nTCRs,breaks=fine))]
      color_tags_color <- c(palette(fine)[1],palette(fine)[50],palette(fine)[100])
    }
  }

  if(!is.null(color_legend_breaks)){
    color_tags <- c(as.character(color_legend_breaks[1]), paste0(as.character(color_legend_breaks[-length(color_legend_breaks)]+1), "-", as.character(color_legend_breaks[-1])))
  }else{
    min_nTCR <- as.character(min(V(igraph.obj)$nTCRs))
    max_nTCR <- as.character(max(V(igraph.obj)$nTCRs))
    medium_nTCR <- as.character(round(as.numeric(max_nTCR)/2))
    color_tags <- c(min_nTCR,medium_nTCR,max_nTCR)
  }

  # Plot.
  tmp.file.name <- paste0(results.path, '/GLIPH2Results_iGraph_', file.suffix, '.pdf')
  pdf(file=tmp.file.name, 10, 7)
  plot(igraph.obj, alpha=0.6, vertex.label=NA, main= paste0("TCR motifs found by GLIPH2 in ", celltype, "+ Tcells")) #, vertex.label.cex = 0.6, vertex.label.dist=2
  a <- legend('bottomleft',legend=unique(legend_breaks),pt.cex=legend_breaks_scaled*10/2000,col='white',
          pch=21, pt.bg='white', title='clone size', bty='n', yjust = 10)
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  y <- y - seq(0, 0.1, length.out = length(y))
  symbols(x,y,circles=legend_breaks_scaled*10/2000,inches=FALSE,add=TRUE,bg='white')
  legend(x='bottomright', legend = color_tags, pch=21, col=color_tags_color, pt.bg=color_tags_color, pt.cex=1, cex=.8, ncol=1, bty='n', title='# clonotypes')
  dev.off()
}

# fname <- paste0("/.GLIPH2Input_TCR-Data_All-TRB_", celltype); if(hinojosa) paste0(fname,"_Hinojosa")
fname <- paste0("/GLIPH2Input_TCR-Data_All-TRB_", celltype, hto_suffix); if(hinojosa) fname <- paste0(fname,"_Hinojosa")
gliph.all.beta.info <- read.csv(paste0(input.path, fname, ".tsv"), sep="\t")
suffix <- celltype; if(hinojosa) suffix <- paste0(suffix,"_Hinojosa"); if(LG) suffix <- paste0(suffix,"_LG"); if(HTO) suffix <- paste0(suffix,"_hto")
# Run
check.gliph2.output(file.suffix=suffix, beta.info=data.table(gliph.all.beta.info), hinojosa = hinojosa, min_dotsize = min_dotsize, max_dotsize = max_dotsize, legend_breaks = legend_breaks, cs_upper_limit = cs_upper_limit, color_legend_breaks = color_legend_breaks)
