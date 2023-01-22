# PEDIATRIC BRAIN TUMORS
------------

Description
------------

Brain tumors are the most common solid tumors in children and outcomes remain dismal for a high proportion of patients. Immunotherapies potentiate anti-tumor responses of T cells, however translation of their clinical benefit to pediatric brain tumors has been hindered by a lack of understanding of immune responses within the brain and the relatively small disease population.   

In order to evaluate anti-tumor immune response in pediatric brain tumors, we performed single-cell RNA sequencing (scRNA-seq) and paired single-cell TCR sequencing (scTCR-seq) of patient-derived brain tumor-infiltrating T cells to couple T cell molecular program with TCR repertoire and clonality. We generated single-cell transcriptomic profiles of FACS-sorted T cells isolated from brain tumors of 24 pediatric patients with varying diagnoses and histological grades, those being: Pilocytic astrocytoma, Ganglioglioma, Diffuse astrocytoma, Craniopharryngioma, Choroid plexus papilloma, High grade glioma, Medulloblastoma and Anaplastic ependymoma.  

We also did scRNA-seq and scTCR-seq of patients (n=10) bearing NSCLC, a tumor type with high tumor mutational burden where ICBs have shown promising outcomes, to compare features of clonally expanded cells and TCR repertoir diversity with what we see in pediatric brain tumors.

This repository contains the scripts used to analyze the aforementioned samples focusing on: CD4+ and CD8+ T cells.  

Requirements
------------

This project was done using the following modules/programs:

* [R](https://cran.r-project.org/) (v3.6.1)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.1.0)
* [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)
* [Seurat](https://satijalab.org/seurat) (v3.1.5)

Raw data
------------
The single-cell RNA-seq raw and processed files can be downloaded through the following GEO accession number: [GSE221776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221776). 

For survival analysis, bulk RNA-se data from the Pediatric Brain Tumor Atlas (PBTA) were downloaded from the [Gabriella Miller Kids First Data Resource Portal](https://portal.kidsfirstdrc.org/login) through the [CAVATICA](https://www.cavatica.org/) cloud-based platform, and clinical data were accessed using [PedcBioPortal](https://pedcbioportal.kidsfirstdrc.org/).



Raw data pre-processing
------------

* To do the 10x demultiplexing and mapping just pull [our in-house pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.
* To do the donor the multiplexing just pull [our in-house pipeline](https://github.com/vijaybioinfo/ab_capture).
* To do the single-cell quality control just pull [our in-house pipeline](https://github.com/vijaybioinfo/quality_control).
* To do the doublet detection use [our in-house pipeline](https://github.com/vijaybioinfo/doublet_detection) using Scrublet. 
* To generate the clustering of single-cell data just pull [our in-house pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.
* To do the aggregation of VDJ libraries just pull [our in-house pipeline](https://github.com/vijaybioinfo/VDJ_aggr).

For more specific information about the data generation and processing, please check the "methods" section within the manuscript.  

> All relevant scripts all located in ./pre-processing  

Figures
------------
> All relevant scripts all located in ./figures

Downstream Analysis
------------
* DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)
* [GLIPH2](http://50.255.35.37:8080/)
> All relevant scripts all located in ./downstream_analysis


Usage & Citation
--------------

If you want to clone this repository run:
```bash
git clone https://github.com/vijaybioinfo/PBT_2023.git
```
Please cite the following manuscript if you are using this repository:


Maintainers
-----------

Current maintainers:
* Kevin Meza-Landeros (kmlanderos@lji.org) 
* Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org)

Vijayanand Lab.  
Division of Vaccine Discovery La Jolla Institute for Immunology La Jolla, CA 92037, USA

Contact
-----------
Please email Kevin Meza-Landeros (kmlanderos@lji.org), Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org) and/or Vijayanand Pandurangan (vijay@lji.org).
