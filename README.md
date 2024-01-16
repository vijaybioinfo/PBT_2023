# Intra-tumoral T cells in pediatric brain tumors display clonal expansion and effector properties
------------

Description
------------

Brain tumors are the most common solid tumors in children and outcomes remain dismal for a high proportion of patients. Immunotherapies potentiate anti-tumor responses of T cells, however, translation of their clinical benefit to pediatric brain tumors has been hindered by a lack of understanding of immune responses within the brain and the relatively small disease population.   

In order to evaluate anti-tumor immune response in pediatric brain tumors (PBT), we performed single-cell RNA sequencing (scRNA-seq) and paired single-cell TCR sequencing (scTCR-seq) of patient-derived brain tumor-infiltrating T cells to couple T cell molecular program with TCR repertoire and clonality. We generated single-cell transcriptomic profiles of FACS-sorted T cells isolated from brain tumors of 38 pediatric patients with varying diagnoses and histological grades, those being: Pilocytic astrocytoma, Ganglioglioma, Diffuse astrocytoma, Craniopharyngioma, Choroid plexus papilloma, High grade glioma, Medulloblastoma, Anaplastic ependymoma, Meningioma and Embryonal Tumor with Multilayered Rosettes.

We proved that tumor-infiltrating lymphocytes (TILs) expressing checkpoint transcripts (PDCD1, LAG3), also express effector molecules (IFNG), display cytotoxic functions (GZMA, GZMB), and are capable of proliferation. Therefore we suggest that these cells might be a good target for immune checkpoint blockade (ICB) and promote their reinvigoration.

In order to prove the importance of neoantigen-reactive T cell gene signatures in survival outcomes, we used whole tumor bulk RNA-seq from pediatric tumors and assessed the survival status of patients with high and low values for the gene signatures across time.

We then compared proportions of CD4+ TREGs and CD4-CTLs among total CD4+ tumor-infiltrating lymphocytes from our pediatric cohort with what is observed in adult brain tumors (n=26), which have been described as non-responsive upon ICB therapy. Finally, we used scRNA-seq data of TILs from patients with baseline NSCLC (n=10) or from PDCD1-responsive patients with basal cell carcinoma (BCC; baseline tumors), to compare features of clonally expanded cells with what we see in pediatric brain tumors.

This repository contains the data and scripts used to analyze the samples mentioned above.

Requirements
------------

This project was done using the following modules/programs:

* [R](https://cran.r-project.org/) (v3.6.1)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.1.0)
* [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)
* [Seurat](https://satijalab.org/seurat) (v3.1.5)
* [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) (v1.2.1)

Raw data
------------
* The single-cell RNA-seq raw and processed files from PBT and NSCLC (our study) can be downloaded through the following GEO accession number: [GSE221776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221776) 

* For survival analysis, previously published bulk RNA-seq data from the Pediatric Brain Tumor Atlas (PBTA) were downloaded from the [Gabriella Miller Kids First Data Resource Portal](https://portal.kidsfirstdrc.org/login) through the [CAVATICA](https://www.cavatica.org/) cloud-based platform, and clinical data were accessed using [PedcBioPortal](https://pedcbioportal.kidsfirstdrc.org/). For more specific information about the data processing, please check the "methods" section within the manuscript.  

*  Previously published scRNA-seq data (Mathewson et.al 2021) of adult high-grade gliomas can be downloaded using the following GEO accession number: [GSE163108](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163108)

* Previously published scRNA-seq data (Yost et.al 2019) derived from patients with basal cell carcinoma (BCC) pre-anti-PD-1 therapy can be downloaded using the following GEO accession number: [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813)


Raw data pre-processing  
------------

### Single-cell
* To do the 10x demultiplexing and mapping pull [our in-house pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.
* To do the donor demultiplexing pull [our in-house pipeline](https://github.com/vijaybioinfo/ab_capture).
* To do the single-cell quality control pull [our in-house pipeline](https://github.com/vijaybioinfo/quality_control).
* To do the doublet detection use [our in-house pipeline](https://github.com/vijaybioinfo/doublet_detection) using Scrublet. 
* To generate the clustering of single-cell data just pull [our in-house pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.
* To do the aggregation of VDJ libraries pull [our in-house pipeline](https://github.com/vijaybioinfo/VDJ_aggr).  

> Relevant scripts are located in: ./pre-processing  

For more specific information about the data generation and processing, please check the "methods" section within the manuscript.  

### Bulk   
Patients diagnosed with High-Grade Gliomas (HGG; n=120) with records of age (<20 years old) and sex were considered for the survival analysis. RNA-seq samples of the donors classified as metastasis or with unavailable TUMOR TYPE information were filtered out; samples derived from cell lines or obtained from the spine region were also excluded. The final list of donor samples can be found in: ./pre-processing 


Figures
------------
> Relevant scripts are located in: ./figures


Downstream Analysis
------------
* DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)
* [GLIPH2](http://50.255.35.37:8080/) - Clustering of CDR3β sequences
* [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) - TCR diversity estimation
> Relevant scripts all located in: ./downstream_analysis


Usage & Citation
--------------

If you want to clone this repository run:
```bash
git clone https://github.com/vijaybioinfo/PBT_2023.git
```
Please cite the following manuscript if you are using this repository:  
Upadhye, A., Meza Landeros, K.E., Ramírez-Suástegui, C. et al. Intra-tumoral T cells in pediatric brain tumors display clonal expansion and effector properties. Nat Cancer (2024). https://doi.org/10.1038/s43018-023-00706-9

Maintainers
-----------

Current maintainers:
* Kevin Meza Landeros (kmlanderos@lji.org) 
* Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org)

Vijayanand Lab.  
Division of Vaccine Discovery La Jolla Institute for Immunology La Jolla, CA 92037, USA


Contact
-----------
Please email Kevin Meza Landeros (kmlanderos@lji.org), Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org) and/or Vijayanand Pandurangan (vijay@lji.org).
