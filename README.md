# PEDIATRIC_BRAIN_TUMORS
------------

Description
------------

This repository contains the scripts used to analyze our single-cell RNA-seq samples coming from pediatric brain tumors; focusing on CD4+ and CD8+ T cells.


REQUIREMENTS
------------

This project was done using the following modules/programs:

* [R](https://cran.r-project.org/) (v3.6.1)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.1.0)
* [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)
* [Seurat](https://satijalab.org/seurat) (v3.1.5)

Data
------------
Raw data files and the metadata can be downloaded though the following GEO accession number: [GSE]()

Data pre-processing
------------

* To do the 10x demultiplexing and mapping just pull [our in-house pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.
* To do the donor the multiplexing just pull [our in-house pipeline](https://github.com/vijaybioinfo/ab_capture).
* To do the single-cell quality control just pull [our in-house pipeline](https://github.com/vijaybioinfo/quality_control).
* To do the doublet detection use [our in-house pipeline](https://github.com/vijaybioinfo/doublet_detection) using Scrublet. 
* To generate the clustering of single-cell data just pull [our pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.
* To to the aggregation of VDJ libraries just pull [our pipeline](https://github.com/vijaybioinfo/VDJ_aggr).

For more specific information about the data generation and processing, please check the methods.

Downstream Analysis
------------
* DGEA - You can follow [our DGEA pipeline](https://github.com/vijaybioinfo/dgea)
* Figures - Files inside `/scripts/figures`


Usage & Citation
--------------

If you want to clone this repository run:
```bash
git clone https://github.com/vijaybioinfo/PEDIATRIC_BRAIN_TUMORS.git
```
Please cite the following manuscript if you are using this repository:


MAINTAINERS
-----------

Current maintainers:
* Kevin Meza-Landeros (kmlanderos@lji.org) 
* Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org)

Vijayanand Lab.  
Division of Vaccine Discovery La Jolla Institute for Immunology La Jolla, CA 92037, USA

Contact
-----------
Please email Kevin Meza-Landeros (kmlanderos@lji.org), Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org) and/or Vijayanand Pandurangan (vijay@lji.org).
