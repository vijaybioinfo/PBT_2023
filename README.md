# PEDIATRIC_BRAIN_TUMORS
------------

This repository contains the scripts used to analyze our samples pediatric brain tumors focused on T cells.
The data contains RNA-seq samples from single-cell.

REQUIREMENTS
------------

This project requires the following modules/programs:

* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.1.0)
* [Seurat](https://satijalab.org/seurat) (v3.1.5)

INSTALLATION
------------

* You can follow the instructions for [our mapping pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.
* For the single-cell quality control just pull [our in-house script](https://github.com/vijaybioinfo/quality_control).
* Clustering of single-cell data with [our pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.

Global description
------------

*Demultiplexing libraries*: Cell Ranger was used to demultiplex the 10x libraries.

*Quality control*: An in-house [script](https://github.com/vijaybioinfo/quality_control)
was used to explore the quality of the data and select the thresholds.

*Clustering*: Seurat was used to cluster the data.

For more specific information about the data generation and processing, please check the methods.

MAINTAINERS
-----------

Current maintainers:
* Kevin Meza-Landeros (kmlanderos@lji.org)

Contact
-----------
Please email Kevin Meza-Landeros (kmlanderos@lji.org) and Vijayanand Pandurangan (vijay@lji.org).
