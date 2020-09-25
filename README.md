# Maize Single Cell Cis Regulatory Atlas

This repository contains scripts used in the analysis of maize scATAC-seq atlas data. Scripts are provided as is. However, there may be an effort in the future to translate this code into an R package for streamlined analysis depending on demand. 


Users should be aware that the majority of this code was designed for our specific use case, and as such, has been written to analyze the data sets detailed in the manuscript. Thus, all users should carefully look over the code to better understand if the default parameters make sense for your experiment/species/conditions/etc. Most of the provided scripts were designed to run from the commandline. For example, 

```
Rscript pseudotime.R <sparse> <motif.deviations> <gene> <TFgene> <meta> <svd> <prefix> <config> <threads>
```

Inputs for sparse, gene, and TFgene are triplet format sparse matrices saved as text files. Examples of input files are provided in the **example_files** directory.

## Provided scripts
* **Regularized quasibinomial logistic regression, UMAP, and Louvain clustering**  
These scripts contain the main functions for clustering barcode x feature (windows/ACRs/peaks) sparse matrices, producing UMAP embeddings and Louvain neighborhood groupings.
	* snATAC_rr.R
	* snATAC_rr.cluster_Utils.R

* **Identification of co-accessible ACRs and estimates of gene accessibility scores**  
These scripts identify co-accessible ACRs using empirical FDR thresholds. Co-accessible ACRs are used tie distal ACRs to genes, which along with TSS localized ACRs, are used to estimate gene accessible scores predictive of gene expression via a graphical LASSO model. Initial identification of co-accessible ACRs and gene accessibility scores leverage the [cicero](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/) scATAC-seq framework. [1]
	* call_coaccessible.R
	* call_coaccessible_UTILs.R

* **Pseudotime analysis**  
Scripts used to order and estimate pseudotime for each cell. Differential accessibility of ACRs, genes, and motifs are determined via F-statistic based hypothesis testing with natural splines linear regression. Pseudotime estimates are derived from the SVD space and rely heavily on code developed by Ryan Corces and Jeffrey Granja from the Greenleaf lab [archR](https://www.archrproject.com). [2]
	* pseudotime_analysis.R

* **Pseudotime alignment**  
R scripts that align homologous (putative 1-1 orthologs) genes between *Zea mays* and *Arabidopsis thaliana* along a pseudotime trajectory via dynamic time warping. 
	* alignTrajectories.R


## Additional scripts
This repository contains the major scripts used in the analysis. Requests for additional code can be made in the issues tab of this repository, or by contacting me by email (marand@uga.edu). 

## Genome browser of cell-type resolved ATAC-seq profiles
In addition to raw and processed data available from NCBI GEO, we also provide coverage tracks of *Zea mays* and *Arabidopsis thaliana* cell-types on our Jbrowse genome browser. Cell-type aggregate chromatin accessibility profiles can be found under the following URLs:
* [*Zea mays*](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=zea_mays_v4&cat=Maize%20Epigenome&loc=8%3A171225459..171227370&tracks=genes&highlight=)
* [*Arabidopsis thaliana*](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=a_thaliana_tair10&loc=chr5%3A19883361..19903660&tracks=genes&highlight=)

Simply load the relevant tracks under the tab labeled **scATAC_celltypes**

## References
1. Pliner et al. (2018). [Cicero Predicts cis-Regulatory DNA interactions from Single-Cell Chromatin Accessibility Data](https://doi.org/10.1016/j.molcel.2018.06.044). *Molecular Cell*
2. Granja, Corces, et al. (2020). [ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis](https://www.biorxiv.org/content/10.1101/2020.04.28.066498v1). *bioRxiv*
