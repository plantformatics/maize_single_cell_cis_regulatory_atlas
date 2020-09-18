# Maize Single Cell Cis Regulatory Atlas

This repository contains scripts used in the analysis of maize scATAC-seq atlas data. Scripts are provided as is. However, there may be an effort in the future to translate this code into an R package for streamlined analysis depending on demand. 


Users should be aware that the majority of this code was designed for our specific use case, and as such, has been written to analyze the data sets detailed in the manuscript. Thus, all users should carefully look over the code to better understand if the default parameters make sense for your experiment/species/conditions/etc. Most of the provided scripts were designed to run from the commandline. For example, 

```
Rscript pseudotime.R <sparse> <motif.deviations> <gene> <TFgene> <meta> <svd> <prefix> <config> <threads>
```

Inputs for sparse, gene, and TFgene are triplet format sparse matrices saved as text files. Examples of input files are provided in the **example_files** directory.

## Provided scripts
* Quasibinomial logistic regression
	* snATAC_rr.R
	* snATAC_rr.cluster_Utils.R

* Identification of co-accessible ACRs and estimates of gene accessibility scores
	* call_coaccessible.R
	* call_coaccessible_UTILs.R

* Pseudotime analysis
	* pseudotime_analysis.R

* Pseudotime alignment
	* alignTrajectories.R


## Additional scripts
This repository contains the major scripts used in the analysis. Requests for additional code can be made in the issues tab of this repository, or by contact me my email (marand@uga.edu). 

## Genome browser of cell-type resolved ATAC-seq profiles
In addition to raw and processed data available from NCBI GEO, we also provide coverage tracks of *Zea mays* and *Arabidopsis thaliana* cell-types on our Jbrowse genome browser. Cell-type aggregate chromatin accessibility profiles can be found under the following URLs:
* [maize](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=zea_mays_v4&cat=Maize%20Epigenome&loc=8%3A171225459..171227370&tracks=genes&highlight=)
* [arabidopsis](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=a_thaliana_tair10&loc=chr5%3A19883361..19903660&tracks=genes&highlight=)

Simply load the relevant tracks under the tab labeled **scATAC_celltypes**
