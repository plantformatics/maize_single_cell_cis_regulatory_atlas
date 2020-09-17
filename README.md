# maize_single_cell_cis_regulatory_atlas

This repository contains scripts used in the analysis of maize scATAC-seq atlas data. Scripts are provided as is. However, there may be an effort in the future to translate this code into an R package for streamlined analysis. 

## Provided scripts
* Quasibinomial logistic regression
	* snATAC_rr.R
	* snATAC_rr.cluster_Utils.R

* Pseudotime analysis
	* pseudotime_analysis.R

* Pseudotime alignment
	* alignTrajectories.R


## Additional scripts
Requests for additional code can be made in the issues tab of this repository. 

## Genome browser of cell-type resolved ATAC-seq profiles
In addition to raw and processed data available from NCBI GEO, we also provide coverage tracks of *Zea mays* and *Arabidopsis thaliana* cell-types on our Jbrowse genome browser. Cell-type aggregate chromatin accessibility profiles can be found under the following URLs:
* [maize](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=zea_mays_v4&cat=Maize%20Epigenome&loc=8%3A171225459..171227370&tracks=genes&highlight=)
* [arabidopsis](http://epigenome.genetics.uga.edu/PlantEpigenome/?data=a_thaliana_tair10&loc=chr5%3A19883361..19903660&tracks=genes&highlight=)

Simply load the relevant tracks under the tab labeled **scATAC_celltypes**
