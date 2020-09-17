###################################################################################################
## Cicero co-accessibility analysis ---------------------------------------------------------------
###################################################################################################

# verbose
message("")
message("==========================================")
message("    running Cicero/trajectory analysis    ")
message("==========================================")
message("")

# libraries 
library(cicero)
library(Matrix)
library(parallel)
library(doSNOW)
library(methods)
library(tcltk)
library(iterators)
library(itertools)

# add annotation
loadMeta <- function(cds, jac){
    
    # svd and raw
    ids <- colnames(cds)
    jac <- jac[ids,]
    jac$UMAP1 <- jac$umapsub_1
    jac$UMAP2 <- jac$umapsub_2
    ids.2 <- rownames(jac)
    cds <- cds[,colnames(cds) %in% ids.2]
    umap.d <- t(jac[,c("UMAP1","UMAP2")])
    
    # UMAP output
    cds@reducedDimA <- umap.d
    cds@reducedDimS <- umap.d
    cds@dim_reduce_type <- "UMAP"
    pData(cds)$tissue_cluster <- jac$tissue_cluster
    pData(cds)$Cluster <- as.factor(jac$tissue_cluster)
    pData(cds)$preclust <- as.factor(pData(cds)$tissue_cluster)
    
    # return data
    return(cds)
}