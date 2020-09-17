###############################################################################
## Cicero trajectories
###############################################################################

# default variables
threads <- 1
dims <- 50

# commandline arguments
args = commandArgs(TRUE)
if(length(args)==4){
    
    # read in commandline arguments
    threads  <- as.numeric(args[1])
    output   <- as.character(args[2])
    input    <- as.character(args[3])
    metafile <- as.character(args[4])
        
}else{
    stop("wrong number of arguments")
}


###############################################################################
## Load data, functions, and wd
###############################################################################

# setwd
setwd(getwd())

# load libraries
source("call_coaccessible_UTILs.R")

# load data
message(" ... Loading data")
a <- read.table(input)
zm <- read.table("/scratch/apm25309/single_cell/ATACseq/v3/bedfiles/resources/Zm.v4.genome")
genes <- read.table("/scratch/apm25309/single_cell/ATACseq/v3/bedfiles/resources/Zm.geneAnnotation.bed", header=T)
meta <- read.table(metafile)


###################################################################################################
## Cluster cells 
###################################################################################################

# create cicero object
message(" ... Creating CDS")
a <- a[as.character(a$V2) %in% rownames(meta),]
a$V1 <- droplevels(a$V1)
a$V2 <- droplevels(a$V2)
shuf <- a
shuf$V1 <- shuf$V1[sample(length(shuf$V1))]
shuf$V2 <- shuf$V2[sample(length(shuf$V2))]
cds <- make_atac_cds(a, binarize=T)
shufcds <- make_atac_cds(shuf, binarize=T)

c.colSums <- Matrix::colSums(exprs(cds))
c.rowSums <- Matrix::rowSums(exprs(cds))
s.colSums <- Matrix::colSums(exprs(shufcds))
s.rowSums <- Matrix::rowSums(exprs(shufcds))

# check shuffled
message("   * orig. matrix = ",nrow(cds), " | ", ncol(cds))
message("   * shuf. matrix = ",nrow(shufcds), " | ", ncol(shufcds))
message("   * orig. matrix colSums = ", paste(c.colSums[1:5], collapse=", "))
message("   * shuf. matrix colSums = ", paste(s.colSums[1:5], collapse=", "))
message("   * orig. matrix rowSums = ", paste(c.rowSums[1:5], collapse=", "))
message("   * shuf. matrix rowSums = ", paste(s.rowSums[1:5], collapse=", "))

# add metadata, filter, and run TFIDF/library regression/batch effect removal
cds <- cds[,colnames(cds) %in% rownames(meta)]
shufcds <- shufcds[,colnames(cds) %in% rownames(meta)]
pData(cds) <- meta[colnames(exprs(cds)),]
pData(shufcds) <- meta[colnames(exprs(shufcds)),]
cds <- cds[Matrix::rowSums(exprs(cds))>0,]
cds <- cds[,Matrix::colSums(exprs(cds))>0]
shufcds <- shufcds[Matrix::rowSums(exprs(shufcds))>0,]
shufcds <- shufcds[,Matrix::colSums(exprs(shufcds))>0]

# process basic
cds <- detectGenes(cds)
shufcds <- detectGenes(shufcds)
cds <- estimateSizeFactors(cds)
shufcds <- estimateSizeFactors(shufcds)

# load results from jaccard
message(" ... Loading Jaccard-based clustering results and reduced dimensions")
cds <- loadMeta(cds, meta)
shufcds <- loadMeta(shufcds, meta)
    
    
###################################################################################################
## Estimate connections, modules and gene activities			       
###################################################################################################

# run cicero to get co-accessible sites and modules BY CLUSTER
meta2 <- pData(cds)
meta2$Cluster <- as.character(meta2$Cluster)
print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
cl <- makeSOCKcluster(threads)
registerDoSNOW(cl)
tasks <- length(clusts)
pb <- txtProgressBar(max = tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
package.labs <- c("cicero", "Matrix")
message(" ... Initializing per cluster cicero run")

# # run in parallel
# gact <- list()
# gact <- foreach(i=clusts, .combine='c', .packages=package.labs, .options.snow=opts) %dopar% {
#     
#     # get umap coordinates and make cicero CDS
#     message("###--- Creating cicero object, cluster",i, " ---###")
#     ids <- rownames(meta2[meta2$Cluster==i,])
#     index.keep <- colnames(exprs(cds)) %in% ids
#     s.cds <- cds[,index.keep]
#     
#     # only consider sites accessible in at least 1% of cells in cluster
#     s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>(ncol(exprs(s.cds))*0.00),]
#     s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
#     print(head(exprs(s.cds)[,1:5]))
#     message(" - number of sites for cluster ", i, " = ", nrow(s.cds))
#     
#     # get UMAP coordinates
#     umap_coords <- t(reducedDimA(s.cds))
#     umap_coords <- umap_coords[colnames(s.cds),]
#     message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
#     rownames(umap_coords) <- colnames(exprs(s.cds))
#     cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=15)
#     
#     # run cicero (connections)
#     message(" ... Running cicero")
#     conns <- run_cicero(cicero_cds, zm, window=500000, sample_num=100)
#     
#     # write results to disk
#     write.table(conns, file=paste("bycluster",i,".",output,".cicero.loops.txt",sep=""),
#                 sep="\t",quote=F, col.names=F, row.names=F)
# }
# close(pb)
# stopCluster(cl)
#   
# 
# ###################################################################################################    
# ## SHUFFLE ----------------------------------------------------------------------------------------
# ###################################################################################################
# 
# # for each cluster
# meta2 <- pData(shufcds)
# meta2$Cluster <- as.character(meta2$Cluster)
# print(table(meta2$Cluster))
# clusts <- unique(meta2$Cluster)
# cell_ids <- c()
# 
# # iterate
# its <- 0
# 
# # foreach parameters
# cl <- makeSOCKcluster(threads)
# registerDoSNOW(cl)
# tasks <- length(clusts)
# pb <- txtProgressBar(max = tasks, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# package.labs <- c("cicero", "Matrix")
# message(" ... Initializing shuffled per cluster cicero run")
# 
# ## shuffled ##
# gact <- list()
# gact <- foreach(i=clusts, .combine='c', .packages=package.labs, .options.snow=opts) %dopar% {
#     
#     # get umap coordinates and make cicero CDS
#     message("###--- Creating cicero object, cluster",i, " ---###")
#     ids <- rownames(meta2[meta2$Cluster==i,])
#     index.keep <- colnames(exprs(shufcds)) %in% ids
#     s.cds <- shufcds[,index.keep]
#     
#     # only consider sites accessible in at least 1% of cells in cluster
#     s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>(ncol(exprs(s.cds))*0.00),]
#     s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
#     print(head(exprs(s.cds)[,1:5]))
#     message(" - number of sites for cluster ", i, " = ", nrow(s.cds))
#     
#     # get UMAP coordinates
#     umap_coords <- t(reducedDimA(s.cds))
#     umap_coords <- umap_coords[colnames(s.cds),]
#     message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
#     rownames(umap_coords) <- colnames(exprs(s.cds))
#     cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=15)
#     
#     # run cicero (connections)
#     message(" ... Running cicero")
#     conns <- run_cicero(cicero_cds, zm, window=500000, sample_num=100)
#     
#     # write results to disk
#     write.table(conns, file=paste("shuffled_bycluster",i,".",output,".cicero.loops.txt",sep=""),
#                 sep="\t",quote=F, col.names=F, row.names=F)
#     
# }
# close(pb)
# stopCluster(cl)


###################################################################################################    
## COMPUTE GENE ACTIVITY --------------------------------------------------------------------------
###################################################################################################

# for each cluster
meta2 <- pData(cds)
meta2$Cluster <- as.character(meta2$Cluster)
print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
message(" ... Initializing per cluster cicero run - GENE ACTIVITY")
gascores <- mclapply(clusts, function(x){

    # load connections
    id.true <- paste("bycluster",x,".",output,".cicero.loops.txt",sep="")
    id.false <- paste("shuffled_bycluster",x,".",output,".cicero.loops.txt",sep="")
    t.conns <- read.table(id.true)
    s.conns <- read.table(id.false)
    
    # filter loops --------------------------------------------------------------------------------
    
    # empty vector
    b.sub <- c()
    lims <- seq(from=0, to=0.99, length.out=100)
    
    # find cut-off
    for(j in lims){
        b.sub <- c(b.sub, nrow(subset(s.conns, s.conns$V3 >= j)))
    }
    fdr <- b.sub/nrow(s.conns)
    threshold <- min(lims[which(fdr < 0.05)])
    message(" - threshold = ", threshold, " | ", id.true)
    
    # filter loops
    a.sub <- subset(t.conns, t.conns$V3 >= threshold)
    id <- gsub("bycluster", "filtered", id.true)
    write.table(a.sub, file=id, quote=F, row.names=F, col.names=F, sep="\t")
    colnames(a.sub) <- c("Peak1", "Peak2", "coaccess")
    
    # get gene activity scores --------------------------------------------------------------------
    message("--- estimating gene activity scores for cluster ",x)
    ids <- rownames(meta2[meta2$Cluster==x,])
    index.keep <- colnames(exprs(cds)) %in% ids
    s.cds <- cds[,index.keep]
    
    # only consider sites accessible in at least 1% of cells in cluster
    s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>0,]
    s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
    print(head(exprs(s.cds)[,1:5]))
    message(" - number of sites for cluster ", x, " = ", nrow(s.cds))
    
    # get UMAP coordinates
    umap_coords <- t(reducedDimA(s.cds))
    umap_coords <- umap_coords[colnames(s.cds),]
    message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
    rownames(umap_coords) <- colnames(exprs(s.cds))
    cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=15)
    
    # estimate gene activity
    message(" ... Estimating gene activity scores")
    pos <- subset(genes, strand == "+")
    pos <- pos[order(pos$start),] 
    pos <- pos[!duplicated(pos$transcript),]
    pos$end <- pos$start
    pos$start <- pos$start - 1000
    neg <- subset(genes, strand == "-")
    neg <- neg[order(neg$start, decreasing = T),] 
    neg <- neg[!duplicated(neg$transcript),] 
    neg$start <- neg$end
    neg$end <- neg$end + 1000
    
    # merge
    gene_ann2 <- rbind(pos, neg)
    gene_ann2 <- gene_ann2[,c(1:3, 8)]
    gene_ann2 <- gene_ann2[order(gene_ann2$start, decreasing=F),]
    colnames(gene_ann2)[4] <- "gene"
    
    # annotate genes
    message("     - annotate genes by peaks ...")
    s.cds <- annotate_cds_by_site(s.cds, gene_ann2, all=F)
    
    # estimate un-normalized activity
    message("     - build gene activity matrix ... ")
    unnorm_ga <- build_gene_activity_matrix(s.cds, a.sub, dist_thresh = 500000, coaccess_cutoff = 0)
    unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]
    
    # gene activity per cluster
    num_genes <- pData(s.cds)$num_genes_expressed
    names(num_genes) <- row.names(pData(s.cds))
    
    # normalize
    cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
    geneact <- as.data.frame(summary(cicero_gene_activities))
    geneact$i <- rownames(cicero_gene_activities)[geneact$i]
    geneact$j <- colnames(cicero_gene_activities)[geneact$j]
    
    # output
    write.table(geneact, file=paste("filtered",x,".",output,".cicero.geneActivity.txt",sep=""),
                sep="\t",quote=F, col.names=F, row.names=F)
    
    # return
    return(unnorm_ga)
}, mc.cores=threads)

