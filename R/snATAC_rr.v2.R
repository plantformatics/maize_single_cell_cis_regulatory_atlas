## cluster scATAC ##

###################################################################################################
# global vars -------------------------------------------------------------------------------------
###################################################################################################
n.pcs       <- 50
k.near      <- 15
min.reads   <- 1e6
m.clst      <- 50
batch       <- NULL
sub.batch   <- NULL
subpeaks    <- NULL
threshold   <- 4
res         <- 0.05
k.nearC     <- 15
type        <- "response"
theta       <- 2
sub.theta   <- 2
lambda      <- 0.5
sub.lambda  <- 0.5
stdize      <- T
center      <- F
cap         <- F
modtype     <- "regModel"
scaleVar    <- F
useTop      <- NULL
useTop.th   <- NULL
doPlot      <- T
weights     <- NULL
clustOB     <- "umap"
link        <- "logit"
method      <- "elasticNet"
alpha       <- 0.5
variates    <- "~log10nSites"
r.variates  <- "~log10nSites"
var.th      <- 0.5
subres      <- 0.1
subthresh   <- 4
subK        <- 15
s.pcs       <- 30
sub.mreads  <- 5e5
sub.quant   <- 0
reclustType <- 2
sub.do.Clean <- T
do.Clean    <- F
clustType   <- "louvain"
stdLSI      <- F
sub.stdLSI  <- F
doL2        <- 1
sub.doL2    <- F
e.thresh    <- 5


###################################################################################################
# arguments ---------------------------------------------------------------------------------------
###################################################################################################

# load arguments
args <- commandArgs(T)
if(length(args)!=5){stop("Rscript snATAC_cluster.R <sparse> <meta> <prefix> <config> <nthreads>")}

# vars
input  <- as.character(args[1])
meta   <- as.character(args[2])
prefix <- as.character(args[3])
config <- as.character(args[4])
ncpus  <-   as.numeric(args[5])


###################################################################################################
# libraries ---------------------------------------------------------------------------------------
###################################################################################################
library(Matrix)
library(irlba)
library(harmony)
library(uwot)
library(tcltk)
library(iterators)
library(itertools)
library(parallel)
library(Seurat)
library(FNN)
library(RColorBrewer)
library(scales)
library(viridis)
library(mclust)
library(gplots)
library(glmnet)
library(dplyr)
library(gtools)
library(densityClust)
library(dbscan)
suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
suppressMessages(library(caret))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(reshape2))
source("/scratch/apm25309/single_cell/ATACseq/v3/step1_clustering/model_based/snATAC_rr.marker_accessibility_Utils.R")
source("/scratch/apm25309/single_cell/ATACseq/v3/step1_clustering/model_based/snATAC_rr.cluster_Utils.v2.R")
source("/scratch/apm25309/single_cell/ATACseq/v3/step1_clustering/model_based/snATAC_rr.celltype_call_Utils.R")


###################################################################################################
# process data ------------------------------------------------------------------------------------
###################################################################################################

# load data
source(config)
if(file.exists(paste0(prefix,".sparseA.rds")) & file.exists(paste0(prefix,".metaB.rds"))){
    a <- readRDS(paste0(prefix,".sparseA.rds"))
    b <- readRDS(paste0(prefix,".metaB.rds"))
}else{
    dat <- loadData(input, meta)
    a <- dat$a
    b <- dat$b

    # filter data
    a <- cleanData(a, b, min.p=min.p, min.t=min.t, min.c=min.c, max.t=0.0001)
    
    # set up vars
    constant      <- 1/(ncol(a))
    b             <- b[colnames(a),]
    b$nSites      <- Matrix::colSums(a)
    b$p_tss       <- b$tss/b$unique
    b$p_ptmt      <- b$PtMt/b$total
    b$p_nuclear   <- b$unique/b$total
    b$log10nSites <- log10(Matrix::colSums(a))
    b$logUNIQUE   <- as.numeric(scale(log(b$unique)))
    b$logTSS      <- as.numeric(scale(log(b$p_tss + constant)))
    b$logPTMT     <- as.numeric(scale(log(b$p_ptmt + constant)))
    b$logNUCLEAR  <- as.numeric(scale(log(b$p_nuclear + constant)))
    b$log10tss    <- as.numeric(scale(log10(b$tss+1)))
    b$log10ptmt   <- as.numeric(scale(log10(b$PtMt+1)))

    # save a/b as RDS
    saveRDS(a, file=paste0(prefix,".sparseA.rds"))
    saveRDS(b, file=paste0(prefix,".metaB.rds"))
}

# get residuals
r.dev <- runModel(a,b, variates, r.variates,
                  subpeaks = subpeaks,
                  chunks = 256,
                  type = type,
                  modtype = modtype,
                  nthreads = ncpus,
                  prefix = prefix,
                  doPlot = doPlot,
                  weights = weights,
                  link = link,
                  method = method,
                  var.th = var.th)

# std dev
if(center==T){
    dev <- stddev(r.dev, threads=ncpus)
    rm(r.dev)
    gc()
}else{
    dev <- r.dev
    rm(r.dev)
    gc()
}

# reduce dimensionality
out.pcs  <- reduceDims(dev, b, 
                       n.pcs = n.pcs,
                       batch = batch, 
                       theta = theta, 
                       lambda = lambda, 
                       doL2 = doL2, 
                       cap = cap,
                       prefix = prefix,
                       scaleVar = scaleVar,
                       stdLSI = stdLSI,
                       raw = a, 
                       center = center,
                       useTop = useTop,
                       useTop.th = useTop.th)
rm(dev)
gc()

# project with UMAP
out.umap <- projectUMAP(out.pcs, m.dist=0.01, k.near=k.nearC, metric="euclidean")

# call clusters
b <- callClusters(a, b, out.pcs, out.umap, 
                  min.reads = min.reads, 
                  m.clst = m.clst, 
                  threshold3 = threshold,
                  k.near = k.nearC,
                  res = res,
                  prefix = prefix,
                  clustOB = clustOB,
                  clustType = clustType,
                  cleanCluster = do.Clean,
                  e.thresh = e.thresh,
                  cl.method = 1)

# write stats to shell
reportResults(b)

# output
outData(out.pcs, b, prefix=prefix, dev=NULL)

# plot UMAP
b$LouvainClusters <- factor(b$LouvainClusters)
plotUMAP(b, prefix=prefix)
plotSTATS(b, prefix=prefix)

# recluster
if(reclustType==1){
    subclusters <- reclust(a, b, variates, r.variates,
                           modtype="regModel",
                           chunks=256,
                           type=type,
                           subpeaks=subpeaks,
                           ncpus=ncpus,
                           stdize=stdize,
                           link="logit",
                           var.th=var.th,
                           scaleVar=scaleVar,
                           doL2=doL2,
                           prefix=prefix,
                           n.pcs=s.pcs,
                           batch=batch,
                           theta=theta,
                           lambda=lambda,
                           min.reads=sub.mreads,
                           m.clst=50,
                           threshold3=subthresh,
                           k.nearC=subK,
                           res=subres,
                           clustOB="svd",
                           nthreads=ncpus,
                           cleanCluster = do.Clean)
}else if(reclustType==2){
    subclusters <- reclust2(a, b, out.pcs,
                            prefix=prefix,
                            n.pcs=s.pcs,
                            batch=sub.batch,
                            theta=sub.theta,
                            lambda=sub.lambda,
                            min.reads=sub.mreads,
                            m.clst=25,
                            threshold3=subthresh,
                            k.nearC=subK,
                            res=subres,
                            clustOB=sub.clustOB,
                            clustType=sub.clustType,
                            nthreads=ncpus,
                            cleanCluster=sub.do.Clean,
                            e.thresh=e.thresh,
                            doSTDize=sub.stdLSI,
                            doL2=sub.doL2)
}else{
    message(" - skipping reclustering")
}

# output
if(reclustType==1 | reclustType==2 | reclustType==3){
    rm(a)
    gc()
    write.table(subclusters, file=paste0(prefix,"_subclusters.Louvain.txt"),
                quote=F, row.names=T, col.names=T, sep="\t")
    plotSTATS(subclusters, prefix=paste0(prefix,".subclusters"))
}else{
    rm(a)
    gc()
    subclusters <- b
}

# ###################################################################################################
# ## call celltypes ---------------------------------------------------------------------------------
# ###################################################################################################
# sparse <- "/scratch/apm25309/single_cell/ATACseq/v3/step3_gene_activity_and_coaccessibility/automated_annotation/custom_annotation/all.sct.exGBcounts.sparse"
# markers <- "/scratch/apm25309/single_cell/ATACseq/v3/step3_gene_activity_and_coaccessibility/automated_annotation/custom_annotation/marker_list.txt"
# obj <- loadDataS(sparse, subclusters, markers)
# a <- obj$a
# b <- obj$b
# m <- obj$m
# gm <- read.table("/scratch/apm25309/single_cell/ATACseq/v3/bedfiles/genes/Zea_mays.AGPv4.36.Allgene.expanded.bed")
# 
# # get gene tests
# all <- runDESeq2(a, b, prefix=prefix, threads=ncpus, reps=5)
# 
# # coarse level celltypes
# celltype.results <- classifyTypes.v1(all, m, b, 
#                                      num.markers=10, 
#                                      clustID="top_tissue_cluster", 
#                                      prefix=prefix, 
#                                      orderVar="statistic",
#                                      threshold=5e-2,
#                                      q.thresh=1,
#                                      adjustM=T,
#                                      scaleStat=F)
# 
# # fine level celltypes
# pred.celltypes <- classifyTypes.v2(a, celltype.results$b, m, 
#                                    celltype.results$topAnn, 
#                                    all, 
#                                    gm,
#                                    threads=ncpus, 
#                                    balance=F, 
#                                    marker.qval=1e-2,
#                                    normalize=T,
#                                    useweights=T,
#                                    enrichMethod=4,
#                                    threshold=log2(5),
#                                    doHard=T,
#                                    preCalcEnrich=NULL,
#                                    prefix=prefix)
# 
# # write output
# write.table(pred.celltypes$meta, file="metadata.step6.txt", quote=F, row.names=T, col.names=T, sep="\t")


# ###################################################################################################
# # check markers -----------------------------------------------------------------------------------
# ###################################################################################################
# 
# # load marker/gene activity
# datt <- loadMData(subclusters, out.pcs, geneact, mark, clustID="tissue_cluster")
# b <- datt$b
# activity <- datt$activity
# h.pcs <- datt$h.pcs
# marker.info <- datt$marker.info
# 
# # normalize per cell activity by cluster average and size factors
# results <- normalize.activity(b, activity, output=prefix, logTransform=F, scaleP=F)
# activity <- results$norm.act
# row.o <- results$row.o
# 
# # impute gene accessibility scores
# impute.activity <- smooth.data(activity, 
#                                k=15, step=2, npcs=ncol(out.pcs), df=NULL,
#                                rds=h.pcs, cleanExp=F, output=prefix)
# 
# # collect marker accessibility
# plot.act.scores(b, acts=activity, 
#                 info=marker.info, 
#                 logT=T,
#                 outname=paste0(prefix,".normalized.known.Markers.pdf"))
# 
# plot.act.scores(b, acts=impute.activity, 
#                 info=marker.info, 
#                 logT=F,
#                 outname=paste0(prefix,".impute.known.Markers.pdf"))
