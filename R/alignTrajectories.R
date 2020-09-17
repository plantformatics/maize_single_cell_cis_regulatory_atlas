## align species trajectories ##

# arguments
args <- commandArgs(T)
if(length(args) != 5){stop("Rscript alignTrajectories.R <zm.trajectory> <at.trajectory> <zm.tests> <at.tests> <prefix>")}
traj.zm <- as.character(args[1])
traj.at <- as.character(args[2])
test.zm <- as.character(args[3])
test.at <- as.character(args[4])
prefix  <- as.character(args[5])

# load libraries
library(cellAlign)
library(scales)
library(reshape2)
library(mgcv)
library(pheatmap)
library(RANN)
library(Matrix)
library(RColorBrewer)
library(parallel)
library(MASS)
library(gplots)

# load functions
smooth.data      <- function(x, k=15, step=3, npcs=30, df=NULL, rds=NULL, verbose=F){
    
    # verbose
    if(verbose){message(" - imputing gene activity ...")}
    
    # input
    data.use <- x
    
    # verbose
    if(!is.null(rds)){
        
        if(!is.null(df)){
            if(verbose){message("   * using UMAP manifold for smoothing ...")}
            pcs <- df[,c("umap1","umap2")]
        }else{
            if(verbose){message("   * using prior PC space as manifold ...")}
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{
        
        # LSI
        if(verbose){message("   * PC manifold set to NULL, running LSI (TFIDF)...")}
        x[x>0] <- 1
        tf.idf <- tfidf(x)
        
        # get PCS
        if(verbose){message("   * PC manifold set to NULL, running LSI ...")}
        pc <- irlba(t(tf.idf), npcs)
        pcs <- pc$u 
        rownames(pcs) <- colnames(x)
        colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
        
        # do l2-norm
        pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
    }
    
    # get KNN
    if(verbose){message("   * finding knn graph ...")}
    knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
    j <- as.numeric(x = t(x = knn.graph))
    i <- ((1:length(x = j)) - 1) %/% k + 1
    edgeList = data.frame(i, j, 1)
    A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
    
    # Smooth graph
    if(verbose){message("   * smoothing graph ...")}
    A = A + t(A)
    A = A / Matrix::rowSums(A)
    step.size = step
    if(step.size > 1){
        for(i in 1:step.size){
            if(verbose){message("     ~ step ",i)}
            A = A %*% A
        }
    }
    
    # smooth data
    if(verbose){message("   * smoothing activity ...")}
    impute.activity <- t(A %*% t(data.use))
    colnames(impute.activity) <- colnames(x)
    rownames(impute.activity) <- rownames(x)
    
    # clean empty rows (if present) and round to two decimal places
    if(verbose){message("   * remove empty rows/columns and scale to per million ...")}
    impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
    impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
    impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
    impute.activity@x <- round(impute.activity@x, digits=2)
    
    # return sparse Matrix
    return(impute.activity)
}
scalePseudotime  <- function(x){
    x <- x[complete.cases(x$trajectory),]
    x$trajectory <- rescale(x$trajectory, c(0,1))
    x <- x[order(x$trajectory, decreasing=F),]
    return(x)
}
processOrthologs <- function(zm, at, zm.meta, at.meta, zm.tests, at.tests, orthologs, 
                             smooth=F, FDR=0.05, doTPM=F, prefix="trajectory"){
    
    # rescale pseudotime
    zm.meta <- scalePseudotime(zm.meta)
    at.meta <- scalePseudotime(at.meta)
    zm <- zm[,rownames(zm.meta)]
    at <- at[,rownames(at.meta)]
    
    # smooth data
    message(" - smoothing accessibility profiles ...")
    if(smooth){
        zm.svd <- read.table("all.SVD.txt")
        at.svd <- read.table("Arabidopsis.root.rawSVD.txt")
        zm.svd <- zm.svd[rownames(zm.meta),]
        at.svd <- at.svd[rownames(at.meta),]
        zm.svd <- t(apply(zm.svd, 1, function(x) x/sum(sqrt(x^2))))
        at.svd <- t(apply(at.svd, 1, function(x) x/sum(sqrt(x^2))))
        zm <- smooth.data(zm, k=15, step=3, npcs=ncol(zm.svd), rds=zm.svd)
        at <- smooth.data(at, k=15, step=3, npcs=ncol(at.svd), rds=at.svd)
    }
    
    # subset by pseudotime tests
    message(" - constrain to pseudotime-dependent regions ...")
    zm.tests <- subset(zm.tests, zm.tests$qval < FDR)
    at.tests <- subset(at.tests, at.tests$qval < FDR)
    zm <- zm[rownames(zm) %in% rownames(zm.tests),]
    at <- at[rownames(at) %in% rownames(at.tests),]
    message(" - Zea mays genes associated with pseudotime (n=",nrow(zm.tests),")")
    message(" - Arabidopsis thaliana genes associated with pseudotime (n=",nrow(at.tests),")")
    
    # subset to orthologs
    orthologs <- orthologs[which(as.character(orthologs$V1) %in% rownames(zm) & as.character(orthologs$V2) %in% rownames(at)),]
    zm.specific <- rownames(zm)[!rownames(zm) %in% as.character(orthologs$V2)]
    at.specific <- rownames(at)[!rownames(at) %in% as.character(orthologs$V1)]
    write.table(zm.specific, file=paste0("Zm.",prefix,".specific.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    write.table(at.specific, file=paste0("At.",prefix,".specific.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    zm <- zm[as.character(orthologs$V1),]
    at <- at[as.character(orthologs$V2),]
    maizeIDs <- rownames(zm)
    arabIDs <- rownames(at)
    rownames(zm) <- rownames(at)
    
    # remove empty rows/columns
    zm <- zm[,Matrix::colSums(zm)>0]
    at <- at[,Matrix::colSums(at)>0]
    acc.genes <- which(Matrix::rowSums(zm) > 0 & Matrix::rowSums(at) > 0)
    zm <- zm[acc.genes,]
    at <- at[acc.genes,]
    #zm <- log1p(zm)
    #at <- log1p(at)
    maizeIDs <- maizeIDs[acc.genes]
    arabIDs <- arabIDs[acc.genes]
    
    # scale to per 100K
    if(doTPM){
        zm <- zm %*% Diagonal(x=1e5/Matrix::colSums(zm))
        at <- at %*% Diagonal(x=1e5/Matrix::colSums(at))
    }
    
    # update meta
    zm.meta <- zm.meta[colnames(zm),]
    at.meta <- at.meta[colnames(at),]
    orthologs <- orthologs[as.character(orthologs$V2) %in% rownames(at),]
    message(" - Retained orthologs (At/Zm: n=",nrow(at),"|",nrow(zm),")")
    
    # plot heatmaps of At/Zm gene by cell
    pdf("Arabidopsis.orthologs.pdf", width=7, height=6)
    pheatmap(as.matrix(at), scale="row", show_rownames = F, show_colnames = F, cluster_cols = F)
    dev.off()
    
    pdf("Maize.orthologs.pdf", width=7, height=6)
    pheatmap(as.matrix(zm), scale="row", show_rownames = F, show_colnames = F, cluster_cols = F)
    dev.off()
    
    # return
    return(list(zm=zm, at=at, zm.meta=zm.meta, at.meta=at.meta, orthologs=orthologs, mIDs=maizeIDs, aIDs=arabIDs))
    
}
processOrtholog2 <- function(zm, at, zm.meta, at.meta, zm.tests, at.tests, orthologs, 
                             smooth=F, FDR=0.05, doTPM=F, prefix="trajectory"){
    
    # rescale pseudotime
    zm.meta <- scalePseudotime(zm.meta)
    at.meta <- scalePseudotime(at.meta)
    zm <- zm[,rownames(zm.meta)]
    at <- at[,rownames(at.meta)]
    
    # smooth data
    message(" - smoothing accessibility profiles ...")
    if(smooth){
        zm.svd <- read.table("all.SVD.txt")
        at.svd <- read.table("Arabidopsis.root.rawSVD.txt")
        zm.svd <- zm.svd[rownames(zm.meta),]
        at.svd <- at.svd[rownames(at.meta),]
        zm.svd <- t(apply(zm.svd, 1, function(x) x/sum(sqrt(x^2))))
        at.svd <- t(apply(at.svd, 1, function(x) x/sum(sqrt(x^2))))
        zm <- smooth.data(zm, k=15, step=3, npcs=ncol(zm.svd), rds=zm.svd)
        at <- smooth.data(at, k=15, step=3, npcs=ncol(at.svd), rds=at.svd)
    }
    
    # subset by pseudotime tests
    message(" - constrain to pseudotime-dependent regions ...")
    zm.tests <- subset(zm.tests, zm.tests$qval < FDR)
    at.tests <- subset(at.tests, at.tests$qval < FDR)
    zm <- zm[rownames(zm) %in% rownames(zm.tests),]
    at <- at[rownames(at) %in% rownames(at.tests),]
    message(" - Zea mays genes associated with pseudotime (n=",nrow(zm.tests),")")
    message(" - Arabidopsis thaliana genes associated with pseudotime (n=",nrow(at.tests),")")
    
    # subset to orthologs
    message(" - select orthologs ...")
    orthologs <- orthologs[which(as.character(orthologs$V2) %in% rownames(zm) & as.character(orthologs$V1) %in% rownames(at)),]
    zm.specific <- rownames(zm)[!rownames(zm) %in% as.character(orthologs$V2)]
    at.specific <- rownames(at)[!rownames(at) %in% as.character(orthologs$V1)]
    write.table(zm.specific, file=paste0("Zm.",prefix,".specific.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    write.table(at.specific, file=paste0("At.",prefix,".specific.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    zm <- zm[as.character(orthologs$V2),]
    at <- at[as.character(orthologs$V1),]
    rownames(zm) <- rownames(at)
    message(" - number of orthologs kept in At|Zm: ",nrow(at),"|",nrow(zm))
    
    # remove empty rows/columns
    zm <- zm[,Matrix::colSums(zm)>0]
    at <- at[,Matrix::colSums(at)>0]
    acc.genes <- which(Matrix::rowSums(zm) > 0 & Matrix::rowSums(at) > 0)
    zm <- zm[acc.genes,]
    at <- at[acc.genes,]
    
    # scale to per 100K
    if(doTPM){
        zm <- zm %*% Diagonal(x=1e5/Matrix::colSums(zm))
        at <- at %*% Diagonal(x=1e5/Matrix::colSums(at))
    }
    
    zm <- log1p(zm)
    at <- log1p(at)
    
    # update meta
    zm.meta <- zm.meta[colnames(zm),]
    at.meta <- at.meta[colnames(at),]
    orthologs <- orthologs[as.character(orthologs$V1) %in% rownames(at),]
    message(" - Retained orthologs (n=",nrow(at),")")
    
    # plot heatmaps of At/Zm gene by cell
    pdf("Arabidopsis.orthologs.pdf", width=7, height=6)
    pheatmap(as.matrix(at), scale="row", show_rownames = F, show_colnames = F, cluster_cols = F)
    dev.off()
    
    pdf("Maize.orthologs.pdf", width=7, height=6)
    pheatmap(as.matrix(zm), scale="row", show_rownames = F, show_colnames = F, cluster_cols = F)
    dev.off()
    
    # return
    return(list(zm=zm, at=at, zm.meta=zm.meta, at.meta=at.meta, orthologs=orthologs))
    
}
findRelatedGenes <- function(zm, at, zm.meta, at.meta, num.pts=250, threshold=0){
    
    # order rows
    zm.meta <- zm.meta[order(zm.meta$trajectory, decreasing=F),]
    at.meta <- at.meta[order(at.meta$trajectory, decreasing=F),]
    zm <- zm[,rownames(zm.meta)]
    at <- at[,rownames(at.meta)]
    
    # interpolate
    pseudo.y <- seq(from=0, to=1, length.out=num.pts)
    
    # iterate over each gene
    outs <- lapply(rownames(zm), function(x){
        
        # get raw pseudo-time sorted accessibility
        zm.df <- data.frame(zm.gene=as.numeric(zm[x,]), 
                            zm.dev=zm.meta$trajectory)
        at.df <- data.frame(at.gene=as.numeric(at[x,]),
                            at.dev=at.meta$trajectory)
    
        # model
        zm.mod <- gam(zm.gene~s(zm.dev, bs="cr"), data=zm.df)
        at.mod <- gam(at.gene~s(at.dev, bs="cr"), data=at.df)
        
        # predict
        zm.pred <- predict(zm.mod, data.frame(zm.dev=pseudo.y))
        at.pred <- predict(at.mod, data.frame(at.dev=pseudo.y))
        
        # correlation
        cor(zm.pred, at.pred, method="spearman")

    })
    names(outs) <- rownames(zm)
    outs <- do.call(c, outs)
    
    # filter genes
    keep <- names(outs[outs > threshold])
    zm <- zm[keep,]
    at <- at[keep,]
    
    # return
    return(list(zm=zm, at=at))
    
}
naiveAlign       <- function(zm.meta, at.meta, pts=250, prefix="test"){
    
    # remove NA
    zm.meta <- zm.meta[complete.cases(zm.meta),]
    at.meta <- at.meta[complete.cases(at.meta),]
    zm.meta$celltype <- as.character(zm.meta$celltype)
    at.meta$celltype <- as.character(at.meta$celltype)
    
    # order
    zm.meta <- zm.meta[order(zm.meta$trajectory, decreasing=F),]
    at.meta <- at.meta[order(at.meta$trajectory, decreasing=F),]
    
    # map celltypes to pseudotimes
    zm.cts <- unique(zm.meta$celltype)
    at.cts <- unique(at.meta$celltype)
    zm.uniform <- lapply(zm.cts, function(x){
        zm.sub <- subset(zm.meta, zm.meta$celltype==x)
        zm.t.range <- quantile(zm.sub$trajectory, c(0.025, 0.975))
        zm.brks <- seq(from=zm.t.range[1], to=zm.t.range[2], length.out=pts)
        names(zm.brks) <- paste0(x,"-", seq(1:pts))
        zm.brks
    })
    zm.uniform <- do.call(c, zm.uniform)
    at.uniform <- lapply(at.cts, function(x){
        at.sub <- subset(at.meta, at.meta$celltype==x)
        at.t.range <- quantile(at.sub$trajectory, c(0.025, 0.975))
        at.brks <- seq(from=at.t.range[1], to=at.t.range[2], length.out=pts)
        names(at.brks) <- paste0(x,"-", seq(1:pts))
        at.brks
    })
    at.uniform <- do.call(c, at.uniform)
    zm.uniform <- zm.uniform[order(zm.uniform, decreasing=F)]
    at.uniform <- at.uniform[order(at.uniform, decreasing=F)]
    df <- data.frame(zm=zm.uniform, at=at.uniform, row.names=names(zm.uniform))
    
    # plot 
    pdf(paste0(prefix, ".naiveTrajectory.pdf"), width=5, height=5)
    plot(df$zm, df$at, type="l", xlab="Zea mays", ylab="Arabidopsis thaliana")
    segments(0, 0, 1, 1, lty=2, col="red")
    types <- as.data.frame(do.call(rbind, strsplit(rownames(df), "-")))
    types <- as.character(types$V1)
    cols <- brewer.pal(length(unique(types)), "Paired")
    points(seq(from=0,to=1,length.out=nrow(df)), rep(0, times=nrow(df)), col=cols[factor(types)], pch=15, cex=0.25)
    points(rep(0, times=nrow(df)), seq(from=0, to=1, length.out=nrow(df)), col=cols[factor(types)], pch=15, cex=0.25)
    dev.off()
    
    # return
    return(df)
    
}
runGBLalign      <- function(zm.meta, at.meta, zm, at, prefix="test", numPts=250){
    
    message(" - weighting cells and interpolating ...")
    zm.traj <- zm.meta$trajectory
    at.traj <- at.meta$trajectory
    names(zm.traj) <- rownames(zm.meta)
    names(at.traj) <- rownames(at.meta)
    interGlobalZm = cellAlign::interWeights(expDataBatch = zm, trajCond = zm.traj,
                                            winSz = 0.1, numPts = numPts)
    interGlobalAt = cellAlign::interWeights(expDataBatch = at, trajCond = at.traj,
                                            winSz = 0.1, numPts = numPts)
    
    # intersect markers
    message(" - aligning trajectories ...")
    # sharedMarkers = intersect(rownames(zm), rownames(at))
    # whichgene = sharedMarkers[1]
    # selectedZm <- interGlobalZm$interpolatedVals[whichgene,]
    # selectedAt <- interGlobalAt$interpolatedVals[whichgene,]
    # 
    # dfZmi = data.frame(traj = interGlobalZm$traj, value=(selectedZm), error=interGlobalZm$error[whichgene,])
    # dfZm = data.frame(traj = zm.traj, t(zm[whichgene,]))
    # dfAti = data.frame(traj = interGlobalAt$traj, value=(selectedAt), error=interGlobalAt$error[whichgene,])
    # dfAt = data.frame(traj = at.traj, t(at[whichgene,]))
    # dfZmM = melt(dfZm, id.vars = 'traj')
    # dfAtM = melt(dfAt, id.vars = 'traj')
    
    # scale gene accessibility
    interScaledGlobalZm = cellAlign::scaleInterpolate(interGlobalZm)
    interScaledGlobalAt = cellAlign::scaleInterpolate(interGlobalAt)
    print(head(interScaledGlobalZm$scaledData))
    
    # trajectory alignment
    alignment = globalAlign(interScaledGlobalZm$scaledData, interScaledGlobalAt$scaledData,
                            scores = list(query = interScaledGlobalAt$traj, 
                                          ref = interScaledGlobalZm$traj),
                            sigCalc = T, numPerm = 100)
    
    pdf(paste0(prefix,".alnTrajectory.pdf"), width=8, height=7)
    plotAlign(alignment)
    dev.off()
    
    # print overall distance
    print(alignment$normalizedDistance)
    
    # save results
    saveRDS(alignment, file=paste0(prefix,".cellAlign.rds"))
    
    
    
}
runLCLalign      <- function(zm.meta, at.meta, zm, at, prefix="test", Thresh=0.2, numPts=250){

    # defaults
    zm.traj <- zm.meta$trajectory
    at.traj <- at.meta$trajectory
    names(zm.traj) <- rownames(zm.meta)
    names(at.traj) <- rownames(at.meta)
    interLocalzm = interWeights(expDataBatch = zm, trajCond = zm.traj, winSz = 0.1, numPts = numPts)
    interLocalat = interWeights(expDataBatch = at, trajCond = at.traj, winSz = 0.1, numPts = numPts)
    interScaledLocalzm = cellAlign::scaleInterpolate(interLocalzm)
    interScaledLocalat = cellAlign::scaleInterpolate(interLocalat)
    
    # dist matrix
    A=calcDistMat(interScaledLocalat$scaledData,interScaledLocalzm$scaledData, dist.method = 'Euclidean')
    A[A > 10*Thresh] <- max(A)
    alignment = localAlign(interScaledLocalat$scaledData,interScaledLocalzm$scaledData, threshPercent = Thresh)
    
    # plot local alignment
    pdf(paste0(prefix,".localAlignment.pdf"), width=8, height=8)
    plotAlign(alignment)
    dev.off()
    
    # return cost matrix
    costMat = t(apply(A,1,function(x){return(as.numeric(x))}))
    linearInd = cellAlign::sub2ind(nrow(A), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
    costMat[linearInd] = NA
    costMat = data.frame(costMat, row.names=1:numPts)
    colnames(costMat) = 1:numPts
    print(head(costMat))
    
    # heatmap
    pdf(paste0(prefix,".costMatrix.pdf"), width=8, height=7)
    pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
             main = 'gated search region',
             show_rownames = F, show_colnames = F)
    dev.off()
    
    # get vals
    BLPS=colMeans(interScaledLocalzm$scaledData)
    BPAM=colMeans(interScaledLocalat$scaledData)
    p=unique(alignment$align[[1]]$index1)
    q=unique(alignment$align[[1]]$index2)
    
    # plot 
    pdf(paste0(prefix,".localTrajectoryAlignment.pdf"), width=6, height=6)
    plot(1:numPts,BPAM,xlab = "pseudotime",ylab = "mean interpolated expression", main = "unaligned mean expression",ylim = c(0,1.1))
    points(p,BPAM[p],col="red")
    points(1:numPts,BLPS,col="grey60")
    points(q,BLPS[q],col="red")
    text(90,1,"Zea mays")
    text(150,1,"Arabidopsis thaliana")
    text(125,.3,"red points are conserved")
    dev.off()
    
}
runGBLalignGene  <- function(zm.meta, at.meta, zm, at, prefix="test", numPts=250, threads=1){
    
    message(" - weighting cells and interpolating ...")
    zm.traj <- zm.meta$trajectory
    at.traj <- at.meta$trajectory
    names(zm.traj) <- rownames(zm.meta)
    names(at.traj) <- rownames(at.meta)
    #rownames(zm) <- paste0(rownames(zm),"_",seq(1:nrow(zm)))
    #rownames(at) <- paste0(rownames(at),"_",seq(1:nrow(at)))
    interGlobalZm = cellAlign::interWeights(expDataBatch = zm, trajCond = zm.traj,
                                            winSz = 0.1, numPts = numPts)
    interGlobalAt = cellAlign::interWeights(expDataBatch = at, trajCond = at.traj,
                                            winSz = 0.1, numPts = numPts)
    
    # intersect markers
    message(" - aligning trajectories ...")
    
    # scale gene accessibility
    interScaledGlobalZm = cellAlign::scaleInterpolate(interGlobalZm)
    interScaledGlobalAt = cellAlign::scaleInterpolate(interGlobalAt)
    
    # cluster by pseudotime shift
    cluster <- pseudotimeClust(x=interScaledGlobalZm$scaledData, y=interScaledGlobalAt$scaledData, k = 2)
    
    # trajectory alignment
    per.gene <- lapply(seq(1:length(rownames(interScaledGlobalZm$scaledData))), function(x){
        alignment = globalAlign(interScaledGlobalZm$scaledData[x,], interScaledGlobalAt$scaledData[x,],
                                scores = list(query = interScaledGlobalAt$traj,
                                              ref = interScaledGlobalZm$traj),
                                sigCalc = T, numPerm = 100)
        if(x==1){
            pdf("test.singlegeneALN.pdf", width=6, height=5.5)
            plotAlign(alignment)
            dev.off()
        }
        return(alignment)
    })
    names(per.gene) <- rownames(interScaledGlobalZm$scaledData)
    return(list(per.gene=per.gene, cluster=cluster, zm.exp=interScaledGlobalZm$scaledData, at.exp=interScaledGlobalAt$scaledData))
    
    
}
oldcluster       <- function(out){
    
    pers <- out$per.gene
    clusts <- out$cluster
    pdf("all.per_gene.pdf", width=5, height=5)
    
    con <- list()
    cl1 <- list()
    cl2 <- list()
    cl3 <- list()
    #ran <- list()
    
    for(i in 1:length(pers)){
        if(pers[[i]]$normalizedDistance > 0.1 & clusts$cluster[i] == 1){
            cl1[[i]] <- as.matrix(sparseMatrix(i=as.numeric(pers[[i]]$align[[1]]$index1),
                                               j=as.numeric(pers[[i]]$align[[1]]$index2),
                                               x=rep(1,length(pers[[i]]$align[[1]]$index2)),
                                               dimnames=list(seq(1:200),seq(1:200))))
        }else if(pers[[i]]$normalizedDistance > 0.1 & clusts$cluster[i] == 3){
            cl3[[i]] <- as.matrix(sparseMatrix(i=as.numeric(pers[[i]]$align[[1]]$index1),
                                               j=as.numeric(pers[[i]]$align[[1]]$index2),
                                               x=rep(1,length(pers[[i]]$align[[1]]$index2)),
                                               dimnames=list(seq(1:200),seq(1:200))))
        }else if(pers[[i]]$normalizedDistance <= 0.1){
            con[[i]] <- as.matrix(sparseMatrix(i=as.numeric(pers[[i]]$align[[1]]$index1),
                                               j=as.numeric(pers[[i]]$align[[1]]$index2),
                                               x=rep(1,length(pers[[i]]$align[[1]]$index2)),
                                               dimnames=list(seq(1:200),seq(1:200))))
        }else if(pers[[i]]$normalizedDistance > 0.1 & clusts$cluster[i] == 2){
            cl2[[i]] <- as.matrix(sparseMatrix(i=as.numeric(pers[[i]]$align[[1]]$index1),
                                               j=as.numeric(pers[[i]]$align[[1]]$index2),
                                               x=rep(1,length(pers[[i]]$align[[1]]$index2)),
                                               dimnames=list(seq(1:200),seq(1:200))))
        }
    }
    cl1 <- cl1[-which(sapply(cl1, is.null))]
    cl2 <- cl2[-which(sapply(cl2, is.null))]
    cl3 <- cl3[-which(sapply(cl3, is.null))]
    con <- con[-which(sapply(con, is.null))]
    
    
    names(con) <- paste0("gene_",seq(1:length(con)))
    names(cl1) <- paste0("gene_",seq(1:length(cl1)))
    names(cl2) <- paste0("gene_",seq(1:length(cl2)))
    names(cl3) <- paste0("gene_",seq(1:length(cl3)))
    n.con <- length(con)
    n.cl1 <- length(cl1)
    n.cl2 <- length(cl2)
    n.cl3 <- length(cl3)
    print(head(cl1[[1]]))
    print(head(cl2[[1]]))
    print(head(cl3[[1]]))
    print(head(con[[1]]))
    
    # average conserved
    con <- Reduce("+", con)/length(con)
    cl1 <- Reduce("+", cl1)/length(cl1)
    cl2 <- Reduce("+", cl2)/length(cl2)
    cl3 <- Reduce("+", cl3)/length(cl2)
    
    # dens
    conxx <- colMeans(con)
    conyy <- rowMeans(con)
    cl1xx <- colMeans(cl1)
    cl1yy <- rowMeans(cl1)
    cl2xx <- colMeans(cl2)
    cl2yy <- rowMeans(cl2)
    cl3xx <- colMeans(cl3)
    cl3yy <- rowMeans(cl3)
    
    plot(conxx, conyy, col="red", type="l")
    lines(cl1xx, cl1yy, col="dodgerblue4", type="l")
    lines(cl2xx, cl2yy, col="dodgerblue1", type="l")
    lines(cl3xx, cl3yy, col="grey75", type="l")
    dev.off()
    
    
    
    
}

# load data
message(" - load data ...")
orthologs <- read.table("Zm_At_orthologs.txt")
zm <- readRDS("Zea_mays.normalizedActivity.rds")
at <- readRDS("Arabidopsis_thaliana.normalizedActivity.rds")
zm.meta <- read.table(traj.zm, comment.char="")
at.meta <- read.table(traj.at, comment.char="")
zm.tests <- read.table(test.zm)
at.tests <- read.table(test.at)

# process orthologs
message(" - processing orthologs ...")
processed <- processOrthologs(zm, at, zm.meta, at.meta, zm.tests, at.tests, orthologs, 
                              smooth=T, FDR=1, doTPM=F, prefix=prefix)
zm <- processed$zm
zm.meta <- processed$zm.meta
at <- processed$at
at.meta <- processed$at.meta
orthologs <- processed$orthologs

# naive alignment
aligned <- naiveAlign(zm.meta, at.meta, pts=300, prefix=prefix)

# global alignment
out <- runGBLalignGene(zm.meta, at.meta, zm, at, prefix=prefix, numPts=200)
saveRDS(out, file=paste0(prefix,".geneGlobalAlignTrajectory.rds"))

# get output
pers <- out$per.gene
clusts <- out$cluster

# split by cluster
conx <- c()
cl1x <- c()
cl2x <- c()
cony <- c()
cl1y <- c()
cl2y <- c()
con_gene <- list()
cl1_gene <- list()
cl2_gene <- list()
n.dists <- c()

# populate lists/vectors
for(i in 1:length(pers)){
    if(pers[[i]]$normalizedDistance > 0.15 & clusts$cluster[i] == 1){
        cl1x <- c(cl1x, as.numeric(pers[[i]]$align[[1]]$index1))
        cl1y <- c(cl1y, as.numeric(pers[[i]]$align[[1]]$index2))
        cl1_gene[names(pers[i])] <- pers[[i]]$normalizedDistance
    }else if(pers[[i]]$normalizedDistance > 0.15 & clusts$cluster[i] == 2){
        cl2x <- c(cl2x, as.numeric(pers[[i]]$align[[1]]$index1))
        cl2y <- c(cl2y, as.numeric(pers[[i]]$align[[1]]$index2))
        cl2_gene[names(pers[i])] <- pers[[i]]$normalizedDistance
    }else if(pers[[i]]$normalizedDistance <= 0.15){
        conx <- c(conx, as.numeric(pers[[i]]$align[[1]]$index1))
        cony <- c(cony, as.numeric(pers[[i]]$align[[1]]$index2))
        con_gene[names(pers[i])] <- pers[[i]]$normalizedDistance
    }
}

# make dfs
cl1_gene <- do.call(c, cl1_gene)
cl2_gene <- do.call(c, cl2_gene)
con_gene <- do.call(c, con_gene)
cl1_top <- names(cl1_gene[order(cl1_gene, decreasing=T)])[1]
cl2_top <- names(cl2_gene[order(cl2_gene, decreasing=T)])[1]
con_top <- names(con_gene[order(con_gene, decreasing=F)])[1]

# save clustering
all <- c(cl1_gene, cl2_gene, con_gene)
id <- rep(c("cl1", "cl2", "con"), c(length(cl1_gene),length(cl2_gene),length(con_gene)))
df <- data.frame(n.dists=all, type=id, genes=names(all))
write.table(df, file=paste0(prefix,".ptShiftClusters.txt"), quote=F, row.names=F, col.names=T, sep="\t")

# print number of genes in each group
message(" cluster 1 = ", length(cl1_gene))
message(" cluster 2 = ", length(cl2_gene))
message(" conserved = ", length(con_gene))

# plot gene from each
cols1 <- brewer.pal(3, "Paired")
cols2 <- brewer.pal(4, "Paired")
pdf("cluster1_geneExample.pdf", width=6, height=6)
layout(matrix(c(1:2), nrow=2))
par(mar=c(3,5,1,1))
plot(zm.meta$trajectory, as.numeric(zm[cl1_top,]), col=cols2[as.factor(zm.meta$celltypeID)], pch=16, xlab="", ylab="Zea mays", main=cl1_top, xaxt="n")
plot(at.meta$trajectory, as.numeric(at[cl1_top,]), col=cols1[as.factor(at.meta$celltypeID)], pch=16, xlab="Pseudotime", ylab="Arabidopsis thaliana")
dev.off()

pdf("cluster2_geneExample.pdf", width=6, height=6)
layout(matrix(c(1:2), nrow=2))
par(mar=c(3,5,1,1))
plot(zm.meta$trajectory, as.numeric(zm[cl2_top,]), col=cols2[as.factor(zm.meta$celltypeID)], pch=16, xlab="", ylab="Zea mays", main=cl2_top, xaxt="n")
plot(at.meta$trajectory, as.numeric(at[cl2_top,]), col=cols1[as.factor(at.meta$celltypeID)], pch=16, xlab="Pseudotime", ylab="Arabidopsis thaliana")
dev.off()

pdf("conserved_geneExample.pdf", width=6, height=6)
layout(matrix(c(1:2), nrow=2))
par(mar=c(3,5,1,1))
plot(zm.meta$trajectory, as.numeric(zm[con_top,]), col=cols2[as.factor(zm.meta$celltypeID)], pch=16, xlab="", ylab="Zea mays", main=con_top, xaxt="n")
plot(at.meta$trajectory, as.numeric(at[con_top,]), col=cols1[as.factor(at.meta$celltypeID)], pch=16, xlab="Pseudotime", ylab="Arabidopsis thaliana")
dev.off()

# density 
cl1den <- kde2d(cl1x, cl1y, n=300)
cl2den <- kde2d(cl2x, cl2y, n=300)
conden <- kde2d(conx, cony, n=300)

# plots
cols <- colorRampPalette(c("white","grey85",brewer.pal(9, "YlGnBu")))(100)

# cl1
pdf("cluster1_ptshift.pdf", width=6, height=5)
pheatmap(cl1den$z, col=cols, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, useRaster=T)
dev.off()

# cl2
pdf("cluster2_ptshift.pdf", width=6, height=5)
pheatmap(cl2den$z, col=cols, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, useRaster=T)
dev.off()

# con
pdf("conserved_ptshift.pdf", width=6, height=5)
pheatmap(conden$z, col=cols, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, useRaster=T)
dev.off()

# plot pt shift clusters
ptshift <- lapply(seq_along(pers), function(x){
    pers[[x]]$ptShift
})
ptshift <- do.call(rbind, ptshift)
rownames(ptshift) <- names(pers)
ptshift <- ptshift[c(names(cl1_gene),names(cl2_gene),names(con_gene)),]
ids <- rep(c("A","B","C"), c(length(cl1_gene), length(cl2_gene), length(con_gene)))
cols <- c("darkorchid", "dodgerblue", "grey75")
pdf("PseudotimeShiftHeatmap.pdf", width=6, height=6)
heatmap.2(as.matrix(ptshift), trace='none', scale='none', dendrogram='none', Colv=F, Rowv=F,
          labRow = F, labCol = F, RowSideColors = cols[factor(ids)],
          col=colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100), useRaster=T)
dev.off()
