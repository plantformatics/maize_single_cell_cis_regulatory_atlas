## pseudotime ##

# load arguments
arg <- commandArgs(T)
print(arg)
if(length(arg)!= 9){stop("pseudotime.R <sparse> <motif.deviations> <gene> <TFgene> <meta> <svd> <prefix> <config> <threads>")}
sparse <- as.character(arg[1])
motif <- as.character(arg[2])
gene <- as.character(arg[3])
tfgene <- as.character(arg[4])
meta <- as.character(arg[5])
svd <- as.character(arg[6])
prefix <- as.character(arg[7])
config <- as.character(arg[8])
threads <- as.numeric(arg[9])

# defaults
column <- "top_tissue_cluster"
LC <- 1
featureMin <- 0

# load config file
source(config)
umap1 <- "umapPT_1"
umap2 <- "umapPT_2"

# load libraries
library(Matrix)
library(nabor)
library(dplyr)
library(viridis)
library(data.table)
library(scales)
library(mgcv)
library(gplots)
library(RColorBrewer)
library(parallel)
library(speedglm)
library(gtools)
library(uwot)
library(splines)
library(RANN)

# load functions
smooth.data<- function(x, k=15, step=3, npcs=30, df=NULL, rds=NULL, verbose=F){
    
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
loadData   <- function(sm, mt, gn, tf, ma, rd, doL2=0, traj, column="top_tissue_cluster"){
    message(" - loading data ...")
    sm <- readRDS(sm)
    mt <- t(read.table(mt))
    gn <- read.table(gn)
    gn <- sparseMatrix(i=as.numeric(gn$V1),
                       j=as.numeric(gn$V2),
                       x=as.numeric(gn$V3),
                       dimnames=list(levels(gn$V1),levels(gn$V2)))
    tf <- read.table(tf)
    tf <- sparseMatrix(i=as.numeric(tf$V1),
                       j=as.numeric(tf$V2),
                       x=as.numeric(tf$V3),
                       dimnames=list(levels(tf$V1),levels(tf$V2)))
    ma <- read.table(ma, comment.char="")
    rd <- read.table(rd)
    ids <- intersect(rownames(ma), colnames(sm))
    sm <- sm[,ids]
    ma <- ma[ids,]
    sm <- sm[Matrix::rowSums(sm)>0,]
    sm <- sm[,Matrix::colSums(sm)>0]
    ma <- ma[colnames(sm),]
    rd <- rd[colnames(sm),]
    
    # subset cells by cluster
    ma <- ma[as.character(ma[,column]) %in% as.character(traj),]
    tf <- tf[,rownames(ma)]
    gn <- gn[,rownames(ma)]
    gn <- gn[Matrix::rowSums(gn)>0,]
    sm <- sm[,rownames(ma)]
    rd <- rd[rownames(ma),]
    #gn <- smooth.data(gn, npcs=ncol(rd), rds=rd)
    #tf <- smooth.data(tf, npcs=ncol(rd), rds=rd)
    
    # L2
    if(doL2==1){
        rd <- t(apply(rd, 1, function(x) x/sqrt(sum(x^2))))
    }else if(doL2==2){
        rd <- apply(rd, 2, function(x) x/sqrt(sum(x^2)))
    }
    
    # re-run UMAP
    new.umap <- umap(rd, min_dist=0.1, n_neighbors=15)
    colnames(new.umap) <- c("umapPT_1", "umapPT_2")
    rownames(new.umap) <- rownames(rd)
    ma <- cbind(ma, new.umap)
    
    # return list
    return(list(a=sm, m=mt, gn=gn, tf=tf, b=ma, d=rd))
}
calcPseudo <- function(b, d, column="top_tissue_cluster", traj=NULL, cell.dist1=0.9, cell.dist2=0.99,
                       dof=250, spar=1, new.column="trajectory"){
    
    
    #### the following code was repurposed from ArchR written by Jeff Granja & Ryan Corces             ####
    #### for the original source, see https://github.com/GreenleafLab/ArchR/blob/master/R/Trajectory.R ####
    
    # initiation checks
    if(is.null(traj)){
        stop(" - Argument (vector) is missing to object: traj ...")
    }
    if(!is.character(traj)){
        stop(" - traj argument must be a character vector ...")
    }
    if(!is.character(b[,c(column)])){
        stop(" - supplied columnID, ",column,", is not a character vector ...")
    }
    
    # hidden functions - taken from archR 
    .getQuantiles <- function(v = NULL, len = length(v)){
        if(length(v) < len){
            v2 <- rep(0, len)
            v2[seq_along(v)] <- v
        }else{
            v2 <- v
        }
        p <- trunc(rank(v2))/length(v2)
        if(length(v) < len){
            p <- p[seq_along(v)]
        }
        return(p)
    }
    
    # get average coordinates for specified clusters
    t.cells <- lapply(traj, function(x){
        b.sub <- b[b[,column]==as.character(x),]
        d.sub <- d[rownames(b.sub),]
        aves <- colMeans(d.sub)
        
        # filter to top 5% of cells
        per.traj <- sqrt(colSums((t(d.sub) - aves)^2))
        idx <- which(per.traj <= quantile(per.traj, cell.dist1))
        ids.keep <- rownames(d.sub)[idx]
        return(ids.keep)
    })
    names(t.cells) <- traj
    
    # get initial trajectory
    raw.pseudo <- unlist(lapply(seq_along(traj), function(x){
        
        # current cluster
        c.i <- d[t.cells[[x]],]
        
        # compute trajectory        
        if(x != length(traj)){
            c.j <- colMeans(d[t.cells[[(x+1)]],])
            dif <- sqrt(colSums((t(c.i) - c.j)^2))
            pt <- (1 - .getQuantiles(dif)) + x
        }else{
            c.j <- colMeans(d[t.cells[[(x-1)]],])
            dif <- sqrt(colSums((t(c.i) - c.j)^2))
            pt <- .getQuantiles(dif) + x
        }
        
        return(pt)
        
    }))
    
    # fit spline
    d.filt <- d[names(raw.pseudo),]
    d.spline <- lapply(seq_len(ncol(d.filt)), function(x){
        
        # fit split
        stats::smooth.spline(x = raw.pseudo, 
                             y = d.filt[,x], 
                             df = dof, 
                             spar = spar)[[2]]
        
        
    }) %>% Reduce("cbind",.) %>% data.frame()
    
    # KNN fit versus actual fit
    knnObj <- nabor::knn(data = d.spline, query = d.filt, k = 3)
    
    # place along trajectory
    knnIdx <- knnObj[[1]]
    knnDist <- knnObj[[2]]
    knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
    knnDistQ <- .getQuantiles(knnDist[,1])
    
    #Filter Outlier Cells to Trajectory for High Resolution
    idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], cell.dist2))
    dfTrajectory <- data.frame(row.names = rownames(d.filt),
                               distance = knnDist[, 1],
                               distanceIdx = knnIdx[, 1] + knnDiff * knnDist)[idxKeep, , drop = FALSE]
    dfTrajectory$trajectory <- 100 * .getQuantiles(dfTrajectory[,2])
    dfTrajectory <- dfTrajectory[rownames(b),]
    b[,c(new.column)] <- dfTrajectory$trajectory
    
    # return
    return(b)
}
plotPT     <- function(meta, subsetCluster=NULL, t.id="trajectory", umap1="umapsub_1", umap2="umapsub_2",
                       smoother=5, addArrow=T, cex=0.5, xlab="umap1",ylab="umap2", bty='o'){
    
    # functions
    .centerRollMean <- function(v = NULL, k = NULL){
        o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
        if(k%%2==0){
            o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
        }else if(k%%2==1){
            o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
        }else{
            stop("Error!")
        }
        o2
    }
    
    # plot test
    if(!is.null(subsetCluster)){
        meta <- meta[meta$LouvainClusters %in% subsetCluster,]
    }
    
    # define pseudotime coloring
    cols <- viridis(100)
    
    # separate by NA
    test.1 <- meta[is.na(meta[,c(t.id)]),]
    test.2 <- meta[!is.na(meta[,c(t.id)]),]
    test.3 <- test.2
    
    # layout
    layout(matrix(c(1:2), nrow=1))
    
    # plot grey first
    plot(test.1[,c(umap1)], test.1[,c(umap2)], xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
         col=NA, pch=16, cex=cex, xlab=xlab, ylab=ylab, bty=bty)
    
    # plot trajectory values
    points(test.2[,c(umap1)], test.2[,c(umap2)], col=cols[cut(test.2[,c(t.id)], breaks=101)], pch=16, cex=cex)
    
    # add arrow
    if(addArrow){
        test.2 <- test.2[c(umap1, umap2, t.id)]
        dfArrow <-  split(test.2, floor(test.2[,c(t.id)] / 1.01)) %>% 
            lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame(row.names = NULL)
        dfArrow[,c(umap1)] <- .centerRollMean(dfArrow[,c(umap1)], smoother)
        dfArrow[,c(umap2)] <- .centerRollMean(dfArrow[,c(umap2)], smoother)
        dfArrow.smooth <- smooth.spline(dfArrow[,umap1],dfArrow[,umap2], spar=1)
        
        # plot
        lines(dfArrow[,umap1], dfArrow[,umap2], lwd=3, col=alpha("black", 0.8))
        dfArrow <- dfArrow[!duplicated(dfArrow[,c(umap1,umap2)]),]
        len.df <- nrow(dfArrow)
        a.x1 <- dfArrow[,umap1][len.df]
        a.y1 <- dfArrow[,umap2][len.df]
        dir1 <- a.x1 - mean(dfArrow[,umap1][(len.df-1):(len.df-1)])
        dir2 <- a.y1 - mean(dfArrow[,umap2][(len.df-1):(len.df-1)])
        a.x2 <- a.x1+dir1
        a.y2 <- a.y1+dir2
        arrows(a.x1, a.y1, a.x2, a.y2, lwd=3, length=0.1, col=alpha("black",0.8))
    }
    
    # plot celltypes
    test.2 <- test.3
    print(head(test.2))
    plot(test.2[,c(umap1)], test.2[,c(umap2)], col=as.character(test.2$sub_cluster_color), pch=16,
         xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
         cex=cex, xlab=xlab, ylab=ylab, bty=bty)
    legend("topright", legend=sort(unique(as.character(test.2$celltypeID))), 
           col=as.character(test.2$subcluster_color)[factor(sort(unique(as.character(test.2$celltypeID))),levels=sort(unique(as.character(test.2$celltypeID))))])
    
    # add arrow
    if(addArrow){
        test.2 <- test.2[c(umap1, umap2, t.id)]
        dfArrow <-  split(test.2, floor(test.2[,c(t.id)] / 1.01)) %>% 
            lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame(row.names = NULL)
        dfArrow[,c(umap1)] <- .centerRollMean(dfArrow[,c(umap1)], smoother)
        dfArrow[,c(umap2)] <- .centerRollMean(dfArrow[,c(umap2)], smoother)
        dfArrow.smooth <- smooth.spline(dfArrow[,umap1],dfArrow[,umap2], spar=1)
        
        # plot
        lines(dfArrow[,umap1], dfArrow[,umap2], lwd=3, col=alpha("black", 0.8))
        dfArrow <- dfArrow[!duplicated(dfArrow[,c(umap1,umap2)]),]
        len.df <- nrow(dfArrow)
        a.x1 <- dfArrow[,umap1][len.df]
        a.y1 <- dfArrow[,umap2][len.df]
        dir1 <- a.x1 - mean(dfArrow[,umap1][(len.df-1):(len.df-1)])
        dir2 <- a.y1 - mean(dfArrow[,umap2][(len.df-1):(len.df-1)])
        a.x2 <- a.x1+dir1
        a.y2 <- a.y1+dir2
        arrows(a.x1, a.y1, a.x2, a.y2, lwd=3, length=0.1, col=alpha("black",0.8))
    }
}
sigPseudo  <- function(obj, meta, n.pseudo.cells=NULL, num.bins=NULL, threads=1, type="ACRs", test="binomial"){
    
    # hidden functions
    .estimate_sf_sparse <- function(counts,
                                    round_exprs=TRUE,
                                    method="mean-geometric-mean-total"){
        if (round_exprs)
            counts <- round(counts)
        
        if(method == 'mean-geometric-mean-total') {
            cell_total <- Matrix::colSums(counts)
            sfs <- cell_total / exp(mean(log(cell_total)))
        }else if(method == 'mean-geometric-mean-log-total') {
            cell_total <- Matrix::colSums(counts)
            sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
        }
        
        sfs[is.na(sfs)] <- 1
        sfs
    }
    
    # use
    if(type=="ACRs"){
        use <- obj$a
    }else if(type=="TFs"){
        use <- obj$tf
    }else if(type=="MTs"){
        use <- obj$m - min(obj$m, na.rm=T)
        print(head(use[,1:5]))
    }else if(type=="Genes"){
        use <- obj$gn
    }
    
    # order cells
    meta <- meta[!is.na(meta$trajectory),]
    meta <- meta[order(meta$trajectory, decreasing=F),]
    if(is.null(n.pseudo.cells) & !is.null(num.bins)){
        num.per <- ceiling(nrow(meta)/num.bins)
        meta$pseudoCells <- ceiling(1:nrow(meta)/num.per)
    }else if(!is.null(n.pseudo.cells)){
        meta$pseudoCells <- ceiling(1:nrow(meta)/n.pseudo.cells)
        if(max(meta$pseudoCells) < 10){
            num.per <- ceiling(nrow(meta)/10)
            meta$pseudoCells <- ceiling(1:nrow(meta)/num.per)
        }
    }else{
        meta$pseudoCells <- ceiling(1:nrow(meta)/200)
    }
    pseudocell <- split(meta, factor(meta$pseudoCells))
    
    # collapse into pseudocells
    mat <- matrix(nrow=nrow(use), ncol=length(unique(meta$pseudoCells)), 
                  dimnames=list(rownames(use), paste0("pseudo_",names(pseudocell))))
    p.cells <- lapply(names(pseudocell), function(x){
        df <- pseudocell[[x]]
        ids <- rownames(df)
        a.sub <- use[,colnames(use) %in% ids]
        if(type == "ACRs"){
            a.sum <- Matrix::rowSums(a.sub)
        }else if(type=="TFs"){
            a.sum <- Matrix::rowMeans(a.sub)
        }else if(type == "MTs"){
            a.sum <- Matrix::rowMeans(a.sub)
        }else if(type== "Genes"){
            a.sum <- Matrix::rowMeans(a.sub)
        }
        p.id <- paste0("pseudo_",x)
        mat[,p.id] <<- a.sum
        mean.site <- sum(a.sum)
        mean.traj <- mean(df$trajectory)
        outs <- c(mean.site, mean.traj, ncol(a.sub))
        names(outs) <- c("nSites", "trajectory", "num_cells")
        return(outs)
    })
    p.cells <- as.data.frame(do.call(rbind, p.cells))
    rownames(p.cells) <- paste0("pseudo_", names(pseudocell))
    mat <- Matrix(mat, sparse=T)
    mat <- mat[Matrix::rowSums(mat)>0,]
    p.cells$size_factors <- .estimate_sf_sparse(mat)

    # run tests for each site
    outs <- mclapply(seq(1:nrow(mat)), function(x){
        if((x %% 1000)==0){message(" - iterated over ",x, " sites ...")}
        df.sub <- mat[x,rownames(p.cells)]
        df.sub <- cbind(p.cells, df.sub)
        colnames(df.sub) <- c("nSites","trajectory", "num_cells","size_factors","accessibility")
        if(test=="binomial"){
            df.dat <- as.matrix(data.frame(succ=df.sub$accessibility, failure=df.sub$num_cells-df.sub$accessibility))
            #df.sub$trajectory <- as.factor(df.sub$trajectory)
            mod <- glm(df.dat ~ trajectory + nSites, 
                       data=df.sub, 
                       family=stats::binomial())
        }else if(type != "ACRs"){
            df.sub$norm_access <- df.sub$accessibility / df.sub$size_factors
            df.sub$norm_access <- round(df.sub$norm_access)
            mod <- glm(norm_access ~ trajectory + nSites, 
                       data=df.sub, 
                       family=stats::quasipoisson())
        }
        mod.sum <- summary(mod)
        traj <- mod.sum$coefficients["trajectory",]
        names(traj) <- c("Estimate","se","Tval","pval")
        return(traj)
    }, mc.cores=threads)
    outs <- as.data.frame(do.call(rbind, outs))
    rownames(outs) <- rownames(mat)
    
    # estimate q-values
    outs1 <- subset(outs, outs$pval == 1)
    outs2 <- subset(outs, outs$pval != 1)
    outs1$qval <- rep(1, nrow(outs1))
    outs2$qval <- p.adjust(outs2$pval, method="fdr")
    outs <- rbind(outs1, outs2)
    outs <- outs[mixedorder(rownames(outs), decreasing=F),]
    # return
    return(list(pseudotests=outs, aggmat=mat, pseudometa=p.cells))
    
}
sigPseudo2 <- function(obj, meta, type="ACRs", threads=1){
    
    # use
    if(type=="ACRs"){
        use <- obj$a
    }else if(type=="TFs"){
        use <- obj$tf
        use <- use %*% Diagonal(x=1e5/Matrix::colSums(use))
    }else if(type=="MTs"){
        use <- obj$m
    }else if(type=="Genes"){
        use <- obj$gn
        use <- use %*% Diagonal(x=1e5/Matrix::colSums(use))
    }
    
    # align cells
    ids <- intersect(colnames(use), rownames(meta))
    use <- use[,ids]
    meta <- meta[ids,]
    meta$celltypeID <- factor(meta$celltypeID)
    meta$celltypeID <- droplevels(meta$celltypeID)
    
    lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
    }
    
    meta$lib_size <- Matrix::colSums(use)
    
    # run anova per gene
    outs <- mclapply(seq(1:nrow(use)), function(x){
        if((x%%1000)==0){message("   ~ iterated over ",x," gene/peak/motifs ...")}
        if(type=="ACRs"){
            peak1 <- as.numeric(residuals(glm(as.numeric(use[x,])~meta$log10nSites, family=binomial()), type="response"))
        }else{
            peak1 <- as.numeric(use[x,])
        }
        mod <- lm(peak1~ns(meta$trajectory, df=6))
        res <- lmp(mod)
        names(res) <- rownames(use)[x]
        return(res)
    }, mc.cores=threads)
    pval <- do.call(c, outs)
    ids <- names(pval)
    qval <- p.adjust(pval, method="fdr")
    df <- data.frame(pval=pval,qval=qval,row.names=ids)
    return(df)
}
plotTrajHM <- function(obj, pt, cluster=1, prefix="temp", threads=1, top=30000, featureMin=0, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt$LouvainClusters %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    binary <- obj$a[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    binary <- binary[Matrix::rowSums(binary)>featureMin,]
    binary <- binary[,Matrix::colSums(binary)>0]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x ACR matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory, lib=pt$log10nSites)
        mod.scores <- glm(acc~lib, data=df)
        df$res.acc <- residuals(mod.scores)
        mod <- gam(res.acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
        return(zscore)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    write.table(fit, file=paste0(prefix,".ACR_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # plot top 50% by variance
    fit <- t(apply(fit, 1, function(x){rescale(x, c(-1,1))}))
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    
    # plot
    message(" - plotting cell trajectory ...")
    cols <- colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
    pdf(paste0(prefix,".trajectoryHM.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T, 
              ylab=paste("ACRs", paste0("(n=",nrow(fit),")"), sep=" "))
    dev.off()
    
    # return
    return(fit)
    
}
plotTrajMT <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt$LouvainClusters %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    motif <- obj$m[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        motif <- motif[rownames(motif) %in% rownames(tests),]
    }
    pt <- pt[colnames(motif),]
    message(" - filtered zero-sum columns/rows of cell x motif matrix : ",ncol(motif)," | ", nrow(motif))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(motif)),function(x){
        df <- data.frame(acc=as.numeric(motif[x,]), p.time=pt$trajectory)
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
        return(zscore)
    }, mc.cores=threads)
    names(fit) <- rownames(motif)
    fit <- do.call(rbind, fit)
    fit <- t(apply(fit, 1, function(x){rescale(x, c(-1,1))}))
    write.table(fit, file=paste0(prefix,".Mt_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    
    # plot
    message(" - plotting cell trajectory ...")
    cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
    pdf(paste0(prefix,".trajectoryMT.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("Motifs", paste0("(n=",nrow(fit),")"), sep=" "))
    dev.off()
    
    # return
    return(fit)
    
}
plotTrajTF <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt$LouvainClusters %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    obj$tf <- smooth.data(obj$tf, npcs=ncol(obj$d), rds=obj$d)
    binary <- obj$tf[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    binary <- binary[,Matrix::colSums(binary)>0]
    binary <- binary[Matrix::rowSums(binary)>0,]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x gene matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory)
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    write.table(fit, file=paste0(prefix,".TF_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # filter by variances
    fit <- t(apply(fit, 1, function(x){
        rescale(x, c(0,1))
    }))
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]

    # plot
    message(" - plotting cell trajectory ...")
    cols <- colorRampPalette(c("grey80", "grey75",brewer.pal(7, "YlGnBu")[2:7]))(100)
    pdf(paste0(prefix,".trajectoryTF.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("TFs", paste0("(n=",nrow(fit),")"), sep=" "))
    dev.off()
    
    # return
    return(fit)
    
}
plotTrajGN <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt$LouvainClusters %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    obj$gn <- smooth.data(obj$gn, npcs=ncol(obj$d), rds=obj$d)
    binary <- obj$gn[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    binary <- binary[,Matrix::colSums(binary)>0]
    binary <- binary[Matrix::rowSums(binary)>0,]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x gene matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory)
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    write.table(fit, file=paste0(prefix,".Gene_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # filter by variances
    fit <- t(apply(fit, 1, function(x){
        rescale(x, c(0,1))
    }))
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    
    # plot
    message(" - plotting cell trajectory ...")
    cols <- viridis(100)
    pdf(paste0(prefix,".trajectoryGene.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("Genes", paste0("(n=",nrow(fit),")"), sep=" "))
    dev.off()
    
    # return
    return(fit)
    
}
plotSites  <- function(acr, mt, tf, x){
    
    # get sites
    acr.1 <- as.numeric(acr[which(rownames(acr) %in% x[1]),])
    mt.1 <- as.numeric(mt[which(rownames(mt) %in% x[2]),])
    tf.1 <- as.numeric(tf[which(rownames(tf) %in% x[3]),])
    xvals <- seq(from=0, to=100, length.out=length(acr.1))
    
    # plot
    plot(xvals, acr.1, type="l", col="darkorchid", lwd=2)
    lines(xvals, mt.1, type="l", col="darkorange", lwd=2)
    lines(xvals, tf.1, type="l", col="dodgerblue", lwd=2)
    legend("topright", legend=x, pch=16, col=c("darkorchid", "darkorange", "dodgerblue"))
}
findCor    <- function(acr, mt, tf, x){
    acr.use <- as.numeric(acr[x,])
    best.mt <- cor(t(mt), acr.use, method="pearson")
    best.tf <- cor(t(tf), acr.use, method="pearson")
    top.mt <- rownames(best.mt)[which.max(best.mt[,1])]
    top.tf <- rownames(best.tf)[which.max(best.tf[,1])]
    out <- c(x, top.mt, top.tf)
    return(out)
}

###################################################################################################
# load and process raw data -----------------------------------------------------------------------
###################################################################################################
obj <- loadData(sparse, motif, gene, tfgene, meta, svd, doL2=0, traj)
#saveRDS(obj, file=paste0(prefix,".pseudotime.RDS"))

# iterate over types of trajectories
obj$b[,column] <- as.character(obj$b[,column])
out <- calcPseudo(obj$b, obj$d, column=column, traj=traj, cell.dist1=0.95, cell.dist2=0.95)

# plot UMAP traj
pdf(paste0(prefix,".trajectory.pdf"), width=10, height=5)
plotPT(out, subsetCluster=LC, smoother=5, t.id="trajectory", umap1=umap1, umap2=umap2)
dev.off()
write.table(out, file=paste0(prefix,".trajectory.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# find significant peaks across pseudotime
message(" - finding significant ACRs across pseudotime ...")
diff.peaks <- sigPseudo2(obj, out, threads=threads, type="ACRs")
write.table(diff.peaks, file=paste0(prefix,".ACR.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# find significant TFs across pseudotime
message(" - finding significant TFs across pseudotime ...")
diff.TFs <- sigPseudo2(obj, out, threads=threads, type="TFs")
write.table(diff.TFs, file=paste0(prefix,".TFs.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# find significant Genes across pseudotime
message(" - finding significant Genes across pseudotime ...")
diff.Gns <- sigPseudo2(obj, out, threads=threads, type="Genes")
write.table(diff.Gns, file=paste0(prefix,".Genes.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# find significant motifs across pseudotime
message(" - finding significant MTs across pseudotime ...")
diff.MTs <- sigPseudo2(obj, out, threads=threads, type="MTs")
write.table(diff.MTs, file=paste0(prefix,".Motifs.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# plot heatmap traj
acr <- plotTrajHM(obj, out, cluster=LC, prefix=prefix, threads=threads, featureMin=featureMin, tests=diff.peaks, FDR=0.05)
mt <- plotTrajMT(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.MTs, FDR=0.05)
tf <- plotTrajTF(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.TFs, FDR=0.05)
gn <- plotTrajGN(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.Gns, FDR=0.05)

