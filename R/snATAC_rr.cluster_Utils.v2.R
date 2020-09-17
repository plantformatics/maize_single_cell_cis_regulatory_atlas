###################################################################################################
# functions ---------------------------------------------------------------------------------------
###################################################################################################

# visualization -----------------------------------------------------------------------------------
plotUMAP      <- function(b, prefix="out", column="LouvainClusters", m1="umap1", m2="umap2"){
    
    # plot space
    pdf(paste0(prefix,".UMAP_clusters.pdf"), width=6, height=6)
    
    # test if column is present
    if(!column %in% colnames(b)){
        stop(" column header: ", column, ", is missing from input dataframe ...")
    }
    
    if(m1 != "umap1" | m2 != "umap2"){
        b$umap1 <- b[,m1]
        b$umap2 <- b[,m2]
    }
    
    # cols
    if(is.factor(b[,c(column)])){
        b <- b[sample(nrow(b)),]
        cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
        colv <- cols[factor(b[,column])]
    }else if(is.character(b[,column])){
        b[,column] <- factor(b[,column])
        b <- b[sample(nrow(b)),]
        cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
        colv <- cols[factor(b[,column])]
    }else if(is.numeric(b[,column])){
        b <- b[order(b[,column], decreasing=F),]
        cols <- viridis(100)
        colv <- cols[cut(b[,column], breaks=101)]
    }
    
    # plot
    plot(b$umap1, b$umap2, pch=16, cex=0.5, col=colv,
         xlim=c(min(b$umap1), max(b$umap1)+(abs(max(b$umap1))*0.5)))
    
    if(is.factor(b[,column])){
        legend("right", legend=sort(unique(b[,column])), 
               fill=cols[sort(unique(b[,column]))])
    }
    
    # device off
    dev.off()
    
}
plotUMAPstats <- function(x, column="LouvainClusters", palette="Paired", m1="umap1", m2="umap2"){
    
    # require viridis
    require(viridis)
    
    # update coordinates
    x$umap11 <- as.numeric(x[,c(m1)])
    x$umap22 <- as.numeric(x[,c(m2)])
    x <- x[,c("umap11","umap22",column)]
    x <- x[complete.cases(x),]

    # set up color palette
    if(is.factor(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        colc <- colorRampPalette(brewer.pal(12,palette)[1:8])(length(unique(x[,column])))
        colv <- colc[factor(x[,column])]
        ids.names <- levels(factor(x[,column]))
        color.names <- colc[factor(ids.names, levels=ids.names)]
        type <- "fac"
    }else if(is.character(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        colc <- colorRampPalette(brewer.pal(12,palette)[1:8])(length(unique(x[,column])))
        colv <- colc[factor(x[,column])]
        ids.names <- levels(factor(x[,column]))
        color.names <- colc[factor(ids.names, levels=ids.names)]
        type <- "fac"
    }else if(is.numeric(x[,column]) | is.integer(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        numeric.nums <- as.numeric(x[,column])
        message("   *  converted to numeric ...")
        colc <- colorRampPalette(c("grey75","darkorchid4"))(100)
        message("   * successfully specified viridis color palette ...")
        colv <- colc[cut(numeric.nums, breaks=101)]
        type <- "num"
    }else{
        message("   * class for selected column is unsupported: ", class(x[,column]))
    }
    
    # plot
    plot(x$umap11, x$umap22, pch=16, cex=0.2, col=colv, 
         bty="n", xaxt='n', yaxt='n',
         xlab="UMAP1", ylab="UMAP2", main=column,
         xlim=c(min(x$umap11),max(x$umap11)+abs(max(x$umap11)*0.5)))
    
    # legend
    if(type=="fac"){
        legend("right", legend=ids.names, fill=color.names, border=NA, cex=0.5)
    }else if(type=="num"){
        cc <- x[,column]
        half <- signif((min(cc)+max(cc))/2, digits=2)
        minn <- signif(min(cc), digits=2)
        maxx <- signif(max(cc), digits=2)
        legend("right", legend=c(minn,half,maxx), fill=c(colc[1],colc[51],colc[100]), border=NA)
    }
    
    # axes
    axis(1, lwd.tick=0, labels=FALSE)
    axis(2, lwd.tick=0, labels=FALSE)
    
}
plotSTATS     <- function(a, prefix="out", m1="umap1", m2="umap2"){
    
    # adjust factors
    a <- as.data.frame(a)
    a$LouvainClusters <- factor(a$LouvainClusters)
    a$LSI_clusters <- factor(a$LSI_clusters)
    a$tissue <- factor(a$tissue)
    a$library <- factor(a$library)
    a$cellcycle_stage <- factor(a$cellcycle_stage)
    a$logUNIQUE <- as.numeric(as.character(a$logUNIQUE))
    a$log10nSites <- as.numeric(as.character(a$log10nSites))
    a$logTSS <- as.numeric(as.character(a$logTSS))
    a$logPTMT <- as.numeric(as.character(a$logPTMT))
    a$logNUCLEAR <- as.numeric(as.character(a$logNUCLEAR))
    
    # plot UMAP
    message(" - plotting embedding stats ...")
    pdf(paste0(prefix,".UMAP_stats.pdf"), width=16, height=6)
    layout(matrix(c(1:10), nrow=2, byrow=T))
    
    plotUMAPstats(a, m1=m1, m2=m2, column="LouvainClusters")
    plotUMAPstats(a, m1=m1, m2=m2, column="LSI_clusters")
    plotUMAPstats(a, m1=m1, m2=m2, column="tissue")
    plotUMAPstats(a, m1=m1, m2=m2, column="library")
    plotUMAPstats(a, m1=m1, m2=m2, column="cellcycle_stage")
    
    plotUMAPstats(a, m1=m1, m2=m2, column="logUNIQUE")
    plotUMAPstats(a, m1=m1, m2=m2, column="log10nSites")
    plotUMAPstats(a, m1=m1, m2=m2, column="logTSS")
    plotUMAPstats(a, m1=m1, m2=m2, column="logPTMT")
    plotUMAPstats(a, m1=m1, m2=m2, column="logNUCLEAR")
    
    dev.off()
    
    # barplot of proportions
    colc <- brewer.pal(length(unique(a$tissue)),"Paired")
    prop <- prop.table(table(a$tissue, a$LouvainClusters),2)
    
    # plot
    pdf(paste0(prefix,".TissueLouvain.barplot.pdf"), width=15, height=4)
    layout(matrix(c(1:2), nrow=1), widths=c(10,1))
    par(mar=c(2,2,1,1))
    barplot(prop, col=colc, border=NA)
    plot.new()
    legend("center",legend=sort(unique(a$tissue)), fill=colc[sort(unique(a$tissue))], border=NA)
    dev.off()
    
    # adjusted rand index
    a$LSI_clusters <- paste(a$tissue,a$LSI_clusters, sep="_")
    t.ari <- adjustedRandIndex(a$LouvainClusters,a$tissue)
    lsi.ari <- adjustedRandIndex(a$LouvainClusters,a$LSI_clusters)
    
    # verbose
    message(" - tissue ARI = ", t.ari)
    message(" - LSI cluster ARI = ", lsi.ari)
}
plotSTATS2    <- function(a, prefix="out", m1="umap1", m2="umap2"){
    
    # adjust factors
    a$LouvainClusters <- factor(a$LouvainClusters)
    a$LSI_clusters <- factor(a$LSI_clusters)
    a$tissue <- factor(a$tissue)
    a$library <- factor(a$library)
    a$cellcycle_stage <- factor(a$cellcycle_stage)
    a$logUNIQUE <- as.numeric(as.character(a$logUNIQUE))
    a$log10nSites <- as.numeric(as.character(a$log10nSites))
    a$logTSS <- as.numeric(as.character(a$logTSS))
    a$logPTMT <- as.numeric(as.character(a$logPTMT))
    a$logNUCLEAR <- as.numeric(as.character(a$logNUCLEAR))

    # plot UMAP
    message(" - plotting embedding stats ...")
    pdf(paste0(prefix,".UMAP_stats.pdf"), width=16, height=6)
    layout(matrix(c(1:10), nrow=2, byrow=T))
    
    plotUMAPstats(a, m1=m1, m2=m2, column="LouvainClusters")
    plotUMAPstats(a, m1=m1, m2=m2, column="LSI_clusters")
    plotUMAPstats(a, m1=m1, m2=m2, column="tissue")
    plotUMAPstats(a, m1=m1, m2=m2, column="library")
    plotUMAPstats(a, m1=m1, m2=m2, column="cellcycle_stage")
    plotUMAPstats(a, m1=m1, m2=m2, column="logUNIQUE")
    plotUMAPstats(a, m1=m1, m2=m2, column="log10nSites")
    plotUMAPstats(a, m1=m1, m2=m2, column="logTSS")
    plotUMAPstats(a, m1=m1, m2=m2, column="logPTMT")
    plotUMAPstats(a, m1=m1, m2=m2, column="LouvainCluster_sub")
    
    dev.off()
    
    # barplot of proportions
    colc <- brewer.pal(length(unique(as.character(a$tissue))),"Paired")
    prop <- prop.table(table(a$tissue, a$LouvainCluster_sub),2)
    
    # plot
    pdf(paste0(prefix,".TissueLouvain.barplot.pdf"), width=15, height=4)
    layout(matrix(c(1:2), nrow=1), widths=c(10,1))
    par(mar=c(2,2,1,1))
    barplot(prop, col=colc, border=NA)
    plot.new()
    legend("center",legend=sort(unique(as.character(a$tissue))), fill=colc[sort(unique(as.character(a$tissue)))], border=NA)
    dev.off()
    
    # adjusted rand index
    a$LSI_clusters <- paste(a$tissue,a$LSI_clusters, sep="_")
    t.ari <- adjustedRandIndex(a$LouvainCluster_sub,a$tissue)
    lsi.ari <- adjustedRandIndex(a$LouvainCluster_sub,a$LSI_clusters)
    
    # verbose
    message(" - tissue ARI = ", t.ari)
    message(" - LSI cluster ARI = ", lsi.ari)
}
plotRes       <- function(res, prefix="out", cap=T, center=T, verbose=F){
    
    # verbose
    if(verbose){message(" - plotting residual distributions ...")}
    
    # plot set up
    hmcols <- colorRampPalette(c("dodgerblue4","deepskyblue","grey75","darkorange","firebrick3"))(100)
    
    # select top peaks
    set.seed(1111)
    res <- res[sample(nrow(res),10000),]
    
    # process
    if(center){
        if(verbose){message("   * centering residuals ...")}
        res <- (res - rowMeans(res))#/apply(res, 1, var)
    }
    if(cap){
        if(verbose){message("   * capping residual scores ...")}
        cap.t <- sqrt(ncol(res)/50)
        res[res > cap.t] <- cap.t
        res[res < -cap.t] <- -cap.t
        
        res <- t(apply(res, 1, function(z){ 
            #lims <- quantile(z, c(0.05, 0.95))
            #z[z < lims[1]] <- lims[1]
            #z[z > lims[2]] <- lims[2]
            z <- exp(z)/(1+exp(z))
            return(z)
        }))
        #cap.t <- sqrt(ncol(res)/30)
        #res[res > cap.t] <- cap.t
        #res[res < -cap.t] <- -cap.t
    }
    
    # run LSA
    if(verbose){message("   * running LSA ...")}
    LLL <- runLSA(res, 11)
    sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
    diag(sk_diag) <- LLL$sk
    sk_diag[1,1] = 0
    
    # cluster cells and peaks
    if(verbose){message("   * clustering cells ...")}
    LLL_d <- t(sk_diag %*% t(LLL$dk))
    hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
    TTT_d <- t(sk_diag %*% t(LLL$tk))
    if(verbose){message("   * clustering peaks ...")}
    hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")
    
    # plot heatmap
    if(verbose){message("   * plotting heatmap ...")}
    png(paste0(prefix,".residualheatmap.png"), width=6, height=6, units="in", res=300, type="cairo")
    heatmap.2(res, 
              col=hmcols,
              Rowv = as.dendrogram(hclust_genes), 
              Colv = as.dendrogram(hclust_cells),
              labRow=FALSE, labCol=FALSE, 
              trace="none", scale="none",
              useRaster=TRUE)
    
    # device off
    dev.off()
    
}

# processing --------------------------------------------------------------------------------------
loadData      <- function(input, meta){
    
    # verbose
    message(" - loading sparse matrix and meta data ... ")
    
    # format sparse
    a <- read.table(input)
    b <- read.table(meta)
    a <- sparseMatrix(i=as.numeric(a$V1),
                      j=as.numeric(a$V2),
                      x=as.numeric(a$V3),
                      dimnames=list(levels(a$V1),levels(a$V2)))
    
    # order
    both <- intersect(rownames(b), colnames(a))
    a <- a[,both]
    b <- b[both,]
    
    # set ploidy to factor
    if("ploidy" %in% colnames(b)){
        b$ploidy <- factor(b$ploidy)
    }
    if("LSI_clusters" %in% colnames(b)){
        b$LSI_clusters <- factor(b$LSI_clusters)
    }
    if("rdnum" %in% colnames(b)){
	    b$rdnum <- factor(b$rdnum)
    }
    
    # make sure nSites is calculated
    b$nSites   <- Matrix::colSums(a)
    b$log10nSites <- log10(b$nSites)
    
    # return
    return(list(a=a,b=b))
}
outData       <- function(pcs, o.clst, prefix, dev=NULL){
    
    # verbose
    message(" - writing data to disk ...")
    
    # write data to disk
    write.table(pcs, file=paste0(prefix,".SVD.txt"),quote=F,row.names=T,col.names=T,sep="\t")
    write.table(o.clst, file=paste0(prefix,".LouvainCluster.txt"),quote=F,row.names=T,col.names=T,sep="\t")
    
    # check if residuals are passed
    if(!is.null(dev)){
        write.table(dev, file=paste0(prefix,".residuals.txt"),quote=F,row.names=T,col.names=T,sep="\t")
    }
}
cleanData     <- function(x, y, min.p=0.01, min.t=0.001, min.c=100, max.t=0.005){
    
    # verbose
    message("   * Input: cells = ", ncol(x), " | peaks = ", nrow(x))
    
    # min peak count from dist
    y <- y[colnames(x),]
    
    # split cells by cluster
    num.clusts <- table(y$LSI_clusters)
    num.clusts <- num.clusts[num.clusts > 10]
    clusts <- names(num.clusts)
    keep <- lapply(clusts, function(i){
        cells <- rownames(subset(y, y$LSI_clusters==i))
        sub.cells <- x[,colnames(x) %in% cells]
        if(length(cells)>1){
            peak.props <- Matrix::rowMeans(sub.cells)
            return(names(peak.props)[which(peak.props >= min.p)])
        }else{
            return(names(sub.cells)[which(sub.cells > 0)])
        }
        
    })
    keep <- do.call(c, keep)
    keep <- unname(keep, force=T)
    keep <- unique(keep)
    x <- x[rownames(x) %in% keep,]
    message("   * Filtered.1: cells = ", ncol(x), " | peaks = ", nrow(x))
    
    # filter peaks by over all frequency
    x <- x[Matrix::rowSums(x)>(ncol(x)*min.t),]
    x <- x[Matrix::rowSums(x)<(quantile(Matrix::rowSums(x), c(1-max.t))),]
    
    # final clean
    if(min.c < 1){
        x <- x[,Matrix::colSums(x) > quantile(Matrix::colSums(x), min.c)]
    }else{
        x <- x[,Matrix::colSums(x)>min.c]
    }
    
    # last
    x <- x[Matrix::rowSums(x)>0,]
    x <- x[,Matrix::colSums(x)>0]
    
    # verbose
    message("   * Filtered.2: cells = ", ncol(x), " | peaks = ", nrow(x))
    
    # return
    return(x)
    
}
downSample    <- function(x, maxSites=5000){
    
    # if not using scalar
    #maxSites <- min(Matrix::colSums(x))
    maxSites <- as.integer(quantile(Matrix::colSums(x), c(0.9)))
    message(" - downsampling cells to ",maxSites, " max sites ...")
    
    # reads per peak
    x <- x[Matrix::rowSums(x) > 0,]
    x <- x[,Matrix::colSums(x) > 0]
    siteRanks <- Matrix::rowSums(x)
    siteRanks <- siteRanks[order(siteRanks, decreasing=T)]
    
    # downsample, removing constituitive sites first
    cellSize <- Matrix::colSums(x)
    dif <- cellSize - maxSites
    ds.cells <- which(dif > 0)
    num.ds <- length(ds.cells)
    message(" - ",num.ds, " cells above threshold ...")
    for(i in as.numeric(ds.cells)){
        cell <- x[order(siteRanks, decreasing=T),i]
        num.0 <- which(cumsum(cell)==dif[i])
        x[1:num.0,i] <- 0
    }
    
    # finished
    print(head(warnings()))
    return(x)
}
filterSingle  <- function(pro, k=50, threshold1=3, type="umap", prefix="out",
                          m1="umap1", m2="umap2", doplot=T){
    
    # get nearest x neighbors
    if(type=="umap"){
        vars <- c(m1, m2)
    }else if(type=="tsne"){
        vars <- c("tsne1","tsne2")
    }
    topk <- get.knn(pro[,vars], k=k)
    cell.dists <- as.matrix(topk$nn.dist)
    rownames(cell.dists) <- rownames(pro)
    colnames(cell.dists) <- paste0("k",seq(1:ncol(cell.dists)))
    aves <- apply(cell.dists, 1, mean)
    zscore <- as.numeric(scale(aves))
    names(zscore) <- rownames(pro)
    
    # thresholds
    p.zscore <- zscore[order(zscore, decreasing=T)]
    num.pass <- length(zscore[zscore < threshold1])
    
    # plot results
    if(doplot){
        pdf(paste0(prefix,".singletonFilter.pdf"), width=6, height=5)
        plot(p.zscore, pch=16, cex=0.8, col="grey75")
        abline(h=threshold1, col="red", lty=2)
        legend("bottomleft", legend=paste("Pass=",num.pass," (",length(zscore),")"))
        dev.off()
    }
    
    # filter
    prop.good <- zscore[zscore < threshold1]
    ids <- names(prop.good)
    out <- pro[rownames(pro) %in% ids,]
    
    return(out)
}
filtDistClst  <- function(b, cname="seurat_clusters", umap1="umap1", umap2="umap2",threshold=0.1, 
                          prefix="out"){
    
    # hidden functions
    .getZ <- function(x){
    
        # get cluster trimmed mean coordinates
        coord1 <- x[,1]
        coord2 <- x[,2]
        quants1 <- quantile(coord1, c(0.25, 0.75))
        quants2 <- quantile(coord2, c(0.25, 0.75))
        mean1 <- median(coord1[coord1 > quants1[1] & coord1 < quants1[2]])
        mean2 <- median(coord2[coord2 > quants2[1] & coord2 < quants2[2]])
        aves <- c(mean1, mean2)
        
        return(aves)
    }    
    
    # housekeeping
    b$refine_clusters <- b[,c(cname)]
    pass.cells <- subset(b, b$refine_clusters != 0)
    pass.cells$refine_clusters <- paste0("cluster_", pass.cells$refine_clusters)
    non.cells <- subset(b, b$refine_clusters == 0)
    
    # get cluster averages
    clust.aves <- lapply(unique(pass.cells$refine_clusters), function(z){
        b.sub <- subset(pass.cells, pass.cells$refine_clusters==z)
        umap.out.sub <- b.sub[,c(umap1, umap2)]
        clust.mean.coord <- .getZ(umap.out.sub)
        return(clust.mean.coord)
    })
    clust.aves <- data.frame(do.call(rbind, clust.aves))
    rownames(clust.aves) <- unique(pass.cells$refine_clusters)
    colnames(clust.aves) <- c(umap1, umap2)

    # estimate distance to cluster trimmed mean
    cell.dists <- lapply(rownames(clust.aves), function(z){
        b.sub <- subset(pass.cells, pass.cells$refine_clusters==z)
        umap.out.sub <- b.sub[,c(umap1, umap2)]
        cell.dist <- sqrt(colSums((t(umap.out.sub)-as.numeric(clust.aves[z,]))^2))
        cell.z <- as.numeric(scale(cell.dist))
        names(cell.z) <- rownames(umap.out.sub)
        return(cell.z)
    })
    names(cell.dists) <- rownames(clust.aves)
    
    # prep plot
    clusters <- unique(pass.cells$refine_clusters)
    n.cols <- 4
    if(length(clusters) < 4){
        n.cols <- length(clusters)
        n.rows <- 1
    }else{
        n.rows <- ceiling(length(clusters)/n.cols)
    }
    n.total <- n.rows * n.cols
    pdf(paste0(prefix,".refineClusters.step1.pdf"), width=n.cols*4, height=n.rows*4)
    layout(matrix(c(1:n.total), ncol=n.cols, nrow=n.rows, byrow=T))
    par(mar=c(3,3,1,1), oma=c(1,1,1,1))
    
    # filter cells
    cell.rm <- list()
    cell.filt <- lapply(rownames(clust.aves), function(z){
        cell.z <- cell.dists[[z]]
        drop.scores <- vector(length=(length(cell.z)-1))
        cell.z <- cell.z[order(cell.z, decreasing=T)]
        for(i in 1:(length(drop.scores))){
            j <- i + 1
            drop.scores[i] <- abs(cell.z[i] - cell.z[j])
        }
        drop.scores <- drop.scores[cell.z > 0]
        if(any(which(drop.scores >= threshold))){
            idx.thresh <- max(which(drop.scores >= threshold, arr.ind=T))
            z.thresh <- (cell.z[idx.thresh] + cell.z[idx.thresh+1])/2
            if(z.thresh > 3){
                z.thresh <- 3
            }else if(z.thresh < 0){
                z.thresh <- 3
            }
            keep <- names(cell.z)[cell.z < z.thresh]
            rm.cell <- names(cell.z)[!names(cell.z) %in% keep]
        }else{
            z.thresh <- 3
            keep <- names(cell.z)[cell.z < z.thresh]
            rm.cell <- names(cell.z)[!names(cell.z) %in% keep]
        }
        cell.rm[[z]] <<- rm.cell
        
        # plot
        plot(cell.z, pch=16, col=ifelse(cell.z < z.thresh, "grey75", "dodgerblue"),
             main=z)
        abline(h=z.thresh, col="red", lty=2)
        legend("topright", legend=paste("Pass=",length(keep)," (",length(cell.z),")"))
        
        # return
        return(keep)
    })
    
    # device off
    dev.off()
    
    # reformat
    names(cell.filt) <- rownames(clust.aves)
    names(cell.rm) <- rownames(clust.aves)
    
    # make filtered cell data frame
    cell.filt <- unique(do.call(c, cell.filt))
    cell.rm <- unique(do.call(c, cell.rm))
    cell.rm <- cell.rm[!is.na(cell.rm)]
    cell.rm <- c(cell.rm, rownames(non.cells))
    b.filt <- as.data.frame(b[cell.filt,])
    b.updated <- as.data.frame(b[cell.rm,])
    
    # plot updated UMAP
    pdf(paste0(prefix,".refineClusters.step2.pdf"), width=5, height=5)
    cols <- brewer.pal(length(unique(b.filt$refine_clusters)), "Paired")
    plot(b.filt[,c(umap1)], b.filt[,c(umap2)], pch=16, cex=0.3, col=cols[factor(b.filt$refine_clusters)],
         xlab="umap1",ylab="umap2", bty='n')
    dev.off()
    
    # update cluster aves
    clust.aves <- lapply(unique(b.filt$refine_clusters), function(z){
        b.sub <- subset(b.filt, b.filt$refine_clusters==z)
        umap.out.sub <- b.sub[,c(umap1, umap2)]
        clust.mean.coord <- .getZ(umap.out.sub)
        return(clust.mean.coord)
    })
    clust.aves <- data.frame(do.call(rbind, clust.aves))
    rownames(clust.aves) <- unique(b.filt$refine_clusters)
    colnames(clust.aves) <- c(umap1, umap2)
    
    # iterate over filtered cells, assign to nearest cluster and redo filtering
    new.cell.clust <- lapply(cell.rm, function(z){
        info <- b.updated[z,]
        cell.coord <- b.updated[z,c(umap1,umap2)]
        cell.dist <- sqrt(colSums(t(clust.aves)-as.numeric(cell.coord))^2)
        names(cell.dist) <- rownames(clust.aves)
        new.clst <- names(cell.dist)[which.min(cell.dist)]
        info$refine_clusters <- new.clst
        
        return(info)
    })
    new.cell.clust <- as.data.frame(do.call(rbind, new.cell.clust))
    rownames(new.cell.clust) <- new.cell.clust$cellID
    
    # merge clean and updated cell clusters
    step.2.b <- rbind(b.filt, new.cell.clust)
    
    # prep plot
    clusters <- unique(step.2.b$refine_clusters)
    n.cols <- 4
    if(length(clusters) < 4){
        n.cols <- length(clusters)
        n.rows <- 1
    }else{
        n.rows <- ceiling(length(clusters)/n.cols)
    }
    n.total <- n.rows * n.cols
    pdf(paste0(prefix,".refineClusters.step3.pdf"), width=n.cols*4, height=n.rows*4)
    layout(matrix(c(1:n.total), ncol=n.cols, nrow=n.rows, byrow=T))
    par(mar=c(3,3,1,1), oma=c(1,1,1,1))
    
    # redo cell distances
    final <- lapply(clusters, function(z){
        
        # get cell distances
        b.sub <- subset(step.2.b, step.2.b$refine_clusters==z)
        b.sub.umap <- b.sub[,c(umap1,umap2)]
        aves <- as.numeric(clust.aves[z,])
        cell.dist <- sqrt(colSums((t(b.sub.umap)-aves)^2))
        cell.z <- as.numeric(scale(cell.dist))
        names(cell.z) <- rownames(b.sub)
        cell.z <- cell.z[!is.na(cell.z)]
        
        # estimate drop scores and filter
        drop.scores <- vector(length=(length(cell.z)-1))
        cell.z <- cell.z[order(cell.z, decreasing=T)]
        for(i in 1:(length(drop.scores))){
            j <- i + 1
            drop.scores[i] <- abs(cell.z[i] - cell.z[j])
        }
        drop.scores <- drop.scores[cell.z > 0]
        if(any(which(drop.scores >= threshold))){
            idx.thresh <- max(which(drop.scores >= threshold, arr.ind=T))
            z.thresh <- (cell.z[idx.thresh] + cell.z[idx.thresh+1])/2
            if(z.thresh > 3){
                z.thresh <- 3
            }else if(z.thresh < 0){
                z.thresh <- 3
            }
            keep <- names(cell.z)[cell.z < z.thresh]
        }else{
            z.thresh <- 3
            keep <- names(cell.z)[cell.z < z.thresh]
        }
        
        # plot
        plot(cell.z, pch=16, col=ifelse(cell.z < z.thresh, "grey75", "dodgerblue"),
             main=z)
        abline(h=z.thresh, col="red", lty=2)
        legend("topright", legend=paste("Pass=",length(keep)," (",length(cell.z),")"))
        
        # return
        return(keep)
    })
    
    # close device
    dev.off()
    
    # filter
    final <- do.call(c, final)
    step.2.b <- step.2.b[final,]
    step.2.b$refine_clusters <- gsub("cluster_","",step.2.b$refine_clusters)
    step.2.b$refine_clusters <- as.numeric(factor(step.2.b$refine_clusters))
    step.2.b[,c(cname)] <- step.2.b$refine_clusters
    step.2.b$refine_clusters <- NULL
    b.filt <- step.2.b[intersect(rownames(step.2.b), rownames(b.filt)),]
    
    # plot updated UMAP
    pdf(paste0(prefix,".refineClusters.step4.pdf"), width=5, height=5)
    cols <- colorRampPalette(brewer.pal(6, "Paired"))(length(unique(step.2.b[,cname])))
    plot(step.2.b[,c(umap1)], step.2.b[,c(umap2)], pch=16, cex=0.3, 
         col=cols[factor(step.2.b[,cname])],
         xlab="umap1",ylab="umap2", bty='n')
    dev.off()
    
    # return
    return(list(bb=step.2.b, bb.filt=b.filt))
}
filtDistClst2 <- function(b, umap1="umap1", umap2="umap2", threshold2=2, prefix="out"){
 
    # iterate over each cluster
    clusts <- unique(b$seurat_clusters)
    out <- lapply(clusts, function(x){
        b.sub <- subset(b, b$seurat_clusters == x)
        b.umap <- b.sub[,c(umap1, umap2)]
        out.b.umap <- filterSingle(b.umap, k=25, threshold1=threshold2, doplot=F, m1=umap1, m2=umap2)
        return(rownames(out.b.umap))
    })
    out <- do.call(c, out)
    b.out <- b[rownames(b) %in% as.character(out),]
    message("   * total number of cells surviving subcluster filtering = ", nrow(b.out))
    return(b.out)
       
}
formatRes     <- function(dev, cap=F, center=F, stdize=F){
    
    # just do stdize if center/std both specified
    if(center==T & stdize==T){
        center <- F
    }
    
    # cap
    if(cap){
        top <- sqrt(ncol(dev)/50)
        dev[dev > top] <- top
        dev[dev < -top] <- -top
    }
    
    # center
    if(center){
        dev <- dev - rowMeans(dev)
    }
    
    # standarize
    if(stdize){
        dev <- t(scale(t(dev)))
        dev <- as.matrix(new)
    }
    return(dev)
}
reportResults <- function(b, column="LouvainClusters"){
    
    # reads per cluster
    print(aggregate(b$unique~b[,c(column)], FUN=sum))
    
    # cells per cluster
    print(table(b[,c(column)]))
    
    # total cells
    print(sum(table(b[,c(column)])))
    
}
stddev        <- function(r.dev, threads=1, center=T, scale=T, bins=1024){
    
    # verbose
    message(" - standardizing deviance scores ...")
    
    # set up bins
    bin_ind <- ceiling(x = 1:nrow(r.dev) / bins)
    max_bin <- max(bin_ind)
    ids <- rownames(r.dev)
    
    # run in bins
    dev <- lapply(seq(1:max_bin), function(x){
        peaks_bin <- rownames(r.dev)[bin_ind == x]
        if(center==T & scale==F){
            t(apply(r.dev[peaks_bin,], 1, function(z){z-mean(z, na.rm=T)}))
        }else{
            t(apply(r.dev[peaks_bin,], 1, function(z){(z-mean(z, na.rm=T))/sd(z, na.rm=T)}))
        }
    })
    rm(r.dev)
    
    # merge and return
    dev <- do.call(rbind, dev)
    dev <- dev[ids,]
    
    return(dev)
}

# reduce dims -------------------------------------------------------------------------------------
reduceDims    <- function(x, df, n.pcs=50, batch=NULL, scaleVar=F, doL2=1, lambda=0.5, theta=2, 
                          stdize=F, cap=F, prefix="out", verbose=F, raw=NULL, center=F, cor.max=0.75,
                          useTop=NULL, useTop.th=1.25, stdLSI=NULL){
    
    # verbose
    message(" - reduce dimensions with SVD ... ")
    if(cap){
        message("   * capping residual scores ...")
        cap.t <- 10 #sqrt(ncol(x)/30)
        x[x > cap.t] <- cap.t
        x[x < -cap.t] <- -cap.t
    }
    
    # run
    message("   * running irlba-SVD ...")
    pcs <- irlba(t(x), nv=(n.pcs+1))
    pc <- pcs$u
    if(scaleVar){
        message("   * scaling PCs by variance ... ")
        pc <- pc %*% diag(pcs$d)
    }
    
    # add names
    rownames(pc) <- colnames(x)
    colnames(pc) <- paste0("PC_", seq(1:ncol(pc)))
    write.table(pc, file=paste0(prefix,".rawSVD.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # plot variance explained
    var.exp <- pcs$d^2/sum(pcs$d^2)
    pdf(paste0(prefix,".perc_variance.pdf"), width=5, height=5)
    plot(var.exp, xlab="PCs", ylab="% variances", pch=16)
    dev.off()
    pc <- pc[,c(1:n.pcs)]
    
    # plot correlation with read-depth
    if(!is.null(raw)){
        raw1 <- raw[,rownames(pc)]
        depth <- Matrix::colSums(raw1)
        cors <- apply(pc, 2, function(u) cor(u,depth,method="spearman"))
        pc <- pc[,abs(cors) < cor.max]
        pdf(paste0(prefix,".PC_depth_correlation.pdf"), width=5, height=5)
        plot(cors, xlab="PCs", ylab="Spearman's rho", type="h", lwd=3, col="grey75")
        abline(h=cor.max, col="red", lty=2)
        dev.off()
        rm(raw1)
    }
    
    # remove batch effects
    if(!is.null(batch)){
        
        # verbose
        message(" - removing batch effects with harmony", batch, " ... ")
        df <- df[rownames(pc),]
        df$rdnum <- factor(ntile(df$unique, 5))
        
        # run
        pc <- HarmonyMatrix(data_mat=pc, 
                            meta_data=df,
                            vars_use=batch,
                            do_pca=F,
                            max.iter.harmony=30,
                            theta=theta,
                            lambda=lambda, 
                            sigma=0.1, 
                            nclust=100, 
                            tau=5)
        
    }
    
    # convert l2 norm
    if(!is.null(doL2)){
        if(doL2 != FALSE){
            
            # verbose
            message(" - L2 norm of reduced dimensions ... ")
            
            if(doL2==1){
                pc <- t(apply(pc, 1, function(x) x/(sqrt(sum(x^2)))))
            }else if(doL2==2){
                pc <- apply(pc, 2, function(x) x/(sqrt(sum(x^2))))
            }else if(doL2==3){
                pc2 <- pc^2
                pc1 <- pc/abs(pc)
                pc <- pc2*pc1
            }
            
            if(!is.null(raw)){
                raw1 <- raw[,rownames(pc)]
                depth <- Matrix::colSums(raw1)
                cors <- apply(pc, 2, function(u) cor(u,depth,method="spearman"))
                pc <- pc[,abs(cors) < cor.max]
                pdf(paste0(prefix,".L2PC_depth_correlation.pdf"), width=5, height=5)
                plot(cors, xlab="PCs", ylab="Spearman's rho", type="h", lwd=3, col="grey75")
                abline(h=cor.max, col="red", lty=2)
                dev.off()
                rm(raw1)
            }
        }
    }
    
    # standarize reduced dimensions per cell
    if(!is.null(stdLSI)){
        if(stdLSI != FALSE){
            
            # verbose
            if(stdLSI==1){
                message(" - standardizing reduced dimensions per cell ... ")
                pc <- t(apply(pc, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
            }else if(stdLSI==2){
                message(" - standardizing reduced dimensions by componenet ... ")
                pc <- apply(pc, 2, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)})
            }else{
                message(" - standardizing reduced dimensions by componenet ... ")
                pc <- apply(pc, 2, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)})
            }
        }
    }
    
    # return
    return(pc)
    
}
projectUMAP   <- function(pcs, m.dist=0.1, k.near=15, metric="cosine"){
    
    # verbose
    message(" - non-linear dimensionality reduction with UMAP ...")
    
    # run UMAP (uwot implementation)
    umap.res <- umap(pcs, verbose=F, min_dist=m.dist, n_neighbors=k.near, metric=metric)
    umap.out <- as.data.frame(umap.res)
    rownames(umap.out) <- rownames(pcs)
    colnames(umap.out) <- c("umap1", "umap2")
    return(umap.out)
}
tfidf         <- function(bmat, frequencies=T, log_scale_tf=T, scale_factor=100000){
        
        # hidden functions
        .safe_tfidf       <- function(tf, idf,  block_size=2000e6){
            result = tryCatch({
                result = tf * idf
                result
            }, error = function(e) {
                options(DelayedArray.block.size=block_size)
                DelayedArray:::set_verbose_block_processing(TRUE)
                
                tf = DelayedArray(tf)
                idf = as.matrix(idf)
                
                result = tf * idf
                result
            })
            return(result)
        }
        
        # Use either raw counts or divide by total counts in each cell
        if (frequencies) {
            # "term frequency" method
            tf = t(t(bmat) / Matrix::colSums(bmat))
        } else {
            # "raw count" method
            tf = bmat
        }
        
        # Either TF method can optionally be log scaled
        if (log_scale_tf) {
            if (frequencies) {
                tf@x = log1p(tf@x * scale_factor)
            } else {
                tf@x = log1p(tf@x * 1)
            }
        }
        
        # IDF w/ "inverse document frequency smooth" method
        idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
        
        # TF-IDF
        tf_idf_counts = .safe_tfidf(tf, idf)
        rownames(tf_idf_counts) = rownames(bmat)
        colnames(tf_idf_counts) = colnames(bmat)
        return(Matrix(tf_idf_counts, sparse=T))
    }

# clustering --------------------------------------------------------------------------------------
callClusters  <- function(a, b, pc, umap.out, res=0.02, k.near=20, clustOB="umap", clustType="louvain",
                          cname="LouvainClusters", min.reads=1e6, m.clst=25, threshold3=2, 
                          prefix="out", umap1="umap1", umap2="umap2", key="UMAP_", e.thresh=2, 
                          dynamic=F, cleanCluster=T, cl.method=1){
    
    # filter umap coordinates
    message(" - filtering outliers in UMAP manifold (z-score e.thresh = ", e.thresh, ") ...")
    umap.out <- filterSingle(umap.out, threshold1=e.thresh, prefix=prefix, m1=umap1, m2=umap2)
    a <- a[,rownames(umap.out)]
    b <- b[rownames(umap.out),]
    pc <-  pc[rownames(umap.out),]
    
    # densityCluster or Louvain clustering
    if(clustType=="densityClust"){
        
        if(clustOB=="svd"){
            use <- pc
        }else if(clustOB=="umap"){
            use <- as.data.frame(umap.out)
        }else{
            use <- as.data.frame(umap.out)
        }
        
        # run DCLUST
        message(" - running density-based clustering ... ")
        dclust <- densityClust(dist(use), gaussian=T)
        dclust <- findClusters(dclust, rho=2, delta=1)
        
        # plot dclust thresholds
        pdf(paste0(prefix,".dclust_thresholds.pdf"), width=5, height=5)
        plot(dclust)
        dev.off()
        
        # filter by cluster size
        agg.reads <- aggregate(b$unique~dclust$clusters, FUN=sum)
        colnames(agg.reads) <- c("clusters","readDepth")
        table(dclust$clusters)
        sro.meta <- b
        sro.meta$clusters <- dclust$clusters
        clust.cnts <- table(sro.meta$clusters)
        agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
        agg.pass <- subset(agg.reads, agg.reads$num_cells>=m.clst & agg.reads$readDepth>=min.reads)
        sro.filt <- sro.meta[sro.meta$clusters %in% agg.pass$clusters,]
        sro.filt$clusters <- as.numeric(as.factor(sro.filt$clusters))
        umap.out <- umap.out[rownames(sro.filt),]
        b <- b[rownames(sro.filt),]
        b[,cname] <- sro.filt[rownames(b),]$clusters
        b[,c(umap1)] <- umap.out[,c(umap1)]
        b[,c(umap2)] <- umap.out[,c(umap2)]
        
        # if clean cluster
        if(cleanCluster){
            b <- filtDistClst(b, cname=cname, umap1=umap1, umap2=umap2, threshold=threshold)
        }
        
    }else{
        
        # verbose
        message(" - running default clustering: louvain ...")
        message(" - creating seurat object for graph-based clustering ...")
        
        # dynamic
        start.cells <- nrow(b)
        finish <- 0
        its <- 0
        while(finish==0){
            
            # update iteration
            its <- its + 1
            message("   * louvain clustering iteration = ", its)
            
            # create Seurat object, find clusters
            sro <- CreateSeuratObject(a, min.cells=0, min.features=0)
            sro[["svd"]] <- CreateDimReducObject(embeddings = pc, key = "PC_", assay = DefaultAssay(sro))
            sro[["umap"]] <- CreateDimReducObject(embeddings=as.matrix(umap.out), key="umapsub_", assay=DefaultAssay(sro))
            sro <- AddMetaData(sro, b)
            sro <- FindNeighbors(sro, dims = 1:ncol(sro[[clustOB]]), reduction=clustOB, nn.eps=0, k.param=k.near)
            sro <- FindClusters(sro, resolution=res, reduction=clustOB, n.start=100, algorithm=cl.method)
            sro.meta <- data.frame(sro@meta.data)
            sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
            
            # remove outliers?
            if(cleanCluster){
                
                # verbose
                message("   * removing low quality clusters ...")
                
                # prep
                sro.umap <- umap.out[rownames(sro.meta),]
                colnames(sro.umap) <- c(umap1, umap2)
                sro.meta <- cbind(sro.meta, sro.umap)
                
                # filter by cluster size and number of cells
                agg.reads <- aggregate(sro.meta$unique~sro.meta$seurat_clusters, FUN=sum)
                colnames(agg.reads) <- c("clusters","readDepth")
                clust.cnts <- table(sro.meta$seurat_clusters)
                agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
                agg.pass <- subset(agg.reads, agg.reads$num_cells >= m.clst & agg.reads$readDepth >= min.reads)
                sro.filt <- sro.meta[as.character(sro.meta$seurat_clusters) %in% as.character(agg.pass$clusters),]
                sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
                sro.filt$seurat_clusters <- as.numeric(factor(sro.filt$seurat_clusters))
                
                # remove outliers in the embedding
                message("   * filtering per-cluster outliers (z-score filtDistClst2 = ", threshold3, ") ...")
                sro.meta <- filtDistClst2(sro.filt, umap1=umap1, umap2=umap2, threshold2=threshold3)
                sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
            }
            
            # filter by cluster size
            agg.reads <- aggregate(sro.meta$unique~sro.meta$seurat_clusters, FUN=sum)
            colnames(agg.reads) <- c("clusters","readDepth")
            clust.cnts <- table(sro.meta$seurat_clusters)
            agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
            agg.pass <- subset(agg.reads, agg.reads$num_cells>=m.clst & agg.reads$readDepth>=min.reads)
            sro.filt <- sro.meta[sro.meta$seurat_clusters %in% agg.pass$clusters,]
            sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
            sro.filt$seurat_clusters <- as.numeric(as.factor(sro.filt$seurat_clusters))
            
            # check
            if(dynamic==F){
                finish <- 1
            }else if(nrow(sro.filt) >= ceiling(0.9*start.cells)){
                finish <- 1
            }else{
                if(res >= 1){
                    res <- res*0.95
                }else{
                    res <- res^2
                }
            }
            
            # exit if too many iterations
            if(res < 1e-2){
                sro.filt <- sro.meta
                sro.filt$seurat_clusters <- 1
                finish <- 1
            }
        }
        
        umap.out <- umap.out[rownames(sro.filt),]
        clusts <- sro.filt$seurat_clusters
        b <- b[rownames(sro.filt),]
        b[,cname] <- clusts
        b[,c(umap1)] <- umap.out[,c(umap1)]
        b[,c(umap2)] <- umap.out[,c(umap2)]
        
    }
    
    
    # return
    return(b)
}

# models ------------------------------------------------------------------------------------------
logisticModel <- function(x, y, r.variates="~log10nSites", variates="~log10nSites", method="glm", 
                          nthreads=1, chunks=256, subpeaks=NULL, logout=T, output="model.log", 
                          type="response", prefix="out", doPlot=F, link="logit", writeRawRes=T,
                          doCenter=F, doSTD=T, cap=F, topSites=25000, var.th=1){
    
    # verbose
    message(" - running logistic regression ...")
    
    # set-up data
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    yo <- y
    form <- as.formula(paste0("z",r.variates))
    message("   * formula: ", form)
    
    # select subsample of peaks
    if(!is.null(subpeaks)){
        set.seed(1111)
        log_peak_mean <- log10(Matrix::rowMeans(x))
        log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
        sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean)$y + .Machine$double.eps)
        peak.sub <- sample(x=rownames(x), size=subpeaks, prob=sample_prob)
        log_peak_mean_sub <- log10(Matrix::rowMeans(x[peak.sub,]))
        x <- x[peak.sub,]
        message("   * subsampling: cells = ", ncol(x), " | peaks = ", nrow(x))
    }
    
    # use method
    do.func <- function(form, df){
        
        # run GLM
        suppressWarnings(mod <- glm(form, 
                                    data=df, 
                                    family=quasibinomial(link = link)))
        
        # return residuals
        return(residuals.glm(mod, type=type))
    }
    
    # parallel
    if(nthreads == 1){
        
        # transpose
        x <- t(x)
        
        # begin iterate
        its <- 0
        res <- lapply(seq(1:ncol(x)), function(i, y, form){
            
            # verbose about progress
            its <<- its + 1
            if(its %% 1000 == 0){message("   * estimated residual devaince for ", its, " peaks ...")}
            
            # create DF
            df <- cbind(y, x[,i])
            colnames(df) <- c(colnames(y), "z")
            
            # run GLM
            do.func(form, df)
            
        }, y=y, form=form)
        
        # reform res
        res <- as.matrix(do.call(rbind, res))
        rownames(res) <- colnames(x)
        colnames(res) <- rownames(x)
        
    }else{
        
        # log file
        if(file.exists(output)){
            unlink(output)
        }
        
        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)
        
        # control looping/chunks
        peaknames <- rownames(x)
        niter <- nrow(x)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }
        
        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        package.labs <- c("Matrix")
        
        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it) 
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }
        
        # run logistic regression for each peak independently - in parallel
        res <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts,
                       .packages=package.labs) %dopar% {
                           
                           # print to log
                           if(logout){
                               ctime <- Sys.time()
                               chunkn <- ceiling(n$i/chunks)
                               lastp <- n$i + n$m - 1
                               cat(paste("SENT: chunk # ",chunkn,"\t","| peaks","\t", n$i, " - ", lastp,
                                         "\t","| current:", peaknames[n$i], "\t","| current - ",
                                         ctime, "\n",sep=""),
                                   file=paste(getwd(),"/", output, sep=""), append=T)
                           }
                           
                           # run chunks
                           its <- c(n$i:(n$i+n$m-1))
                           dfout <- lapply(its, function(j){
                               
                               # create DF
                               df <- cbind(y, x[j,])
                               colnames(df) <- c(colnames(y), "z")
                               
                               # run GLM
                               do.func(form, df)
                               
                           })
                           
                           # reformat
                           dfout <- do.call(rbind, dfout)
                           rownames(dfout) <- rownames(x)[n$i:(n$i+n$m-1)]
                           return(dfout)
                           
                       }
        
        # close connections
        close(pb)
        stopCluster(cl)
        
        # return results
        colnames(res) <- colnames(x)
        
        # plot
        if(doPlot){
            res.range <- range(res)
            message("   * residual range: ",res.range[1], " - ", res.range[2])
            plotRes(res, prefix=paste0(prefix,".uncorrected"))
        }
        
        
        ###########################################################################################
        # variables to regress from residuals -----------------------------------------------------
        ###########################################################################################
        if(!is.null(variates)){
            
            # verbose
            message(" - removing confounding variables ...")
            
            # model
            form.nr <- as.formula(paste0("z",variates))
            message("   * formula: ", form.nr)
            
            # parallel parameters
            cl <- makeSOCKcluster(nthreads)
            registerDoSNOW(cl)
            
            # control looping/chunks
            chunks <- 1000
            peaknames <- rownames(x)
            niter <- nrow(x)
            tasks <- ceiling(niter/chunks)
            if(tasks < nthreads){
                message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
            }
            
            # track vars
            pb <- txtProgressBar(max = tasks, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            
            # functions
            idivix <- function(n, chunkSize) {
                i <- 1
                it <- idiv(n, chunkSize=chunkSize)
                nextEl <- function() {
                    m <- nextElem(it) 
                    value <- list(i=i, m=m)
                    i <<- i + m
                    value
                }
                obj <- list(nextElem=nextEl)
                class(obj) <- c('abstractiter', 'iter')
                obj
            }
            
            # calc first QR
            y1 <- yo[colnames(res),]
            regression.mat <- cbind(y1, res[1,])
            colnames(regression.mat) <- c(colnames(y1),"z")
            qr <- lm(form.nr, data = regression.mat, qr = TRUE)$qr
            rm(regression.mat)
            r.names <- rownames(res)
            
            # run logistic regression for each peak independently - in parallel
            res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts) %dopar% {
                
                # run chunks
                its <- c(n$i:(n$i+n$m-1))
                dfout <- lapply(its, function(j){
                    
                    # run GLM
                    qr.resid(qr = qr, y = res[j,])
                    
                })
                
                # reformat
                dfout <- do.call(rbind, dfout)
                rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
                return(dfout)
                
            }
            
            # close connections
            close(pb)
            stopCluster(cl)
            
            # clean
            res <- res2
            rm(res2)
            res <- res[r.names,]
            
            # plot
            if(doPlot){
                res.range <- range(res)
                message("   * residual range: ",res.range[1], " - ", res.range[2])
                plotRes(res, prefix=paste0(prefix,".corrected"))
            }
            
        }
        
        
        ###########################################################################################
        # return highly variable peaks ------------------------------------------------------------
        ###########################################################################################
        p.vars <- apply(res, 1, var)
        p.aves <- Matrix::rowMeans(x)
        p.vars <- p.vars[order(p.vars, decreasing=T)]
        p.aves <- p.aves[order(p.vars, decreasing=T)]
        
        # keep peaks
        keepers <- names(p.vars[p.vars>var.th])
        if(length(keepers) < topSites){
            keepers <- names(p.vars[1:topSites])
        }
        res <- res[rownames(res) %in% keepers, ]
        num.keep <- nrow(res)
        ids.c <- colnames(res)
        message(" - kept top ",num.keep, " variable peaks ...")
        if(doSTD==T){
            res <- t(apply(t(res), 2, function(r) (r-mean(r, na.rm=T))/sd(r, na.rm=T) ))
        }
        colnames(res) <- ids.c
        if(cap){
            res[res > 10] <- 10
            res[res < -10] <- -10
        }
        
        # plot peak usage vs. residual variance
        pdf(paste0(prefix,".residualVariance.pdf"), width=5, height=5)
        plot(p.aves, p.vars, pch=16, cex=0.5, col=ifelse(p.vars < var.th, "black", "red"),
             xlab="Mean", ylab="Residual variance", ylim=range(p.vars), 
             xlim=range(p.aves))
        dev.off()
        
        # return
        return(res)   
    }
}
regModel      <- function(x, y, r.variates='~log10nSites', variates="~log10nSites", subpeaks=5000, 
                          bins=256, bw_adjust=10, type="response", verbose=T, plotCoef=T, 
                          link="logit", nthreads=1, doPlot=F, weights=NULL, prefix="oh", doSTD=T,
                          doCenter=F, method="elasticNet", alpha=0.5, writeRawRes=T, var.th=1,
                          topSites=25000){
    
    library(doSNOW)
    
    # fix type if necessary
    if(type!="deviance" & type !="pearson" & type!="response" & type !="working"){
        message(" - *type* incorrect, setting residual to default (pearson) ...")
        type <- "pearson"
    }
    
    # functions
    is_outlier          <- function(y, x, th=10) {
        bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
        eps <- .Machine$double.eps * 10
        breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
        breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
        score1 <- robust_scale_binned(y, x, breaks1)
        score2 <- robust_scale_binned(y, x, breaks2)
        return(pmin(abs(score1), abs(score2)) > th)
    }
    robust_scale_binned <- function(y, x, breaks) {
        bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
        tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
        score <- rep(0, length(x))
        o <- order(bins)
        if (inherits(x = tmp$x, what = 'list')) {
            score[o] <- unlist(tmp$x)
        } else {
            score[o] <- as.numeric(t(tmp$x))
        }
        return(score)
    }
    deviance_residual   <- function(y, mu, wt){
        d.res <- sqrt(pmax((binomial()$dev.resid)(y, mu, wt),0))
        d.res <- ifelse(y > mu, d.res, -d.res)
        return(d.res)
    }
    pearson_residual    <- function(y, mu, theta){(y-mu)/sqrt((mu*(1-mu))*theta)}
    robust_scale        <- function(x){return((x - median(x)) / (mad(x) + .Machine$double.eps))}
    ChunkPoints         <- function(dsize, csize) {
        return(vapply(
            X = 1L:ceiling(x = dsize / csize),
            FUN = function(i) {
                return(c(
                    start = (csize * (i - 1L)) + 1L,
                    end = min(csize * i, dsize)
                ))
            },
            FUN.VALUE = numeric(length = 2L)
        ))
    }
    invprobit           <- function(x){
        thresh <- -qnorm(.Machine$double.eps)
        x <- pmin(pmax(x, -thresh), thresh)
        pnorm(x)
    }
    
    
    ###############################################################################################
    # start ---------------------------------------------------------------------------------------
    ###############################################################################################
    message(" - regularizing logistic model parameters ...")
    
    # set-up data
    yo <- y
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    form <- as.formula(paste0("z",r.variates))
    message("   * formula: ", form)
    
    #select subsample of peaks
    set.seed(1111)
    con <- 1
    if(verbose){message("   * estimating geometric mean ...")}
    log_peak_mean <- log(Matrix::rowMeans(x))
    if(verbose){message("   * density sampling on peak space ...")}
    log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
    sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean)$y + .Machine$double.eps)
    peak.sub <- sample(x=rownames(x), size=subpeaks, prob=sample_prob)
    x.sub <- x[peak.sub,]
    log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))
    
    # fit model to subset
    message("   * fitting parameters to ",nrow(x.sub)," peaks ...")
    pars <- lapply(seq(1:nrow(x.sub)), function(j){
        df <- cbind(y, x.sub[j,])
        colnames(df) <- c(colnames(y), "z")
        suppressWarnings(mod <- glm(form, 
                                    data=df, 
                                    family=quasibinomial(link = link)))
        theta <- summary(mod)$dispersion
        names(theta) <- "theta"
        return(c(theta, mod$coefficients))
    })
    
    # merge
    pars <- do.call(rbind, pars)
    rownames(pars) <- rownames(x.sub)
    
    
    ###############################################################################################
    ### regularize --------------------------------------------------------------------------------
    ###############################################################################################
    
    # outliers
    message("   * finding outliers ...")
    peaks <- names(log_peak_mean)
    outliers <- apply(pars, 2, function(k) is_outlier(k, log_peak_mean_sub))
    outliers <- apply(outliers, 1, any)
    if (sum(outliers) > 0){
        message("   * found ", sum(outliers), " outliers ...")
        pars <- pars[!outliers, ]
        peaks_step1 <- rownames(pars)
        log_peak_mean_sub <- log_peak_mean_sub[!outliers]
    }
    
    # select bw
    bw <- bw.SJ(log_peak_mean) * bw_adjust
    
    # parameter for predictions
    x_points <- pmax(log_peak_mean, min(log_peak_mean))
    x_points <- pmin(x_points, max(log_peak_mean))
    
    # take results from step 1 and fit/predict parameters to all genes
    message("   * regularizing all coefficients ...")
    o <- order(x_points)
    pars_fit <- matrix(NA_real_, length(peaks), ncol(pars),
                       dimnames = list(peaks, colnames(pars)))
    
    # global fit / regularization for all coefficients
    for (i in 1:ncol(pars)){
        pars_fit[o, i] <- ksmooth(x = log_peak_mean_sub, y = pars[, i],
                                  x.points = x_points, bandwidth = bw, kernel='normal')$y
        #suppressWarnings(fit <- smooth.spline(log_peak_mean_sub, pars[,i], cv=T))
        #pars_fit[o, i] <- predict(fit, x=x_points[o])$y
    }
    
    # plot coefficients
    if(plotCoef){
        
        # verbose
        message("   * plotting coefficients before and after regularization ...")
        
        # set up plot
        n.cols <- ncol(pars_fit)
        n.rows <- 1
        pdf(paste0(prefix,".coefficient.pdf"), width=(5*n.cols), height=(n.rows*5))
        layout(matrix(c(1:(n.cols*n.rows)), nrow=n.rows, byrow=T))
        raw.coeff <- as.data.frame(pars)
        raw.coeff$logpeakmeans <- log_peak_mean_sub
        write.table(raw.coeff, file=paste0(prefix,".coefficients.sub.txt"), quote=F, row.names=T, col.names=T, sep="\t")
        
        # iterate over coefficients
        for(i in 1:ncol(pars)){
            plot(log_peak_mean_sub, pars[,i], pch=19, cex=0.7, col=alpha("black", 0.5), 
                 main=colnames(pars)[i])
            df <- data.frame(x=x_points, y=pars_fit[,i])
            df <- df[order(df$x, decreasing=F),]
            lines(df$x, df$y, col=alpha("red", 0.5), lwd=2)
        }
        
        # device off
        dev.off()
        
    }
    
    
    ###############################################################################################
    # fit data ------------------------------------------------------------------------------------
    ###############################################################################################
    
    # affine transform
    message("   * estimating deviance residuals on the full data set ...")
    regressor_data <- model.matrix(as.formula(r.variates), data=y)
    
    # split peaks into bins
    bin_ind <- ceiling(x = 1:length(x = peaks) / bins)
    max_bin <- max(bin_ind)
    
    # prepare residual  matrix
    res <- matrix(NA_real_, length(peaks), nrow(regressor_data), 
                  dimnames = list(peaks, rownames(regressor_data)))
    
    # iterate
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    for (i in 1:max_bin){
        peaks_bin <- peaks[bin_ind == i]
        if(link=="logit"){
            mu <- exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)) / 
                (1+exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }else if(link == "probit"){
            mu <- invprobit(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data))
        }else if(link =="cloglog"){
            mu <- 1-exp(-exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }
        y <- as.matrix(x[peaks_bin, , drop=FALSE])
        if(type=="deviance"){
            res[peaks_bin, ] <- deviance_residual(y, mu, 1)
        }else if(type=="pearson"){
            res[peaks_bin,] <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }else if(type=="response"){
            res[peaks_bin,] <- (y - mu)
        }else{
            res[peaks_bin,] <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # remove na
    res[is.na(res)] <- 0
    res <- res[rowSums(res == 0) != ncol(res),]
    
    # report intial residuals
    res.range <- range(res)
    message("   * residual range: ",res.range[1], " - ", res.range[2])
    
    ###########################################################################################
    # return highly variable peaks ------------------------------------------------------------
    ###########################################################################################
    if(var.th > 0){
        p.vars <- apply(res, 1, var)
        p.vars[is.na(p.vars)] <- 0
        p.aves <- Matrix::rowMeans(x)
        p.order <- order(p.vars, decreasing=T)
        p.vars <- p.vars[p.order]
        p.aves <- p.aves[p.order]
        
        # keep peaks
        keepers <- names(p.vars[p.vars>var.th])
        if(length(keepers) < topSites & topSites < nrow(res)){
            keepers <- names(p.vars[1:topSites])
        }
        res <- res[rownames(res) %in% keepers, ]
        num.keep <- nrow(res)
        ids.c <- colnames(res)
        message(" - kept top ",num.keep, " variable peaks ...")
        
        # plot peak usage vs. residual variance
        pdf(paste0(prefix,".residualVariance.pdf"), width=5, height=5)
        plot(p.aves, p.vars, pch=16, cex=0.5, col=ifelse(p.vars < var.th, "black", "red"),
             xlab="Mean", ylab="Residual variance", ylim=range(p.vars), 
             xlim=range(p.aves))
        dev.off()
        
        # var matrix
        df.vars <- data.frame(peaks=names(p.vars), nsites=p.aves, vars=p.vars)
        write.table(df.vars, file=paste0(prefix,".residualVariance.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    }
    
    # do stnd?
    if(doSTD == T & is.null(variates)){
        res <- as.matrix(t(scale(t(res), scale=T, center=T)))
    }else if(doCenter == T & is.null(variates)){
        res <- as.matrix(t(scale(t(res), scale=F, center=T)))
    }
    
    # plot residuals
    if(doPlot){
        
        # heatmap
        keepers <- names(p.vars[1:topSites])
        plotRes(res[rownames(res) %in% keepers,], prefix=paste0(prefix,".uncorrected"))
    }
    
    ###############################################################################################
    # variables to regress from residuals ---------------------------------------------------------
    ###############################################################################################
    
    if(!is.null(variates) && method=="lm"){
        
        # verbose
        message(" - removing effects from confounding variables ...")
        
        # model
        form.nr <- as.formula(paste0("z",variates))
        message("   * formula: ", form.nr)
        
        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)
        
        # control looping/chunks
        chunks <- 1000
        peaknames <- rownames(x)
        niter <- nrow(x)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }
        
        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        
        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it) 
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }
        
        # calc first QR
        y1 <- yo[colnames(res),]
        regression.mat <- cbind(y1, res[1,])
        colnames(regression.mat) <- c(colnames(y1),"z")
        qr <- lm(form.nr, data = regression.mat, qr = TRUE)$qr
        rm(regression.mat)
        r.names <- rownames(res)
        
        # run logistic regression for each peak independently - in parallel
        res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts) %dopar% {
            
            # run chunks
            its <- c(n$i:(n$i+n$m-1))
            dfout <- lapply(its, function(j){
                
                # run GLM
                val <- qr.resid(qr = qr, y = res[j,])
                if(doSTD==T){
                    val <- (val - mean(val, na.rm=T))/sd(val, na.rm=T)
                }else if(doCenter==T){
                    val <- (val - mean(val, na.rm=T))
                }
                return(val)
                
            })
            
            # reformat
            dfout <- do.call(rbind, dfout)
            rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
            return(dfout)
            
        }
        
        # close connections
        close(pb)
        stopCluster(cl)
        
        # clean
        res <- res2
        rm(res2)
        res <- res[r.names,]
        
        # plot
        if(doPlot){
            res.range <- range(res)
            message("   * residual range: ",res.range[1], " - ", res.range[2])
            
            # res distribution
            pdf(paste0(prefix,".corrected.res.dist.pdf"), width=5, height=5)
            hist(as.numeric(res))
            dev.off
            
            # heatmap
            plotRes(res, prefix=paste0(prefix,".corrected"))
        }
    }else if(!is.null(variates) & method=="elasticNet"){
        
        # verbose
        message(" - removing effects from confounding variables with elasticNet ...")
        
        # prep variates
        y1 <- yo[colnames(res),]
        message("   * number of cells in meta data = ", nrow(y1))
        ind.vars <- gsub("~","",variates)
        ind.vars <- gsub(" ","",ind.vars)
        ind.vars <- strsplit(ind.vars, "\\+")
        ind.vars <- do.call(c, ind.vars)
        n.vars <- c()
        f.vars <- c()
        for(i in ind.vars){
            if(is.numeric(y1[,i])){
                n.vars <- c(n.vars, i)
            }else if(is.factor(y1[,i])){
                f.vars <- c(f.vars, i)
            }else if(is.character(y1[,i])){
                f.vars <- c(f.vars, i)
                y1[,i] <- factor(y1[,i])
            }
        }
        sqr <- paste("I(",n.vars,"^0.5)",sep="")
        sqr <- paste(sqr, collapse=" + ")
        sq <- paste("I(",n.vars,"^2)",sep="")
        sq <- paste(sqr, collapse=" + ")
        if(length(n.vars) > 1){
            combs <- combn(ind.vars, 2)
            combs <- paste(combs[1,], combs[2,], sep=":")
            combs <- paste(combs, collapse=" + ")
            final <- as.formula(paste(variates, "+", combs, "+", sqr, "+", sq, sep=" "))
        }else if(is.null(n.vars) & length(f.vars) > 0){
            combs <- combn(ind.vars, 2)
            combs <- paste(combs[1,], combs[2,], sep=":")
            combs <- paste(combs, collapse=" + ")
            final <- as.formula(paste(variates,"+",combs,sep=" "))
        }else{
            final <- as.formula(variates)
            #final <- as.formula(paste(variates, "+", sqr, "+", sq, sep=" "))
        }
        
        # prep train/full
        if(length(f.vars) > 1){
            message("   * factor variables > 1")
            xvars <- sparse.model.matrix(final, y1, 
                                         contrasts.arg=lapply(y1[,colnames(y1) %in% f.vars], contrasts, contrasts=F))
        }else if(length(f.vars)==1){
            message("   * factor variables = 1")
            f.list <- list(f.vars=contrasts(y1[,f.vars], contrasts=F))
            names(f.list) <- f.vars
            xvars <- sparse.model.matrix(final, y1, contrasts.arg=f.list)
        }else{
            message("   * factor variables = 0")
            xvars <- sparse.model.matrix(final, y1)
        }
        
        # save peak ids
        r.names <- rownames(res)
        
        # get initial estimates of lambda
        message("   * getting regularized value of lambda from 1000 random peaks ...")
        lambda <- c()
        for(i in 1:1000){
            rand <- sample(nrow(res), 1)
            yvar <- res[rand,]
            cv.lasso <- cv.glmnet(xvars, yvar, alpha=alpha, family="gaussian")
            lambda <- c(lambda, cv.lasso$lambda.min)
        }
        ave.lambda <- median(lambda[lambda > 0], na.rm=T)
        message("   * average lambda = ", ave.lambda)
        
        # plot lambda permutation
        pdf(paste0(prefix,".lambdaPerm.pdf"), width=5, height=5)
        plot(density(lambda))
        abline(v=ave.lambda, col="red")
        dev.off()
        
        # regularize glmnet function
        regularizedR  <- function(y, xvars, ave.lambda){
            
            # elastic-net
            #yvar <- y[split]
            fit <- glmnet(xvars, y, family="gaussian", alpha=alpha, lambda=ave.lambda)
            return(as.numeric(y-predict(fit, s=ave.lambda, newx=xvars, type="response")[,1]))
        }
        
        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)
        
        # control looping/chunks
        chunks <- bins
        peaknames <- rownames(res)
        niter <- nrow(res)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }
        
        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        package.list <- c("glmnet")
        
        # output raw residuals
        res.out <- paste0(prefix, ".rawresiduals.txt")
        if(file.exists(res.out)){
            unlink(res.out)
        }
        fakedf <- matrix(0, ncol=ncol(res), nrow=0)
        colnames(fakedf) <- colnames(res)
        suppressWarnings(write.table(fakedf, file=res.out, append=T, sep="\t", col.names=T, row.names=F, quote=F))
        
        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it) 
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }
        
        # run logistic regression for each peak independently - in parallel
        message("   * fitting elasticNet regressions and extracting residuals ...")
        res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts,
                        .packages=package.list) %dopar% {
                            
                            # run chunks
                            its <- c(n$i:(n$i+n$m-1))
                            dfout <- lapply(its, function(j){
                                
                                # run GLMNET
                                return(regularizedR(res[j,], xvars, ave.lambda))
                                
                            })
                            
                            
                            # reformat
                            dfout <- do.call(rbind, dfout)
                            
                            # write raw residuals to disk
                            if(writeRawRes){
                                rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
                                colnames(dfout) <- colnames(res)
                                
                                # write raw residuals to disk
                                write.table(dfout, file=res.out, append=T, sep="\t", col.names=F, row.names=T, quote=F)
                            }
                            
                            # standarize
                            if(doSTD){
                                dfout <- t(apply(t(dfout), 2, function(w){
                                    as.numeric(scale(w))
                                }))
                            }else if(doCenter){
                                dfout <- t(apply(t(dfout), 2, function(w){
                                    as.numeric(scale(w, center=T, scale=F))
                                }))
                            }
                            rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
                            colnames(dfout) <- colnames(res)
                            
                            # return
                            return(dfout)
                            
                        }
        
        # close connections
        close(pb)
        stopCluster(cl)
        
        # clean
        res <- res2
        rm(res2)
        res <- res[r.names,]
        
        # plot
        if(doPlot){
            
            # heatmap
            plotRes(res, prefix=paste0(prefix,".corrected"))
        }
    }
    
    detach("package:doSNOW", unload=TRUE)
    
    # return
    return(res)
}
regModel2     <- function(x, y, r.variates='~log10nSites', subpeaks=5000, subcells=1000, propC=0.1,
                          cellNames=NULL, bins=256, bw_adjust=10, type="response", verbose=T, 
                          plotCoef=T, link="logit", nthreads=1, doPlot=F, weights=NULL, 
                          prefix="out", var.th=1, topSites=50000, do.sample=F, stdDev=F){
    
    library(doSNOW)
    
    # fix type if necessary
    if(type!="deviance" & type !="pearson" & type!="response" & type !="working"){
        message(" - *type* incorrect, setting residual to default (pearson) ...")
        type <- "pearson"
    }
    
    # functions
    is_outlier          <- function(y, x, th = 10) {
        bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
        eps <- .Machine$double.eps * 10
        breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
        breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
        score1 <- robust_scale_binned(y, x, breaks1)
        score2 <- robust_scale_binned(y, x, breaks2)
        return(pmin(abs(score1), abs(score2)) > th)
    }
    robust_scale_binned <- function(y, x, breaks) {
        bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
        tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
        score <- rep(0, length(x))
        o <- order(bins)
        if (inherits(x = tmp$x, what = 'list')) {
            score[o] <- unlist(tmp$x)
        } else {
            score[o] <- as.numeric(t(tmp$x))
        }
        return(score)
    }
    deviance_residual   <- function(y, mu, wt){
        d.res <- sqrt(pmax((binomial()$dev.resid)(y, mu, wt),0))
        d.res <- ifelse(y > mu, d.res, -d.res)
        return(d.res)
    }
    pearson_residual    <- function(y, mu, theta){(y-mu)/sqrt((mu*(1-mu))*theta)}
    robust_scale        <- function(x){return((x - median(x)) / (mad(x) + .Machine$double.eps))}
    ChunkPoints         <- function(dsize, csize) {
        return(vapply(
            X = 1L:ceiling(x = dsize / csize),
            FUN = function(i) {
                return(c(
                    start = (csize * (i - 1L)) + 1L,
                    end = min(csize * i, dsize)
                ))
            },
            FUN.VALUE = numeric(length = 2L)
        ))
    }
    invprobit           <- function(x){
        thresh <- -qnorm(.Machine$double.eps)
        x <- pmin(pmax(x, -thresh), thresh)
        pnorm(x)
    }
    .sampleCells        <- function(x, min.cells.per=1000, propLib=0.05){
        
        # number of cells 
        x$library <- as.character(x$library)
        libs <- unique(x$library)
        out <- lapply(libs, function(z){
            
            # subset by library
            y <- x[x$library == z,]
            
            # num cells
            num.cells <- ceiling(nrow(y)*propLib)
            if(num.cells < min.cells.per){
                num.cells <- min.cells.per
            }
            
            # sample
            rownames(y)[sample(seq(1:nrow(y)), size=num.cells)]
        })
        out <- do.call(c, out)
        message("   * sampled a total of ",length(out)," cells ...")
        return(out)
        
    }
    stddev              <- function(r.dev, bins=256, threads=1){
        
        # verbose
        message(" - standardizing deviance scores ...")
        
        # set up bins
        bin_ind <- ceiling(x = 1:nrow(r.dev) / bins)
        max_bin <- max(bin_ind)
        ids <- rownames(r.dev)
        
        # run in bins
        dev <- lapply(seq(1:max_bin), function(x){
            peaks_bin <- rownames(r.dev)[bin_ind == x]
            t(as.matrix(scale(t(r.dev[peaks_bin,]))))
        })#, mc.cores=threads)
        rm(r.dev)
        
        # merge and return
        dev <- do.call(rbind, dev)
        dev <- dev[ids,]
        
        return(dev)
    }
    
    
    ###############################################################################################
    # start ---------------------------------------------------------------------------------------
    ###############################################################################################
    message(" - regularizing logistic model parameters ...")
    
    # set-up data
    yo <- y
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    form <- as.formula(paste0("z",r.variates))
    message("   * formula: ", form)
    
    # sample cells
    if(do.sample){
        message("   * sampling cells ...")
        sub.cells <- .sampleCells(y, min.cells.per=subcells, propLib=propC)
        x.sub1 <- x[,sub.cells]
        x.sub1 <- x.sub1[Matrix::rowSums(x.sub1)>0,]
        x.sub1 <- x.sub1[,Matrix::colSums(x.sub1)>0]
        y.sub1 <- y[colnames(x.sub1),]
        
        #select subsample of peaks
        set.seed(1111)
        con <- 1
        if(verbose){message("   * estimating geometric mean ...")}
        log_peak_mean <- log(Matrix::rowMeans(x))
        log_peak_mean_red <- log(Matrix::rowMeans(x.sub1))
        r.lpm <- range(log_peak_mean)
        log_peak_mean_red <- log_peak_mean_red[log_peak_mean_red > r.lpm[1] & log_peak_mean_red < r.lpm[2]]
        
        if(verbose){message("   * density sampling on peak space ...")}
        log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
        sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean_red)$y + .Machine$double.eps)
        peak.sub <- sample(x=names(log_peak_mean_red), size=subpeaks, prob=sample_prob)
        x.sub <- x.sub1[peak.sub,]
        log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))
        
        # fit model to subset
        message("   * fitting parameters to ",nrow(x.sub)," peaks ...")
        pars <- lapply(seq(1:nrow(x.sub)), function(j){
            df <- cbind(y.sub1, x.sub[j,])
            colnames(df) <- c(colnames(y), "z")
            suppressWarnings(mod <- glm(form, 
                                        data=df, 
                                        family=quasibinomial(link = link)))
            theta <- summary(mod)$dispersion
            names(theta) <- "theta"
            return(c(theta, mod$coefficients))
        })
        
    }else{
        
        #select subsample of peaks
        set.seed(1111)
        con <- 1
        if(verbose){message("   * estimating geometric mean ...")}
        log_peak_mean <- log(Matrix::rowMeans(x))
        
        if(verbose){message("   * density sampling on peak space ...")}
        log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
        sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean)$y + .Machine$double.eps)
        peak.sub <- sample(x=rownames(x), size=subpeaks, prob=sample_prob)
        x.sub <- x[peak.sub,]
        log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))
        
        # fit model to subset
        message("   * fitting parameters to ",nrow(x.sub)," peaks ...")
        pars <- lapply(seq(1:nrow(x.sub)), function(j){
            df <- cbind(y, x.sub[j,])
            colnames(df) <- c(colnames(y), "z")
            suppressWarnings(mod <- glm(form, 
                                        data=df, 
                                        family=quasibinomial(link = link)))
            theta <- summary(mod)$dispersion
            names(theta) <- "theta"
            return(c(theta, mod$coefficients))
        })
        
    }
    
    # merge
    pars <- do.call(rbind, pars)
    rownames(pars) <- rownames(x.sub)
    
    
    ###############################################################################################
    ### regularize --------------------------------------------------------------------------------
    ###############################################################################################
    
    # outliers
    message("   * finding outliers ...")
    peaks <- names(log_peak_mean)
    outliers <- apply(pars, 2, function(k) is_outlier(k, log_peak_mean_sub))
    outliers <- apply(outliers, 1, any)
    if (sum(outliers) > 0){
        message("   * found ", sum(outliers), " outliers ...")
        pars <- pars[!outliers, ]
        peaks_step1 <- rownames(pars)
        log_peak_mean_sub <- log_peak_mean_sub[!outliers]
    }
    
    # select bw
    bw <- bw.SJ(log_peak_mean) * bw_adjust
    
    # parameter for predictions
    x_points <- pmax(log_peak_mean, min(log_peak_mean))
    x_points <- pmin(x_points, max(log_peak_mean))
    
    # take results from step 1 and fit/predict parameters to all genes
    message("   * regularizing all coefficients ...")
    o <- order(x_points)
    pars_fit <- matrix(NA_real_, length(peaks), ncol(pars),
                       dimnames = list(peaks, colnames(pars)))
    
    # global fit / regularization for all coefficients
    for (i in 1:ncol(pars)){
        pars_fit[o, i] <- ksmooth(x = log_peak_mean_sub, y = pars[, i],
                                  x.points = x_points, bandwidth = bw, kernel='normal')$y
    }
    
    # plot coefficients
    if(plotCoef){
        
        # verbose
        message("   * plotting coefficients before and after regularization ...")
        
        # set up plot
        n.cols <- ncol(pars_fit)
        n.rows <- 1
        pdf(paste0(prefix,".coefficient.pdf"), width=(5*n.cols), height=(n.rows*5))
        layout(matrix(c(1:(n.cols*n.rows)), nrow=n.rows, byrow=T))
        raw.coeff <- as.data.frame(pars)
        raw.coeff$logpeakmeans <- log_peak_mean_sub
        write.table(raw.coeff, file=paste0(prefix,".coefficients.sub.txt"), quote=F, row.names=T, col.names=T, sep="\t")
        
        # iterate over coefficients
        for(i in 1:ncol(pars)){
            plot(log_peak_mean_sub, pars[,i], pch=19, cex=0.7, col=alpha("black", 0.5), 
                 main=colnames(pars)[i])
            df <- data.frame(x=x_points, y=pars_fit[,i])
            df <- df[order(df$x, decreasing=F),]
            lines(df$x, df$y, col=alpha("red", 0.5), lwd=2)
        }
        
        # device off
        dev.off()
        
    }
    
    
    ###############################################################################################
    # fit data ------------------------------------------------------------------------------------
    ###############################################################################################
    
    # affine transform
    message("   * estimating deviance residuals on the full data set ...")
    regressor_data <- model.matrix(as.formula(r.variates), data=y)
    
    # split peaks into bins
    bin_ind <- ceiling(x = 1:length(x = peaks) / bins)
    max_bin <- max(bin_ind)
    
    # iterate
    res <- mclapply(seq(1:max_bin), function(i){
        peaks_bin <- peaks[bin_ind == i]
        if(link=="logit"){
            mu <- exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)) / 
                (1+exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }else if(link == "probit"){
            mu <- invprobit(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data))
        }else if(link =="cloglog"){
            mu <- 1-exp(-exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }
        y <- as.matrix(x[peaks_bin, , drop=FALSE])
        if(type=="deviance"){
            out.1 <- deviance_residual(y, mu, 1)
        }else if(type=="pearson"){
            out.1 <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }else if(type=="response"){
            out.1 <- (y - mu)
        }else{
            out.1 <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }
        colnames(out.1) <- rownames(regressor_data)
        rownames(out.1) <- peaks_bin
        return(out.1)
    }, mc.cores=nthreads)
    
    # merge
    res <- do.call(rbind, res)
    res <- res[peaks, ]
    
    # remove na
    res[is.na(res)] <- 0
    res <- res[rowSums(res)!=0,]
    
    # report intial residuals
    res.range <- range(res)
    message("   * residual range: ",res.range[1], " - ", res.range[2])
    
    ###########################################################################################
    # return highly variable peaks ------------------------------------------------------------
    ###########################################################################################
    if(var.th > 0){
        p.vars <- apply(res, 1, var)
        p.vars[is.na(p.vars)] <- 0
        p.aves <- Matrix::rowMeans(x)
        p.order <- order(p.vars, decreasing=T)
        p.vars <- p.vars[p.order]
        p.aves <- p.aves[p.order]
        
        # keep peaks
        keepers <- names(p.vars[p.vars>var.th])
        if(length(keepers) < topSites & topSites < nrow(res)){
            keepers <- names(p.vars[1:topSites])
        }
        res <- res[rownames(res) %in% keepers, ]
        num.keep <- nrow(res)
        ids.c <- colnames(res)
        message(" - kept top ",num.keep, " variable peaks ...")
        
        # plot peak usage vs. residual variance
        pdf(paste0(prefix,".residualVariance.pdf"), width=5, height=5)
        plot(p.aves, p.vars, pch=16, cex=0.5, col=ifelse(p.vars < var.th, "black", "red"),
             xlab="Mean", ylab="Residual variance", ylim=range(p.vars), 
             xlim=range(p.aves))
        dev.off()
        
        # var matrix
        df.vars <- data.frame(peaks=names(p.vars), nsites=p.aves, vars=p.vars)
        write.table(df.vars, file=paste0(prefix,".residualVariance.txt"), quote=F, row.names=F, col.names=F, sep="\t")
    }
    
    # detach doSNOW
    detach("package:doSNOW", unload=TRUE)
    
    # if stddev
    if(stdDev==T){
        res <- stddev(res)
    }
    
    # return
    return(res)
}
runModel      <- function(a, b, variates, r.variates, modtype="regModel", chunks=1000, 
                          type="response", subpeaks=1000, nthreads=1, weights=NULL, doPlot=F, 
                          prefix="out", stdize=F, doCenter=F, link="logit", method="elasticNet", 
                          alpha=0.5, var.th=1, ...){
    
    # choose model
    if(modtype=="regModel"){
        regModel(a,b,variates=variates,r.variates=r.variates,
                 subpeaks=subpeaks,
                 type=type, 
                 weights=weights,
                 doPlot=doPlot, 
                 prefix=prefix, 
                 doSTD=stdize, 
                 doCenter=doCenter,
                 link=link,
                 var.th=var.th,
                 ...)
    }else if(modtype=="regModel2"){
        regModel2(a,b,r.variates=r.variates,
                  subcells=1000, propC=0.1,
                  do.sample=T,
                  stdDev=F,
                  subpeaks=subpeaks,
                  type=type, 
                  weights=weights,
                  doPlot=doPlot, 
                  prefix=prefix, 
                  link=link,
                  var.th=var.th,
                  ...)
    }else if(modtype=="logisticModel"){
        logisticModel(a,b,variates=variates,r.variates=r.variates,
                      nthreads=nthreads,
                      subpeaks=subpeaks,
                      chunks=chunks, 
                      doPlot=doPlot,
                      type=type,
                      prefix=prefix,
                      link=link,
                      doSTD=stdize, 
                      doCenter=doCenter,
                      var.th=var.th,
                      ...)
    }else if(modtype=="tfidf"){
        message(" - normalizing peak x matrix with TFIDF ...")
        tfidf(a, frequencies=T, log_scale_tf=T, scale_factor=1e5)
    }else{
        regModel(a,b,variates=variates,r.variates=r.variates,
                 subpeaks=subpeaks,
                 type=type,
                 weights=weights,
                 doPlot=doPlot, 
                 prefix=prefix, 
                 doSTD=stdize, 
                 link=link, 
                 var.th=var.th,
                 ...)
    }
    
}

# recluster ---------------------------------------------------------------------------------------
reclust       <- function(a, b, variates, r.variates, modtype="regModel", chunks=1000, 
                          type="response", subpeaks=1000, ncpus=1, doPlot=F, prefix="out", 
                          stdize=T, doCenter=F, link="logit", method="elasticNet", alpha=0.5, 
                          var.th=1, scaleVar=T, doL2=1, cap=F, useTop=NULL, n.pcs=50,
                          batch=NULL, theta=1, lambda=1, useTop.th=0.75, min.reads=5e5,
                          m.clst=25, threshold= -0.25, k.nearC=15, res=0.03, clustOB="umap", 
                          key="UMAP", ...){
    
    # get clusters to iterate over
    b$LouvainClusters <- factor(b$LouvainClusters)
    LC <- mixedsort(unique(b$LouvainClusters))
    
    # iterate
    finalresults <- mclapply(LC, function(i){
        
        # update
        prefix1 <- paste0(prefix,"_",i)
        clustname <- "LouvainCluster_sub"
        
        # verbose
        message(" ------------ ", i, " ------------")
        
        # clust
        clust <- rownames(subset(b, b$LouvainClusters==i))
        aa <- a[,clust]
        bb <- b[colnames(aa),]
        aa <- aa[Matrix::rowSums(aa)>=(ncol(aa)*0.005),]
        aa <- aa[,Matrix::colSums(aa)>50]
        aa <- aa[Matrix::rowSums(aa)>0,]
        aa <- aa[,Matrix::colSums(aa)>0]
        
        if(nrow(bb) < 200 | sum(bb$unique) < 2e6){
            message(" - skipping Louvain cluster ", i, ": too few cells or reads")
            bb$umapsub_1 <- bb$umap1
            bb$umapsub_2 <- bb$umap2
            bb[,c(clustname)] <- factor(bb$LouvainClusters)
            return(bb)
        }
        
        # rerun pipeline
        dev <- runModel(aa, bb, variates, r.variates,
                        subpeaks = subpeaks,
                        chunks = 256,
                        type = type,
                        modtype = modtype,
                        nthreads = ncpus,
                        prefix = prefix1,
                        doPlot = doPlot,
                        stdize = stdize,
                        link = link,
                        doCenter = doCenter,
                        method = method,
                        alpha = alpha,
                        var.th = var.th)
        
        # std dev
        dev <- stddev(dev)
        
        # reduce dims
        out.pcs  <- reduceDims(dev, bb, 
                               n.pcs = n.pcs,
                               batch = batch, 
                               theta = theta, 
                               lambda = lambda, 
                               doL2 = doL2, 
                               cap = cap,
                               prefix = prefix1,
                               scaleVar = scaleVar,
                               raw = aa, 
                               center = doCenter,
                               useTop = useTop,
                               useTop.th = useTop.th)
        
        # project with UMAP
        out.umap <- projectUMAP(out.pcs)
        out.umap$umapsub_1 <- out.umap$umap1
        out.umap$umapsub_2 <- out.umap$umap2
        out.umap$umap1 <- NULL
        out.umap$umap2 <- NULL
        
        # call clusters
        bb <- callClusters(aa, bb, out.pcs, out.umap, 
                           umap1="umapsub_1", umap2="umapsub_2",
                           min.reads = min.reads, 
                           m.clst = m.clst, 
                           threshold = threshold,
                           k.near = k.nearC,
                           res = res,
                           prefix = prefix1,
                           clustOB = clustOB,
                           cname = clustname,
                           key = "umapsub_")
        
        # make column factor
        bb[,c(clustname)] <- factor(bb[,c(clustname)])
        bb$sub_louvain <- paste(bb$LouvainClusters, bb[,c(clustname)], sep="_")
        
        # write stats to shell
        reportResults(bb, column=clustname)
        
        # output
        outData(out.pcs, bb, prefix=prefix1, dev=NULL)
        
        # plot UMAP
        plotUMAP(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1, column=clustname)
        plotSTATS2(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1)
        
        # update list
        finalresults[[i]] <- bb
        
        # check number of sub-clusters
        if(length(unique(bb$sub_louvain))>1){
            
            # load marker/gene activity
            bb$umap1 <- bb$umapsub_1
            bb$umap2 <- bb$umapsub_2
            datt <- loadMData(bb, out.pcs, geneact, mark, clustID="sub_louvain")
            bdat <- datt$b
            activity <- datt$activity
            h.pcs <- datt$h.pcs
            marker.info <- datt$marker.info
            print(head(bdat))
            print(head(activity[,1:5]))
            print(head(h.pcs))
            
            # normalize per cell activity by cluster average and size factors
            results <- normalize.activity(bdat, activity, output=prefix1, logTransform=F, scaleP=F)
            activity <- results$norm.act
            row.o <- results$row.o
            
            # impute gene accessibility scores
            impute.activity <- smooth.data(activity, k=15, step=3, npcs=s.pcs, df=NULL,
                                           rds=h.pcs, cleanExp=F, output=prefix1)
            
            # collect marker accessibility
            plot.act.scores(bdat, acts=activity, 
                            info=marker.info, 
                            logT=T,
                            outname=paste0(prefix1,".normalized.known.Markers.pdf"))
            
            plot.act.scores(bdat, acts=impute.activity, 
                            info=marker.info, 
                            logT=F,
                            outname=paste0(prefix1,".impute.known.Markers.pdf"))
            
        }
    }, mc.cores=nthreads)
    
    # condense list
    df <- do.call(rbind, finalresults)
    rownames(df) <- df$cellID
    df$sub_louvain <- paste(df$LouvainClusters, df$LouvainCluster_sub, sep="_")
    df$tissue_cluster <- as.numeric(factor(df$sub_louvain, levels=mixedsort(unique(as.character(df$sub_louvain)))))
    return(df)
    
}
reclust2      <- function(a, b, pcs, prefix="out", n.pcs=50, batch=NULL, theta=1, nthreads=1,
                          lambda=1, min.reads=5e5, m.clst=25, threshold3=0.25, doL2=F,
                          clustType="densityClust", k.nearC=15, res=0.03, clustOB="umap", 
                          cleanCluster=T, e.thresh=2, doSTDize=F){
    
    # get clusters to iterate over
    LC <- mixedsort(unique(as.character(b$LouvainClusters)))
    
    # load data for marker analysis
    mark <- "/scratch/apm25309/single_cell/ATACseq/v3/step1_clustering/model_based/markers.bed"
    geneact <- "/scratch/apm25309/single_cell/ATACseq/v3/sparse/genes/gene_activity/all.ex.GBaccessibility.sparse"
    
    # resolutions
    #res <- rep(res, length(LC))
    message(" - initializing sub-clustering ... ")
    print(LC)
    
    # adjust number of threads
    if(length(LC) < nthreads){
        nthreads <- length(LC)
    }
    
    # iterate
    it <- 0
    finalresults <- mclapply(LC, function(i){
        
        # update
        prefix1 <- paste0(prefix,"_",i)
        clustname <- "LouvainCluster_sub"
        
        # verbose
        message("---------------------------------------------")
        message("------------------ ", i, " ------------------")
        message("---------------------------------------------")
        
        # clust
        clust <- rownames(subset(b, as.character(b$LouvainClusters)==i))
        aa <- a[,colnames(a) %in% clust]
        bb <- b[colnames(aa),]
        
        # use minimum among the number of available PCs and requested components
        if(ncol(pcs) < n.pcs){
            n.pcs <- ncol(pcs)
        }
        out.pcs1 <- pcs[colnames(aa),1:n.pcs]
        
        # standarize embeddings
        if(doSTDize == T){
            out.pcs1 <- t(apply(out.pcs1, 1, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)))
        }
        
        # L2 norm embeddings
        if(doL2 == T){
            out.pcs1 <- t(apply(out.pcs1, 1, function(x) x/sqrt(sum(x^2))))
        }
        
        # verbose
        message(" - loaded ", nrow(out.pcs1), " cells ...")
        
        # do harmony batch correciton
        if(!is.null(batch)){
            out.pcs1 <- HarmonyMatrix(out.pcs1, bb, do_pca=F, vars_use=batch,
                                      theta=2, lambda=1, nclust=100, tau=0, 
                                      sigma=0.1, block.size=0.01)
        }
        
        # project with UMAP
        umap.met <- "euclidean"
        out.umap <- as.data.frame(projectUMAP(out.pcs1, m.dist=0.1, k.near=k.nearC, metric=umap.met))
        colnames(out.umap) <- c("umap1","umap2")
        message(" - finished UMAP ...")
        out.umap$umapsub_1 <- out.umap$umap1
        out.umap$umapsub_2 <- out.umap$umap2
        out.umap$umap1 <- NULL
        out.umap$umap2 <- NULL
        
        # call clusters
        message(" - begin clustering ...")
        bb <- callClusters(aa, bb, out.pcs1, out.umap, 
                           umap1="umapsub_1", umap2="umapsub_2",
                           min.reads = min.reads, 
                           m.clst = m.clst, 
                           threshold3 = threshold3,
                           k.near = k.nearC,
                           res = res,
                           prefix = prefix1,
                           clustOB = clustOB,
                           cname = clustname,
                           clustType = clustType,
                           dynamic = F,
                           e.thresh = e.thresh,
                           cleanCluster = cleanCluster,
                           cl.method = 2)
        
        # make column factor
        bb$sub_louvain <- factor(paste(bb$LouvainClusters, bb[,c(clustname)], sep="_"))
        
        # write stats to shell
        reportResults(bb, column=clustname)
        
        # plot UMAP
        bb$LouvainClusters <- factor(bb$LouvainClusters)
        bb[,clustname] <- factor(bb[,clustname])
        plotUMAP(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1, column=clustname)
        plotSTATS2(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1)

        # update list
        dd <- bb
        message(" - clusters contain ", nrow(bb), " cells ... (out of = ",nrow(out.pcs1),")")
        
        # check number of sub-clusters
        if(length(unique(bb$sub_louvain))>1){
            
            # load marker/gene activity
            bb$umap1 <- bb$umapsub_1
            bb$umap2 <- bb$umapsub_2
            datt <- loadMData(bb, out.pcs1, geneact, mark, clustID="sub_louvain")
            bdat <- datt$b
            activity1 <- datt$activity
            h.pcs1 <- datt$h.pcs
            marker.info1 <- datt$marker.info
            print(head(bdat))
            print(head(activity1[,1:5]))
            print(head(h.pcs1))
            
            # normalize per cell activity by cluster average and size factors
            results <- normalize.activity(bdat, activity1, output=prefix1, logTransform=F, scaleP=F)
            activity1 <- results$norm.act
            row.o <- results$row.o
            
            # impute gene accessibility scores
            impute.activity1 <- smooth.data(activity1, k=15, step=3, npcs=ncol(h.pcs1), df=NULL,
                                            rds=h.pcs1, cleanExp=F, output=prefix1)
            
            # collect marker accessibility
            plot.act.scores(bdat, acts=activity1, 
                            info=marker.info1, 
                            logT=T,
                            outname=paste0(prefix1,".normalized.known.Markers.pdf"))
            
            plot.act.scores(bdat, acts=impute.activity1, 
                            info=marker.info1, 
                            logT=F,
                            outname=paste0(prefix1,".impute.known.Markers.pdf"))
            
        }
        
        # return results
        return(dd)
        
    }, mc.cores=nthreads)
    
    # condense list
    message(" - combining child processes ... ")
    df <- as.data.frame(do.call(rbind, finalresults))
    df <- df[!duplicated(df$cellID),]
    df <- df[!is.na(df$cellID),]
    rownames(df) <- df$cellID
    df$tissue_cluster <- as.numeric(factor(df$sub_louvain, levels=mixedsort(unique(as.character(df$sub_louvain)))))
    tissue.props <- prop.table(table(df$tissue_cluster,df$tissue),1)
    tissue.calls <- apply(tissue.props, 1, function(z) names(z)[which.max(z)] )
    df$top_tissue_cluster <- paste(tissue.calls[df$tissue_cluster], names(tissue.calls[df$tissue_cluster]), sep="_")
    print(warnings())
    return(df)
    
}




