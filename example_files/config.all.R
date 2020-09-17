# config file

# process features/cells
min.p <- 0
min.t <- 5e-3
min.c <- 50

# function parameters #

# harmony
batch <- NULL #"library"
sub.batch <- NULL
theta <- 0
sub.theta <- 0
lambda <- 10
sub.lambda <- 4

# regModel
modtype <- "regModel"
subpeaks <- 5000

# reclustering
reclustType <- 2
subres <- 0.75
subthresh <- 1.5
subK <- 20
sub.mreads <- 5e5
sub.clustOB <- "svd"
sub.clustType <- "louvain"
sub.do.Clean <- T
sub.doL2 <- F
sub.stdLSI <- T
e.thresh <- 3

# cluster/SVD
threshold <- 3
res <- 0.02
clustOB <- "umap"
clustType <- "louvain"
k.nearC <- 50
doL2 <- NULL
scaleVar <- T
n.pcs <- 25
s.pcs <- 20
center <- T
stdize <- F
stdLSI <- 1
useTop <- F
var.th <- 0
do.Clean <- T

# cluster definitions
min.reads <- 1.5e6
m.clst <- 50

# model parameters
r.variates <- "~log10nSites"
type <- "pearson"
link <- "logit"
