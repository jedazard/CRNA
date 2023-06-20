##########################################################################################################################################
# E2 study
##########################################################################################################################################

#=========================================================================================#
# Set some global options
#=========================================================================================#
options(contrasts=c("contr.treatment", "contr.poly"))
options(max.print=1e6)
options(digits=12)
options(echo=TRUE)
#options(warn=-1)
#options(error=utils::recover)

#=========================================================================================#
# Loading Package Libraries
#=========================================================================================#
library("parallel")
library("limma")
library("Matrix")
library("abind")
library("igraph")
library("glmnet")
library("testit")
library("pheatmap")
library("grid")
library("gridExtra")

#library("foreach")
#library("doParallel")
library("magrittr")

#=========================================================================================#
# Loading Additional Libraries for Simulations
#=========================================================================================#
library("Biobase")
library("ROCR")
library("CausalR")

#==========================================================================================#
# Set and check the working directory
#==========================================================================================#
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/REAL", fsep = .Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

#==========================================================================================#
# Set some path shortcuts
#==========================================================================================#
CODE.path <- "CODES"
DATA.path <- "DATASETS/METHODS/CRNA/REAL/E2_Study"
RESULT.path <- "RESULTS/METHODS/CRNA/REAL"
NATM.path <- "PUBLICATIONS/PAPERS_IN_PREP/NATM_2021"

#==========================================================================================#
# Saving the workspace
#==========================================================================================#
#save.image(file=file.path(HOME.path, RESULT.path, "CRNA_REAL_E2.RData", fsep = .Platform$file.sep))

#==========================================================================================#
# Load the workspace
#==========================================================================================#
load(file=file.path(HOME.path, RESULT.path, "CRNA_REAL_E2.RData", fsep = .Platform$file.sep))

#==========================================================================================#
# Set and check the working directory
#==========================================================================================#
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/REAL", fsep = .Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

#==========================================================================================#
# Set some path shortcuts
#==========================================================================================#
CODE.path <- "CODES"
DATA.path <- "DATASETS/METHODS/CRNA/REAL/E2_Study"
RESULT.path <- "RESULTS/METHODS/CRNA/REAL"
NATM.path <- "PUBLICATIONS/PAPERS_IN_PREP/NATM_2021"

#==========================================================================================#
# Source some R procedure files
#==========================================================================================#
source(file=file.path(HOME.path, CODE.path, "R/METHODS/CRNA/SIMS/CRNA_procedures.r", fsep=.Platform$file.sep))

#==========================================================================================#
# List of global environment covariates in the workspace
#==========================================================================================#
ls(.GlobalEnv) #same as ls() or objects()

#==========================================================================================#
# Close all graphical devices
#==========================================================================================#
if (!is.null(dev.list())) graphics.off()

#==========================================================================================#
# Erase the random seed if it exists and set it up to the default one
#==========================================================================================#
if (exists(".Random.seed")) rm(.Random.seed)
RNGkind(kind="L'Ecuyer-CMRG")

#==========================================================================================#
# Restore graphical parameters upon exit
#==========================================================================================#
oldpar <- par()
oldpar$new <- FALSE
on.exit(par(oldpar))
dev.off()

#==========================================================================================#
# Load datasets
#==========================================================================================#

#==========================================================================================#
# Date and time
#==========================================================================================#
date()

#==========================================================================================#
# Functions, constants, and parameters
#==========================================================================================#
fdrmin <- 0.01
fdropt <- 0.05
fdrmax <- 0.10

#=========================================================================================#
# Cluster configuration
#=========================================================================================#
if (require("parallel")) {
  print("'parallel' is attached correctly \n")
} else {
  stop("'parallel' must be attached first \n")
}

if (.Platform$OS.type == "windows") {
  # Windows OS
  cpus <- detectCores(logical = TRUE)
  conf <- list("spec"=rep("localhost", cpus),
               "type"="SOCKET",
               "homo"=TRUE,
               "verbose"=TRUE,
               "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
} else if (.Platform$OS.type == "unix") {
  # Linux or Mac OS
  argv <- commandArgs(trailingOnly=TRUE)
  if (is.empty(argv)) {
    # Mac OS : No argument "argv" in this case
    cpus <- detectCores(logical = TRUE)
    conf <- list("spec"=rep("localhost", cpus),
                 "type"="SOCKET",
                 "homo"=TRUE,
                 "verbose"=TRUE,
                 "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
  } else {
    # Linux OS : Retrieving arguments "type" and "cpus" from the SLURM script
    type <- as.character(argv[1])
    cpus <- as.numeric(Sys.getenv("SLURM_NTASKS"))
    if (type == "SOCKET") {
      conf <- list("spec"=rep("localhost", cpus),
                   "type"="SOCKET",
                   "homo"=TRUE,
                   "verbose"=TRUE,
                   "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
    } else if (type == "MPI") {
      if (require("Rmpi")) {
        print("Rmpi is loaded correctly \n")
      } else {
        stop("Rmpi must be installed first \n")
      }
      conf <- list("spec"=cpus,
                   "type"="MPI",
                   "homo"=TRUE,
                   "verbose"=TRUE,
                   "outfile"=file.path(getwd(), "output_MPI.txt", fsep=.Platform$file.sep))
    } else {
      stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type\n")
    }
  }
} else {
  stop("Unrecognized platform: you must specify a \"windows\" or \"unix\" platform type\n")
}

cat("Cluster configuration:\n")
print(conf)

#==========================================================================================#
# Save the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_REAL_E2.RData", fsep = .Platform$file.sep))

#==========================================================================================#
# Importing Interaction data
#==========================================================================================#
interactions9606 <- OmnipathR::import_AllInteractions(select_organism = 9606)
# Downloaded 113179 interactions

#==========================================================================================#
# Importing expression data
#==========================================================================================#
# Import the expression data and convert it to a matrix
EXP <- readRDS(file.path(HOME.path, DATA.path, "GSE11352.RDS", fsep = .Platform$file.sep))
X <- Biobase::exprs(EXP)

# Import the graph data from the data frame and convert it to an 'igraph' object
PPI <- readRDS(file.path(HOME.path, DATA.path, "PPI_graph.RDS", fsep = .Platform$file.sep))
PPI.ig <- igraph::graph_from_data_frame(d=PPI, directed=TRUE) %>% simplify(edge.attr.comb = "sum")
E(PPI.ig)$sign <- sign(E(PPI.ig)$sign)

igraph::vcount(PPI.ig)  # |V|=3743 initial nodes
igraph::ecount(PPI.ig)  # |E|=12086 initial edges

#==========================================================================================#
# Expression Data
# Build expanded matrix accounting for multiple gene names
#==========================================================================================#
nodes_names <- Reduce(union, sapply(X=as.character(rownames(X)),
                                    simplify=FALSE,
                                    FUN=function(x){trimws(unlist(strsplit(x, split='///', fixed=TRUE)))}))

id <- pmatch(x=nodes_names, table=rownames(X))
w1 <- which(!is.na(id-1))
w2 <- which(is.na(id))-1
w <- intersect(w1,w2)
Y <- X[id,,drop=FALSE]
for (i in 1:length(w)) {
  j <- w[i]+1
  while ((j <= length(id)) && (is.na(id[j]))) {
    Y[j,] <- Y[w[i],,drop=TRUE]
    j <- j+1
  }
}
rownames(Y) <- nodes_names
X <- Y

dim(X)      # 21333 x 18

rm(w1, w2, w, id, i, j, Y, nodes_names)

#==========================================================================================#
# Intersection of expression matrix genes and graph nodes
# Final Expression Data matrix
# Final Knowledge Base graph
#==========================================================================================#
#----------------------- Updating final Matrix dimensions ---------------------#
nodes_names <- igraph::V(PPI.ig)$name
inter_names <- intersect(x=rownames(X), y=nodes_names)
X <- X[inter_names,]
n <- ncol(X)  # n=18
p <- nrow(X)  # p=3635

#--------------------- Updating final Graph characteristics -------------------#
PPI.ig <- delete_vertices(graph=PPI.ig, v=setdiff(x=nodes_names, y=inter_names))
e <- igraph::ecount(PPI.ig)  # e=11741 final edges
p <- igraph::vcount(PPI.ig)  # p=3635 final nodes

dist <- igraph::distances(
  graph = PPI.ig,
  v = igraph::V(PPI.ig)$name,
  to = igraph::V(PPI.ig)$name,
  mode = "out",
  weights = NULL)[1,]
dist <- dist[!is.infinite(dist)]
for (f in 1:length(fdr.set)) {
  cat("fdr: ", fdr.set[f], "\n", sep="")
  cat("d/p: ", d[[f]]/p, "\n", sep="")
}
cat("Median degree (\\nu): ", median(igraph::degree(PPI.ig)), "\n", sep="")
cat("Median path length (\\lambda): ", median(dist), "\n", sep="")
cat("e: ", e, "\n", sep="")
cat("p: ", p, "\n", sep="")
cat("e/p: ", e/p, "\n", sep="")

rm(f, inter_names, nodes_names)

fdr: 1e-07
d/p: 0.0107290233838
fdr: 1e-06
d/p: 0.0211829436039
fdr: 1e-05
d/p: 0.042090784044
Median degree (\nu): 3
Median path length (\lambda): 3.658
e: 11741
p: 3635
e/p: 3.22998624484

#==========================================================================================#
# Definition of the factors:
# Factor #1: Treatment Fator (TF)
# Factor #2: Duration Factor (DF)
#==========================================================================================#
TF <- factor(c(rep(1,3), rep(2,3),
               rep(1,3), rep(2,3),
               rep(1,3), rep(2,3)), ordered=FALSE, labels=c("C", "T"))
TF.lev <- levels(TF)
TF.ng <- nlevels(TF)
TF.def <- vector(mode="list", length=TF.ng)
TF.tab <- numeric(TF.ng)
for (g in 1:TF.ng) {
  TF.def[[g]] <- which(TF == TF.lev[g])
  TF.tab[g] <- length(TF.def[[g]])
}

DF <- factor(c(rep(1,6), 
               rep(2,6), 
               rep(3,6)), ordered=FALSE, labels=c("12h", "24h", "48h"))
DF.lev <- levels(DF)
DF.ng <- nlevels(DF)
DF.def <- vector(mode="list", length=DF.ng)
DF.tab <- numeric(DF.ng)
for (g in 1:DF.ng) {
  DF.def[[g]] <- which(DF == DF.lev[g])
  DF.tab[g] <- length(DF.def[[g]])
}

DT <- factor(c(rep(1,3), rep(2,3),
               rep(3,3), rep(4,3),
               rep(5,3), rep(6,3)), ordered=FALSE, labels=as.character(outer(X=TF.lev, Y=DF.lev, FUN=function(x,y) {paste(x,y, sep=".")})))
DT.lev <- levels(DT)
DT.ng <- nlevels(DT)
DT.def <- vector(mode="list", length=DT.ng)
DT.tab <- numeric(DT.ng)
for (g in 1:DT.ng) {
  DT.def[[g]] <- which(DT == DT.lev[g])
  DT.tab[g] <- length(DT.def[[g]])
}

rm(g)

#==========================================================================================#
# Sample Names
#==========================================================================================#
sample.names <- c(paste0(DT.lev[1], ".", DT.def[[1]]),
                  paste0(DT.lev[2], ".", DT.def[[1]]),
                  paste0(DT.lev[3], ".", DT.def[[2]]),
                  paste0(DT.lev[4], ".", DT.def[[2]]),
                  paste0(DT.lev[5], ".", DT.def[[3]]),
                  paste0(DT.lev[6], ".", DT.def[[3]]))

#==========================================================================================#
# Drop outs
#==========================================================================================#
TF <- factor(c(rep(1,3), rep(2,3),
               rep(1,3), rep(2,3)), ordered=FALSE, labels=c("C", "T"))
TF.lev <- levels(TF)
TF.ng <- nlevels(TF)
TF.def <- vector(mode="list", length=TF.ng)
TF.tab <- numeric(TF.ng)
for (g in 1:TF.ng) {
  TF.def[[g]] <- which(TF == TF.lev[g])
  TF.tab[g] <- length(TF.def[[g]])
}

DF <- factor(c(rep(2,6), rep(3,6)), ordered=FALSE, labels=c("24h", "48h"))
DF.lev <- levels(DF)
DF.ng <- nlevels(DF)
DF.def <- vector(mode="list", length=DF.ng)
DF.tab <- numeric(DF.ng)
for (g in 1:DF.ng) {
  DF.def[[g]] <- which(DF == DF.lev[g])
  DF.tab[g] <- length(DF.def[[g]])
}

DT <- factor(c(rep(3,3), rep(4,3),
               rep(5,3), rep(6,3)), ordered=FALSE, labels=as.character(outer(X=TF.lev, Y=DF.lev, FUN=function(x,y) {paste(x,y, sep=".")})))
DT.lev <- levels(DT)
DT.ng <- nlevels(DT)
DT.def <- vector(mode="list", length=DT.ng)
DT.tab <- numeric(DT.ng)
for (g in 1:DT.ng) {
  DT.def[[g]] <- which(DT == DT.lev[g])
  DT.tab[g] <- length(DT.def[[g]])
}

drop.out <- c(1:6)
sample.names <- sample.names[-drop.out]
X <- X[,-drop.out]

rm(g)

#==========================================================================================#
# Updating final Matrix dimensions
#==========================================================================================#
n <- ncol(X)  # n=12 samples
p <- nrow(X)  # p=3635 nodes/variables

#==========================================================================================#
# Renaming of expression matrix sample names after definition of groups and drop outs  
#==========================================================================================#
colnames(X) <- sample.names

#==========================================================================================#
# Initializations - Constants - Parameters
#==========================================================================================#
seed <- 777
B <- 50
R <- 100
P <- 100
FR <- 0.1
FC <- 1

alpha.set <- seq(0.10, 0.90, by=0.10)
size.set <- ceiling(seq(from=p/1000, to=p, by=5))
quant.set <- rev(1-seq(from=0.001, to=0.1, by=0.0015))
fdr.set <- 10^seq(-7, -5, by=1)

#==========================================================================================#
# Estimating Simulated DEGs for a range of DE FDR values
#==========================================================================================#
degfdr.obj <- deg.fdr(fit=deg(express=X, factor=TF), fdr=fdr.set)
DEGFDR <- vector(mode="list", length=length(fdr.set))
DEGTOPTABLE <- vector(mode="list", length=length(fdr.set))
DEGs <- vector(mode="list", length=length(fdr.set))
d <- vector(mode="list", length=length(fdr.set))
for (f in 1:length(fdr.set)) {
  cat("DE FDR: ", fdr.set[f], "\n")
  names(degfdr.obj$d[[f]]) <- paste("DE FDR=", format(x=fdr.set[f], digits=3, scientific=TRUE), sep="")
  DEGFDR[[f]] <- fdr.set[f]
  DEGTOPTABLE[[f]] <- degfdr.obj$degstoptable[[f]]
  DEGs[[f]] <- degfdr.obj$degs[[f]]
  d[[f]] <- degfdr.obj$d[[f]]
  cat("DEGs: ", DEGs[[f]], "\n")
  cat("d: ", d[[f]], "\n")
}

rm(f, degfdr.obj)

#==========================================================================================#
# Tuning of RWR parameters q (model size) and alpha (restart probability)
# For fixed DE FDR
#==========================================================================================#

#==============================================================================#
# Computation of OOB PE as a function of RWR parameters for fixed DE FDR
crna.fit.sl <- crna.tuning.sl(graph=PPI.ig,
                              express=X,
                              factors=list(TF),
                              B=B,
                              alpha=alpha.set,
                              size=size.set,
                              fdr=fdr.set,
                              lag=2, span=0.40, degree=2, family="gaussian",
                              conf=conf,
                              parallel=TRUE,
                              seed=seed)
PE.mu.sl <- crna.fit.sl$PE.mu
PE.se.sl <- crna.fit.sl$PE.se
alpha.cv.sl <- crna.fit.sl$alpha.cv
size.cv.sl <- crna.fit.sl$size.cv

#==============================================================================#
# RWR parameters minimizing OOB PE for each DE FDR
f <- 1   # fixed DE FDR = 1e-7
q <- 1   # fixed size qm of maximum variance

amin <- vector(mode="list", length=length(fdr.set))
amin.id <- vector(mode="list", length=length(fdr.set))
a1se <- vector(mode="list", length=length(fdr.set))
a1se.id <- vector(mode="list", length=length(fdr.set))
qmin <- vector(mode="list", length=length(fdr.set))
qmin.id <- vector(mode="list", length=length(fdr.set))
q1se <- vector(mode="list", length=length(fdr.set))
q1se.id <- vector(mode="list", length=length(fdr.set))
for (f in 1:length(fdr.set)) {
  cat("DE FDR: ", fdr.set[f], "\n")
  amin[[f]] <- alpha.cv$amin[q,f]
  amin.id[[f]] <- alpha.cv$amin.id[q,f]
  a1se[[f]] <- alpha.cv$a1se[q,f]
  a1se.id[[f]] <- alpha.cv$a1se.id[q,f]
  qmin[[f]] <- size.cv$qmin[amin.id[[f]],f]
  qmin.id[[f]] <- size.cv$qmin.id[amin.id[[f]],f]
  q1se[[f]] <- size.cv$q1se[amin.id[[f]],f]
  q1se.id[[f]] <- size.cv$q1se.id[amin.id[[f]],f]
}

# amin : 0.9
# amin.id : 9
# a1se : 0.6
# a1se.id : 6

#==============================================================================#
# PE curve as a function of restart probability alpha, for fixed size and DE FDR
# Selection of optimal parameters minimizing the PE surface, for fixed DE FDR
f <- 1   # fixed DE FDR = 1e-7
q <- 1   # fixed size qm of maximum variance

my.file <- paste("/CV PE profile as a function of restart probability (alpha), for fixed size (q=", size.set[q], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - E2 Study.ps", sep="")
postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=6)

plot.profile.sl(pe.mu=PE.mu.sl[[f]], 
                pe.se=PE.se.sl[[f]], 
                alpha=NULL,
                size=size.set[q], 
                fdr=fdr.set[f], 
                onese=FALSE, 
                xlab=expression(alpha), 
                ylab="PE",
                main=paste("CV PE profile as a function of restart probability (alpha), \n for fixed size (q = ", size.set[q], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n E2 Study", sep=""),
                pch=16, 
                cex=0.5, 
                lty=1, 
                lwd=1, 
                lag=2, 
                span=0.40, 
                degree=2, 
                family="gaussian")

dev.off()

#==============================================================================#
# PE curve as a function of size, for fixed alpha and DE FDR
# Selection of optimal parameters minimizing the PE surface, for fixed DE FDR
f <- 1   # fixed DE FDR = 1e-7

my.file <- paste("/CV PE profile as a function of size (q), for fixed restart probability (alpha=", amin[[f]], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - E2 Study.ps", sep="")
postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=6)

plot.profile.sl(pe.mu=PE.mu.sl[[f]], 
                pe.se=PE.se.sl[[f]], 
                alpha=amin[[f]], 
                size=NULL,
                fdr=fdr.set[f], 
                onese=FALSE, 
                xlab="size (q)", 
                ylab="PE",
                main=paste("CV PE profile as a function of size (q), \n for fixed restart probability (alpha = ", amin[[f]], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n E2 Study", sep=""),
                pch=16, 
                cex=0.5, 
                lty=1, 
                lwd=1, 
                lag=2, 
                span=0.40, 
                degree=2, 
                family="gaussian")

dev.off()

#==============================================================================#
# PE surface as a function of size and restart probability alpha, for fixed DE FDR
f <- 1   # fixed DE FDR = 1e-7

my.file <- paste("/CV PE surface as a function of size (q), for fixed restart probability (alpha=", amin[[f]], ") and DE FDR=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), " - E2 Study.ps", sep="")
postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=FALSE, width=5.5, height=6)

plot.surface.sl(pe.mu=PE.mu.sl[[f]],
                alpha=NULL, 
                size=size.cv.sl[[f]],
                fdr=fdr.set[f], 
                proj=c(1,2),
                onese=FALSE,
                xlab="size (q)", 
                ylab="alpha",
                zlab="PE",
                theta=20, 
                phi=10,
                scale=TRUE, 
                d=5, 
                r=1, 
                shade=0.1, 
                expand=0.5,
                main=paste("CV PE surface as a function of size (q), \n for fixed restart probability (alpha = ", amin[[f]], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n E2 Study", sep=""),
                cex.axis=0.5, 
                lty=1, 
                lwd=0.5)

dev.off()

rm(q, f, crna.fit.sl, my.file)

#==========================================================================================#
# Tuning of RWR parameters q (quantile threshold) 
# For fixed alpha (restart probability)
# For fixed DE FDR
#==========================================================================================#

#==============================================================================#
# Computation of OOB PE as a function of RWR parameters for fixed DE FDR
f <- 1   # fixed DE FDR = 1e-7
a <- 1   # fixed alpha = amin = 0.9

crna.fit.as <- crna.tuning.as(graph=PPI.ig,
                              express=X,
                              factors=list(TF),
                              B=B,
                              alpha=amin[[f]],
                              quant=quant.set,
                              fdr=fdr.set,
                              lag=2, span=0.40, degree=2, family="gaussian",
                              conf=conf,
                              parallel=TRUE,
                              seed=seed)
PE.mu.as <- crna.fit.as$PE.mu
PE.se.as <- crna.fit.as$PE.se
quant.cv.as <- crna.fit.as$quant.cv

#==============================================================================#
# RWR parameters minimizing OOB PE for each DE FDR
f <- 1   # fixed DE FDR = 1e-7
a <- 1   # fixed alpha = amin = 0.9

qmin <- quant.cv.as$qmin[a,f]
qmin.id <- quant.cv.as$qmin.id[a,f]
q1se <- quant.cv.as$q1se[a,f]
q1se.id <- quant.cv.as$q1se.id[a,f]

# fixed DE FDR = 1e-7
# qmin : 0.9855
# qmin.id : 58
# q1se : 0.972
# q1se.id : 49

#==============================================================================#
# PE curve as a function of quant, for fixed alpha and DE FDR
# Selection of optimal parameters minimizing the PE surface, for fixed DE FDR
f <- 1   # fixed DE FDR = 1e-7

my.file <- paste("/CV PE profile as a function of quant (1-q), for fixed restart probability (alpha=", amin[[f]], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - E2 Study.ps", sep="")
postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=6)

plot.profile.as(pe.mu=PE.mu.as[[f]], 
                pe.se=PE.se.as[[f]], 
                alpha=amin[[f]], 
                quant=NULL,
                fdr=fdr.set[f], 
                onese=FALSE, 
                xlab="quant (1-q)", 
                ylab="PE",
                main=paste("CV PE profile as a function of quant (1-q), \n for fixed restart probability (alpha = ", amin[[f]], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n E2 Study", sep=""),
                pch=16, 
                cex=0.5, 
                lty=1, 
                lwd=1, 
                lag=2, 
                span=0.40, 
                degree=2, 
                family="gaussian")

dev.off()

rm(q, a, f, crna.fit.as, my.file)

#==========================================================================================#
# CRNA
# For fixed quantile threshold q, fixed restart probability alpha, and fixed DE FDR
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

CRNA.output <- crna.as(graph=PPI.ig,
                       express=X,
                       factors=list(TF),
                       alpha=amin[[f]],
                       quant=qmin[[f]],
                       FC=FC,
                       FR=FR,
                       R=R,
                       fdr=fdr.set[f],
                       conf=conf,
                       parallel=TRUE,
                       seed=seed)

# Visualization of inferred network graph - Exporting lists to Cytoscape
write.table(CRNA.output$edges, 
            file.path(HOME.path, RESULT.path, paste("CRNA_E2_reduced_edgelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(CRNA.output$states, 
            file.path(HOME.path, RESULT.path, paste("CRNA_E2_reduced_nodelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")

#==========================================================================================#
# IPA
# Fixed DE FDR
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

IPA.output <- ipa(graph=PPI.ig,
                  express=X,
                  factors=list(TF),
                  FR=FR,
                  R=R,
                  fdr=fdr.set[f],
                  conf=conf,
                  parallel=TRUE,
                  seed=seed)

# Visualization of inferred network graph - Exporting lists to Cytoscape
write.table(IPA.output$edges, 
            file.path(HOME.path, RESULT.path, paste("IPA_E2_reduced_edgelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(IPA.output$nodes, 
            file.path(HOME.path, RESULT.path, paste("IPA_E2_reduced_nodelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")

#==========================================================================================#
# CausalR
# Fixed DE FDR
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

CAUSALR.output <- causal(graph=PPI.ig,
                         express=X,
                         factors=list(TF),
                         FR=FR,
                         R=R,
                         fdr=fdr.set[f],
                         conf=conf,
                         parallel=TRUE,
                         seed=seed)

# Visualization of inferred network graph - Exporting lists to Cytoscape
write.table(CAUSALR.output$edges, 
            file.path(HOME.path, RESULT.path, paste("CAUSALR_E2_reduced_edgelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(CAUSALR.output$nodes, 
            file.path(HOME.path, RESULT.path, paste("CAUSALR_E2_reduced_nodelist,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")

#==========================================================================================#
# Save the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_REAL_E2.RData", fsep = .Platform$file.sep))


#===============================================================================================================================#
# MANUSCRIPT FIGURES
#===============================================================================================================================#

#==========================================================================================#
# FIGURE 13:
# PE curve as a function of restart probability alpha, for fixed size and DE FDR
# Selection of optimal parameters minimizing the PE surface, for fixed DE FDR
#==========================================================================================#
# PE curve as a function of quant, for fixed alpha and DE FDR
# Selection of optimal parameters minimizing the PE surface, for fixed DE FDR
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7
q <- 1   # fixed size qm of maximum variance

my.file <- file.path(HOME.path, NATM.path, "Figure13.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=4, paper="special")
par(mfrow=c(1, 2), oma=c(0, 0, 0, 0), mar=c(3, 3, 2, 2), mgp=c(2, 1, 0), xpd=TRUE)

plot.profile.sl(pe.mu=PE.mu.sl[[f]], 
                pe.se=PE.se.sl[[f]], 
                alpha=NULL,
                size=size.set[q], 
                fdr=fdr.set[f], 
                onese=FALSE, 
                xlab=expression(alpha), 
                ylab="PE",
                pch=16, 
                cex=0.5, 
                lty=1, 
                lwd=1, 
                lag=2, 
                span=0.40, 
                degree=2, 
                family="gaussian")

plot.profile.as(pe.mu=PE.mu.as[[f]], 
                pe.se=PE.se.as[[f]], 
                alpha=amin[[f]], 
                quant=NULL,
                fdr=fdr.set[f], 
                onese=FALSE, 
                xlab="quant (1-q)", 
                ylab="PE",
                pch=16, 
                cex=0.5, 
                lty=1, 
                lwd=1, 
                lag=2, 
                span=0.40, 
                degree=2, 
                family="gaussian")

dev.off()

rm(q, my.file)

#==========================================================================================#
# FIGURE 14:
# Heatmap plots
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

my.file <- file.path(HOME.path, NATM.path, "Figure14.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")
par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(3.5, 2.5, 2.5, 2.5), mgp=c(1.5, 0.25, 0), xpd=TRUE)

  plot.heatmap(crna=CRNA.output,
               degs=DEGs[[f]],
               score.colors=c("royalblue", "white", "orangered"),
               express.colors=c("white", "orangered"),
               score.paletteLength=101,
               express.paletteLength=100,
               score.transformation="zscore",
               express.transformation="tanh",
               score.truncation=0.5,
               express.truncation=1.5,
               fontsize.row=5,
               fontsize.col=5,
               lwd=0.05,
               main=list("Activity Scores Heatmap \n E2 Study",
                         "Expression Data Heatmap \n E2 Study"))

dev.off()

rm(my.file)

#==========================================================================================#
# FIGURE 15:
# Visualization of CRNA inferred network graph - Exporting lists to Cytoscape
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

write.table(CRNA.output$edges, 
            file.path(HOME.path, NATM.path, paste("Figure15_edgelist_CRNA_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(CRNA.output$states, 
            file.path(HOME.path, NATM.path, paste("Figure15_nodelist_CRNA_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")

#==========================================================================================#
# FIGURE 19:
# Visualization of IPA & CausalR inferred network graph - Exporting lists to Cytoscape
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

write.table(IPA.output$edges, 
            file.path(HOME.path, NATM.path, paste("Figure19_edgelist_IPA_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(IPA.output$nodes, 
            file.path(HOME.path, NATM.path, paste("Figure19_nodelist_IPA_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")

#==========================================================================================#
# FIGURE 19:
# Visualization of IPA & CausalR inferred network graph - Exporting lists to Cytoscape
#==========================================================================================#
f <- 1   # fixed DE FDR = 1e-7

write.table(CAUSALR.output$edges, 
            file.path(HOME.path, NATM.path, paste("Figure21_edgelist_CAUSALR_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")
write.table(CAUSALR.output$nodes, 
            file.path(HOME.path, NATM.path, paste("Figure21_nodelist_CAUSALR_E2,fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
            quote=FALSE, 
            row.names=FALSE, 
            sep="\t")


