##########################################################################################################################################
# Graph Simulations
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
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/SIMS", fsep=.Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

#==========================================================================================#
# Set some path shortcuts
#==========================================================================================#
CODE.path <- "CODES"
DATA.path <- "DATASETS/METHODS/CRNA/SIMS"
RESULT.path <- "RESULTS/METHODS/CRNA/SIMS"
NATM.path <- "PUBLICATIONS/PAPERS_IN_PREP/NATM_2021"

#==========================================================================================#
# Saving the workspace
#==========================================================================================#
#save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# Load the workspace
#==========================================================================================#
load(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# Set and check the working directory
#==========================================================================================#
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/SIMS", fsep=.Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

#==========================================================================================#
# Set some path shortcuts
#==========================================================================================#
CODE.path <- "CODES"
DATA.path <- "DATASETS/METHODS/CRNA/SIMS"
RESULT.path <- "RESULTS/METHODS/CRNA/SIMS"
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

#=========================================================================================#
# Dynamic loading: we assume that the shared object/dynamic library has already been created
# Linux shared object     : crna.so
# Windows dynamic library : crna.dll
#=========================================================================================#
#dyn.load(x=paste(Sys.getenv("HOME"), "/CODES/C/CRNA/crna", .Platform$dynlib.ext, sep=""),
#         local=TRUE,
#         now=TRUE)

#=========================================================================================#
# Checking if the symbol name of the C/C++ function is loaded and searchable
# and hence available for use in '.C'
#=========================================================================================#
#is.loaded(symbol="name_of_C_function", PACKAGE="crna")

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
          stop("Rmpi must be installed first to use MPI\n")
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
# Saving the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# Initializations - Constants - Parameters
#==========================================================================================#
seed <- 777

nnodes <- 5000   # |V| initial nodes
ndepths <- 4
ntdegs <- 30
depth.levels <- 1:ndepths

#==========================================================================================#
# Generation of Knowledge Base Graph
#==========================================================================================#
set.seed(seed)

#---------------------------- Generating Initial Graph ------------------------#
# Generating scale-free PPI graph according to the Barabasi-Albert model with |V| nodes 
PPI.seed <- igraph::sample_pa(n=nnodes, directed=TRUE, m=8, zero.appeal=.3) %>%
  igraph::rewire(igraph::each_edge(p=.5, loops=FALSE, multiple=TRUE))

igraph::vcount(PPI.seed)  # |V|=5000 initial nodes
igraph::ecount(PPI.seed)  # |E|=39964 initial edges

# randomly assign weights to edges
igraph::E(PPI.seed)$weight <- sample(x=100, size=igraph::ecount(PPI.seed), replace=TRUE)/100

# Filtering edges that have low weights to reduce connectivity
PPI.seed.ig <- PPI.seed - igraph::edge(igraph::E(PPI.seed)[which(igraph::E(PPI.seed)$weight < 0.7)])

# Edge sign for directed graph
igraph::E(PPI.seed.ig)$sign <- 1

# Vertex name
igraph::V(PPI.seed.ig)$name <- as.character(1:igraph::vcount(PPI.seed.ig))

# Removing singleton nodes
PPI.seed.ig <- igraph::graph_from_data_frame(igraph::as_data_frame(PPI.seed.ig))

#-------------------------- Generating Regulatory Graph -----------------------#
# Randomly selecting true DEGs from previous graph
tdegs <- as.character(sample(igraph::V(PPI.seed.ig)$name, ntdegs, FALSE))

# Rebuilding the graph to include all the true DEGS to true DEGs paths including non true DEGs in between
edges.set <- igraph::shortest_paths(PPI.seed.ig, from=tdegs[1], to=tdegs[-1], output="epath", mode="out")$epath[[1]]
for (j in 1:length(tdegs)) {
  tdegs.sp <- igraph::shortest_paths(PPI.seed.ig, from=tdegs[j], to=tdegs[-j], output="epath", mode="out")$epath
  for (k in 1:length(tdegs.sp)) {
    edges.set <- c(edges.set, tdegs.sp[[k]])
  }
}
PPI.seed.ig <- igraph::subgraph.edges(PPI.seed.ig, eids=unique(edges.set), delete.vertices=TRUE)

# Updating true DEGs
tdegs <- as.character(intersect(igraph::V(PPI.seed.ig)$name, tdegs))

rm(j, k, PPI.seed, tdegs.sp, edges.set)

#==========================================================================================#
# Definition of the factors:
# Factor #1: Treatment Fator
#==========================================================================================#
TF <- factor(c(rep(1,50), rep(2,50)), ordered=FALSE, labels=c("C", "T"))
TF.lev <- levels(TF)
TF.ng <- nlevels(TF)
TF.def <- vector(mode="list", length=TF.ng)
TF.tab <- numeric(TF.ng)
for (g in 1:TF.ng) {
  TF.def[[g]] <- which(TF == TF.lev[g])
  TF.tab[g] <- length(TF.def[[g]])
}

rm(g)

#==========================================================================================#
# Sample Names
#==========================================================================================#
sample.names <- c(paste(TF.lev[1], 1:50, sep=""), 
                  paste(TF.lev[2], 1:50, sep=""))

#==========================================================================================#
# Simulated Regulatory Network and Expression values at each depth
# For a given simulation design (ascending, descending, none)
#==========================================================================================#
mu <- c(0.5, 1, 1.5, 2)
sigma2 <- c(0.5, 1, 1.5, 2)
ln.param(mu, sigma2, slow="mu")

sim.obj <- sim(depth.levels=depth.levels,
               tdegs=tdegs,
               ppi.seed.ig=PPI.seed.ig,
               meanlog=mu[1],
               sdlog=sigma2[3],
               factor=TF,
               descent=FALSE,
               seed=seed)

X <- vector(mode="list", length=ndepths)
TREGs <- vector(mode="list", length=ndepths)
TDEGs <- vector(mode="list", length=ndepths)
PPI.ig <- vector(mode="list", length=ndepths)
PPI.ud <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  X[[depth]] <- sim.obj$express[[depth]]
  TREGs[[depth]] <- sim.obj$tregs[[depth]]
  TDEGs[[depth]] <- sim.obj$tdegs[[depth]]
  PPI.ig[[depth]] <- sim.obj$ppi.ig[[depth]]
}

rm(depth, sim.obj)

#==========================================================================================#
# Intersection of expression matrix genes and graph nodes at each depth
# Final Expression Data matrix
# Final Knowledge Base graph
#==========================================================================================#
#----------------------- Updating final Matrix dimensions ---------------------#
for (depth in depth.levels) {
  nodes_names <- igraph::V(PPI.ig[[depth]])$name
  inter_names <- intersect(x=rownames(X[[depth]]), y=nodes_names)
  X[[depth]] <- X[[depth]][inter_names,]
}
n <- ncol(X[[1]])                 # n=100 Common sample size of the data at each depth
p <- nrow(X[[1]])                 # p=421 Common graph size and dimensionality of the data at each depth

#--------------------- Updating final Graph characteristics -------------------#
for (depth in depth.levels) {
  PPI.ig[[depth]] <- delete_vertices(graph=PPI.ig[[depth]], v=setdiff(x=nodes_names, y=inter_names))
}
e <- igraph::ecount(PPI.ig[[1]])  # e=553 edges
p <- igraph::vcount(PPI.ig[[1]])  # p=421 nodes/variables

for (depth in depth.levels) {
  cat("\n Depth: ", depth, "\n", sep="")
  dist <- igraph::distances(
    graph = PPI.ig[[depth]],
    v = igraph::V(PPI.ig[[depth]])$name,
    to = igraph::V(PPI.ig[[depth]])$name,
    mode = "out",
    weights = NULL)[1,]
  dist <- dist[!is.infinite(dist)]
  for (f in 1:length(fdr.set)) {
    cat("fdr: ", fdr.set[f], "\n", sep="")
    cat("d/p: ", d[[depth]][[f]]/p, "\n", sep="")
  }
  cat("Median degree (\\nu): ", median(igraph::degree(PPI.ig[[depth]])), "\n", sep="")
  cat("Median path length (\\lambda): ", median(dist), "\n", sep="")
  cat("e/p: ", e/p, "\n", sep="")
}

rm(f, depth, inter_names, nodes_names)

Depth: 1
fdr: 0.01
d/p: 0.0380047505938
fdr: 0.05
d/p: 0.0546318289786
fdr: 0.1
d/p: 0.0641330166271
Median degree (\nu): 2
Median path length (\lambda): 11.215
e/p: 1.3135391924

Depth: 2
fdr: 0.01
d/p: 0.0736342042755
fdr: 0.05
d/p: 0.114014251781
fdr: 0.1
d/p: 0.118764845606
Median degree (\nu): 2
Median path length (\lambda): 11.215
e/p: 1.3135391924

Depth: 3
fdr: 0.01
d/p: 0.140142517815
fdr: 0.05
d/p: 0.182897862233
fdr: 0.1
d/p: 0.192399049881
Median degree (\nu): 2
Median path length (\lambda): 11.215
e/p: 1.3135391924

Depth: 4
fdr: 0.01
d/p: 0.20190023753
fdr: 0.05
d/p: 0.24703087886
fdr: 0.1
d/p: 0.268408551069
Median degree (\nu): 2
Median path length (\lambda): 11.215
e/p: 1.3135391924

#==========================================================================================#
# Renaming of expression matrix sample names after definition of groups and drop outs  
#==========================================================================================#
for (depth in depth.levels) {
  colnames(X[[depth]]) <- sample.names
}

rm(depth)

#==========================================================================================#
# Initializations - Constants - Parameters
#==========================================================================================#
seed <- 777
R <- 100
B <- 50
FR <- 0.1
FC <- 2

alpha.set <- seq(0.10, 0.90, by=0.10)
size.set <- ceiling(seq(from=p/100, to=p, length=100))
fdr.set <- sort(c(10^seq(-2, -1, by=1), 5* 10^(-2)))

fc.set <- seq(from=1, to=3, by=0.5)
fr.set <- seq(from=0.0, to=0.2, by=0.02)
rep.set <- seq(from=10, to=100, by=10)

#==========================================================================================#
# Estimated simulated DEGs at each depth
# For a given Simulation design (ascending, descending, fixed)
# For a range of DE FDR values 
#==========================================================================================#
degfdr.obj <- vector(mode="list", length=ndepths)
DEGFDR <- vector(mode="list", length=ndepths)
DEGTOPTABLE <- vector(mode="list", length=ndepths)
DEGs <- vector(mode="list", length=ndepths)
d <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth=", depth, "\n")
  degfdr.obj[[depth]] <- deg.fdr(fit=deg(express=X[[depth]], factor=TF), fdr=fdr.set)
  DEGFDR[[depth]] <- vector(mode="list", length=length(fdr.set))
  DEGTOPTABLE[[depth]] <- vector(mode="list", length=length(fdr.set))
  DEGs[[depth]] <- vector(mode="list", length=length(fdr.set))
  d[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    names(degfdr.obj[[depth]]$d[[f]]) <- paste("DE FDR=", format(x=fdr.set[f], digits=3, scientific=TRUE), sep="")
    DEGFDR[[depth]][[f]] <- fdr.set[f]
    DEGTOPTABLE[[depth]][[f]] <- degfdr.obj[[depth]]$degstoptable[[f]]
    DEGs[[depth]][[f]] <- degfdr.obj[[depth]]$degs[[f]]
    d[[depth]][[f]] <- degfdr.obj[[depth]]$d[[f]]
    cat("DEGs: ", DEGs[[depth]][[f]], "\n")
    cat("d: ", d[[depth]][[f]], "\n")
  }
}

rm(f, depth, degfdr.obj)

#==========================================================================================#
# Tuning of RWR model parameters q (model size) and alpha (restart probability)
# For each DE FDR
# For each depth 
# For a given Simulation design (ascending, descending, fixed)
#==========================================================================================#

#==============================================================================#
# Computation of OOB PE as a function of RWR parameters for each DE FDR
alpha.cv <- vector(mode="list", length=ndepths)
size.cv <- vector(mode="list", length=ndepths)
PE.mu <- vector(mode="list", length=ndepths)
PE.se <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  crna.fit <- crna.tuning.sl(graph=PPI.ig[[depth]],
                             express=X[[depth]],
                             factors=list(TF),
                             B=B,
                             alpha=alpha.set,
                             size=size.set,
                             fdr=fdr.set,
                             lag=2, span=0.40, degree=2, family="gaussian",
                             conf=conf,
                             parallel=TRUE,
                             seed=seed)
  PE.mu[[depth]] <- crna.fit$PE.mu
  PE.se[[depth]] <- crna.fit$PE.se
  alpha.cv[[depth]] <- crna.fit$alpha.cv
  size.cv[[depth]] <- crna.fit$size.cv
}

#==============================================================================#
# RWR parameters minimizing OOB PE for each DE FDR
q <- 1   # fixed size qm of maximum variance

amin <- vector(mode="list", length=ndepths)
amin.id <- vector(mode="list", length=ndepths)
a1se <- vector(mode="list", length=ndepths)
a1se.id <- vector(mode="list", length=ndepths)
qmin <- vector(mode="list", length=ndepths)
qmin.id <- vector(mode="list", length=ndepths)
q1se <- vector(mode="list", length=ndepths)
q1se.id <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  amin[[depth]] <- vector(mode="list", length=length(fdr.set))
  amin.id[[depth]] <- vector(mode="list", length=length(fdr.set))
  a1se[[depth]] <- vector(mode="list", length=length(fdr.set))
  a1se.id[[depth]] <- vector(mode="list", length=length(fdr.set))
  qmin[[depth]] <- vector(mode="list", length=length(fdr.set))
  qmin.id[[depth]] <- vector(mode="list", length=length(fdr.set))
  q1se[[depth]] <- vector(mode="list", length=length(fdr.set))
  q1se.id[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    amin[[depth]][[f]] <- alpha.cv[[depth]]$amin[q,f]
    amin.id[[depth]][[f]] <- alpha.cv[[depth]]$amin.id[q,f]
    a1se[[depth]][[f]] <- alpha.cv[[depth]]$a1se[q,f]
    a1se.id[[depth]][[f]] <- alpha.cv[[depth]]$a1se.id[q,f]
    qmin[[depth]][[f]] <- size.cv[[depth]]$qmin[amin.id[[depth]][[f]],f]
    qmin.id[[depth]][[f]] <- size.cv[[depth]]$qmin.id[amin.id[[depth]][[f]],f]
    q1se[[depth]][[f]] <- size.cv[[depth]]$q1se[amin.id[[depth]][[f]],f]
    q1se.id[[depth]][[f]] <- size.cv[[depth]]$q1se.id[amin.id[[depth]][[f]],f]
  }
}

#==========================================================================================#
# Tuning Plots of RWR model parameters q (model size) and alpha (restart probability) 
# For each DE FDR, at each depth
#==========================================================================================#

#==============================================================================#
# PE curve as a function of restart probability alpha, for fixed size and DE FDR, at each depth
# Selection of optimal parameters minimizing the PE surface, for each DE FDR, at each depth
f <- 1   # fixed DE FDR = 1e-2
q <- 1   # fixed size qm of maximum variance

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  
  my.file <- paste("/CV PE profile as a function of restart probability (alpha), for fixed size (q=", size.set[q], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - Depth=", depth, ".ps", sep="")
  postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=6)
  
  plot.profile.sl(pe.mu=PE.mu[[depth]], 
                  pe.se=PE.se[[depth]], 
                  alpha=NULL,
                  size=size.set[q], 
                  fdr=fdr.set[f], 
                  onese=FALSE, 
                  xlab=expression(alpha), 
                  ylab="PE",
                  main=paste("CV PE profile as a function of restart probability (alpha), \n for fixed size (q = ", size.set[q], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n Depth = ", depth, sep=""),
                  pch=16, 
                  cex=0.5, 
                  lty=1, 
                  lwd=1, 
                  lag=2, 
                  span=0.40, 
                  degree=2, 
                  family="gaussian")
  
  dev.off()
}

#==============================================================================#
# PE curve as a function of size, for fixed alpha and DE FDR, at each depth
# Selection of optimal parameters minimizing the PE surface, for each DE FDR, at each depth
f <- 1     # fixed DE FDR = 1e-2

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  my.file <- paste("/CV PE profile as a function of size (q), for fixed restart probability (alpha=", amin[[depth]][[f]], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - Depth=", depth, ".ps", sep="")
  postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=6)
  
  plot.profile.sl(pe.mu=PE.mu[[depth]], 
                  pe.se=PE.se[[depth]], 
                  alpha=amin[[depth]][[f]], 
                  size=NULL,
                  fdr=fdr.set[f], 
                  onese=FALSE, 
                  xlab="size (q)", 
                  ylab="PE",
                  main=paste("CV PE profile as a function of size (q), \n for fixed restart probability (alpha = ", amin[[depth]][[f]], ") and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n Depth = ", depth, sep=""),
                  pch=16, 
                  cex=0.5, 
                  lty=1, 
                  lwd=1, 
                  lag=2, 
                  span=0.40, 
                  degree=2, 
                  family="gaussian")
  
  dev.off()
}

#==============================================================================#
# PE surface as a function of size and restart probability alpha, for each DE FDR, at each depth
f <- 1     # fixed DE FDR = 1e-2

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  
  my.file <- paste("/CV PE surface as a function of size (size), for fixed restart probability (alpha=", amin[[depth]][[f]], ") and DE FDR (fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") - Depth=", depth, ".ps", sep="")
  postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=FALSE, width=5.5, height=6)
  
  plot.surface.sl(pe.mu=PE.mu[[depth]],
                  alpha=NULL, 
                  size=size.cv[[depth]],
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
                  main=paste("CV PE surface as a function of size (size), \n for fixed restart probability (alpha = ", amin[[depth]][[f]], ") \n and DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ") \n Depth = ", depth, sep=""),
                  cex.axis=0.5, 
                  lty=1, 
                  lwd=0.5)
  
  dev.off()
}

rm(q, f, depth, my.file, crna.fit)

#==========================================================================================#
# Empirical vs. Optimal RWR Model Prediction Performances for each DE FDR, at each depth
#==========================================================================================#

#==============================================================================#
# Optimal Performance Metrics as a function of restart probability alpha, for fixed size q and each DE FDR
f <- 1     # fixed DE FDR = 1e-2
q <- 1     # fixed size qm of maximum variance

AUC.size.bo <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  AUC.size.bo[[depth]] <- matrix(data=NA, nrow=R, ncol=length(alpha.set))
}

for (r in 1:R) {
  cat("Replication:", r, "\n")
  set.seed(seed+r)
  bo <- c(sample(x=1:(n/2), size=n/2, replace=TRUE), sample(x=(1+n/2):n, size=n/2, replace=TRUE))
  TF.bo <- TF[bo]
  sim.obj.bo <- sim(depth.levels=depth.levels,
                    tdegs=tdegs,
                    ppi.seed.ig=PPI.seed.ig,
                    meanlog=mu[1],
                    sdlog=sigma2[3],
                    factor=TF.bo,
                    descent=FALSE,
                    seed=seed+r)
  for (depth in depth.levels) {
    X.bo <- sim.obj.bo$express[[depth]]
    tregs.bo <- sim.obj.bo$tregs[[depth]]
    tdegs.bo <- sim.obj.bo$tdegs[[depth]]
    degs.bo <- deg.fdr(fit=deg(express=X.bo, factor=TF.bo), fdr=fdr.set[f])$degs[[f]]
    PPI.ig.bo <- sim.obj.bo$ppi.ig[[depth]]
    PPI.ud.bo <- igraph::as_adjacency_matrix(igraph::as.undirected(PPI.ig.bo), attr="weight")
    for (a in 1:length(alpha.set)) {
      as.bo <- rwr(graph=PPI.ud.bo, 
                   seedSet=degs.bo, 
                   r=alpha.set[a], 
                   random.seed=FALSE, 
                   n.thread=NULL)
      as.bo <- sort(as.bo, decreasing=TRUE)
      preds.bo <- names(as.bo[1:size.set[q]])
      obj <- obspred(X=X.bo, tdegs=tdegs.bo, tregs=tregs.bo, preds=preds.bo)
      perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      AUC.size.bo[[depth]][r,a] <- perfo$auc
    }
  }
}

AUC.size.mu <- matrix(data=NA, nrow=length(alpha.set), ncol=ndepths)
AUC.size.se <- matrix(data=NA, nrow=length(alpha.set), ncol=ndepths)
aucmax.size.id <- numeric(ndepths)
for (depth in depth.levels) {
  AUC.size.mu[,depth] <- apply(X=AUC.size.bo[[depth]], 2, mean, na.rm=TRUE)
  AUC.size.se[,depth] <- apply(X=AUC.size.bo[[depth]], 2, sd, na.rm=TRUE)
  aucmax.size.id[depth] <- which.max(AUC.size.mu[,depth])
}

#==============================================================================#
# Empirical PEs vs. Optimal AUCs as a function of restart probability alpha, for fixed size q and each DE FDR
f <- 1     # fixed DE FDR = 1e-2
q <- 1     # fixed size qm of maximum variance

for (depth in depth.levels) {
  my.file <- paste("/Empirical vs. Optimal AUCs as a function of restart probability (alpha), for fixed size=", size.set[q], " and DE FDR=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), " (Depth=", depth, ").ps", sep="")
  postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=4)
  par(mfrow=c(1,1), oma=c(1, 1, 3, 1), mar=c(2.5, 2.5, 3, 2.5), mgp=c(1.5, 0.25, 0))
  
  plot(x=1:length(alpha.set), y=AUC.size.mu[,depth], type="b", axes=FALSE, xlab="", ylab="", col=1, cex=0.5, ylim=range(AUC.size.mu[,depth]))
  arrows(1:length(alpha.set), AUC.size.mu[,depth], 1:length(alpha.set), AUC.size.mu[,depth]-AUC.size.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  arrows(1:length(alpha.set), AUC.size.mu[,depth], 1:length(alpha.set), AUC.size.mu[,depth]+AUC.size.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  abline(h=AUC.size.mu[aucmax.size.id[depth],depth] - AUC.size.se[aucmax.size.id[depth],depth], col=1, lty=2, lwd=0.5)
  segments(x0=amin.id[[depth]][[f]], y0=0, x1=amin.id[[depth]][[f]], y1=AUC.size.mu[amin.id[[depth]][[f]],depth], col=1, lty=2, lwd=0.5)
  axis(side=1, at=1:length(alpha.set), labels=alpha.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=2, at=pretty(x=AUC.size.mu[,depth]), col=1, col.axis=1, cex.axis=0.5)
  mtext(text=expression(alpha), line=1, side=1, col=1, outer=FALSE)  
  mtext(text="AUC", line=1, side=2, col=1, outer=FALSE)  
  
  par(new=TRUE)
  plot(x=1:length(alpha.set), y=PE.mu[[depth]][q,,f], type="b", axes=FALSE, xlab="", ylab="", col=2, cex=0.5, ylim=range(PE.mu[[depth]][q,,f]))
  arrows(1:length(alpha.set), PE.mu[[depth]][q,,f], 1:length(alpha.set), PE.mu[[depth]][q,,f]-PE.se[[depth]][q,,f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  arrows(1:length(alpha.set), PE.mu[[depth]][q,,f], 1:length(alpha.set), PE.mu[[depth]][q,,f]+PE.se[[depth]][q,,f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  abline(h=PE.mu[[depth]][q,amin.id[[depth]][[f]],f] + PE.se[[depth]][q,amin.id[[depth]][[f]],f], col=2, lty=2, lwd=0.5)
  axis(side=1, at=1:length(alpha.set), labels=alpha.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=4, at=pretty(x=PE.mu[[depth]][q,,f]), col=2, col.axis=2, cex.axis=0.5)
  mtext(text=expression(alpha), line=1, side=1, col=1, outer=FALSE)  
  mtext(text="PE", line=1, side=4, col=2, outer=FALSE)  
  legend(x="top", inset=-0.15, legend=paste("depth = ", depth, sep=""), col=1, lty=1, lwd=1, cex=0.7, xpd=TRUE)
  
  mtext(text=paste("Empirical vs. Optimal AUCs as a function of restart probability (alpha), \n for fixed size (q = ", size.set[q], ") and fixed DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ")", sep=""), 
        line=-1, side=3, outer=TRUE)  
  dev.off()
}

#==============================================================================#
# Optimal Performance Metrics as a function of size q, for fixed restart probability alpha and each DE FDR
f <- 1     # fixed DE FDR = 1e-2

AUC.alpha.bo <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  AUC.alpha.bo[[depth]] <- matrix(data=NA, nrow=R, ncol=length(size.set))
}

for (r in 1:R) {
  cat("Replication:", r, "\n")
  set.seed(seed+r)
  bo <- c(sample(x=1:(n/2), size=n/2, replace=TRUE), sample(x=(1+n/2):n, size=n/2, replace=TRUE))
  TF.bo <- TF[bo]
  sim.obj.bo <- sim(depth.levels=depth.levels,
                    tdegs=tdegs,
                    ppi.seed.ig=PPI.seed.ig,
                    meanlog=mu[1],
                    sdlog=sigma2[3],
                    factor=TF.bo,
                    descent=FALSE,
                    seed=seed+r)
  for (depth in depth.levels) {
    X.bo <- sim.obj.bo$express[[depth]]
    tregs.bo <- sim.obj.bo$tregs[[depth]]
    tdegs.bo <- sim.obj.bo$tdegs[[depth]]
    degs.bo <- deg.fdr(fit=deg(express=X.bo, factor=TF.bo), fdr=fdr.set[f])$degs[[f]]
    PPI.ig.bo <- sim.obj.bo$ppi.ig[[depth]]
    PPI.ud.bo <- igraph::as_adjacency_matrix(igraph::as.undirected(PPI.ig.bo), attr="weight")
    as.bo <- rwr(graph=PPI.ud.bo, 
                 seedSet=degs.bo, 
                 r=amin[[depth]][[f]], 
                 random.seed=FALSE, 
                 n.thread=NULL)
    for (q in 1:length(size.set)) {
      as.bo <- sort(as.bo, decreasing=TRUE)
      preds.bo <- names(as.bo[1:size.set[q]])
      obj <- obspred(X=X.bo, tdegs=tdegs.bo, tregs=tregs.bo, preds=preds.bo)
      perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      AUC.alpha.bo[[depth]][r,q] <- perfo$auc
    }
  }
}

AUC.alpha.mu <- matrix(data=NA, nrow=length(size.set), ncol=ndepths)
AUC.alpha.se <- matrix(data=NA, nrow=length(size.set), ncol=ndepths)
aucmax.alpha.id <- numeric(ndepths)
for (depth in depth.levels) {
  AUC.alpha.mu[,depth] <- apply(X=AUC.alpha.bo[[depth]], 2, mean, na.rm=TRUE)
  AUC.alpha.se[,depth] <- apply(X=AUC.alpha.bo[[depth]], 2, sd, na.rm=TRUE)
  aucmax.alpha.id[depth] <- which.max(AUC.alpha.mu[,depth])
}

#==============================================================================#
# Empirical PEs vs. Optimal AUCs as a function of size q, for fixed restart probability alpha and each DE FDR 
f <- 1     # fixed DE FDR = 1e-2

for (depth in depth.levels) {
  my.file <- paste("/Empirical vs. Optimal AUCs as a function of size (q), for fixed alpha=", amin[[depth]][[f]], " and DE FDR=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), " (Depth=", depth, ").ps", sep="")
  postscript(file=paste(HOME.path, RESULT.path, my.file, sep="/"), horizontal=F, width=8.5, height=4)
  par(mfrow=c(1,1), oma=c(1, 1, 3, 1), mar=c(2.5, 2.5, 3, 2.5), mgp=c(1.5, 0.25, 0))
  
  plot(x=1:length(size.set), y=AUC.alpha.mu[,depth], type="b", axes=FALSE, xlab="", ylab="", col=1, cex=0.5, ylim=range(AUC.alpha.mu[,depth]))
  arrows(1:length(size.set), AUC.alpha.mu[,depth], 1:length(size.set), AUC.alpha.mu[,depth]-AUC.alpha.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  arrows(1:length(size.set), AUC.alpha.mu[,depth], 1:length(size.set), AUC.alpha.mu[,depth]+AUC.alpha.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  abline(h=AUC.alpha.mu[aucmax.alpha.id[depth],depth] - AUC.alpha.se[aucmax.alpha.id[depth],depth], col=1, lty=2, lwd=0.5)
  segments(x0=qmin.id[[depth]][[f]], y0=0, x1=qmin.id[[depth]][[f]], y1=AUC.alpha.mu[qmin.id[[depth]][[f]],depth], col=4, lty=2, lwd=0.5)
  axis(side=1, at=1:length(size.set), labels=size.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=2, at=pretty(x=AUC.alpha.mu[,depth]), col=1, col.axis=1, cex.axis=0.5)
  mtext(text="size (q)", line=1, side=1, col=1, outer=FALSE)  
  mtext(text="AUC", line=1, side=2, col=1, outer=FALSE)  
  
  par(new=TRUE)
  plot(x=1:length(size.set), y=PE.mu[[depth]][,amin.id[[depth]][[f]],f], type="b", axes=FALSE, xlab="", ylab="", col=2, cex=0.5, ylim=range(PE.mu[[depth]][,amin.id[[depth]][[f]],f]))
  arrows(1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f], 1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f]-PE.se[[depth]][,amin.id[[depth]][[f]],f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  arrows(1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f], 1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f]+PE.se[[depth]][,amin.id[[depth]][[f]],f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  segments(x0=qmin.id[[depth]][[f]], y0=0, x1=qmin.id[[depth]][[f]], y1=PE.mu[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f], col=1, lty=2, lwd=0.5)
  abline(h=PE.mu[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f] + PE.se[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f], col=2, lty=2, lwd=0.5)
  fit <- zeroslope(y=PE.mu[[depth]][,amin.id[[depth]][[f]],f], x=1:length(size.set), lag=2, span=0.40, degree=2, family="gaussian", minimum=TRUE)
  lines(x=1:length(size.set), y=fit$loess, type="l", col=4, lty=2, lwd=0.5)
  axis(side=1, at=1:length(size.set), labels=size.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=4, at=pretty(x=PE.mu[[depth]][,amin.id[[depth]][[f]],f]), col=2, col.axis=2, cex.axis=0.5)
  mtext(text="size (q)", line=1, side=1, col=1, outer=FALSE)  
  mtext(text="PE", line=1, side=4, col=2, outer=FALSE)  
  legend(x="top", inset=-0.15, legend=paste("depth = ", depth, sep=""), col=1, lty=1, lwd=1, cex=0.7, xpd=TRUE)
  
  mtext(text=paste("Empirical vs. Optimal AUCs as a function of size (q), \n for fixed restart probability (alpha = ", amin[[depth]][[f]], ") and fixed DE FDR (fdr = ", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ")", sep=""), 
        line=-1, side=3, outer=TRUE)  
  dev.off()
}

rm(a, f, q, r, bo, depth, my.file, fit, perfo,
   sim.obj.bo, TF.bo, X.bo, PPI.ud.bo, PPI.ig.bo, 
   preds.bo, as.bo, degs.bo, tregs.bo, tdegs.bo,
   AUC.alpha.bo, AUC.size.bo)

#==========================================================================================#
# Empirical RWR Prediction Performances, for each DE FDR, at each depth
#==========================================================================================#

#==============================================================================#
# Predicted REGs and DEGs, for fixed alpha, for each DE FDR, at each depth
PREDs <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  PPI.ud[[depth]] <- igraph::as_adjacency_matrix(igraph::as.undirected(PPI.ig[[depth]]), attr="weight")
  PREDs[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    as <- rwr(graph=PPI.ud[[depth]], 
              seedSet=DEGs[[depth]][[f]], 
              r=amin[[depth]][[f]], 
              random.seed=FALSE, 
              n.thread=NULL)
    as <- sort(as, decreasing=TRUE)
    PREDs[[depth]][[f]] <- names(as[1:qmin[[depth]][[f]]])
    cat("Predicted REGs and DEGs: ", PREDs[[depth]][[f]], "\n")
    cat("alpha: ", amin[[depth]][[f]], "\n")
    cat("size: ", qmin[[depth]][[f]], "\n")
  }
}

#==============================================================================#
# Empirical performance metrics with tuning, for each DE FDR, at each depth
perfo <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  perfo[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=PREDs[[depth]][[f]])
    perfo[[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
    cat("RWR TPR: ", perfo[[depth]][[f]]$`tpr`, "\n")
    cat("RWR TNR: ", perfo[[depth]][[f]]$`tnr`, "\n")
    cat("RWR AUC: ", perfo[[depth]][[f]]$`auc`, "\n")
    cat("RWR FDR: ", perfo[[depth]][[f]]$`fdr`, "\n")
    cat("RWR FOR: ", perfo[[depth]][[f]]$`for`, "\n")
  }
}

rm(as, depth, obj)

#==========================================================================================#
# CRNA
# Effect of Affinity Scores Fold Change
# For fixed graph size q, and fixed restart probability alpha
# For fixed DE FDR, fixed frequency, and fixed number of replications
#==========================================================================================#
CRNA.fc <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.fc[[depth]] <- vector(mode="list", length=length(fc.set))
  for (fc in 1:length(fc.set)) {
    cat("fc: ", fc.set[fc], "\n")
    CRNA.fc[[depth]][[fc]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      CRNA.fc[[depth]][[fc]][[f]] <- crna.sl(graph=PPI.ig[[depth]],
                                             express=X[[depth]],
                                             factors=list(TF),
                                             alpha=amin[[depth]][[f]],
                                             size=qmin[[depth]][[f]],
                                             FC=fc.set[fc],
                                             FR=FR,
                                             R=R,
                                             fdr=fdr.set[f],
                                             conf=conf,
                                             parallel=TRUE,
                                             seed=seed)
    }
  }
}

rm(depth, fc, f)

#==========================================================================================#
# CRNA
# Effect of Affinity Scores Fold Change on Method Prediction Performances
# For fixed graph size q, and fixed restart probability alpha
# For fixed DE FDR, fixed frequency, and fixed number of replications
#==========================================================================================#
CRNA.fc.perfo <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.fc.perfo[[depth]] <- vector(mode="list", length=length(fc.set))
  for (fc in 1:length(fc.set)) {
    cat("fc: ", fc.set[fc], "\n")
    CRNA.fc.perfo[[depth]][[fc]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CRNA.fc[[depth]][[fc]][[f]]$nodes)
      CRNA.fc.perfo[[depth]][[fc]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("CRNA TPR: ", CRNA.fc.perfo[[depth]][[fc]][[f]]$`tpr`, "\n")
      cat("CRNA TNR: ", CRNA.fc.perfo[[depth]][[fc]][[f]]$`tnr`, "\n")
      cat("CRNA AUC: ", CRNA.fc.perfo[[depth]][[fc]][[f]]$`auc`, "\n")
      cat("CRNA FDR: ", CRNA.fc.perfo[[depth]][[fc]][[f]]$`fdr`, "\n")
      cat("CRNA FOR: ", CRNA.fc.perfo[[depth]][[fc]][[f]]$`for`, "\n")
    }
  }
}

TPR.fc <- vector(mode="list", length=ndepths)
TNR.fc <- vector(mode="list", length=ndepths)
AUC.fc <- vector(mode="list", length=ndepths)
FDR.fc <- vector(mode="list", length=ndepths)
FOR.fc <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  TPR.fc[[depth]] <- matrix(data=NA, nrow=length(fc.set), ncol=length(fdr.set), dimnames=list(fc.set,fdr.set))
  TNR.fc[[depth]] <- matrix(data=NA, nrow=length(fc.set), ncol=length(fdr.set), dimnames=list(fc.set,fdr.set))
  AUC.fc[[depth]] <- matrix(data=NA, nrow=length(fc.set), ncol=length(fdr.set), dimnames=list(fc.set,fdr.set))
  FDR.fc[[depth]] <- matrix(data=NA, nrow=length(fc.set), ncol=length(fdr.set), dimnames=list(fc.set,fdr.set))
  FOR.fc[[depth]] <- matrix(data=NA, nrow=length(fc.set), ncol=length(fdr.set), dimnames=list(fc.set,fdr.set))
  for (fc in 1:length(fc.set)) {
    cat("fc: ", fc.set[fc], "\n")
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      TPR.fc[[depth]][fc,f] <- c(CRNA.fc.perfo[[depth]][[fc]][[f]]$`tpr`)
      TNR.fc[[depth]][fc,f] <- c(CRNA.fc.perfo[[depth]][[fc]][[f]]$`tnr`)
      AUC.fc[[depth]][fc,f] <- c(CRNA.fc.perfo[[depth]][[fc]][[f]]$`auc`)
      FDR.fc[[depth]][fc,f] <- c(CRNA.fc.perfo[[depth]][[fc]][[f]]$`fdr`)
      FOR.fc[[depth]][fc,f] <- c(CRNA.fc.perfo[[depth]][[fc]][[f]]$`for`)
    }
  }
}

par(mfrow=c(ndepths,5), oma=c(1, 1, 3, 1), mar=c(2.5, 2.5, 3, 2.5), mgp=c(1.5, 0.25, 0))

for (depth in depth.levels) {
  matplot(x=TPR.fc[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TPR.fc[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(TNR.fc[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TNR.fc[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.fc[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(AUC.fc[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.fc[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FDR.fc[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.fc[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FOR.fc[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
  
  legend(x="topright", inset=-0.2, legend=fdr.set, col=c(2,3,4), lty=1, lwd=1, cex=0.3, xpd=TRUE)
  mtext(text=paste("Comparative Effect of Affinity Scores Fold Change on CRNA Prediction Performances by DE FDR", sep=""), line=-1, side=3, outer=TRUE)  
}

rm(depth, fc, f, obj)

#==========================================================================================#
# CRNA
# Effect of Frequency of Occurrence in Replications 
# For fixed graph size q, fixed restart probability alpha, and fixed fold change
# For fixed DE FDR, and fixed number of replications
#==========================================================================================#
CRNA.fr <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.fr[[depth]] <- vector(mode="list", length=length(fr.set))
  for (fr in 1:length(fr.set)) {
    cat("fr: ", fr.set[fr], "\n")
    CRNA.fr[[depth]][[fr]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      CRNA.fr[[depth]][[fr]][[f]] <- crna.sl(graph=PPI.ig[[depth]],
                                             express=X[[depth]],
                                             factors=list(TF),
                                             alpha=amin[[depth]][[f]],
                                             size=qmin[[depth]][[f]],
                                             FC=FC,
                                             FR=fr.set[fr],
                                             R=R,
                                             fdr=fdr.set[f],
                                             conf=conf,
                                             parallel=TRUE,
                                             seed=seed)
    }
  }
}

rm(depth, fr, f)

#==========================================================================================#
# CRNA
# Effect of Frequency of Occurrence in Replications on Method Prediction Performances
# For fixed graph size q, fixed restart probability alpha, and fixed fold change
# For fixed DE FDR, and fixed number of replications
#==========================================================================================#
CRNA.fr.perfo <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.fr.perfo[[depth]] <- vector(mode="list", length=length(fr.set))
  for (fr in 1:length(fr.set)) {
    cat("fr: ", fr.set[fr], "\n")
    CRNA.fr.perfo[[depth]][[fr]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CRNA.fr[[depth]][[fr]][[f]]$nodes)
      CRNA.fr.perfo[[depth]][[fr]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("CRNA TPR: ", CRNA.fr.perfo[[depth]][[fr]][[f]]$`tpr`, "\n")
      cat("CRNA TNR: ", CRNA.fr.perfo[[depth]][[fr]][[f]]$`tnr`, "\n")
      cat("CRNA AUC: ", CRNA.fr.perfo[[depth]][[fr]][[f]]$`auc`, "\n")
      cat("CRNA FDR: ", CRNA.fr.perfo[[depth]][[fr]][[f]]$`fdr`, "\n")
      cat("CRNA FOR: ", CRNA.fr.perfo[[depth]][[fr]][[f]]$`for`, "\n")
    }
  }
}

TPR.fr <- vector(mode="list", length=ndepths)
TNR.fr <- vector(mode="list", length=ndepths)
AUC.fr <- vector(mode="list", length=ndepths)
FDR.fr <- vector(mode="list", length=ndepths)
FOR.fr <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  TPR.fr[[depth]] <- matrix(data=NA, nrow=length(fr.set), ncol=length(fdr.set), dimnames=list(fr.set,fdr.set))
  TNR.fr[[depth]] <- matrix(data=NA, nrow=length(fr.set), ncol=length(fdr.set), dimnames=list(fr.set,fdr.set))
  AUC.fr[[depth]] <- matrix(data=NA, nrow=length(fr.set), ncol=length(fdr.set), dimnames=list(fr.set,fdr.set))
  FDR.fr[[depth]] <- matrix(data=NA, nrow=length(fr.set), ncol=length(fdr.set), dimnames=list(fr.set,fdr.set))
  FOR.fr[[depth]] <- matrix(data=NA, nrow=length(fr.set), ncol=length(fdr.set), dimnames=list(fr.set,fdr.set))
  for (fr in 1:length(fr.set)) {
    cat("fr: ", fr.set[fr], "\n")
    for (f in 1:length(fdr.set)) {
      cat("fdr: ", fdr.set[f], "\n")
      TPR.fr[[depth]][fr,f] <- c(CRNA.fr.perfo[[depth]][[fr]][[f]]$`tpr`)
      TNR.fr[[depth]][fr,f] <- c(CRNA.fr.perfo[[depth]][[fr]][[f]]$`tnr`)
      AUC.fr[[depth]][fr,f] <- c(CRNA.fr.perfo[[depth]][[fr]][[f]]$`auc`)
      FDR.fr[[depth]][fr,f] <- c(CRNA.fr.perfo[[depth]][[fr]][[f]]$`fdr`)
      FOR.fr[[depth]][fr,f] <- c(CRNA.fr.perfo[[depth]][[fr]][[f]]$`for`)
    }
  }
}

par(mfrow=c(ndepths,5), oma=c(1, 1, 3, 1), mar=c(2.5, 2.5, 3, 2.5), mgp=c(1.5, 0.25, 0))

for (depth in depth.levels) {
  matplot(x=TPR.fr[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TPR.fr[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=TNR.fr[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TNR.fr[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.fr[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(AUC.fr[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.fr[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FDR.fr[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.fr[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FOR.fr[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
  
  legend(x="topright", inset=-0.2, legend=fdr.set, col=c(2,3,4), lty=1, lwd=1, cex=0.3, xpd=TRUE)
  mtext(text=paste("Comparative Effect of Frequencies on CRNA Prediction Performances by DE FDR", sep=""), line=-1, side=3, outer=TRUE)  
}

rm(depth, fr, f, obj)

#==========================================================================================#
# CRNA
# Effect of Procedural Replications 
# For fixed graph size q, fixed restart probability alpha, and fixed fold change
# For fixed DE FDR, and fixed frequency
#==========================================================================================#
CRNA.rep <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.rep[[depth]] <- vector(mode="list", length=length(rep.set))
  for (r in 1:length(rep.set)) {
    cat("r: ", rep.set[r], "\n")
    CRNA.rep[[depth]][[r]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      CRNA.rep[[depth]][[r]][[f]] <- crna.sl(graph=PPI.ig[[depth]],
                                             express=X[[depth]],
                                             factors=list(TF),
                                             alpha=amin[[depth]][[f]],
                                             size=qmin[[depth]][[f]],
                                             FC=FC,
                                             FR=FR,
                                             R=rep.set[r],
                                             fdr=fdr.set[f],
                                             conf=conf,
                                             parallel=TRUE,
                                             seed=seed)
    }
  }
}

rm(depth, r, f)

#==========================================================================================#
# Comparative Effect of Procedural Replications on Method Prediction Performances
# For fixed DE FDR, and fixed frequency
#==========================================================================================#
CRNA.rep.perfo <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.rep.perfo[[depth]] <- vector(mode="list", length=length(rep.set))
  for (r in 1:length(rep.set)) {
    cat("r: ", rep.set[r], "\n")
    CRNA.rep.perfo[[depth]][[r]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CRNA.rep[[depth]][[r]][[f]]$nodes)
      CRNA.rep.perfo[[depth]][[r]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("CRNA TPR: ", CRNA.rep.perfo[[depth]][[r]][[f]]$`tpr`, "\n")
      cat("CRNA TNR: ", CRNA.rep.perfo[[depth]][[r]][[f]]$`tnr`, "\n")
      cat("CRNA AUC: ", CRNA.rep.perfo[[depth]][[r]][[f]]$`auc`, "\n")
      cat("CRNA FDR: ", CRNA.rep.perfo[[depth]][[r]][[f]]$`fdr`, "\n")
      cat("CRNA FOR: ", CRNA.rep.perfo[[depth]][[r]][[f]]$`for`, "\n")
    }
  }
}

TPR.rep <- vector(mode="list", length=ndepths)
TNR.rep <- vector(mode="list", length=ndepths)
AUC.rep <- vector(mode="list", length=ndepths)
FDR.rep <- vector(mode="list", length=ndepths)
FOR.rep <- vector(mode="list", length=ndepths)
for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  TPR.rep[[depth]] <- matrix(data=NA, nrow=length(rep.set), ncol=length(fdr.set), dimnames=list(rep.set,fdr.set))
  TNR.rep[[depth]] <- matrix(data=NA, nrow=length(rep.set), ncol=length(fdr.set), dimnames=list(rep.set,fdr.set))
  AUC.rep[[depth]] <- matrix(data=NA, nrow=length(rep.set), ncol=length(fdr.set), dimnames=list(rep.set,fdr.set))
  FDR.rep[[depth]] <- matrix(data=NA, nrow=length(rep.set), ncol=length(fdr.set), dimnames=list(rep.set,fdr.set))
  FOR.rep[[depth]] <- matrix(data=NA, nrow=length(rep.set), ncol=length(fdr.set), dimnames=list(rep.set,fdr.set))
  for (r in 1:length(rep.set)) {
    cat("rep: ", rep.set[r], "\n")
    for (f in 1:length(fdr.set)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      TPR.rep[[depth]][r,f] <- c(CRNA.rep.perfo[[depth]][[r]][[f]]$`tpr`)
      TNR.rep[[depth]][r,f] <- c(CRNA.rep.perfo[[depth]][[r]][[f]]$`tnr`)
      AUC.rep[[depth]][r,f] <- c(CRNA.rep.perfo[[depth]][[r]][[f]]$`auc`)
      FDR.rep[[depth]][r,f] <- c(CRNA.rep.perfo[[depth]][[r]][[f]]$`fdr`)
      FOR.rep[[depth]][r,f] <- c(CRNA.rep.perfo[[depth]][[r]][[f]]$`for`)
    }
  }
}

par(mfrow=c(ndepths,5), oma=c(1, 1, 3, 1), mar=c(2.5, 2.5, 3, 2.5), mgp=c(1.5, 0.25, 0))

for (depth in depth.levels) {
  matplot(x=TPR.rep[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TPR.rep[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=TNR.rep[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(TNR.rep[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.rep[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(AUC.rep[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.rep[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FDR.rep[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.rep[[depth]], type="l", axes=FALSE, xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(FOR.rep[[depth]], na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
  
  legend(x="topright", inset=-0.2, legend=fdr.set, col=c(2,3,4), lty=1, lwd=1, cex=0.3, xpd=TRUE)
  mtext(text=paste("Comparative Effect of Procedural Replications on CRNA Prediction Performances by DE FDR", sep=""), line=-1, side=3, outer=TRUE)  
}

rm(depth, r, f, obj)

#==========================================================================================#
# CRNA
# For fixed graph size q, fixed restart probability alpha, and fixed fold change
# For fixed number of replications, and fixed frequency
#==========================================================================================#
CRNA.output <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.output[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    CRNA.output[[depth]][[f]]  <- crna.sl(graph=PPI.ig[[depth]],
                                          express=X[[depth]],
                                          factors=list(TF),
                                          alpha=amin[[depth]][[f]],
                                          size=qmin[[depth]][[f]],
                                          FC=FC,
                                          FR=FR,
                                          R=R,
                                          fdr=fdr.set[f],
                                          conf=conf,
                                          parallel=TRUE,
                                          seed=seed)
  }
}

# Exporting lists to Cytoscape
for (depth in depth.levels) {
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    write.table(CRNA.output[[depth]][[f]]$edges, 
                file.path(HOME.path, RESULT.path, paste("CRNA_SIMS_ascending_edgelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
    write.table(CRNA.output[[depth]][[f]]$states, 
                file.path(HOME.path, RESULT.path, paste("CRNA_SIMS_ascending_nodelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
  }
}

rm(depth, f)

#==========================================================================================#
# IPA
# For fixed number of replications, and fixed frequency
#==========================================================================================#
IPA.output <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  IPA.output[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    IPA.output[[depth]][[f]] <- ipa(graph=PPI.ig[[depth]],
                                    express=X[[depth]],
                                    factors=list(TF),
                                    FR=FR,
                                    R=R,
                                    fdr=fdr.set[f],
                                    conf=conf,
                                    parallel=TRUE,
                                    seed=seed)
  }
}

# Exporting lists to Cytoscape
for (depth in depth.levels) {
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    write.table(IPA.output[[depth]][[f]]$edges, 
                file.path(HOME.path, RESULT.path, paste("IPA_SIMS_ascending_edgelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
    write.table(IPA.output[[depth]][[f]]$nodes, 
                file.path(HOME.path, RESULT.path, paste("IPA_SIMS_ascending_nodelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
  }
}

rm(depth, f)

#==========================================================================================#
# CAUSALR
# For fixed number of replications, and fixed frequency
#==========================================================================================#
CAUSALR.output <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CAUSALR.output[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    CAUSALR.output[[depth]][[f]] <- causal(graph=PPI.ig[[depth]],
                                           express=X[[depth]],
                                           factors=list(TF),
                                           FR=FR,
                                           R=R,
                                           fdr=fdr.set[f],
                                           conf=conf,
                                           parallel=TRUE,
                                           seed=seed)
  }
}

# Exporting lists to Cytoscape
for (depth in depth.levels) {
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    write.table(CAUSALR.output[[depth]][[f]]$edges, 
                file.path(HOME.path, RESULT.path, paste("CAUSALR_SIMS_ascending_edgelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
    write.table(CAUSALR.output[[depth]][[f]]$nodes, 
                file.path(HOME.path, RESULT.path, paste("CAUSALR_SIMS_ascending_nodelist - depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
  }
}

rm(depth, f)

#==========================================================================================#
# Comparative Method Prediction Performances
# For fixed DE FDR, fixed number of replications, and fixed frequency
#==========================================================================================#
CRNA.output.perfo <- vector(mode="list", length=ndepths)
IPA.output.perfo <- vector(mode="list", length=ndepths)
CAUSALR.output.perfo <- vector(mode="list", length=ndepths)

for (depth in depth.levels) {
  cat("Depth: ", depth, "\n")
  CRNA.output.perfo[[depth]] <- vector(mode="list", length=length(fdr.set))
  IPA.output.perfo[[depth]] <- vector(mode="list", length=length(fdr.set))
  CAUSALR.output.perfo[[depth]] <- vector(mode="list", length=length(fdr.set))
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CRNA.output[[depth]][[f]]$nodes)
    CRNA.output.perfo[[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
    cat("CRNA TPR: ", CRNA.output.perfo[[depth]][[f]]$`tpr`, "\n")
    cat("CRNA TNR: ", CRNA.output.perfo[[depth]][[f]]$`tnr`, "\n")
    cat("CRNA AUC: ", CRNA.output.perfo[[depth]][[f]]$`auc`, "\n")
    cat("CRNA FDR: ", CRNA.output.perfo[[depth]][[f]]$`fdr`, "\n")
    cat("CRNA FOR: ", CRNA.output.perfo[[depth]][[f]]$`for`, "\n")
    obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=IPA.output[[depth]][[f]]$nodes)
    IPA.output.perfo[[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
    cat("IPA TPR: ", IPA.output.perfo[[depth]][[f]]$`tpr`, "\n")
    cat("IPA TNR: ", IPA.output.perfo[[depth]][[f]]$`tnr`, "\n")
    cat("IPA AUC: ", IPA.output.perfo[[depth]][[f]]$`auc`, "\n")
    cat("IPA FDR: ", IPA.output.perfo[[depth]][[f]]$`fdr`, "\n")
    cat("IPA FOR: ", IPA.output.perfo[[depth]][[f]]$`for`, "\n")
    obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CAUSALR.output[[depth]][[f]]$nodes)
    CAUSALR.output.perfo[[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
    cat("CAUSALR TPR: ", CAUSALR.output.perfo[[depth]][[f]]$`tpr`, "\n")
    cat("CAUSALR TNR: ", CAUSALR.output.perfo[[depth]][[f]]$`tnr`, "\n")
    cat("CAUSALR AUC: ", CAUSALR.output.perfo[[depth]][[f]]$`auc`, "\n")
    cat("CAUSALR FDR: ", CAUSALR.output.perfo[[depth]][[f]]$`fdr`, "\n")
    cat("CAUSALR FOR: ", CAUSALR.output.perfo[[depth]][[f]]$`for`, "\n\n")
  }
}

for (depth in depth.levels) {
  cat("\n Depth: ", depth, "\n")
  for (f in 1:length(fdr.set)) {
    cat("DE FDR: ", fdr.set[f], "\n")
    tpr.m <- rep(NA, R)
    tnr.m <- rep(NA, R)
    auc.m <- rep(NA, R)
    fdr.m <- rep(NA, R)
    for.m <- rep(NA, R)
    for (r in 1:R) {
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CRNA.output[[depth]][[f]]$nodes.rep[[r]])
      obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      tpr.m[r] <- obj.perfo$`tpr`
      tnr.m[r] <- obj.perfo$`tnr`
      auc.m[r] <- obj.perfo$`auc`
      fdr.m[r] <- obj.perfo$`fdr`
      for.m[r] <- obj.perfo$`for`
    }
    tpr.se <- sd(tpr.m, na.rm=TRUE)
    tnr.se <- sd(tnr.m, na.rm=TRUE)
    auc.se <- sd(auc.m, na.rm=TRUE)
    fdr.se <- sd(fdr.m, na.rm=TRUE)
    for.se <- sd(for.m, na.rm=TRUE)
    cat("CRNA TPR se: ", round(tpr.se,2), "\n")
    cat("CRNA TNR se: ", round(tnr.se,2), "\n")
    cat("CRNA AUC se: ", round(auc.se,2), "\n")
    cat("CRNA FDR se: ", round(fdr.se,2), "\n")
    cat("CRNA FOR se: ", round(for.se,2), "\n\n")
  
    tpr.m <- rep(NA, R)
    tnr.m <- rep(NA, R)
    auc.m <- rep(NA, R)
    fdr.m <- rep(NA, R)
    for.m <- rep(NA, R)
    for (r in 1:R) {
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=IPA.output[[depth]][[f]]$nodes.rep[[r]])
      obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      tpr.m[r] <- obj.perfo$`tpr`
      tnr.m[r] <- obj.perfo$`tnr`
      auc.m[r] <- obj.perfo$`auc`
      fdr.m[r] <- obj.perfo$`fdr`
      for.m[r] <- obj.perfo$`for`
    }
    tpr.se <- sd(tpr.m, na.rm=TRUE)
    tnr.se <- sd(tnr.m, na.rm=TRUE)
    auc.se <- sd(auc.m, na.rm=TRUE)
    fdr.se <- sd(fdr.m, na.rm=TRUE)
    for.se <- sd(for.m, na.rm=TRUE)
    cat("IPA TPR se: ", round(tpr.se,2), "\n")
    cat("IPA TNR se: ", round(tnr.se,2), "\n")
    cat("IPA AUC se: ", round(auc.se,2), "\n")
    cat("IPA FDR se: ", round(fdr.se,2), "\n")
    cat("IPA FOR se: ", round(for.se,2), "\n\n")
    
    tpr.m <- rep(NA, R)
    tnr.m <- rep(NA, R)
    auc.m <- rep(NA, R)
    fdr.m <- rep(NA, R)
    for.m <- rep(NA, R)
    for (r in 1:R) {
      obj <- obspred(X=X[[depth]], tdegs=TDEGs[[depth]], tregs=TREGs[[depth]], preds=CAUSALR.output[[depth]][[f]]$nodes.rep[[r]])
      obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      tpr.m[r] <- obj.perfo$`tpr`
      tnr.m[r] <- obj.perfo$`tnr`
      auc.m[r] <- obj.perfo$`auc`
      fdr.m[r] <- obj.perfo$`fdr`
      for.m[r] <- obj.perfo$`for`
    }
    tpr.se <- sd(tpr.m, na.rm=TRUE)
    tnr.se <- sd(tnr.m, na.rm=TRUE)
    auc.se <- sd(auc.m, na.rm=TRUE)
    fdr.se <- sd(fdr.m, na.rm=TRUE)
    for.se <- sd(for.m, na.rm=TRUE)
    cat("CAUSALR TPR se: ", round(tpr.se,2), "\n")
    cat("CAUSALR TNR se: ", round(tnr.se,2), "\n")
    cat("CAUSALR AUC se: ", round(auc.se,2), "\n")
    cat("CAUSALR FDR se: ", round(fdr.se,2), "\n")
    cat("CAUSALR FOR se: ", round(for.se,2), "\n\n")
    
  }
}

rm(depth, f, r, obj, obj.perfo, 
   tpr.m, tnr.m, auc.m, fdr.m, for.m,
   tpr.se, tnr.se, auc.se, fdr.se, for.se)

#===============================================================================================================================#
# Unloading of the shared object/dynamic library
#===============================================================================================================================#
# dyn.unload(x=paste(Sys.getenv("HOME"), "/CODES/C/CRNA/crna", .Platform$dynlib.ext, sep=""))

#==========================================================================================#
# Save the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending.RData", fsep=.Platform$file.sep))



#===============================================================================================================================#
# MANUSCRIPT FIGURES
#===============================================================================================================================#

#==========================================================================================#
# FIGURE 01:
# Visualization of a knowledge base graphs at all depths (1-4), in any case of expression levels
# Exporting lists to Cytoscape
#==========================================================================================#
for (depth in depth.levels) {
  write.table(igraph::as_data_frame(PPI.ig[[depth]]), 
              file.path(HOME.path, NATM.path, paste0("Figure01_edgelist,depth=", depth, ".txt"), fsep=.Platform$file.sep),
              quote=FALSE, 
              row.names=FALSE, 
              sep="\t")
  nodelist <- data.frame("symbol"=igraph::V(PPI.ig[[depth]])$name, 
                         "node_color"="b", 
                         stringsAsFactors=FALSE)
  nodelist[which(nodelist$symbol %in% TDEGs[[depth]]),2] <- "d"
  nodelist[which(nodelist$symbol %in% TREGs[[depth]]),2] <- "r"
  write.table(nodelist, 
              file.path(HOME.path, NATM.path, paste0("Figure01_nodelist,depth=", depth, ".txt"), fsep=.Platform$file.sep), 
              quote=FALSE, 
              row.names=FALSE, 
              sep="\t")
}

rm(depth, nodelist)

#==========================================================================================#
# FIGURE 05:
# PE surface as a function of size and restart probability alpha, for fixed DE FDR, at each depth
#==========================================================================================#
f <- 1     # fixed DE FDR = 1e-2

my.file <- file.path(HOME.path, NATM.path, "Figure05.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=6, paper="special")

par(mfrow=c(2, 2), oma=c(0, 0, 0, 0), mar=c(0, 2, 0, 2), mgp=c(2, 1, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot.surface.sl(pe.mu=PE.mu[[depth]],
                  alpha=NULL, 
                  size=size.cv[[depth]],
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
                  cex.axis=0.5, 
                  lty=1, 
                  lwd=0.5)
  mtext(text=paste("depth = ", depth, sep=""), line=-3, side=3, outer=FALSE)
}

dev.off()

rm(depth, my.file)

#==========================================================================================#
# FIGURE 06:
# Empirical PEs vs. Optimal AUCs as a function of restart probability alpha, for fixed size q and DE FDR
#==========================================================================================#
f <- 1     # fixed DE FDR = 1e-2
q <- 1     # fixed size qm of maximum variance

my.file <- file.path(HOME.path, NATM.path, "Figure06.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(2, 2), oma=c(0, 0, 0, 0), mar=c(3.5, 2.5, 2.5, 2.5), mgp=c(1.5, 0.25, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot(x=1:length(alpha.set), y=PE.mu[[depth]][q,,f], type="b", axes=FALSE, xlab="", ylab="", col=2, cex=0.5, ylim=range(PE.mu[[depth]][q,,f]))
  arrows(1:length(alpha.set), PE.mu[[depth]][q,,f], 1:length(alpha.set), PE.mu[[depth]][q,,f]-PE.se[[depth]][q,,f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  arrows(1:length(alpha.set), PE.mu[[depth]][q,,f], 1:length(alpha.set), PE.mu[[depth]][q,,f]+PE.se[[depth]][q,,f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  abline(h=PE.mu[[depth]][q,amin.id[[depth]][[f]],f] + PE.se[[depth]][q,amin.id[[depth]][[f]],f], col=2, lty=2, lwd=0.5)
  axis(side=1, at=1:length(alpha.set), labels=alpha.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=4, at=pretty(x=PE.mu[[depth]][q,,f]), col=2, col.axis=2, cex.axis=0.5)
  mtext(text=expression(alpha), line=1, side=1, col=1, outer=FALSE)  
  mtext(text="PE", line=1, side=4, col=2, outer=FALSE)  
  
  par(new=TRUE)
  plot(x=1:length(alpha.set), y=AUC.size.mu[,depth], type="b", axes=FALSE, xlab="", ylab="", col=1, cex=0.5, ylim=range(AUC.size.mu[,depth]))
  arrows(1:length(alpha.set), AUC.size.mu[,depth], 1:length(alpha.set), AUC.size.mu[,depth]-AUC.size.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  arrows(1:length(alpha.set), AUC.size.mu[,depth], 1:length(alpha.set), AUC.size.mu[,depth]+AUC.size.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  abline(h=AUC.size.mu[aucmax.size.id[depth],depth] - AUC.size.se[aucmax.size.id[depth],depth], col=1, lty=2, lwd=0.5)
  segments(x0=amin.id[[depth]][[f]], y0=0, x1=amin.id[[depth]][[f]], y1=AUC.size.mu[amin.id[[depth]][[f]],depth], col=1, lty=2, lwd=0.5)
  axis(side=1, at=1:length(alpha.set), labels=alpha.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=2, at=pretty(x=AUC.size.mu[,depth]), col=1, col.axis=1, cex.axis=0.5)
  mtext(text=expression(alpha), line=1, side=1, col=1, outer=FALSE)  
  mtext(text="AUC", line=1, side=2, col=1, outer=FALSE)  
  mtext(text=paste("depth = ", depth, sep=""), line=0, side=3, outer=FALSE)
}

dev.off()

rm(f, q, depth, my.file)

#==========================================================================================#
# FIGURE 07:
# Empirical PEs vs. Optimal AUCs as a function of size q, for fixed restart probability alpha and DE FDR 
#==========================================================================================#
f <- 1     # fixed DE FDR = 1e-2

my.file <- file.path(HOME.path, NATM.path, "Figure07.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(2, 2), oma=c(0, 0, 0, 0), mar=c(3.5, 2.5, 2.5, 2.5), mgp=c(1.5, 0.25, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot(x=1:length(size.set), y=PE.mu[[depth]][,amin.id[[depth]][[f]],f], type="b", axes=FALSE, xlab="", ylab="", col=2, cex=0.5, ylim=range(PE.mu[[depth]][,amin.id[[depth]][[f]],f]))
  arrows(1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f], 1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f]-PE.se[[depth]][,amin.id[[depth]][[f]],f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  arrows(1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f], 1:length(size.set), PE.mu[[depth]][,amin.id[[depth]][[f]],f]+PE.se[[depth]][,amin.id[[depth]][[f]],f],
         length=0.03, angle=90, code=2, col=2, lwd=0.3)
  segments(x0=qmin.id[[depth]][[f]], y0=0, x1=qmin.id[[depth]][[f]], y1=PE.mu[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f], col=1, lty=2, lwd=0.5)
  abline(h=PE.mu[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f] + PE.se[[depth]][qmin.id[[depth]][[f]],amin.id[[depth]][[f]],f], col=2, lty=2, lwd=0.5)
  fit <- zeroslope(y=PE.mu[[depth]][,amin.id[[depth]][[f]],f], x=1:length(size.set), lag=2, span=0.40, degree=2, family="gaussian", minimum=TRUE)
  lines(x=1:length(size.set), y=fit$loess, type="l", col=4, lty=2, lwd=0.5)
  axis(side=1, at=1:length(size.set), labels=size.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=4, at=pretty(x=PE.mu[[depth]][,amin.id[[depth]][[f]],f]), col=2, col.axis=2, cex.axis=0.5)
  mtext(text="size (q)", line=1, side=1, col=1, outer=FALSE)  
  mtext(text="PE", line=1, side=4, col=2, outer=FALSE)  
  
  par(new=TRUE)
  plot(x=1:length(size.set), y=AUC.alpha.mu[,depth], type="b", axes=FALSE, xlab="", ylab="", col=1, cex=0.5, ylim=range(AUC.alpha.mu[,depth]))
  arrows(1:length(size.set), AUC.alpha.mu[,depth], 1:length(size.set), AUC.alpha.mu[,depth]-AUC.alpha.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  arrows(1:length(size.set), AUC.alpha.mu[,depth], 1:length(size.set), AUC.alpha.mu[,depth]+AUC.alpha.se[,depth],
         length=0.03, angle=90, code=2, col=1, lwd=0.3)
  abline(h=AUC.alpha.mu[aucmax.alpha.id[depth],depth] - AUC.alpha.se[aucmax.alpha.id[depth],depth], col=1, lty=2, lwd=0.5)
  segments(x0=qmin.id[[depth]][[f]], y0=0, x1=qmin.id[[depth]][[f]], y1=AUC.alpha.mu[qmin.id[[depth]][[f]],depth], col=4, lty=2, lwd=0.5)
  axis(side=1, at=1:length(size.set), labels=size.set, col=1, col.axis=1, cex.axis=0.5)
  axis(side=2, at=pretty(x=AUC.alpha.mu[,depth]), col=1, col.axis=1, cex.axis=0.5)
  mtext(text="size (q)", line=1, side=1, col=1, outer=FALSE)  
  mtext(text="AUC", line=1, side=2, col=1, outer=FALSE)  
  mtext(text=paste("depth = ", depth, sep=""), line=0, side=3, outer=FALSE)
}

dev.off()

rm(f, depth, my.file)

#==========================================================================================#
# FIGURE 08:
# Effect of Affinity Scores Fold Change
#==========================================================================================#
my.file <- file.path(HOME.path, NATM.path, "Figure08.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(ndepths,6), oma=c(0, 0, 0, 0), mar=c(2.5, 1.5, 1.5, 1.5), mgp=c(1, 0.5, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot(x=0, y=0, type="n", axes=FALSE, xlab="", ylab="")
  text(x=0, y=-0.5, label=paste("depth: ", depth, sep=""), cex=1)
  legend(x="top", inset=0, legend=paste("DE FDR ", fdr.set, sep=""), col=c(2,3,4), lty=1, lwd=1, pt.cex = 0.5, cex=0.5, xpd=TRUE)
  
  matplot(x=TPR.fc[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=TNR.fc[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.fc[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.fc[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.fc[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fc.set), labels=fc.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Fold Change", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
}

dev.off()

rm(depth)

#==========================================================================================#
# FIGURE 09:
# Effect of Frequency of Occurrence in Replications 
#==========================================================================================#
my.file <- file.path(HOME.path, NATM.path, "Figure09.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(ndepths,6), oma=c(0, 0, 0, 0), mar=c(2.5, 1.5, 1.5, 1.5), mgp=c(1, 0.5, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot(x=0, y=0, type="n", axes=FALSE, xlab="", ylab="")
  text(x=0, y=-0.5, label=paste("depth: ", depth, sep=""), cex=1)
  legend(x="top", inset=0, legend=paste("DE FDR ", fdr.set, sep=""), col=c(2,3,4), lty=1, lwd=1, pt.cex = 0.5, cex=0.5, xpd=TRUE)
  
  matplot(x=TPR.fr[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=TNR.fr[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.fr[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.fr[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.fr[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(fr.set), labels=fr.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Frequencies", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
}

dev.off()

rm(depth)

#==========================================================================================#
# FIGURE 10:
# Effect of Procedural Replications
#==========================================================================================#
my.file <- file.path(HOME.path, NATM.path, "Figure10.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(ndepths,6), oma=c(0, 0, 0, 0), mar=c(2.5, 1.5, 1.5, 1.5), mgp=c(1, 0.5, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot(x=0, y=0, type="n", axes=FALSE, xlab="", ylab="")
  text(x=0, y=-0.5, label=paste("depth: ", depth, sep=""), cex=1)
  legend(x="top", inset=0, legend=paste("DE FDR ", fdr.set, sep=""), col=c(2,3,4), lty=1, lwd=1, pt.cex = 0.5, cex=0.5, xpd=TRUE)
  
  matplot(x=TPR.rep[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, cex=0.5, side=1, col=1, outer=FALSE)  
  mtext(text="TPR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=TNR.rep[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="TNR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=AUC.rep[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="AUC", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FDR.rep[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FDR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)  
  
  matplot(x=FOR.rep[[depth]], type="l", axes=FALSE, ylim=c(0,1), xlab="", ylab="", col=2:4, lty=1, lwd=1)
  axis(side=1, at=1:length(rep.set), labels=rep.set, col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  axis(side=2, at=pretty(x=range(0, 1, na.rm=TRUE, finite=TRUE)), col=1, col.axis=1, cex.axis=0.5, cex.lab=0.5)
  mtext(text="Replications", line=1.5, side=1, cex=0.5, col=1, outer=FALSE)  
  mtext(text="FOR", line=1.5, side=2, cex=0.5, col=1, outer=FALSE)
}

dev.off()

rm(depth)

#==========================================================================================#
# FIGURE 11:
# Heatmap plots
#==========================================================================================#
f <- 1     # fixed DE FDR = 1e-2

my.file <- file.path(HOME.path, NATM.path, "Figure11.ps", fsep=.Platform$file.sep)
postscript(file=my.file, horizontal=FALSE, width=8.5, height=5, paper="special")

par(mfrow=c(2, 2), oma=c(0, 0, 0, 0), mar=c(3.5, 2.5, 2.5, 2.5), mgp=c(1.5, 0.25, 0), xpd=TRUE)
for (depth in depth.levels) {
  plot.heatmap(crna=CRNA.output[[depth]][[f]],
               degs=DEGs[[depth]][[f]],
               score.colors=c("royalblue", "white", "orangered"),
               express.colors=c("white", "orangered"),
               score.paletteLength=rep(51,ndepths)[depth],
               express.paletteLength=rep(100,ndepths)[depth],
               score.transformation="tanh",
               express.transformation="tanh",
               score.truncation=1.5,
               express.truncation=1.5,
               fontsize.row=3,
               fontsize.col=5,
               lwd=0.05,
               main=list(paste("Activity Scores Heatmap \n Depth = ", depth, sep=""),
                         paste("Expression Data Heatmap \n Depth = ", depth, sep="")))
}

dev.off()

rm(depth, f, my.file)

#==========================================================================================#
# FIGURE 12:
# Visualization of CRNA inferred network graphs at all depths (1-4), in the Ascending case of expression levels
# Exporting lists to Cytoscape
#==========================================================================================#
f <- 1     # fixed DE FDR = 1e-2

for (depth in depth.levels) {
  write.table(CRNA.output[[depth]][[f]]$edges, 
              file.path(HOME.path, NATM.path, paste("Figure12_edgelist,depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
              quote=FALSE, 
              row.names=FALSE, 
              sep="\t")
  nodelist <- data.frame("symbol"=CRNA.output[[depth]][[f]]$states$symbol, 
                         "node_color"="b", stringsAsFactors=FALSE)
  nodelist[which(nodelist$symbol %in% TDEGs[[depth]]),2] <- "d"
  nodelist[which(nodelist$symbol %in% TREGs[[depth]]),2] <- "r"
  write.table(nodelist, 
              file.path(HOME.path, NATM.path, paste("Figure12_nodelist,depth=", depth, ",fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
              quote=FALSE, 
              row.names=FALSE, 
              sep="\t")
}

rm(depth, f, nodelist)
