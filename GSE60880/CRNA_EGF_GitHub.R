library(igraph)
library(dplyr)
library(data.table)
library(dnet)
library(pheatmap)
library(predictionet)

#########################################################################################
# DFS in RWR
#########################################################################################

RWR_dfs <- function(Root,Network,DEG_List){
  passedNodes <- c(Root)
  passedEdges <- c()
  regulatedDEG_count <- 0
  regulatedDEG <- c()
  ### Find r1
  possible_r1 <- unique(neighbors(Network,Root, "out")$name)
  possible_r1 <- possible_r1[!possible_r1 %in% passedNodes]      ### That's r1
  
  ###if Root has childs(r1)
  if (length(possible_r1)>0) {  ### Root has r1
    for (i in 1:length(possible_r1)) {
      r1 <- possible_r1[i]
      passedNodes <- c(passedNodes,r1)
      if (r1 %in% rownames(DEG_List)) {
        #    passedEdges <- c(passedEdges, get.edge.ids(Network,c(Root, r1)))
        regulatedDEG_count <- regulatedDEG_count + 1
        regulatedDEG <- c(regulatedDEG,r1)
      }else{
        possible_r2 <- unique(neighbors(Network,r1,"out")$name)
        possible_r2 <- possible_r2[!possible_r2 %in% passedNodes]      ### That's r2
        ###if r1 has childs(r2)
        if (length(possible_r2)>0) {  ### r1 has r2
          for (i in 1:length(possible_r2)) {
            r2 <- possible_r2[i]
            passedNodes <- c(passedNodes,r2)
            if (r2 %in% rownames(DEG_List)) {
              passedEdges <- c(passedEdges, get.edge.ids(Network,c(Root, r1))
                               #,get.edge.ids(Network,c(r1, r2))
              )
              # passedEdges <- c(passedEdges, get.edge.ids(Network,c(Root, r1)))
              regulatedDEG_count <- regulatedDEG_count + 1
              regulatedDEG <- c(regulatedDEG,r2)
            }else{
              possible_r3 <- unique(neighbors(Network,r2,"out")$name)
              possible_r3 <- possible_r3[!possible_r3 %in% passedNodes]      ### That's r3
              ###if r2 has childs(r3)
              if (length(possible_r3)>0) {  ### r2 has r3
                for (i in 1:length(possible_r3)) {
                  r3 <- possible_r3[i]
                  passedNodes <- c(passedNodes,r3)
                  if (r3 %in% rownames(DEG_List)) {
                    passedEdges <- c(passedEdges,get.edge.ids(Network,c(Root, r1)),
                                     get.edge.ids(Network,c(r1, r2)))
                    regulatedDEG_count <- regulatedDEG_count + 1
                    regulatedDEG <- c(regulatedDEG,r3)
                  }else{
                    possible_r4 <- unique(neighbors(Network,r3,"out")$name)
                    possible_r4 <- possible_r4[!possible_r4 %in% passedNodes]      ### That's r3
                    ###if r2 has childs(r3)
                    if (length(possible_r4)>0) {  ### r2 has r3
                      for (i in 1:length(possible_r4)) {
                        r4 <- possible_r4[i]
                        passedNodes <- c(passedNodes,r4)
                        if (r4 %in% rownames(DEG_List)) {
                          passedEdges <- c(passedEdges,get.edge.ids(Network,c(Root, r1)),
                                           get.edge.ids(Network,c(r1, r2)),
                                           get.edge.ids(Network,c(r2, r3)))
                          regulatedDEG_count <- regulatedDEG_count + 1
                          regulatedDEG <- c(regulatedDEG,r4)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } 
    }    
  }   
  message(paste(Root, "regulates", regulatedDEG_count,"DEGs. \t"))
  return(subgraph.edges(Network,eids = passedEdges))
}






#########################################################################################
# Loading data
#########################################################################################

# scaled_GSE60880_exprs <- read.table("data_for_markdown/scaled_exprs_limma.tab",sep = "\t")
omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
DEGList <- read.csv("data_for_markdown/GSE60880_DEGList.csv",row.names = 1)
GSE60880_exprs_for_JTREE <- read.csv("data_for_markdown/GSE60880_exprs_for_JTREE.csv",row.names = 1)
DEGList <- unique(tidyr::separate_rows(DEGList,SYMBOL,sep = " /// "))
DEGList <- DEGList %>% group_by(SYMBOL) %>%
  summarize(logFC = mean(logFC),
            AveExpr = mean(AveExpr),
            t = mean(t),
            P.Value = mean(P.Value),
            adj.P.Val = mean(adj.P.Val),
            B = mean(B))

#########################################################################################
# Building graph
#########################################################################################

omniNet <- graph_from_data_frame(d = omnipath,directed = T)

#########################################################################################
# RWR
#########################################################################################

RWRnet <- as.undirected(omniNet)
miarray <- as.data.frame(DEGList[which(abs(DEGList$logFC) >= log2(2)),])
rownames(miarray) <- miarray$SYMBOL
noChanged <- as.data.frame(DEGList[which(abs(DEGList$logFC) < log2(1.2)),])
rownames(noChanged) <- noChanged$SYMBOL

non_changed_data <- GSE60880_exprs_for_JTREE[which(rownames(GSE60880_exprs_for_JTREE) %in% rownames(noChanged)),]


seedSet <-  as.data.frame(matrix(0, nrow = length(V(RWRnet)$name),ncol = 1))
rownames(seedSet) <- V(RWRnet)$name
seedSet[which(rownames(seedSet) %in% rownames(miarray)),1] <- 1

RWRresult <- dRWR(g = RWRnet,normalise = "row",setSeeds = seedSet,restart = 0.75,verbose = F)

FC2_undirected_RWR0.65_nodelist <- seedSet
FC2_undirected_RWR0.65_nodelist[,2] <-  RWRresult[,1]
hist(FC2_undirected_RWR0.65_nodelist$V2,breaks = 80)
FC2_undirected_RWR0.65_selected_node <- rownames(FC2_undirected_RWR0.65_nodelist[FC2_undirected_RWR0.65_nodelist$V2 >quantile(FC2_undirected_RWR0.65_nodelist$V2,0.975),])
FC2_undirected_RWR0.65_DFS_candidates <-  FC2_undirected_RWR0.65_selected_node[which(!FC2_undirected_RWR0.65_selected_node %in% miarray$SYMBOL)]
significant_genes_graph <- induced.subgraph(omniNet,c(FC2_undirected_RWR0.65_selected_node,as.character(miarray$SYMBOL[which(miarray$SYMBOL %in% V(omniNet)$name)])))
significant_genes_graph <- igraph::simplify(graph = significant_genes_graph,remove.multiple = T,
                                            remove.loops = T,edge.attr.comb = "max")

print(paste("Here are", length(FC2_undirected_RWR0.65_selected_node), "significant genes in undirected RWR with 0.65 restart probablility.", length(which(FC2_undirected_RWR0.65_selected_node %in% rownames(miarray))),"of them are DEGs"))

#########################  DFS, then get FC2_MST_V8    ################################################################


FC2_undirected_RWR0.65_JTREE_candidates <- c()
FC2_undirected_RWR0.65_JTREE_edgeList <- as.data.frame(matrix(NA,nrow = 0,ncol = 4),row.names = c("from","to","sign","weight"))
for (i in 1:length(FC2_undirected_RWR0.65_DFS_candidates)) {
  r1 <- FC2_undirected_RWR0.65_DFS_candidates[i]
  r1_RWR0.65_DFS_graph <- RWR_dfs(r1,significant_genes_graph,miarray)
  if(length(E(r1_RWR0.65_DFS_graph)) > 0){
    FC2_undirected_RWR0.65_JTREE_edgeList <- rbind(FC2_undirected_RWR0.65_JTREE_edgeList,
                                                      igraph::as_data_frame(r1_RWR0.65_DFS_graph))
    FC2_undirected_RWR0.65_JTREE_candidates <- c(FC2_undirected_RWR0.65_JTREE_candidates,r1)
    }
}


RWR_FC2_DAG_graph <- graph_from_data_frame(unique(FC2_undirected_RWR0.65_JTREE_edgeList),directed = T)
RWR_FC2_DAG_graph <- igraph::simplify(graph = RWR_FC2_DAG_graph,remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "max")
RWR_FC2_DAG_edgeList <- igraph::as_data_frame(RWR_FC2_DAG_graph)

plot(RWR_FC2_DAG_graph, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(RWR_FC2_DAG_graph)$name, vertex.label.color="black",
     
     vertex.label.cex=.5,layout = layout_nicely(RWR_FC2_DAG_graph),vertex.size = 9)

DAG_test <- RWR_FC2_DAG_graph

is.dag(DAG_test)
adjmax <- as.matrix(as_adjacency_matrix(DAG_test,names = F))
removed_adjmax <- adj.remove.cycles(adjmat = adjmax ,maxlength = 10)$adjmat.acyclic

rownames(removed_adjmax) <- V(DAG_test)$name
colnames(removed_adjmax) <- V(DAG_test)$name

DAG_after <- graph_from_adjacency_matrix(removed_adjmax)

is.dag(DAG_after)

DAG_after_edgeList <- igraph::as_data_frame(DAG_after)

DAG_after_edgeList <- merge(x = DAG_after_edgeList,y = RWR_FC2_DAG_edgeList, by = c("from","to"))

RWR_FC2_DAG_graph <- graph_from_data_frame(unique(DAG_after_edgeList),directed = T)

RWR_FC2_DAG_graph <- igraph::simplify(graph = RWR_FC2_DAG_graph,remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "max")

is.dag(RWR_FC2_DAG_graph)
GSE60880_DAG_graph <- (RWR_FC2_DAG_graph)
GSE60880_DAG_edgeList <- igraph::as_data_frame(GSE60880_DAG_graph)

GSE60880_DAG_nodeList <- as.data.frame(as.character(V(GSE60880_DAG_graph)$name))
GSE60880_DAG_nodeList$DEGs <- 0
colnames(GSE60880_DAG_nodeList) <- c("id","DEGs")
GSE60880_DAG_nodeList[which(as.character(GSE60880_DAG_nodeList$id) %in% as.character(miarray$SYMBOL)),2] <- 1


indegree <- as.data.frame(table(GSE60880_DAG_edgeList$to))
indegree <- indegree[which(indegree$Freq > 6 ),]
row.names(indegree) <- indegree$Var1
delete_edge <- c()

if (nrow(indegree)>0){
  for (i in 1:nrow(indegree)) {
    issue <- as.character(indegree$Var1[i])
    issue_edges <- GSE60880_DAG_edgeList[which(GSE60880_DAG_edgeList$to == issue),]
    issue_edges <- issue_edges[which(issue_edges$from %in% as.character(names(which(degree(GSE60880_DAG_graph,neighbors(GSE60880_DAG_graph,v =issue,mode = "in")$name,mode = "in") == 0 )))),]
    delete_edge <- c(delete_edge,as.numeric(rownames(issue_edges)))
    if (indegree[which(indegree$Var1 == issue),"Freq"] -nrow(issue_edges) > 8 ) {
      issue_edges <- GSE60880_DAG_edgeList[which(GSE60880_DAG_edgeList$to == issue),]
      issue_edges <- issue_edges[order(issue_edges$weight,decreasing = T),]
      delete_edge <- c(delete_edge,as.numeric(rownames(issue_edges[-(1:8),])))
    }
  }
  GSE60880_DAG_edgeList <- GSE60880_DAG_edgeList[-unique(delete_edge),]
}


indegree <- as.data.frame(table(GSE60880_DAG_edgeList$to))

outdegree <- as.data.frame(table(GSE60880_DAG_edgeList$from))
outdegree <- outdegree[which(outdegree$Freq >6),]
row.names(outdegree) <- outdegree$Var1
delete_edge <- c()
if (nrow(outdegree)>0) {
  for (i in 1:nrow(outdegree)) {
    issue <- as.character(outdegree$Var1[i])
    issue_edges <- GSE60880_DAG_edgeList[which(GSE60880_DAG_edgeList$from == issue),]
    issue_edges <- issue_edges[which(issue_edges$to %in% as.character(names(which(degree(GSE60880_DAG_graph,neighbors(GSE60880_DAG_graph,v =issue,mode = "out")$name,mode = "out") == 0 )))),]
    delete_edge <- c(delete_edge,as.numeric(rownames(issue_edges)))
    if (outdegree[which(outdegree$Var1 == issue),"Freq"] -nrow(issue_edges) > 6) {
      issue_edges <- GSE60880_DAG_edgeList[which(GSE60880_DAG_edgeList$from == issue),]
      issue_edges <- issue_edges[order(issue_edges$weight,decreasing = T),]
      delete_edge <- c(delete_edge,as.numeric(rownames(issue_edges[-(1:6),])))
    }
  }
  GSE60880_DAG_edgeList <- GSE60880_DAG_edgeList[-unique(delete_edge),]
  
}

outdegree <- as.data.frame(table(GSE60880_DAG_edgeList$from))


# input for JTREE

GSE60880_DAG_node_part <- as.data.frame(as.character(unique(c(GSE60880_DAG_edgeList$from,GSE60880_DAG_edgeList$to))))
GSE60880_DAG_node_part[,2] <- as.character(NA)
colnames(GSE60880_DAG_node_part) <- c("id","type")
GSE60880_DAG_node_part[is.na(GSE60880_DAG_node_part$type),2] <- "protein"
GSE60880_DAG_node_part <- GSE60880_DAG_node_part[,c(2,1)]

### edge part

GSE60880_DAG_edge_part <- GSE60880_DAG_edgeList[,c(1,2,3)]
GSE60880_DAG_edge_part[which(GSE60880_DAG_edge_part$sign == 1),3] <- "-a>"
GSE60880_DAG_edge_part[which(GSE60880_DAG_edge_part$sign == -1),3] <- "-a|"

### expression part
GSE60880_DAG_GSE60880_exprs_for_JTREE <- GSE60880_exprs_for_JTREE[as.character(GSE60880_DAG_node_part$id),]
rownames(GSE60880_DAG_GSE60880_exprs_for_JTREE) <- as.character(GSE60880_DAG_node_part$id)
GSE60880_DAG_GSE60880_for_JTREE_mRNA <- as.data.frame(t(GSE60880_DAG_GSE60880_exprs_for_JTREE))
GSE60880_DAG_GSE60880_for_JTREE_genome <- GSE60880_DAG_GSE60880_for_JTREE_mRNA
GSE60880_DAG_GSE60880_for_JTREE_genome[,] <- NA


GSE60880_DAG_permuted_exprs_for_JTREE <- as.data.frame(matrix(NA,ncol = 0,nrow = nrow(GSE60880_DAG_GSE60880_exprs_for_JTREE)))
rownames(GSE60880_DAG_permuted_exprs_for_JTREE) <- rownames(GSE60880_DAG_GSE60880_exprs_for_JTREE)
for (i in 1:ncol(GSE60880_DAG_GSE60880_exprs_for_JTREE)) {
  for (j in 1:1000) {
    GSE60880_DAG_permuted_exprs_for_JTREE <- cbind(GSE60880_DAG_permuted_exprs_for_JTREE,
                                                      GSE60880_DAG_GSE60880_exprs_for_JTREE[sample(1:nrow(GSE60880_DAG_GSE60880_exprs_for_JTREE), 
                                                                                                      size = nrow(GSE60880_DAG_GSE60880_exprs_for_JTREE), replace = T),
                                                                                               i])
  }    
}

GSE60880_DAG_permuted_exprs_for_JTREE <- cbind(id = rownames(GSE60880_DAG_GSE60880_exprs_for_JTREE),GSE60880_DAG_permuted_exprs_for_JTREE)

rownames(GSE60880_DAG_permuted_exprs_for_JTREE) <- as.character(GSE60880_DAG_permuted_exprs_for_JTREE$id)
colnames(GSE60880_DAG_permuted_exprs_for_JTREE) <- c("id",1:(ncol(GSE60880_exprs_for_JTREE)*1000))

GSE60880_DAG_permuted_for_JTREE_mRNA <- as.data.frame(t(GSE60880_DAG_permuted_exprs_for_JTREE))

GSE60880_DAG_permuted_for_JTREE_genome <- GSE60880_DAG_permuted_for_JTREE_mRNA

GSE60880_DAG_permuted_for_JTREE_genome[-1,] <- NA



############################################################################################################

# Read in Junction Tree results

############################################################################################################

RWR_GSE60880 <- read.table(paste0("data_for_markdown/Jan20_GSE60880_DAG_GSE60880.txt"),comment.char = ">")
RWR_permuted <- read.table(paste0("data_for_markdown/Jan20_GSE60880_DAG_permuted.txt"),comment.char = ">")
heatmap_candidates <- c()


# Load JTREE results
permuted_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_permuted)/(ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1),ncol = (ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1)))
rownames(permuted_activity_score) <- RWR_permuted$V1[1:(nrow(RWR_permuted)/(ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1))]
for (l in 1:(ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1)) {
  tails <- l*(nrow(RWR_permuted)/(ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1))
  heads <- tails-(nrow(RWR_permuted)/(ncol(GSE60880_DAG_permuted_exprs_for_JTREE)-1)-1)
  permuted_activity_score[,l] <- as.numeric(RWR_permuted[heads:tails,2])
}


GSE60880_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_GSE60880)/(ncol(GSE60880_DAG_GSE60880_exprs_for_JTREE)),
                                                ncol = (ncol(GSE60880_DAG_GSE60880_exprs_for_JTREE))))
rownames(GSE60880_activity_score) <- RWR_GSE60880$V1[1:(nrow(RWR_GSE60880)/(ncol(GSE60880_DAG_GSE60880_exprs_for_JTREE)))]
for (m in 1:(ncol(GSE60880_exprs_for_JTREE))) {
  tails <- m*nrow(RWR_GSE60880)/(ncol(GSE60880_exprs_for_JTREE))
  heads <- tails-(nrow(RWR_GSE60880)/(ncol(GSE60880_exprs_for_JTREE))-1)
  GSE60880_activity_score[,m] <- as.numeric(RWR_GSE60880[heads:tails,2])
}
colnames(GSE60880_activity_score) <- colnames(GSE60880_DAG_GSE60880_exprs_for_JTREE)

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(GSE60880_activity_score), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(GSE60880_activity_score)/paletteLength, max(GSE60880_activity_score), length.out=floor(paletteLength/2)))

annotation <- data.frame(Gene_Type = factor(rownames(GSE60880_activity_score) %in% rownames(miarray)))
rownames(annotation) <- rownames(GSE60880_activity_score) # check out the row names of annotation
pheatmap(as.matrix(GSE60880_activity_score) ,fontsize_row=5,main = paste("RWR results"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

GSE60880_activity_quantile <- as.data.frame(matrix(0,nrow = nrow(GSE60880_activity_score),ncol = ncol(GSE60880_activity_score)))


for (i in 1:nrow(GSE60880_activity_quantile)) {
  for (j in 1:ncol(GSE60880_activity_quantile)) {
    null_distribution <- ecdf(as.numeric(permuted_activity_score[i,(1+(j-1)*1000):(1000*j)]))
    GSE60880_activity_quantile[i,j] <- null_distribution(GSE60880_activity_score[i,j])
  }
  
}

length(which(GSE60880_activity_quantile<=.2|GSE60880_activity_quantile>=.8))

rownames(GSE60880_activity_quantile) <- rownames(GSE60880_activity_score)
colnames(GSE60880_activity_quantile) <- colnames(GSE60880_activity_score)

GSE60880_activity_significant <- GSE60880_activity_quantile

for (i in 1:nrow(GSE60880_activity_significant)) {
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] >= 0.85)] <- 1
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] <= 0.15)] <- -1
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] > 0.15 & GSE60880_activity_quantile[i,]<0.85)] <- 0
}

total_zero_genes <- rownames(GSE60880_activity_significant)[which(rowSums(abs(GSE60880_activity_significant)) == 0)]

GSE60880_activity_significant_only <- GSE60880_activity_significant

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(GSE60880_activity_significant_only), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(GSE60880_activity_significant_only)/paletteLength, max(GSE60880_activity_significant_only), length.out=floor(paletteLength/2)))

pheatmap(as.matrix(GSE60880_activity_significant_only) ,fontsize_row=5,main = paste("Sample Level Gene Activity Score for Causal Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)


phenotypes <- as.data.frame(matrix(0,nrow =nrow(GSE60880_activity_significant),ncol = 8 ),row.names = rownames(GSE60880_activity_significant))

for (i in 1:8) {
  if (i < 4) {
    phenotypes[,i] <- sign(rowSums(GSE60880_activity_significant[,(6*i-5):(6*i)]))
  }else if (i == 4) {
    phenotypes[,i] <- sign(rowSums(GSE60880_activity_significant[,(6*i-5):(6*i-1)]))
  }else{
    phenotypes[,i] <- sign(rowSums(GSE60880_activity_significant[,(6*i-6):(6*i-1)]))
  }
}

colnames(phenotypes) <- c("Ctrl-0.5hrs","Ctrl-1hrs","Ctrl-2hrs","Ctrl-8hrs","Tx-0.5hrs","Tx-1hrs","Tx-2hrs","Tx-8hrs")

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(phenotypes), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(phenotypes)/paletteLength, max(phenotypes), length.out=floor(paletteLength/2)))
pheatmap(as.matrix(phenotypes) ,fontsize_row=12,main = paste("8 hours treated group Gene Activity State for Causal Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

pheatmap(as.matrix(phenotypes[which(phenotypes$`Tx-8hrs` != phenotypes$`Ctrl-8hrs` ),c(4,8)]) ,fontsize_row=12,main = paste("8 hours treated group Gene Activity State for Causal Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

Tx48_majority_vote <- as.data.frame(sign(phenotypes[,8]-phenotypes[,4]) ,row.names = rownames(phenotypes))
colnames(Tx48_majority_vote) <- c("node_sign")
Tx48_majority_vote_graph <- induced.subgraph(GSE60880_DAG_graph,as.character(rownames(Tx48_majority_vote)[which(Tx48_majority_vote$node_sign != 0)]))
Tx48_majority_vote_graph <- induced.subgraph(Tx48_majority_vote_graph,names(which(degree(Tx48_majority_vote_graph,V(Tx48_majority_vote_graph)) >0)))

plot(Tx48_majority_vote_graph, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(Tx48_majority_vote_graph)$name, vertex.label.color="black",
     
     vertex.label.cex=.5,layout = layout_nicely(Tx48_majority_vote_graph),vertex.size = 9)


majority_vote <- cbind(Ctrl = as.data.frame(sign(rowSums(GSE60880_activity_significant[,1:23]))),Tx = as.data.frame(sign(rowSums(GSE60880_activity_significant[,24:47]))))

colnames(majority_vote) <- c("Ctrl","Tx")

majority_vote <- majority_vote[which(majority_vote$Tx-majority_vote$Ctrl !=0),]

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(majority_vote), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(majority_vote)/paletteLength, max(majority_vote), length.out=floor(paletteLength/2)))

pheatmap(as.matrix(majority_vote) ,fontsize_row=12,main = paste("Treatment Level Gene Activity State for Causal Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,labels_col = c("Control","Treatment"))

majority_Vote_graph <- induced.subgraph(GSE60880_DAG_graph,as.character(rownames(majority_vote)))
majority_Vote_graph <- induced.subgraph(majority_Vote_graph,names(which(degree(majority_Vote_graph,V(majority_Vote_graph)) >0)))

majority_Vote_edgeList <- igraph::as_data_frame(majority_Vote_graph)

majority_Vote_nodeList <- cbind(rownames(majority_vote),(as.data.frame(sign(majority_vote$Tx - majority_vote$Ctrl),row.names = rownames(majority_vote))))
majority_Vote_nodeList <- majority_Vote_nodeList[majority_Vote_nodeList$`sign(majority_vote$Tx - majority_vote$Ctrl)` != 0,]
colnames(majority_Vote_nodeList) <- c("id","node_sign")

