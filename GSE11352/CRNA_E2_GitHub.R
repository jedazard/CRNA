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

omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
DEGList <- read.csv("data_for_markdown/GSE11352_DEGList.csv",row.names = 1)
GSE11352_exprs_for_JTREE <- read.csv("data_for_markdown/GSE11352_scaled_exprs.csv",row.names = 1)[,-1]
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

seedSet <-  as.data.frame(matrix(0, nrow = length(V(RWRnet)$name),ncol = 1))
rownames(seedSet) <- V(RWRnet)$name
seedSet[which(rownames(seedSet) %in% rownames(miarray)),1] <- 1

RWRresult <- dRWR(g = RWRnet,normalise = "row",setSeeds = seedSet,restart = 0.75,verbose = F)

FC1.5_undirected_RWR0.75_nodelist <- seedSet
FC1.5_undirected_RWR0.75_nodelist[,2] <-  RWRresult[,1]
hist(FC1.5_undirected_RWR0.75_nodelist$V2,breaks = 80)

FC1.5_undirected_RWR0.75_selected_node <- rownames(FC1.5_undirected_RWR0.75_nodelist[FC1.5_undirected_RWR0.75_nodelist$V2 >quantile(FC1.5_undirected_RWR0.75_nodelist$V2,0.975),])
FC1.5_undirected_RWR0.75_DFS_candidates <-  FC1.5_undirected_RWR0.75_selected_node[which(!FC1.5_undirected_RWR0.75_selected_node %in% miarray$SYMBOL)]


significant_genes_graph <- induced.subgraph(omniNet,c(FC1.5_undirected_RWR0.75_selected_node,as.character(miarray$SYMBOL[which(miarray$SYMBOL %in% V(omniNet)$name)])))
significant_genes_graph <- igraph::simplify(graph = significant_genes_graph,remove.multiple = T,
                                            remove.loops = T,edge.attr.comb = "max")

print(paste("Here are", length(FC1.5_undirected_RWR0.75_selected_node), "significant genes in undirected RWR with 0.75 restart probablility.", length(which(FC1.5_undirected_RWR0.75_selected_node %in% rownames(miarray))),"of them are DEGs"))

#########################  DFS, then get GSE11352_DAG    ################################################################


FC1.5_undirected_RWR0.75_JTREE_candidates <- c()
FC1.5_undirected_RWR0.75_JTREE_edgeList <- as.data.frame(matrix(NA,nrow = 0,ncol = 4),row.names = c("from","to","sign","weight"))
for (i in 1:length(FC1.5_undirected_RWR0.75_DFS_candidates)) {
  r1 <- FC1.5_undirected_RWR0.75_DFS_candidates[i]
  r1_RWR0.75_DFS_graph <- RWR_dfs(r1,significant_genes_graph,miarray)
  if(length(E(r1_RWR0.75_DFS_graph)) > 0){
    FC1.5_undirected_RWR0.75_JTREE_edgeList <- rbind(FC1.5_undirected_RWR0.75_JTREE_edgeList,
                                                        igraph::as_data_frame(r1_RWR0.75_DFS_graph))
    FC1.5_undirected_RWR0.75_JTREE_candidates <- c(FC1.5_undirected_RWR0.75_JTREE_candidates,r1)
  }
}


RWR_FC2_DAG_graph <- graph_from_data_frame(unique(FC1.5_undirected_RWR0.75_JTREE_edgeList),directed = T)
RWR_FC2_DAG_graph <- igraph::simplify(graph = RWR_FC2_DAG_graph,remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "max")
RWR_FC2_DAG_edgeList <- igraph::as_data_frame(RWR_FC2_DAG_graph)

plot(RWR_FC2_DAG_graph, edge.arrow.size=.2, edge.curved=0,

     vertex.color="orange", vertex.frame.color="#555555",

     vertex.label=V(RWR_FC2_DAG_graph)$name, vertex.label.color="black",

     vertex.label.cex=.5,layout = layout_as_star(RWR_FC2_DAG_graph),vertex.size = 9)

DAG_test <- RWR_FC2_DAG_graph

if(!is.dag(DAG_test)){
adjmax <- as.matrix(as_adjacency_matrix(DAG_test,names = F))
removed_adjmax <- adj.remove.cycles(adjmat = adjmax ,maxlength = 10)$adjmat.acyclic

rownames(removed_adjmax) <- V(DAG_test)$name
colnames(removed_adjmax) <- V(DAG_test)$name

DAG_after <- graph_from_adjacency_matrix(removed_adjmax)

DAG_after_edgeList <- igraph::as_data_frame(DAG_after)

DAG_after_edgeList <- merge(x = DAG_after_edgeList,y = RWR_FC2_DAG_edgeList, by = c("from","to"))

RWR_FC2_DAG_graph <- graph_from_data_frame(unique(DAG_after_edgeList),directed = T)

RWR_FC2_DAG_graph <- igraph::simplify(graph = RWR_FC2_DAG_graph,remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "max")

GSE11352_DAG_graph <- (RWR_FC2_DAG_graph)
}


GSE11352_DAG_edgeList <- igraph::as_data_frame(GSE11352_DAG_graph)

GSE11352_DAG_nodeList <- as.data.frame(as.character(V(GSE11352_DAG_graph)$name))
GSE11352_DAG_nodeList$DEGs <- 0
colnames(GSE11352_DAG_nodeList) <- c("id","DEGs")
GSE11352_DAG_nodeList[which(as.character(GSE11352_DAG_nodeList$id) %in% as.character(miarray$SYMBOL)),2] <- 1

# input for JTREE

GSE11352_DAG_node_part <- as.data.frame(as.character(unique(c(GSE11352_DAG_edgeList$from,GSE11352_DAG_edgeList$to))))

GSE11352_DAG_node_part[,2] <- as.character(NA)

colnames(GSE11352_DAG_node_part) <- c("id","type")

GSE11352_DAG_node_part[is.na(GSE11352_DAG_node_part$type),2] <- "protein"

GSE11352_DAG_node_part <- GSE11352_DAG_node_part[,c(2,1)]

### edge part

GSE11352_DAG_edge_part <- GSE11352_DAG_edgeList[,c(1,2,3)]

GSE11352_DAG_edge_part[which(GSE11352_DAG_edge_part$sign == 1),3] <- "-a>"

GSE11352_DAG_edge_part[which(GSE11352_DAG_edge_part$sign == -1),3] <- "-a|"

### expression part
GSE11352_DAG_GSE11352_exprs_for_JTREE <- GSE11352_exprs_for_JTREE[as.character(GSE11352_DAG_node_part$id),]
rownames(GSE11352_DAG_GSE11352_exprs_for_JTREE) <- as.character(GSE11352_DAG_node_part$id)
GSE11352_DAG_GSE11352_for_JTREE_mRNA <- as.data.frame(t(GSE11352_DAG_GSE11352_exprs_for_JTREE))


GSE11352_DAG_permuted_exprs_for_JTREE <- as.data.frame(matrix(NA,ncol = 0,nrow = nrow(GSE11352_DAG_GSE11352_exprs_for_JTREE)))
rownames(GSE11352_DAG_permuted_exprs_for_JTREE) <- rownames(GSE11352_DAG_GSE11352_exprs_for_JTREE)
for (i in 1:ncol(GSE11352_DAG_GSE11352_exprs_for_JTREE)) {
  for (j in 1:1000) {
    GSE11352_DAG_permuted_exprs_for_JTREE <- cbind(GSE11352_DAG_permuted_exprs_for_JTREE,
                                                      GSE11352_DAG_GSE11352_exprs_for_JTREE[sample(1:nrow(GSE11352_DAG_GSE11352_exprs_for_JTREE), 
                                                                                                      size = nrow(GSE11352_DAG_GSE11352_exprs_for_JTREE), replace = T),
                                                                                               i])
  }    
}

GSE11352_DAG_permuted_exprs_for_JTREE <- cbind(id = rownames(GSE11352_DAG_GSE11352_exprs_for_JTREE),GSE11352_DAG_permuted_exprs_for_JTREE)

rownames(GSE11352_DAG_permuted_exprs_for_JTREE) <- as.character(GSE11352_DAG_permuted_exprs_for_JTREE$id)
colnames(GSE11352_DAG_permuted_exprs_for_JTREE) <- c("id",1:(ncol(GSE11352_exprs_for_JTREE)*1000))

GSE11352_DAG_permuted_for_JTREE_mRNA <- as.data.frame(t(GSE11352_DAG_permuted_exprs_for_JTREE))





##################################################################################################################

# Readin Junction Tree result

##################################################################################################################



RWR_GSE11352 <- read.table(paste0("data_for_markdown/Jan20_GSE11352_DAG_GSE11352.txt"),comment.char = ">")
RWR_permuted <- read.table(paste0("data_for_markdown/Jan20_GSE11352_DAG_permuted.txt"),comment.char = ">")
heatmap_candidates <- c()


# Load JTREE results
permuted_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_permuted)/(ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1),ncol = (ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1)))
rownames(permuted_activity_score) <- RWR_permuted$V1[1:(nrow(RWR_permuted)/(ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1))]
for (l in 1:(ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1)) {
  tails <- l*(nrow(RWR_permuted)/(ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1))
  heads <- tails-(nrow(RWR_permuted)/(ncol(GSE11352_DAG_permuted_exprs_for_JTREE)-1)-1)
  permuted_activity_score[,l] <- as.numeric(RWR_permuted[heads:tails,2])
}


GSE11352_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_GSE11352)/(ncol(GSE11352_DAG_GSE11352_exprs_for_JTREE)),
                                                ncol = (ncol(GSE11352_DAG_GSE11352_exprs_for_JTREE))))
rownames(GSE11352_activity_score) <- RWR_GSE11352$V1[1:(nrow(RWR_GSE11352)/(ncol(GSE11352_DAG_GSE11352_exprs_for_JTREE)))]
for (m in 1:(ncol(GSE11352_exprs_for_JTREE))) {
  tails <- m*nrow(RWR_GSE11352)/(ncol(GSE11352_exprs_for_JTREE))
  heads <- tails-(nrow(RWR_GSE11352)/(ncol(GSE11352_exprs_for_JTREE))-1)
  GSE11352_activity_score[,m] <- as.numeric(RWR_GSE11352[heads:tails,2])
}

colnames(GSE11352_activity_score) <- c("12hrs_Veh1","12hrs_Veh2","12hrs_Veh3","12hrs_Tx1","12hrs_Tx2","12hrs_Tx3",
                                       "24hrs_Veh1","24hrs_Veh2","24hrs_Veh3","24hrs_Tx1","24hrs_Tx2","24hrs_Tx3",
                                       "48hrs_Veh1","48hrs_Veh2","48hrs_Veh3","48hrs_Tx1","48hrs_Tx2","48hrs_Tx3")
GSE11352_activity_score <- GSE11352_activity_score[,c(1:3,7:9,13:15,4:6,10:12,16:18)]
permuted_activity_score <- permuted_activity_score[,c(1:3000,6001:9000,12001:15000,3001:6000,9001:12000,15001:18000)]


paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(GSE11352_activity_score), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(GSE11352_activity_score)/paletteLength, max(GSE11352_activity_score), length.out=floor(paletteLength/2)))

pheatmap(as.matrix(GSE11352_activity_score) ,fontsize_row=12,main = paste("RWR results"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)


hist(as.numeric(GSE11352_activity_score["ESR1",]))
hist(as.numeric(c(permuted_activity_score["ESR1",15000:18000])),breaks = 80)
ggplot(data = as.data.frame(t(permuted_activity_score[,15000:18000])), aes(x =ESR1)) +
  geom_histogram(color = "white", bins = 20,fill = "royalblue") +
  labs(title = "ESR1 null distribution",
       x = "activity index",y = "Frequency")

GSE11352_activity_quantile <- as.data.frame(matrix(0,nrow = nrow(GSE11352_activity_score),ncol = ncol(GSE11352_activity_score)))


for (i in 1:nrow(GSE11352_activity_quantile)) {
  for (j in 1:ncol(GSE11352_activity_quantile)) {
    null_distribution <- ecdf(as.numeric(permuted_activity_score[i,(1+(j-1)*1000):(1000*j)]))
    GSE11352_activity_quantile[i,j] <- null_distribution(GSE11352_activity_score[i,j])
  }
  
}

rownames(GSE11352_activity_quantile) <- rownames(GSE11352_activity_score)
colnames(GSE11352_activity_quantile) <- colnames(GSE11352_activity_score)

GSE11352_activity_significant <- GSE11352_activity_quantile

for (i in 1:nrow(GSE11352_activity_significant)) {
  GSE11352_activity_significant[i,which(GSE11352_activity_quantile[i,] >= 0.9)] <- 1
  GSE11352_activity_significant[i,which(GSE11352_activity_quantile[i,] <= 0.1)] <- -1
  GSE11352_activity_significant[i,which(GSE11352_activity_quantile[i,] > 0.1& GSE11352_activity_quantile[i,]<0.9)] <- 0
}

total_zero_genes <- rownames(GSE11352_activity_significant)[which(rowSums(abs(GSE11352_activity_significant)) == 0)]

GSE11352_activity_significant_only <- GSE11352_activity_significant

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(GSE11352_activity_significant_only), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(GSE11352_activity_significant_only)/paletteLength, max(GSE11352_activity_significant_only), length.out=floor(paletteLength/2)))

pheatmap(as.matrix(GSE11352_activity_significant_only) ,fontsize_row=12,main = paste("Sample Level Gene Activity Score for Causal Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)


phenotypes <- as.data.frame(matrix(0,nrow =nrow(GSE11352_activity_significant),ncol = 6 ),row.names = rownames(GSE11352_activity_significant))

for (i in 1:6) {
  phenotypes[,i] <- sign(rowSums(GSE11352_activity_significant[,(3*i-2):(3*i)]))
}

colnames(phenotypes) <- c("Veh-12hrs","Veh-24hrs","Veh-48hrs","Tx-12hrs","Tx-24hrs","Tx-48hrs")
# pdf("../Draft/figures/Jan20_E2_phenotype_heatmap.pdf",width = 4,height = 6)
paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(phenotypes), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(phenotypes)/paletteLength, max(phenotypes), length.out=floor(paletteLength/2)))
pheatmap(as.matrix(phenotypes) ,fontsize_row=12,breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)
# dev.off()




Tx24_majority_vote <- as.data.frame(sign(phenotypes[,4]-phenotypes[,2]) ,row.names = rownames(phenotypes))
colnames(Tx24_majority_vote) <- c("node_sign")
Tx24_majority_vote_graph <- induced.subgraph(GSE11352_DAG_graph,as.character(rownames(Tx24_majority_vote)[which(Tx24_majority_vote$node_sign != 0)]))
Tx24_majority_vote_graph <- induced.subgraph(Tx24_majority_vote_graph,names(which(degree(Tx24_majority_vote_graph,V(Tx24_majority_vote_graph)) >0)))

plot(Tx24_majority_vote_graph, edge.arrow.size=.2, edge.curved=0,

vertex.color="orange", vertex.frame.color="#555555",

vertex.label=V(Tx24_majority_vote_graph)$name, vertex.label.color="black",

vertex.label.cex=.5,layout = layout_nicely(Tx24_majority_vote_graph),vertex.size = 9)


Tx24_majority_vote_edgeList <- igraph::as_data_frame(Tx24_majority_vote_graph)

Tx24_majority_vote_nodeList <- Tx24_majority_vote
colnames(Tx24_majority_vote_nodeList) <- c("node_sign")

majority_vote <- as.data.frame(sign(rowSums(GSE11352_activity_significant[,10:18])-rowSums(GSE11352_activity_significant[,1:9])) ,row.names = rownames(GSE11352_activity_significant))
colnames(majority_vote) <- c("node_sign")
majority_vote_graph <- induced.subgraph(GSE11352_DAG_graph,as.character(rownames(majority_vote)[which(majority_vote$node_sign != 0)]))
majority_vote_graph <- induced.subgraph(majority_vote_graph,names(which(degree(majority_vote_graph,V(majority_vote_graph)) >0)))

plot(majority_vote_graph, edge.arrow.size=.2, edge.curved=0,

     vertex.color="orange", vertex.frame.color="#555555",

     vertex.label=V(majority_vote_graph)$name, vertex.label.color="black",

     vertex.label.cex=.5,layout = layout_nicely(majority_vote_graph),vertex.size = 9)


majority_vote_edgeList <- igraph::as_data_frame(majority_vote_graph)

majority_vote_nodeList <- majority_vote
colnames(majority_vote_nodeList) <- c("node_sign")

