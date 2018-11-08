#########################################################################################
# Loading Library
#########################################################################################
library(igraph)
library(dplyr)
library(data.table)
library(dnet)
library(tidyverse)
library(pheatmap)
library(predictionet)

#########################################################################################
# Function
#########################################################################################
RWR_dfs <- function(r1,graph,DEG){
  passedNodes <- c(r1)
  passedEdges <- c()
  regulatedDEG <- 0
  ### Find r2
  possible_r2 <- unique(neighbors(graph,r1, "out")$name)
  possible_r2 <- possible_r2[!possible_r2 %in% passedNodes]      ### That's r2
  
  ###if r1 has childs(r2)
  if (length(possible_r2)>0) {  ### r1 has r2
    for (i in 1:length(possible_r2)) {
      r2 <- possible_r2[i]
      passedNodes <- c(passedNodes,r2)
      if (r2 %in% rownames(DEG)) {
        passedEdges <- c(passedEdges, get.edge.ids(graph,c(r1, r2)))
        regulatedDEG <- regulatedDEG + 1
      }else{
        possible_r3 <- unique(neighbors(graph,r2,"out")$name)
        possible_r3 <- possible_r3[!possible_r3 %in% passedNodes]      ### That's r3
        ###if r2 has childs(r3)
        if (length(possible_r3)>0) {  ### r2 has r3
          for (i in 1:length(possible_r3)) {
            r3 <- possible_r3[i]
            passedNodes <- c(passedNodes,r3)
            if (r3 %in% rownames(DEG)) {
              passedEdges <- c(passedEdges, get.edge.ids(graph,c(r1, r2)),
                               get.edge.ids(graph,c(r2, r3)))
              regulatedDEG <- regulatedDEG + 1
            }else{
              possible_r4 <- unique(neighbors(graph,r3,"out")$name)
              possible_r4 <- possible_r4[!possible_r4 %in% passedNodes]      ### That's r4
              ###if r3 has childs(r4)
              if (length(possible_r4)>0) {  ### r3 has r4
                for (i in 1:length(possible_r4)) {
                  r4 <- possible_r4[i]
                  passedNodes <- c(passedNodes,r4)
                  if (r4 %in% rownames(DEG)) {
                    regulatedDEG <- regulatedDEG + 1
                  }
                  passedEdges <- c(passedEdges,get.edge.ids(graph,c(r1, r2)),
                                   get.edge.ids(graph,c(r2, r3)),
                                   get.edge.ids(graph,c(r3, r4)))
                }
              }
            }
          }
        }
      } 
    }    
  }   
  message(paste(r1, "regulates", regulatedDEG,"DEGs"))
  return(subgraph.edges(graph,eids = passedEdges))
}

#########################################################################################
# Loading data
#########################################################################################
# scaled_GSE11352_exprs <- read.table("data_for_markdown/scaled_exprs_limma.tab",sep = "\t")
omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
GSE11352_DEGList <- read.csv("data_for_markdown/GSE11352_GSE11352_DEGList.csv",row.names = 1)
GSE11352_exprs_for_JTREE <- read.csv("data_for_markdown/GSE11352_scaled_exprs.csv",row.names = 1)[,-1]
omniNet <- graph_from_data_frame(d = omnipath,directed = T)
miarray <- as.data.frame(DEGList[which(abs(DEGList$logFC) >= log2(2)),])

#########################################################################################
# RWR
########################################################################################
RWRnet <- as.undirected(omniNet)
rownames(miarray) <- miarray$SYMBOL
seedSet <-  as.data.frame(matrix(0, nrow = length(V(RWRnet)$name),ncol = 1))
rownames(seedSet) <- V(RWRnet)$name
seedSet[which(rownames(seedSet) %in% rownames(miarray)),1] <- 1
RWR_probability_data <-  as.data.frame(matrix(0, nrow = length(V(RWRnet)$name),ncol = 1000))
rownames(RWR_probability_data) <- V(RWRnet)$name
for (i in 1:1000) {
  RWRresult <- dRWR(g = RWRnet,normalise = "row",setSeeds = seedSet,restart = 0.75,verbose = F)
  RWR_probability_data[,i] <-  RWRresult[,1]
}

FC2_undirected_RWR0.75_nodelist <- seedSet
FC2_undirected_RWR0.75_nodelist[,2] <- rowMeans(RWR_probability_data)
FC2_undirected_RWR0.75_selected_node <- rownames(FC2_undirected_RWR0.75_nodelist[FC2_undirected_RWR0.75_nodelist$V2 >quantile(FC2_undirected_RWR0.75_nodelist$V2,0.975),])
FC2_undirected_RWR0.75_DFS_candidates <-  FC2_undirected_RWR0.75_selected_node[which(!FC2_undirected_RWR0.75_selected_node %in% miarray$SYMBOL)]

ggplot(data = FC2_undirected_RWR0.75_nodelist, aes(x =V8)) +
  geom_histogram(color = "white", bins = 20,fill = "orangered") +
  labs(title = "Distribution of RWR probability",
       x = "probability",y = "Frequency")

significant_genes_graph <- igraph::simplify(graph = induced.subgraph(omniNet,FC2_undirected_RWR0.75_selected_node),remove.multiple = T,
                                            remove.loops = T,edge.attr.comb = "mean")
print(paste("Here are", length(FC2_undirected_RWR0.75_selected_node), "significant genes in undirected RWR with 0.75 restart probablility.", length(which(FC2_undirected_RWR0.75_selected_node %in% rownames(miarray))),"of them are DEGs"))

#########################  DFS    ################################################################
FC2_undirected_RWR0.75_JTREE_candidates <- c()
FC2_undirected_RWR0.75_JTREE_edgeList <- as.data.frame(matrix(NA,nrow = 0,ncol = 4),row.names = c("from","to","sign","weight"))
for (i in 1:length(FC2_undirected_RWR0.75_DFS_candidates)) {
  r1 <- FC2_undirected_RWR0.75_DFS_candidates[i]
  r1_RWR0.75_DFS_graph <- RWR_dfs(r1,significant_genes_graph,miarray)
  if(length(E(r1_RWR0.75_DFS_graph)) > 0){
    FC2_undirected_RWR0.75_JTREE_edgeList <- unique(rbind(FC2_undirected_RWR0.75_JTREE_edgeList,
                                                    igraph::as_data_frame(r1_RWR0.75_DFS_graph)))
    FC2_undirected_RWR0.75_JTREE_candidates <- c(FC2_undirected_RWR0.75_JTREE_candidates,r1)
  }
}

RWR_FC2_DFS_graph <- igraph::simplify(graph = graph_from_data_frame(unique(FC2_undirected_RWR0.75_JTREE_edgeList),directed = T),remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "mean")

######################### then get DAG    ################################################################
DAG_test <- RWR_FC2_DFS_graph
if (is.dag(DAG_test)) {
  acyclic_adjmax <- adj.remove.cycles(adjmat = as.matrix(as_adjacency_matrix(DAG_test,names = F)) ,maxlength = 10)$adjmat.acyclic
  rownames(acyclic_adjmax) <- V(DAG_test)$name
  colnames(acyclic_adjmax) <- V(DAG_test)$name
  DAG_after <- graph_from_adjacency_matrix(acyclic_adjmax)
}  
is.dag(DAG_after)
DAG_after_edgeList <- merge(x = igraph::as_data_frame(DAG_after),y = RWR_FC2_DAG_edgeList, by = c("from","to"))
RWR_FC2_DAG_graph <- igraph::simplify(graph = graph_from_data_frame(unique(DAG_after_edgeList),directed = T),remove.multiple = T,
                                      remove.loops = T,edge.attr.comb = "mean")
# RWR_FC2_DAG_edgeList <- igraph::as_data_frame(RWR_FC2_DAG_graph)
# write.csv(RWR_FC2_DAG_edgeList,"../Draft/materials/V23_DAG_edgeList.csv",quote = F)
# RWR_FC2_DAG_nodeList <- as.data.frame(as.character(V(RWR_FC2_DAG_graph)$name))
# RWR_FC2_DAG_nodeList$DEGs <- 0
# colnames(RWR_FC2_DAG_nodeList) <- c("id","DEGs")
# RWR_FC2_DAG_nodeList[which(as.character(RWR_FC2_DAG_nodeList$id) %in% as.character(miarray$SYMBOL)),2] <- 1
# write.csv(RWR_FC2_DAG_nodeList,"../Draft/materials/V23_DAG_nodeList.csv",quote = F)

# input for JTREE
#Nodepart
RWR_FC2_DAG_node_part <- as.data.frame(as.character(unique(c(RWR_FC2_DAG_edgeList$from,RWR_FC2_DAG_edgeList$to))))
RWR_FC2_DAG_node_part[,2] <- "protein"
colnames(RWR_FC2_DAG_node_part) <- c("id","type")
RWR_FC2_DAG_node_part <- RWR_FC2_DAG_node_part[,c(2,1)]
RWR_FC2_DAG_edge_part <- RWR_FC2_DAG_edgeList[,c(1,2,3)]

# Edge Part
RWR_FC2_DAG_edge_part[which(RWR_FC2_DAG_edge_part$sign == 1),3] <- "-t>"
RWR_FC2_DAG_edge_part[which(RWR_FC2_DAG_edge_part$sign == -1),3] <- "-t|"
#write.table(RWR_FC2_DAG_node_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t")
#write.table(RWR_FC2_DAG_edge_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t",append = T)

### expression part
RWR_FC2_DAG_GSE11352_exprs_for_JTREE <- GSE11352_exprs_for_JTREE[as.character(RWR_FC2_DAG_node_part$id),]
rownames(RWR_FC2_DAG_GSE11352_exprs_for_JTREE) <- as.character(RWR_FC2_DAG_node_part$id)
RWR_FC2_DAG_GSE11352_for_JTREE_mRNA <- as.data.frame(t(RWR_FC2_DAG_GSE11352_exprs_for_JTREE))
RWR_FC2_DAG_GSE11352_for_JTREE_genome <- RWR_FC2_DAG_GSE11352_for_JTREE_mRNA
RWR_FC2_DAG_GSE11352_for_JTREE_genome <- NA
#write.table(rbind(id = colnames(RWR_FC2_DAG_GSE11352_for_JTREE_genome),RWR_FC2_DAG_GSE11352_for_JTREE_genome),"output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_GSE11352_genome.tab",sep = "\t",quote = FALSE,col.names = FALSE)
#write.table(rbind(id = colnames(RWR_FC2_DAG_GSE11352_for_JTREE_mRNA),RWR_FC2_DAG_GSE11352_for_JTREE_mRNA),"output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_GSE11352_mRNA.tab",sep = "\t",quote = FALSE,col.names = FALSE)

# Permutating samples
RWR_FC2_DAG_permuted_exprs_for_JTREE <- as.data.frame(matrix(NA,ncol = 0,nrow = nrow(RWR_FC2_DAG_GSE11352_exprs_for_JTREE)))
rownames(RWR_FC2_DAG_permuted_exprs_for_JTREE) <- rownames(RWR_FC2_DAG_GSE11352_exprs_for_JTREE)
for (i in 1:ncol(RWR_FC2_DAG_GSE11352_exprs_for_JTREE)) {
  for (j in 1:1000) {
    RWR_FC2_DAG_permuted_exprs_for_JTREE <- cbind(RWR_FC2_DAG_permuted_exprs_for_JTREE,
                                                  RWR_FC2_DAG_GSE11352_exprs_for_JTREE[sample(1:nrow(RWR_FC2_DAG_GSE11352_exprs_for_JTREE), 
                                                                                                          size = nrow(RWR_FC2_DAG_GSE11352_exprs_for_JTREE), replace = T),i])
  }    
}
RWR_FC2_DAG_permuted_exprs_for_JTREE <- cbind(id = rownames(RWR_FC2_DAG_GSE11352_exprs_for_JTREE),RWR_FC2_DAG_permuted_exprs_for_JTREE)
rownames(RWR_FC2_DAG_permuted_exprs_for_JTREE) <- as.character(RWR_FC2_DAG_permuted_exprs_for_JTREE$id)
colnames(RWR_FC2_DAG_permuted_exprs_for_JTREE) <- c("id",1:(ncol(GSE11352_exprs_for_JTREE)*1000))
RWR_FC2_DAG_permuted_for_JTREE_mRNA <- as.data.frame(t(RWR_FC2_DAG_permuted_exprs_for_JTREE))
RWR_FC2_DAG_permuted_for_JTREE_genome <- RWR_FC2_DAG_permuted_for_JTREE_mRNA
RWR_FC2_DAG_permuted_for_JTREE_genome[-1,] <- NA
#write.table(RWR_FC2_DAG_permuted_for_JTREE_genome,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_permuted_genome.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
#write.table(RWR_FC2_DAG_permuted_for_JTREE_mRNA, paste0("output_for_significant/RWR_for_JTREE/JTREE_input/RWR_FC2_DAG_permuted_mRNA.tab"),sep = "\t",quote = FALSE,col.names = FALSE)



# Load JTREE results
RWR_GSE11352 <- read.table("output_for_significant/RWR_for_JTREE/JTREE_output/RWR_FC2_DAG_GSE11352.txt",comment.char = ">")
RWR_permuted <- read.table("output_for_significant/RWR_for_JTREE/JTREE_output/RWR_FC2_DAG_permuted.txt",comment.char = ">")

permuted_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_permuted)/(ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1),ncol = (ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1)))
rownames(permuted_activity_score) <- RWR_permuted$V1[1:(nrow(RWR_permuted)/(ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1))]
for (l in 1:(ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1)) {
  tails <- l*(nrow(RWR_permuted)/(ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1))
  heads <- tails-(nrow(RWR_permuted)/(ncol(RWR_FC2_DAG_permuted_exprs_for_JTREE)-1)-1)
  permuted_activity_score[,l] <- as.numeric(RWR_permuted[heads:tails,2])
}

GSE11352_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_GSE11352)/(ncol(RWR_FC2_DAG_GSE11352_exprs_for_JTREE)),
                                                ncol = (ncol(RWR_FC2_DAG_GSE11352_exprs_for_JTREE))))
rownames(GSE11352_activity_score) <- RWR_GSE11352$V1[1:(nrow(RWR_GSE11352)/(ncol(RWR_FC2_DAG_GSE11352_exprs_for_JTREE)))]
for (m in 1:(ncol(GSE11352_exprs_for_JTREE))) {
  tails <- m*nrow(RWR_GSE11352)/(ncol(GSE11352_exprs_for_JTREE))
  heads <- tails-(nrow(RWR_GSE11352)/(ncol(GSE11352_exprs_for_JTREE))-1)
  GSE11352_activity_score[,m] <- as.numeric(RWR_GSE11352[heads:tails,2])
}

#Annotating and rearranging
colnames(GSE11352_activity_score) <- c("12hrs_Ctrl1","12hrs_Ctrl2","12hrs_Ctrl3","12hrs_Tx1","12hrs_Tx2","12hrs_Tx3",
                                       "24hrs_Ctrl1","24hrs_Ctrl2","24hrs_Ctrl3","24hrs_Tx1","24hrs_Tx2","24hrs_Tx3",
                                       "48hrs_Ctrl1","48hrs_Ctrl2","48hrs_Ctrl3","48hrs_Tx1","48hrs_Tx2","48hrs_Tx3")
GSE11352_activity_score <- GSE11352_activity_score[,c(1:3,7:9,13:15,4:6,10:12,16:18)]
permuted_activity_score <- permuted_activity_score[,c(1:3000,6001:9000,12001:15000,3001:6000,9001:12000,15001:18000)]

## null distribution
ggplot(data = as.data.frame(t(permuted_activity_score[,15000:18000])), aes(x =ESR1)) +
  geom_histogram(color = "white", bins = 20,fill = "royalblue") +
  labs(title = "ESR1 null distribution",
       x = "Activity Score",y = "Frequency")

# Significance assessment
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

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(GSE11352_activity_significant), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(GSE11352_activity_significant)/paletteLength, max(GSE11352_activity_significant), length.out=floor(paletteLength/2)))

annotation <- data.frame(Gene_Type = factor(rownames(GSE11352_activity_significant) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(GSE11352_activity_significant) # check out the row names of annotation
pheatmap(as.matrix(GSE11352_activity_significant) ,fontsize_row=5,main = paste("Sample Level Gene Activity Score for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,annotation_row = annotation,annotation_names_row = F)

# Phenotype level Activity state
phenotypes <- as.data.frame(matrix(0,nrow =nrow(GSE11352_activity_significant),ncol = 6 ),row.names = rownames(GSE11352_activity_significant))
for (i in 1:6) {
  phenotypes[,i] <- sign(rowSums(GSE11352_activity_significant[,(3*i-2):(3*i)]))
}
colnames(phenotypes) <- c("Ctrl-12hrs","Ctrl-24hrs","Ctrl-48hrs","Tx-12hrs","Tx-24hrs","Tx-48hrs")
paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(phenotypes), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(phenotypes)/paletteLength, max(phenotypes), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(phenotypes) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(phenotypes) # check out the row names of annotation
pheatmap(as.matrix(phenotypes) ,fontsize_row=5,breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,annotation_row = annotation,annotation_names_row = F)

# For 48 Hours treatment group

Tx48_majority_vote <- as.data.frame(sign(phenotypes[,6]-phenotypes[,3]) ,row.names = rownames(phenotypes))
colnames(Tx48_majority_vote) <- c("node_sign")
Tx48_majority_vote_graph <- induced.subgraph(RWR_FC2_DAG_graph,as.character(rownames(Tx48_majority_vote)[which(Tx48_majority_vote$node_sign != 0)]))
Tx48_majority_vote_graph <- induced.subgraph(Tx48_majority_vote_graph,names(which(degree(Tx48_majority_vote_graph,V(Tx48_majority_vote_graph)) >0)))

plot(Tx48_majority_vote_graph, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(Tx48_majority_vote_graph)$name, vertex.label.color="black",
     
     vertex.label.cex=.5,layout = layout_nicely(Tx48_majority_vote_graph),vertex.size = 9)


regulated_genes <- omnipath[which(omnipath$source_name %in% rownames(Tx48_majority_vote_nodeList)),]
regulated_DEGs <- regulated_genes[which(regulated_genes$target_name %in% miarray$SYMBOL),]
print(paste("48 hours regulatory network directly regulates", length(unique(regulated_genes$target_name)), "DEGs"))

## sample level
majority_vote <- cbind(Ctrl = as.data.frame(sign(rowSums(GSE11352_activity_significant[,1:9]))),Tx = as.data.frame(sign(rowSums(GSE11352_activity_significant[,10:18]))))
colnames(majority_vote) <- c("Ctrl","Tx")
majority_vote <- majority_vote[which(majority_vote$Tx-majority_vote$Ctrl !=0),]

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(majority_vote), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(majority_vote)/paletteLength, max(majority_vote), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(majority_vote) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(majority_vote) # check out the row names of annotation
pheatmap(as.matrix(majority_vote) ,fontsize_row=5,main = paste("Treatment Level Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

majority_Vote_graph <- induced.subgraph(RWR_FC2_MST_V21_graph,as.character(rownames(majority_vote)),names(which(degree(majority_Vote_graph,V(majority_Vote_graph)) >0)))

plot(majority_vote_graph, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(majority_vote_graph)$name, vertex.label.color="black",
     vertex.label.cex=.5,layout = layout_nicely(majority_vote_graph),vertex.size = 9) 

