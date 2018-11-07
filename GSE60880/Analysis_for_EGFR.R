library(igraph)
library(dplyr)
library(data.table)
library(dnet)
library(tidyverse)
library(pheatmap)
library(predictionet)
#########################################################################################
#Calculte p-value
#########################################################################################


#########################################################################################
# DFS in RWR
#########################################################################################

RWR_dfs <- function(r1,graph,miarray){
  passedEdges <- c()
  regulatedDEG <- 0
  ### Find r2
  possible_r2 <- unique(neighbors(graph,r1, "out")$name)
  
  ###if r1 has childs(r2)
  if (length(possible_r2)>0) {  ### r1 has r2
    passedNodes <- c(r1)
    for (i in 1:length(possible_r2)) {
      r2 <- possible_r2[i]
      if (r2 %in% rownames(miarray)) {
        # passedEdges <- c(passedEdges, get.edge.ids(graph,c(r1, r2)))
        regulatedDEG <- regulatedDEG + 1
      }
      else{
        passedNodes <- c(passedNodes,r2)
        possible_r3 <- unique(neighbors(graph,r2,"out")$name)
        possible_r3 <- possible_r3[!possible_r3 %in% passedNodes]      ### That's r3
        ###if r2 has childs(r3)
        if (length(possible_r3)>0) {  ### r2 has r3
          for (i in 1:length(possible_r3)) {
            r3 <- possible_r3[i]
            if (r3 %in% rownames(miarray)) {
              passedEdges <- c(passedEdges, get.edge.ids(graph,c(r1, r2)),
                               get.edge.ids(graph,c(r2, r3)))
              passedEdges <- c(passedEdges, get.edge.ids(graph,c(r1, r2)))
              regulatedDEG <- regulatedDEG + 1
            }
            else{
              passedNodes <- c(passedNodes,r3)
              possible_r4 <- unique(neighbors(graph,r3,"out")$name)
              possible_r4 <- possible_r4[!possible_r4 %in% passedNodes]      ### That's r4
              ###if r3 has childs(r4)
              if (length(possible_r4)>0) {  ### r3 has r4
                for (i in 1:length(possible_r4)) {
                  r4 <- possible_r4[i]
                  passedNodes <- c(passedNodes,r4)
                  if (r4 %in% rownames(miarray)) {
                    regulatedDEG <- regulatedDEG + 1
                    passedEdges <- c(passedEdges,get.edge.ids(graph,c(r1, r2)),
                                     get.edge.ids(graph,c(r2, r3)))
                  }
                  else{passedEdges <- c(passedEdges,get.edge.ids(graph,c(r1, r2)),
                                        get.edge.ids(graph,c(r2, r3)),
                                        get.edge.ids(graph,c(r3, r4)))}
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

# scaled_GSE60880_exprs <- read.table("data_for_markdown/scaled_exprs_limma.tab",sep = "\t")
omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
DEGList <- read.csv("data_for_markdown/GSE60880_DEGList.csv",row.names = 1)
GSE60880_exprs_for_JTREE <- read.csv("data_for_markdown/GSE60880_scaled_exprs.csv",row.names = 1)
omniNet <- graph_from_data_frame(d = omnipath,directed = T)

#########################################################################################
# RWR
#########################################################################################
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

#########################  DFS   ################################################################
FC1.5_undirected_RWR0.65_JTREE_candidates <- c()
FC1.5_undirected_RWR0.65_JTREE_edgeList <- as.data.frame(matrix(NA,nrow = 0,ncol = 4),row.names = c("from","to","sign","weight"))

r1_RWR0.65_DFS_graph <- RWR_dfs("EGFR",significant_genes_graph,miarray)
FC1.5_undirected_RWR0.65_JTREE_edgeList <- rbind(FC1.5_undirected_RWR0.65_JTREE_edgeList,
                                                    igraph::as_data_frame(r1_RWR0.65_DFS_graph))

GSE60880_EGFR_graph <- graph_from_data_frame(unique(FC1.5_undirected_RWR0.65_JTREE_edgeList),directed = T)

GSE60880_EGFR_edgeList <- igraph::as_data_frame(GSE60880_EGFR_graph)

#########################################################################################
# Input for Inference
#########################################################################################

GSE60880_EGFR_node_part <- as.data.frame(as.character(unique(c(GSE60880_EGFR_edgeList$from,GSE60880_EGFR_edgeList$to))))
GSE60880_EGFR_node_part[,2] <- as.character(NA)
colnames(GSE60880_EGFR_node_part) <- c("id","type")
GSE60880_EGFR_node_part[is.na(GSE60880_EGFR_node_part$type),2] <- "protein"
GSE60880_EGFR_node_part <- GSE60880_EGFR_node_part[,c(2,1)]
### edge part
GSE60880_EGFR_edge_part <- GSE60880_EGFR_edgeList[,c(1,2,3)]
GSE60880_EGFR_edge_part[which(GSE60880_EGFR_edge_part$sign == 1),3] <- "-t>"
GSE60880_EGFR_edge_part[which(GSE60880_EGFR_edge_part$sign == -1),3] <- "-t|"
# write.table(GSE60880_EGFR_node_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t")
# write.table(GSE60880_EGFR_edge_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t",append = T)

### expression part
GSE60880_EGFR_GSE60880_exprs_for_JTREE <- GSE60880_exprs_for_JTREE[as.character(GSE60880_EGFR_node_part$id),]
rownames(GSE60880_EGFR_GSE60880_exprs_for_JTREE) <- as.character(GSE60880_EGFR_node_part$id)
GSE60880_EGFR_GSE60880_for_JTREE_mRNA <- as.data.frame(t(GSE60880_EGFR_GSE60880_exprs_for_JTREE))
GSE60880_EGFR_GSE60880_for_JTREE_genome <- GSE60880_EGFR_GSE60880_for_JTREE_mRNA
GSE60880_EGFR_GSE60880_for_JTREE_genome[,] <- NA

# write.table(rbind(id = colnames(GSE60880_EGFR_GSE60880_for_JTREE_genome),GSE60880_EGFR_GSE60880_for_JTREE_genome),paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_GSE60880_genome.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(rbind(id = colnames(GSE60880_EGFR_GSE60880_for_JTREE_mRNA),GSE60880_EGFR_GSE60880_for_JTREE_mRNA), paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_GSE60880_mRNA.tab"),sep = "\t",quote = FALSE,col.names = FALSE)

# Permuting Sample
GSE60880_EGFR_permuted_exprs_for_JTREE <- as.data.frame(matrix(NA,ncol = 0,nrow = nrow(GSE60880_EGFR_GSE60880_exprs_for_JTREE)))
rownames(GSE60880_EGFR_permuted_exprs_for_JTREE) <- rownames(GSE60880_EGFR_GSE60880_exprs_for_JTREE)
for (i in 1:ncol(GSE60880_EGFR_GSE60880_exprs_for_JTREE)) {
  for (j in 1:1000) {
    GSE60880_EGFR_permuted_exprs_for_JTREE <- cbind(GSE60880_EGFR_permuted_exprs_for_JTREE,
                                                       GSE60880_EGFR_GSE60880_exprs_for_JTREE[sample(1:nrow(GSE60880_EGFR_GSE60880_exprs_for_JTREE), 
                                                                                                        size = nrow(GSE60880_EGFR_GSE60880_exprs_for_JTREE), replace = T),i])
  }    
}
GSE60880_EGFR_permuted_exprs_for_JTREE <- cbind(id = rownames(GSE60880_EGFR_GSE60880_exprs_for_JTREE),GSE60880_EGFR_permuted_exprs_for_JTREE)
rownames(GSE60880_EGFR_permuted_exprs_for_JTREE) <- as.character(GSE60880_EGFR_permuted_exprs_for_JTREE$id)
colnames(GSE60880_EGFR_permuted_exprs_for_JTREE) <- c("id",1:(ncol(GSE60880_exprs_for_JTREE)*1000))
GSE60880_EGFR_permuted_for_JTREE_mRNA <- as.data.frame(t(GSE60880_EGFR_permuted_exprs_for_JTREE))
GSE60880_EGFR_permuted_for_JTREE_genome <- GSE60880_EGFR_permuted_for_JTREE_mRNA
GSE60880_EGFR_permuted_for_JTREE_genome[-1,] <- NA
# write.table(GSE60880_EGFR_permuted_for_JTREE_genome,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_permuted_genome.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(GSE60880_EGFR_permuted_for_JTREE_mRNA, paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_permuted_mRNA.tab"),sep = "\t",quote = FALSE,col.names = FALSE)


# 
# system(paste0("chmod +x output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_JTREE.sh"))
# system(paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_JTREE.sh"))
# system(paste0("rm -r output_for_significant/RWR_for_JTREE/JTREE_input/GSE60880_EGFR_JTREE.sh") )


#########################################################################################
# Assessing activity state
#########################################################################################
RWR_GSE60880 <- read.table(paste0("output_for_significant/RWR_for_JTREE/JTREE_output/GSE60880_EGFR_GSE60880.txt"),comment.char = ">")
RWR_permuted <- read.table(paste0("output_for_significant/RWR_for_JTREE/JTREE_output/GSE60880_EGFR_permuted.txt"),comment.char = ">")

# Load JTREE results ###################################################################
permuted_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_permuted)/(ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1),ncol = (ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1)))
rownames(permuted_activity_score) <- RWR_permuted$V1[1:(nrow(RWR_permuted)/(ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1))]
for (l in 1:(ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1)) {
  tails <- l*(nrow(RWR_permuted)/(ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1))
  heads <- tails-(nrow(RWR_permuted)/(ncol(GSE60880_EGFR_permuted_exprs_for_JTREE)-1)-1)
  permuted_activity_score[,l] <- as.numeric(RWR_permuted[heads:tails,2])
}


GSE60880_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_GSE60880)/(ncol(GSE60880_EGFR_GSE60880_exprs_for_JTREE)),
                                                ncol = (ncol(GSE60880_EGFR_GSE60880_exprs_for_JTREE))))
rownames(GSE60880_activity_score) <- RWR_GSE60880$V1[1:(nrow(RWR_GSE60880)/(ncol(GSE60880_EGFR_GSE60880_exprs_for_JTREE)))]
for (m in 1:(ncol(GSE60880_exprs_for_JTREE))) {
  tails <- m*nrow(RWR_GSE60880)/(ncol(GSE60880_exprs_for_JTREE))
  heads <- tails-(nrow(RWR_GSE60880)/(ncol(GSE60880_exprs_for_JTREE))-1)
  GSE60880_activity_score[,m] <- as.numeric(RWR_GSE60880[heads:tails,2])
}
colnames(GSE60880_activity_score) <- colnames(GSE60880_EGFR_GSE60880_exprs_for_JTREE)

#plow raw result

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(GSE60880_activity_score), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(GSE60880_activity_score)/paletteLength, max(GSE60880_activity_score), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(GSE60880_activity_score) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(GSE60880_activity_score) # check out the row names of annotation
pheatmap(as.matrix(GSE60880_activity_score) ,fontsize_row=5,main = paste("RWR results"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

# GFER null distribution
hist(as.numeric(GSE60880_activity_score["EGFR",]))
hist(as.numeric(c(permuted_activity_score["EGFR",41001:47000])),breaks = 80)
ggplot(data = as.data.frame(t(permuted_activity_score[,])), aes(x =EGFR)) +
  geom_histogram(color = "white", bins = 20,fill = "royalblue") +
  labs(title = "EGFR null distribution",
       x = "Activity Index",y = "Frequency")

#assessing siginificance ###################################################################

GSE60880_activity_quantile <- as.data.frame(matrix(0,nrow = nrow(GSE60880_activity_score),ncol = ncol(GSE60880_activity_score)))
for (i in 1:nrow(GSE60880_activity_quantile)) {
  for (j in 1:ncol(GSE60880_activity_quantile)) {
    null_distribution <- ecdf(as.numeric(permuted_activity_score[i,(1+(j-1)*1000):(1000*j)]))
    GSE60880_activity_quantile[i,j] <- null_distribution(GSE60880_activity_score[i,j])
  }
}

rownames(GSE60880_activity_quantile) <- rownames(GSE60880_activity_score)
colnames(GSE60880_activity_quantile) <- colnames(GSE60880_activity_score)

GSE60880_activity_significant <- GSE60880_activity_quantile
for (i in 1:nrow(GSE60880_activity_significant)) {
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] >= 0.9)] <- 1
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] <= 0.1)] <- -1
  GSE60880_activity_significant[i,which(GSE60880_activity_quantile[i,] > 0.1 & GSE60880_activity_quantile[i,]<0.9)] <- 0
}

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(GSE60880_activity_significant), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(GSE60880_activity_significant)/paletteLength, max(GSE60880_activity_significant), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(GSE60880_activity_significant) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(GSE60880_activity_significant) # check out the row names of annotation
pheatmap(as.matrix(GSE60880_activity_significant) ,fontsize_row=5,main = paste("Sample Level Gene Activity Score for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,annotation_row = annotation,annotation_names_row = F)

# calcualting phenotype signficance ###################################################################
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

#Heatmap

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(phenotypes), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(phenotypes)/paletteLength, max(phenotypes), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(phenotypes) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(phenotypes) # check out the row names of annotation
pheatmap(as.matrix(phenotypes) ,fontsize_row=5,main = paste("8 hours treated group Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,annotation_row = annotation,annotation_names_row = F)

# Heatmap for 8hr Tx
pheatmap(as.matrix(phenotypes[which(phenotypes$`Tx-8hrs` != phenotypes$`Ctrl-8hrs` ),c(4,8)]) ,fontsize_row=5,main = paste("8 hours treated group Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE,annotation_row = annotation,annotation_names_row = F)

#Building regulatory network

Tx48_majority_vote <- as.data.frame(sign(phenotypes[,8]-phenotypes[,4]) ,row.names = rownames(phenotypes))
colnames(Tx48_majority_vote) <- c("node_sign")
Tx48_majority_vote_graph <- induced.subgraph(GSE60880_EGFR_graph,as.character(rownames(Tx48_majority_vote)[which(Tx48_majority_vote$node_sign != 0)]))
Tx48_majority_vote_graph <- induced.subgraph(Tx48_majority_vote_graph,names(which(degree(Tx48_majority_vote_graph,V(Tx48_majority_vote_graph)) >0)))

plot(Tx48_majority_vote_graph, edge.arrow.size=.2, edge.curved=0,

     vertex.color="orange", vertex.frame.color="#555555",

     vertex.label=V(Tx48_majority_vote_graph)$name, vertex.label.color="black",

     vertex.label.cex=.5,layout = layout_nicely(Tx48_majority_vote_graph),vertex.size = 9)

Tx48_majority_vote_edgeList <- igraph::as_data_frame(Tx48_majority_vote_graph)


# Treatment level ###################################################################

majority_vote <- cbind(Ctrl = as.data.frame(sign(rowSums(GSE60880_activity_significant[,1:23]))),Tx = as.data.frame(sign(rowSums(GSE60880_activity_significant[,24:47]))))
colnames(majority_vote) <- c("Ctrl","Tx")
majority_vote <- majority_vote[which(majority_vote$Tx-majority_vote$Ctrl !=0),]

# Heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(majority_vote), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(majority_vote)/paletteLength, max(majority_vote), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(majority_vote) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(majority_vote) # check out the row names of annotation
pheatmap(as.matrix(majority_vote) ,fontsize_row=5,main = paste("Treatment Level Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

# Regulatory network
majority_Vote_graph <- induced.subgraph(GSE60880_EGFR_graph,as.character(rownames(majority_vote)))
majority_Vote_graph <- induced.subgraph(majority_Vote_graph,names(which(degree(majority_Vote_graph,V(majority_Vote_graph)) >0)))
plot(majority_Vote_graph, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(Tx48_majority_vote_graph)$name, vertex.label.color="black",
     
     vertex.label.cex=.5,layout = layout_nicely(Tx48_majority_vote_graph),vertex.size = 9)
majority_Vote_edgeList <- igraph::as_data_frame(majority_Vote_graph)


# Phenotype vote? just for try ###################################################################
majority_vote <- cbind(Ctrl = as.data.frame(sign(rowSums(phenotypes[,1:4]))),Tx = as.data.frame(sign(rowSums(phenotypes[,5:8]))))
colnames(majority_vote) <- c("Ctrl","Tx")
majority_vote <- majority_vote[which(majority_vote$Tx-majority_vote$Ctrl !=0),]

#Heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(majority_vote), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(majority_vote)/paletteLength, max(majority_vote), length.out=floor(paletteLength/2)))
annotation <- data.frame(Gene_Type = factor(rownames(majority_vote) %in% rownames(miarray), labels = c("URs","DEGs" )))
rownames(annotation) <- rownames(majority_vote) # check out the row names of annotation
pheatmap(as.matrix(majority_vote) ,fontsize_row=5,main = paste("Treatment Level Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

# Treatment level Regulatory network
majority_Vote_graph <- induced.subgraph(GSE60880_EGFR_graph,as.character(rownames(majority_vote)))
majority_Vote_graph <- induced.subgraph(majority_Vote_graph,names(which(degree(majority_Vote_graph,V(majority_Vote_graph)) >0)))
plot(majority_Vote_graph, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(Tx48_majority_vote_graph)$name, vertex.label.color="black",
     
     vertex.label.cex=.5,layout = layout_nicely(Tx48_majority_vote_graph),vertex.size = 9)




