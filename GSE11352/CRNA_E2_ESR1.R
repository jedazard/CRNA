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

# scaled_GSE11352_exprs <- read.table("data_for_markdown/scaled_exprs_limma.tab",sep = "\t")
omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
DEGList <- read.csv("data_for_markdown/GSE11352_DEGList.csv",row.names = 1)
GSE11352_exprs_for_JTREE <- read.csv("data_for_markdown/GSE11352_scaled_exprs.csv",row.names = 1)
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

significant_genes_graph <- igraph::simplify(graph =  induced.subgraph(omniNet,c(FC2_undirected_RWR0.75_selected_node,as.character(miarray$SYMBOL[which(miarray$SYMBOL %in% V(omniNet)$name)]))),remove.multiple = T,
                                            remove.loops = T,edge.attr.comb = "mean")
print(paste("Here are", length(FC2_undirected_RWR0.75_selected_node), "significant genes in undirected RWR with 0.75 restart probablility.", length(which(FC2_undirected_RWR0.75_selected_node %in% rownames(miarray))),"of them are DEGs"))

#########################  DFS   ################################################################
FC1.5_undirected_RWR0.75_JTREE_candidates <- c()
FC1.5_undirected_RWR0.75_JTREE_edgeList <- as.data.frame(matrix(NA,nrow = 0,ncol = 4),row.names = c("from","to","sign","weight"))

r1_RWR0.75_DFS_graph <- RWR_dfs("ESR1",significant_genes_graph,miarray)
FC1.5_undirected_RWR0.75_JTREE_edgeList <- rbind(FC1.5_undirected_RWR0.75_JTREE_edgeList,
                                                    igraph::as_data_frame(r1_RWR0.75_DFS_graph))

GSE11352_ESR1_graph <- graph_from_data_frame(unique(FC1.5_undirected_RWR0.75_JTREE_edgeList),directed = T)

GSE11352_ESR1_edgeList <- igraph::as_data_frame(GSE11352_ESR1_graph)

#########################################################################################
# Input for Inference
#########################################################################################

GSE11352_ESR1_node_part <- as.data.frame(as.character(unique(c(GSE11352_ESR1_edgeList$from,GSE11352_ESR1_edgeList$to))))
GSE11352_ESR1_node_part[,2] <- as.character(NA)
colnames(GSE11352_ESR1_node_part) <- c("id","type")
GSE11352_ESR1_node_part[is.na(GSE11352_ESR1_node_part$type),2] <- "protein"
GSE11352_ESR1_node_part <- GSE11352_ESR1_node_part[,c(2,1)]
### edge part
GSE11352_ESR1_edge_part <- GSE11352_ESR1_edgeList[,c(1,2,3)]
GSE11352_ESR1_edge_part[which(GSE11352_ESR1_edge_part$sign == 1),3] <- "-t>"
GSE11352_ESR1_edge_part[which(GSE11352_ESR1_edge_part$sign == -1),3] <- "-t|"
# write.table(GSE11352_ESR1_node_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t")
# write.table(GSE11352_ESR1_edge_part,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_pathway.tab"),row.names = F,col.names = F,quote = F,sep = "\t",append = T)

### expression part
GSE11352_ESR1_GSE11352_exprs_for_JTREE <- GSE11352_exprs_for_JTREE[as.character(GSE11352_ESR1_node_part$id),]
rownames(GSE11352_ESR1_GSE11352_exprs_for_JTREE) <- as.character(GSE11352_ESR1_node_part$id)
GSE11352_ESR1_GSE11352_for_JTREE_mRNA <- as.data.frame(t(GSE11352_ESR1_GSE11352_exprs_for_JTREE))
GSE11352_ESR1_GSE11352_for_JTREE_genome <- GSE11352_ESR1_GSE11352_for_JTREE_mRNA
GSE11352_ESR1_GSE11352_for_JTREE_genome[,] <- NA

# write.table(rbind(id = colnames(GSE11352_ESR1_GSE11352_for_JTREE_genome),GSE11352_ESR1_GSE11352_for_JTREE_genome),paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_GSE11352_genome.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(rbind(id = colnames(GSE11352_ESR1_GSE11352_for_JTREE_mRNA),GSE11352_ESR1_GSE11352_for_JTREE_mRNA), paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_GSE11352_mRNA.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(paste0("/usr/bin/time /home/zaza/sbenz-JTREE-c5fc0f9/JTREE -c  /home/zaza/sbenz-JTREE-c5fc0f9/em_simple.cfg -p Jan25_GSE11352_ESR1_pathway.tab -b Jan25_GSE11352_ESR1_GSE11352 > ../JTREE_output/Jan25_GSE11352_ESR1_GSE11352.txt\n\n"),paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_JTREE.sh"),quote = F,row.names = F,col.names = F)

# Permuting Sample
GSE11352_ESR1_permuted_exprs_for_JTREE <- as.data.frame(matrix(NA,ncol = 0,nrow = nrow(GSE11352_ESR1_GSE11352_exprs_for_JTREE)))
rownames(GSE11352_ESR1_permuted_exprs_for_JTREE) <- rownames(GSE11352_ESR1_GSE11352_exprs_for_JTREE)
for (i in 1:ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE)) {
  for (j in 1:1000) {
    GSE11352_ESR1_permuted_exprs_for_JTREE <- cbind(GSE11352_ESR1_permuted_exprs_for_JTREE,
                                                       GSE11352_ESR1_GSE11352_exprs_for_JTREE[sample(1:nrow(GSE11352_ESR1_GSE11352_exprs_for_JTREE), 
                                                                                                        size = nrow(GSE11352_ESR1_GSE11352_exprs_for_JTREE), replace = T),i])
  }    
}
GSE11352_ESR1_permuted_exprs_for_JTREE <- cbind(id = rownames(GSE11352_ESR1_GSE11352_exprs_for_JTREE),GSE11352_ESR1_permuted_exprs_for_JTREE)
rownames(GSE11352_ESR1_permuted_exprs_for_JTREE) <- as.character(GSE11352_ESR1_permuted_exprs_for_JTREE$id)
colnames(GSE11352_ESR1_permuted_exprs_for_JTREE) <- c("id",1:(ncol(GSE11352_exprs_for_JTREE)*1000))
GSE11352_ESR1_permuted_for_JTREE_mRNA <- as.data.frame(t(GSE11352_ESR1_permuted_exprs_for_JTREE))
GSE11352_ESR1_permuted_for_JTREE_genome <- GSE11352_ESR1_permuted_for_JTREE_mRNA
GSE11352_ESR1_permuted_for_JTREE_genome[-1,] <- NA
# write.table(GSE11352_ESR1_permuted_for_JTREE_genome,paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_permuted_genome.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(GSE11352_ESR1_permuted_for_JTREE_mRNA, paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_permuted_mRNA.tab"),sep = "\t",quote = FALSE,col.names = FALSE)
# write.table(paste0("/usr/bin/time /home/zaza/sbenz-JTREE-c5fc0f9/JTREE -c  /home/zaza/sbenz-JTREE-c5fc0f9/em_simple.cfg -p Jan25_GSE11352_ESR1_pathway.tab -b Jan25_GSE11352_ESR1_permuted > ../JTREE_output/Jan25_GSE11352_ESR1_permuted.txt\n\n"),paste0("output_for_significant/RWR_for_JTREE/JTREE_input/Jan25_GSE11352_ESR1_JTREE.sh"),quote = F,append = T,row.names = F,col.names = F)
# 


# 
# system(paste0("chmod +x output_for_significant/RWR_for_JTREE/JTREE_input/GSE11352_ESR1_JTREE.sh"))
# system(paste0("output_for_significant/RWR_for_JTREE/JTREE_input/GSE11352_ESR1_JTREE.sh"))
# system(paste0("rm -r output_for_significant/RWR_for_JTREE/JTREE_input/GSE11352_ESR1_JTREE.sh") )


#########################################################################################
# Assessing activity state
#########################################################################################
RWR_GSE11352 <- read.table(paste0("data_for_markdown/Jan25_GSE11352_ESR1_GSE11352.txt"),comment.char = ">")
RWR_permuted <- read.table(paste0("data_for_markdown/Jan25_GSE11352_ESR1_permuted.txt"),comment.char = ">")

permuted_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_permuted)/(ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1),ncol = (ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1)))
rownames(permuted_activity_score) <- RWR_permuted$V1[1:(nrow(RWR_permuted)/(ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1))]
for (l in 1:(ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1)) {
  tails <- l*(nrow(RWR_permuted)/(ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1))
  heads <- tails-(nrow(RWR_permuted)/(ncol(GSE11352_ESR1_permuted_exprs_for_JTREE)-1)-1)
  permuted_activity_score[,l] <- as.numeric(RWR_permuted[heads:tails,2])
}

GSE11352_activity_score <- as.data.frame(matrix(NA,nrow = nrow(RWR_GSE11352)/(ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE)),
                                                ncol = (ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE))))
rownames(GSE11352_activity_score) <- RWR_GSE11352$V1[1:(nrow(RWR_GSE11352)/(ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE)))]
for (m in 1:(ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE))) {
  tails <- m*nrow(RWR_GSE11352)/(ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE))
  heads <- tails-(nrow(RWR_GSE11352)/(ncol(GSE11352_ESR1_GSE11352_exprs_for_JTREE))-1)
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

# annotation <- data.frame(Gene_Type = factor(rownames(GSE11352_activity_significant) %in% rownames(miarray), labels = c("URs","DEGs" )))
# rownames(annotation) <- rownames(GSE11352_activity_significant) # check out the row names of annotation
pheatmap(as.matrix(GSE11352_activity_significant) ,fontsize_row=5,main = paste("Sample Level Gene Activity Score for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

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
pheatmap(as.matrix(phenotypes) ,fontsize_row=10,breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

## sample level
majority_vote <- cbind(Ctrl = as.data.frame(sign(rowSums(GSE11352_activity_significant[,1:9]))),Tx = as.data.frame(sign(rowSums(GSE11352_activity_significant[,10:18]))))
colnames(majority_vote) <- c("Ctrl","Tx")
majority_vote <- majority_vote[which(majority_vote$Tx-majority_vote$Ctrl !=0),]

paletteLength <- 50
myColor <- colorRampPalette(c("royalblue", "white", "orangered"))(paletteLength)
pheatmap(as.matrix(majority_vote) ,fontsize_row=5,main = paste("Treatment Level Gene Activity State for Generalized Regulatory Network"),breaks = myBreaks,
         color = myColor,cluster_cols = FALSE,scale = "none",cluster_rows = TRUE)

majority_vote_graph <- induced.subgraph(GSE11352_ESR1_graph,as.character(rownames(majority_vote)))

plot(majority_vote_graph, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(majority_vote_graph)$name, vertex.label.color="black",
     vertex.label.cex=.5,layout = layout_nicely(majority_vote_graph),vertex.size = 9) 

