
#rm(list = ls()) 

library(igraph)
library(spectralGraphTopology)
library(quadprog)
library(pals)
library(ggplot2)
library(kernlab)


network_1<-function(g){
    if (is.igraph(g)==T){
      g=g
      }
  ifelse(is.data.frame(g)==T || is.list(g)==T,g=graph_from_data_frame(g,directed = F)
                ,g=graph_from_incidence_matrix(g,directed = F))
  } 

network_1<-function(g){
  
  if (is.igraph(g)==T){
    g=g
  }
  else if(is.data.frame(g)==T || is.list(g)==T){
    g=graph_from_data_frame(g,directed = F)
  }
  else{
    g=graph_from_incidence_matrix(g,directed = F)
  }
  
#---------Basic plot----------------
#components(g, mode = c("weak", "strong"))
  g=delete.vertices(g , which(degree(g)==0))#deletes isolates
  
  plot(g,edge.arrow.size=0.25,
       vertex.color="grey")
  
#---------Graph Matrices---------------
  
  matrix_adjacency=g[]
  
  
  matrix_laplacian=laplacian_matrix(g,sparse = F)
  
  print(matrix_laplacian)
  
  
#------Some graph features--------
Geodesic_distance<-distances(g,weights=NA)#ignore edge weights
Mean_distance<-mean(Geodesic_distance);print(Mean_distance)#average path lentgths  
Diameter<-diameter(g,weights = NA);print(Diameter) 

Degree_graph<-degree(g, mode="all")
print(summary(Degree_graph))#summary of the degrees
hist(Degree_graph, breaks=1:vcount(g)-1, main="Histogram of node degree")


cliques(g) #list of cliques (connected subgraphs)    
clique_num(g)#maximal clique in an undirected network
sapply(cliques(g), length) # clique sizes

#---------Eigen decomposition-------------
  
adjacency_spectrum<-eigen(matrix_adjacency)  
laplacian_spectrum<-eigen(matrix_laplacian)

print(laplacian_spectrum)

#---------Fiedler pair---------------------
Fiedler_val<-laplacian_spectrum$values[length(laplacian_spectrum$value)-1]
Fiedler_vec=laplacian_spectrum$vectors[,ncol(laplacian_spectrum$vectors)-1]
print(Fiedler_val);plot(Fiedler_vec)

#--------Spectral membership and clustering-------------------
laplac_vec1<-Fiedler_vec
laplac_vec2<-laplacian_spectrum$vectors[,ncol(laplacian_spectrum$vectors)-2]
round(laplac_vec2,digits = 2);round(laplac_vec1,digits = 2)
coords<-as.matrix(cbind(laplac_vec1,laplac_vec2))

clust_g<- kmeans(Fiedler_vec,2)#fiedler vector for clustering
clust_g[1]#cluster membership for nodes
plot(1:length(V(g)), rev(laplacian_spectrum$values)[1:length(V(g))], log="y", main = "cluster confirmation")
abline(v=2.25, col="red", lty=2)

#spectral clustering
s<-specc(coords, centers=2)
plot(coords, col=s, pch=4)#estimated classes(2)
points(coords, col=laplacian_spectrum$values, pch=20)#true classes


#---------Spectral drawing (Hall)-------------------
par(mfrow=c(1,2))
plot(laplac_vec1);plot(laplac_vec2)
g$layout<-coords#second and third eigenvectors as cordinate for the spectral drawing of the graph
plot(g,
     vertex.color="grey")
}
