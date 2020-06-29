require(knitr)
require(EpiModel)
opts_chunk$set(comment = NA, message = FALSE, tidy = FALSE)

library(intergraph)
library(network)
library(igraph)
library(ggplot2)
#library(spectralGraphTopology)
library(kernlab)
library(quadprog)
library(pals)
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
  
  #---------cleaning the igraph object----------------
  #components(g, mode = c("weak", "strong"))
  g1=delete.vertices(g , which(degree(g)==0))#deletes isolates
  g1=simplify(g,remove.multiple = T,remove.loops = T)
  n=gorder(g1)
  mean_degree=mean(degree(g1))
  
  #------converting to a network object----------------------------------
  
  network_graph=asNetwork(g1)
  detach("package:igraph")
  network_graph=set.vertex.attribute(network_graph,"risk",rbinom(n,1,0.5))
  
  Edges<-(mean_degree*n)/2
  count_risk_groups<-table(c(get.vertex.attribute(network_graph,"risk")))
  #Node_Factor<-0.75*count_risk_groups[[2]]
  Node_Match<-as.numeric(0.90*Edges)
  #Concurrent_deg<-125#length(c(degree(g))[c(degree(g))>=2])
  
  
  
  
  
  # Set formation formula (standart network model)
  formation <- ~edges + nodematch("risk")# ~edges + nodefactor("risk") + nodematch("risk") + concurrent
  
  # Set target statistics for formation
  target.stats <-c(Edges,Node_Match)#c(Edges,Node_Factor,Node_Match,Concurrent_deg)
  
  #coef.form <- -Inf
  # Obtain the offset coefficients
  # coef.diss <- dissolution_coefs(~offset(edges), 10, 0,0)#dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  # Estimate the STERGM using the edges dissolution approximation
  #obtaining the per capital size invarian of the coefficient of the dyadic inde terms by offset 
  
  
  #---------One mode closed network------------------  
  
  #estimated_network<- netest(network_graph, formation, target.stats, coef.diss,edapprox = TRUE)
  #estimated_network<-netest(network_graph, formation =  ~edges + offset(nodefactor("risk")) + offset(nodematch("risk")),
  #                            target.stats = 15, coef.form = -Inf,
  #                            coef.diss = dissolution_coefs(~offset(edges), 10, 0,0),
  #                            verbose = FALSE)
  
  est <- netest(network_graph, formation, target.stats, coef.diss, edapprox = TRUE)
  
  #set.control.ergm = control.ergm(MCMC.burnin = 1e5,MCMC.interval = 10))
  
  ## ----Timeseries of dynamic net--------------------------------------------
  dx <- netdx(est, nsims = 5, nsteps = 15)
  dx
  #  plot(dx, stats = "edges")
  
  
  ## ----plotNetdx----------------------------------------------------------
  par(mfrow = c(1, 2))#, mgp = c(2,1,0))
  plot(dx, type = "edges",plots.joined = FALSE)
  plot(dx, type = "nodematchrisk",plots.joined = FALSE)
  
  ## ----plotNetdx2----------------------------------------------------------
  par(mfrow = c(1, 2))
  plot(dx, type = "duration")
  plot(dx, type = "dissolution")
  
  # ----netParam1-----------------------------------------------------------
  param <- param.net(inf.prob = 0.1, act.rate = 5, rec.rate = 0.02)
  
  ## ----netInit1------------------------------------------------------------
  status.vector <- c(rbinom(n, 1, 0.1))
  status.vector <- ifelse(status.vector == 1, "i", "s")
  init <- init.net(status.vector = status.vector)
  
  ## ----netControl1, results = "hide"---------------------------------------
  control <- control.net(type = "SIS", nsteps = 15, nsims = 5, epi.by = "risk")
  
  ## ----netsim1, results = "hide", cache = TRUE-----------------------------
  sim1 <- netsim(est, param, init, control)
  
  
  plot(sim1)
  plot(sim1, type = "formation",plots.joined = FALSE)
  #plot(sim1, type = "network")
  ## ----printSim1-----------------------------------------------------------
  sim1
  
  ## ----summSim1------------------------------------------------------------
  summary(sim1, at = 5)
  
  ## ----netAdf1-------------------------------------------------------------
  head(as.data.frame(sim1), 1)
  
  ## ----netAdf2-------------------------------------------------------------
  head(as.data.frame(sim1, out = "mean"), 1)
  
  ## ----netGnet1------------------------------------------------------------
  nn <- get_network(sim1, sim = 5)
  nn
  
  ## ----netGetTm------------------------------------------------------------
  head(get_transmat(sim1, sim = 1), 5)
  
  ## ----netPlot1------------------------------------------------------------
  par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
  plot(sim1)
  
  ## ----netPlot2------------------------------------------------------------
  plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)
  
  ## ----netPlot3------------------------------------------------------------
  plot(sim1, y = c("si.flow", "is.flow"), qnts = 1, legend = TRUE)
  
  ## ----netPlot4------------------------------------------------------------
  plot(sim1, y = c("i.num.risk0", "i.num.risk1"), legend = TRUE)
  
  ## ----net19---------------------------------------------------------------
  par(mfrow = c(1,2), mar = c(0,0,1,0))
  plot(sim1, type = "network", at = 1, col.status = TRUE,
       main = "Prevalence at t1")
  final_plot= plot(sim1, type = "network", at = 15, col.status = TRUE,
                   main = "Prevalence at t15")
  final_plot
  
  library(igraph)
  #final_plot_graph=asIgraph(final_plot)
  
  #Once the number of edges is adjusted to preserve the mean degree, Krivitsky et al. show that all of the dyad independent terms are properly scaled
  
  # plot(g1,edge.arrow.size=0.25,
  #     vertex.color="grey")
  
  #---------Graph Matrices---------------
  
  
  matrix_adjacency=g[]
  
  
  matrix_laplacian=laplacian_matrix(g,sparse = F)
  
  #print(matrix_laplacian)
  
  
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
  print(adjacency_spectrum$values)
  #largest eigen value of the adjacency
  adj_val<-adjacency_spectrum$values[1]
  laplacian_spectrum<-eigen(matrix_laplacian)
  
  #print(laplacian_spectrum)
  
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
  
  #modularity
  m1<- cluster_walktrap(g) # This function tries to find densely connected subgraphs(nick), 
  m2<- cluster_louvain(g) #this is the method used inSah et al 2017 - allows comparison? 
  #also called communities in a graph via random walks.
  
  modularity_g <- modularity(g, membership(m2))
  
  
  #---------Spectral drawing (Hall)-------------------
  par(mfrow=c(1,2))
  plot(laplac_vec1);plot(laplac_vec2)
  g$layout<-coords#second and third eigenvectors as cordinate for the spectral drawing of the graph
  plot(g,edge.arrow.size=0.25,
       vertex.color="grey")
  
  Summary=cbind( modularity_g,Fiedler_val,adj_val,mean_degree,Diameter)
  Summary
  #colnames(Summary)<-c("modularity","fiedler_value","adjacency_value","mean_degree","Diameter")
  #rownames(Summary)<-"graph_object"
  #d=mean(sim1$i.num)
  #v=as.data.frame(sim1, out = "mean");v
}






g<- graph(edges=c(1,20,1,5, 1,19, 1,15, 1,21, 1,16, 1,23, 1,22, 1,8, 2,19, 2,10, 2,26, 2,20, 2,12, 2,8, 3,19, 3,25, 3,16, 3,4,
                  3,18, 3,7, 3,14, 3,10, 3,11, 3,24, 3,9, 4,13, 4,9, 4,18, 4,16, 4,7, 4,14, 4,25,
                  5,19,  5,10, 5,8, 5,16, 6,12, 7,25, 7,19, 7,9, 7,18, 7,23, 7,12, 7,26, 7,22, 7,16,
                  7,8, 7,21, 7,11, 7,17, 7,14, 7,10, 7,24,7,20, 7,15, 8,21 ,8,10, 8,23, 8,17, 8,19, 8,20, 8,16,
                  9,18, 9,13, 9,25, 9,21, 9,16, 9,10, 9,15, 10,17, 10,25, 10,20, 10,15, 10,23, 11,14, 11,25, 11,16,
                  11,20, 12,22, 12,26, 12,17, 12,16, 12,20, 13,24, 13,16, 13,14, 13,20, 13,25, 14,25, 14,24,
                  14,20, 15,17, 15,18, 15,16, 15,20, 16,23, 17,23, 17,20, 17,26, 18,21, 18,19, 18,20, 19,20, 19,25,
                  19,23, 20,25, 20,23, 20,26, 20,22, 21,25, 22,23, 23,26, 26,27),  n=27, directed=F)



