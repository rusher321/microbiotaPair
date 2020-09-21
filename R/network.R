#' @import igraph
#' @import Matrix
#' @import SpiecEasi
#' sparcc network
#'
#' @param dat data.frame or matrix, row is sample ID
#' @param cutoff numberic, the sparcc coefficient of species-species indicate
#' the cutoff of eage
#' @param main character, the network title
#' @param count logistic, indicated the microbiota data is count of relative abundance
#' @param layout character, the paremeter from the igraph indicated the network layout
#' @param ...
#'
#' @return figure
#' @export
#'
#' @usage sparccNet(dat = microbiota, cutoff = 0.3, main = "test", count = F)
#' @examples
#'
#' library(dplyr)
#' data("physeq_data")
#' physeq <- physeq_data
#' microbitota <- otu_table(phyloseq)
#' microbitoFilter <- filterPer(microbiota, 2, 0.5)
#' sparccNet(dat = microbiotaFilter, cutoff = 0.3, main = "test", count = "F)
sparccNet <- function(dat, cutoff, main, count ,layout = "layout.circle",...){

  # library(igraph)
  # library(Matrix)
  # library(SpiecEasi)

  #
  renorm <- function(dat){
    # normlization function
    trans <- function(x){
      y <- x/(sum(x))
      return(y)
    }
    # tran the dataframe
    dat2 <- t(dat)
    dat3 <- t(apply(dat2, 2, trans))

    return(dat3)
  }

  # the function is count data or relative abundan
  if(count){
    datRenorm <- dat
  }else{
    datRenorm <- round(renorm(dat)*10^7)
  }
  basesparcc <- sparcc(datRenorm)
  graph <- basesparcc$Cor
  num.row <- nrow(graph)

  # tran the corelation to adj matrox

  for(i in 1:num.row){
    for(j in 1:num.row){
      a <- graph[i, j]
      graph[i,j] <- ifelse(abs(a) >= cutoff, a, 0)
    }
  }

  diag(graph) <- 0
  igraph <- adj2igraph(Matrix(graph, sparse=TRUE))

  # set edge color，postive correlation
  # postive correlation is red, negative correlation is blue
  igraph.weight = E(igraph)$weight
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#fc9272", ifelse(E.color<0, "#31a354","grey"))
  E(igraph)$color = as.character(E.color)

  # add the edge width
  E(igraph)$width = abs(igraph.weight)*4
  # add the node color
  V(igraph)$color <- "#636363"
  plot(igraph,
       layout=layout,
       vertex.label=colnames(dat),
       main = main ,...)

}


#' pairNet
#' Pair Networt for a specific species on 2 time point
#'
#' @param microbiota a data.frame object, species profile, row is sample
#' @param metadata a data.frame object, metadata table, row is sample
#' @param cutoff numberic, the sparcc coefficient of species-species indicate
#' the cutoff of eage
#' @param species charater, species names, eg "Vibrio_harveyi"
#' @param name charater, provida a name for the specific species of treatment
#'
#' @return matrix
#' @export
#'
#' @examples
pairNet <- function(microbiota, metadata, cutoff, species, name){

  # to compute the sparcc
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  matchdat <- matchdat[order(matchdat[, time_varname]), ]
  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]
  matchmicrobiota <- filterPer(matchmicrobiota, row = 2, percent = 0.2)
  g1microbiota <- round(matchmicrobiota[1:number, ]*10^7)
  g2microbiota <- round(matchmicrobiota[(number+1):(2*number), ]*10^7)

  Cor1 <- sparcc(g1microbiota)$Cor
  Cor2 <- sparcc(g2microbiota)$Cor
  colnames(Cor1) <- rownames(Cor1) <- rownames(Cor2) <- colnames(Cor2) <- colnames(microbiota)
  # generate adj matrix on species

  bacindex <- which(rownames(Cor1) == species)
  retainbac <- names(Cor1[,bacindex][abs(Cor1[,bacindex])>=cutoff])
  retainnet <- Cor1[retainbac, retainbac]
  which(retainbac== species) -> speindex
  retainnet[-speindex, -speindex] <- 0
  retainnet[speindex, speindex] <- 0

  bacindex2 <- which(rownames(Cor2) == species)
  retainbac2 <- names(Cor2[,bacindex][abs(Cor2[,bacindex])>=cutoff])
  retainnet2 <- Cor2[retainbac2, retainbac2]
  which(retainbac2 == species) -> speindex2
  retainnet2[-speindex2, -speindex2] <- 0
  retainnet2[speindex2, speindex2] <- 0

  rownames(retainnet2)[speindex2] <- colnames(retainnet2)[speindex2] <- paste0(species,"_C2")
  retainbac2[speindex2] <- paste0(species, name)
  # combine the two network

  comBac <- unique(c(retainbac, retainbac2))
  out <- matrix(0, nrow=length(comBac), ncol=length(comBac))
  rownames(out) <- colnames(out) <- comBac

  # input the corresponde value
  pmatch(retainbac, comBac) -> netindex1
  netindex1[speindex] -> netspindex1
  out[netindex1, netspindex1] <- retainnet[, speindex]
  out[netspindex1, netindex1] <- retainnet[speindex, ]

  pmatch(retainbac2, comBac) -> netindex2
  netindex2[speindex2] -> netspindex2
  out[netindex2, netspindex2] <- retainnet2[, speindex2]
  out[netspindex2, netindex2] <- retainnet2[speindex2, ]

  return(out)

}

#' @import igraph
#' @import Matrix
#' simpleNet
#' a simple network
#' @param adjmatrix the pairNet result
#' @param main the network title
#'
#' @return figure
#' @export
#'
#' @examples
#'
simplenet <- function(adjmatrix, main,...){
  #library(igraph)
  #library(Matrix)

  graph <- adjmatrix
  num.row <- nrow(graph)
  igraph <- adj2igraph(Matrix(graph, sparse=TRUE))

  # set edge color
  igraph.weight = E(igraph)$weight
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#fc9272",ifelse(E.color<0, "#31a354","grey"))
  E(igraph)$color = as.character(E.color)
  # set node size
  index <- grep(main, colnames(adjmatrix))
  nodesize <- rep(7, num.row)
  nodesize[index] <- 15
  V(igraph)$size <- nodesize

  # add the edge width
  E(igraph)$width = abs(igraph.weight)*4
  # add the node color
  #V(igraph)$color <- "#636363"
  plot(igraph,
       vertex.label=colnames(adjmatrix),
       vertex.label.cex = .6,
       vertex.label.color = "black",
       vertex.color="lightsteelblue2",
       vertex.frame.color="gray",
       main = main ,...)

}

simulateNetwork <- function(dat, repeatN , sampleN){

  # sample the data
  out <- matrix(0, nrow = ncol(dat), ncol = ncol(dat))
  for(i in 1:repeatN){
	print(i)
    index <- sample(c(1:nrow(dat)), sampleN, replace = T)
    dat2 <- dat[index, ]
    out <- out+sparcc(dat2)$Cor
  }

  out2 <- out/repeatN
  return(out2)

}


sparCCnetwork <- function(microdata , rank , phemeta, group, group_var){

  # generate the result
  library(phyloseq)
  library(SpiecEasi)
  library(rmeta)

  phe.cln <- sample_data(phemeta[phemeta[, group_var]==group, ])
  otu.cln <- microdata[,rownames(phe.cln)]
  # simulate the count data
  otuN <- t(round(renorm(t(otu.cln))*10^5))
  otuN <- otu_table(otuN, taxa_are_rows = T)

  pylores <- phyloseq(phe.cln, otu.cln, rank)

  # repeat sparCC
  sparccres <- simulateNetwork(dat = t(otuN), repeatN = 20, sampleN = 30)
  sparcc.graph <- abs(sparccres) >= 0.3
  diag(sparcc.graph) <- 0
  sparcc.graph2 <- Matrix(sparcc.graph, sparse=TRUE)
  ig.sparcc <- adj2igraph(sparcc.graph2, vertex.attr=list(name=taxa_names(pylores)))

  # plot
  fig <- plot_network(ig.sparcc, pylores, type='taxa', color="ta2")

  return(list(fig, sparcc.graph))

}


####


network_sta <- function(cor_matrix, cutoff = 0.3){

  #cutoff <- 0.3
  library(igraph)

  diag(cor_matrix) <- 0
  g0_trans <- ifelse(cor_matrix < -cutoff & cor_matrix < 0 , -1,
                     ifelse(cor_matrix > cutoff & cor_matrix > 0, 1, 0))
  #g0_net <- graph_from_incidence_matrix(g0_trans)
  g0_net <- adj2igraph(Matrix(g0_trans, sparse=TRUE))
  # output
  range_g0 <- quantile(cor_matrix[lower.tri(cor_matrix)])
  edge_g0 <-  c(sum(g0_trans[lower.tri(g0_trans)]==1), sum(g0_trans[lower.tri(g0_trans)]==-1))
  names(edge_g0) <- c("pos", "neg")
  degree_g0 <- degree(g0_net)
  betwenness_g0 <- betweenness(g0_net)
  closeness_g0 <- closeness(g0_net)

  out <- list(edge_g0, range_g0, degree_g0, betwenness_g0, closeness_g0)
  names(out) <- c("edge", "range", "degree", "betwenness", "closeness")

  return(out)

}



# plot
network_tran <- function(adj, abun, rankinf){

  ### edge inf
  n <- ncol(adj)
  node <- colnames(adj)
  from <- c()
  to <- c()
  edge.width <- edge.colour <- c()
  for(i in 1:n){
    for(j in 2:n){
      from <- c(from, node[i])
      to <- c(to, node[j])
      edge.width <- c(edge.width, adj[i,j])
      edge.colour <- c(edge.colour, adj[i,j])
    }
  }

  edge_inf <- data.frame(from, to, edge.width, edge.colour)

  ### node inf
  node_inf <- as.data.frame(apply(abun, 2, mean))
  colnames(node_inf) <- "node.size"
  node_inf$node <- rownames(node_inf)
  node_inf$colour <- rankinf[rownames(node_inf),]

  out <- list(edge_inf, node_inf)

  return(out)

}







