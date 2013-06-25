## No copyright, no warranty
## Dominic John Bennett & Jun Ying Lim
## Last update: 25/06/2013
## Functions for PD analysis BIG project 2013

## Libraries
library(picante)
library(caper)
library(lme4)

## Dom's Functions:
PD <- function(phylo, taxa, type = 1, prop = FALSE,
               display = FALSE, show.tip.label = FALSE) {
  # Calculate Faith's Phylogenetic Diversity and plot phylogeny
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  taxa: vector of taxon names for which to calcualte PD
  #  type: specify the way in which PD is calcualted
  #     1 -- sum of branch lengths of phylogeny consisting solely of the taxa,
  #       default
  #     2 -- sum of branch lengths of taxa to terminal node
  #     3 -- sum of branch lengths represented uniquely by the taxa
  #  prop: if TRUE, return the proportion of the branch not the absolute
  #   value, default FALSE
  #  display: if TRUE, plot phylogeny colouring branches which count towards
  #   PD, default FALSE
  #  show.tip.label: if TRUE, plot tip labels, default FALSE
  #
  # Return:
  #  numeric
  if(!type %in% c(1,2,3)) {
    stop("Type must be an integer: 1, 2 or 3.")
  }
  if(length(taxa) == length(phylo$tip.label)){
    if(display == TRUE){
      plot.phylo(phylo, show.tip.label = show.tip.label)
    }
    if (prop) {
      return(1)
    } else {
      return(sum(phylo$edge.length))
    }
  }
  if(type == 1 & length(taxa) == 1){
    stop("length(taxa) == 1 :
         Cannot calculate PD for a single taxon for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  edges <- match(match(taxa, phylo$tip.label), phylo$edge[,2])
  end.nodes <- phylo$edge[edges, 1]
  term.node <- length(phylo$tip.label) + 1
  if (type == 3){
    while(any(duplicated(end.nodes))){
      start.node <- end.nodes[duplicated(end.nodes)][1]
      if(sum(phylo$edge[,1] %in% start.node) == sum(end.nodes %in% start.node)){
        edge <- match(start.node, phylo$edge[,2])
        end.node <- phylo$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
      }else{
        end.nodes <- end.nodes[end.nodes != start.node]
      }
    }
  } else {
    while(TRUE){
      end.nodes <- sort(end.nodes, TRUE)
      start.node <- end.nodes[1]
      edge <- match(start.node, phylo$edge[,2])
      end.node <- phylo$edge[edge,1]
      edges <- c(edges, edge)
      end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
      if(type == 2){
        if(sum(end.nodes[1] == term.node) == length(end.nodes)){
          break
        }
      } else {
        if(sum(end.nodes[1] == end.nodes) == length(end.nodes)){
          break
        }
      }
    }
  }
  if(display){
    tip.cols <- ifelse(phylo$tip.label %in% taxa, "black", "grey")
    edge.lties <- ifelse(1:nrow(phylo$edge) %in% edges, 1, 3)
    plot.phylo(phylo, edge.lty = edge.lties, tip.color = tip.cols,
               show.tip.label = show.tip.label)
  }
  if (prop) {
    return(sum(phylo$edge.length[edges]) / sum(phylo$edge.length))
  } else {
    return(sum(phylo$edge.length[edges]))
  }
  }


plotComm <- function(comm.data, phylo, groups = rep(1, nrow(comm.data)),
                     no.margin = TRUE){
  # Plot community data on community phylogeny to visualise
  #  differences in community structure. Use colours to distinguish
  #  site groups and alpha to distinguish abundances (works like rainbow())
  #
  # Args:
  #  comm.data: community data matrix (cols taxa, rows sites)
  #  phylo: community phylogeny
  #  groups: site groups
  #
  # Return:
  #  a matrix of community data
  # plot phylogeny, allow space for points
  variable.max <- 5 + (nrow(comm.data)/4)
  #variable.max <- ifelse(variable.max > 200, variable.max, 200) #min is 200
  plot(phylo, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max))
  
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- comm.data/max(comm.data)
  
  # loop init
  ntips <- length(phylo$tip.label)
  spacing <- 0.75
  group <- groups[1]
  j <- 1
  
  # loop through sites and plot points for present species
  for(i in 1:nrow(comm.data)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(comm.data[i,])
    taxa <- phylo$tip.label[pull]
    abunds <- alphas[i, pull]
    tiplabels(tip = match(taxa, phylo$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + 0.25
    group <- groups[i]
  }
}

randCommData <- function(phylo, nsites, nspp, pa = TRUE, lam = 0.1){
  # Generate random community data for testing community phylogenetic
  #  methods
  #
  # Args:
  #   phylo: phylogeny
  #   nsites: number of sites
  #   nspp: number of species for each site
  #   pa: presence/absence data?
  #   lam: lambda for poisson distribution if pa is False
  #
  # Return:
  #   a matrix of community data
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  if (pa) {
    for (i in 1:nsites) {
      output[i, ] <- sample(c(rep(1, nspp), rep(0, ntips - nspp)))
    }
  } else {
    for (i in 1:nsites) {
      output[i, ] <- sample(rpois(ntips, lam))
    }
  }
  return(output)
}

split0 <- function(r = c(1, 10) , n = 2) {
  # evenly split a range into n
  #
  # args:
  #   r: range to be split, must be numerical, can be vector.
  #   n: times to split r
  #
  # returns a vector of splits
  output <- rep(NA, n)
  output[1] <- 0
  division <- 1/(n-1)
  for (i in 2:n) {
    output[i] <- output[i-1] + division
  }
  return(round((output * (max(r) - min(r))) + min(r)))
}

genCommData <- function(phylo, focal, nsites, nspp, fact = 1) {
  # Generate clustered/overdispersed data for testing community
  #  phylogenetic analyses.
  #
  # Args:
  #  phylo - phylo.object, on which the community data will be based
  #  focal - numeric index, indicating which tip to perform cluster/dispersion
  #  nsites - the number of sites to create
  #  nspp - the number of species to occur at each site
  #  clust - if TRUE, data will be clustered else overdispersed
  #  fact - scaling coefficient: larger the number the more pronounced the
  #   effect by a power law (hence 0 not allowed)
  #
  # Return:
  #  a matrix of community data
  invertVector <- function(dists) {
    u <- sort(unique(dists), TRUE)
    s <- sort(u)
    probs <- rep(NA, length(dists))
    for (i in 1:length(u)) {
      probs[u[i] == dists] <- s[i]
    }
    return (probs)
  }
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  pd.dists <- cophenetic.phylo(phylo)
  focal.dists <- pd.dists[ , focal] + 1 # avoid 0
  if (fact > 0) focal.dists <- invertVector(focal.dists)
  probs <- focal.dists^abs(fact)
  probs <- probs/sum(probs)
  for (i in 1:nsites) {
    output[i, ] <- ifelse(1:ntips %in% sample(1:ntips, nspp,
                                              prob = probs), 1, 0)
  }
  return (output)
}

evenCommData <- function(phylo, nsites, nspp) {
  # Generate evenly distributed community phylogenetic data based on a
  #  phylogeny. Takes the phylogeny, calculates cummulative PD from random
  #  taxon, cuts at even intervals based on nspp to maximise distance.
  #
  # Args:
  #  phylo - phylo.object, on which community will be based
  #  nsites - number of sites
  #  nspp - number of species
  #
  # Depends:
  #  split0() - for evenly cutting up the cumsum.pd.dists
  #
  # Return:
  #  a matrix of community data 
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(0, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  pd.dists <- cophenetic.phylo(phylo)
  for (i in 1:nsites) {
    focal.pd.dists <- pd.dists[ , sample(1:ntips, 1)]
    # select focal taxon at random
    splits <- split0(1:sum(focal.pd.dists), nspp)
    cumsum.pd.dists <- cumsum(focal.pd.dists)
    for (j in 1:nspp) {
      pull <- which(abs(cumsum.pd.dists - splits[j])
                    == min(abs(cumsum.pd.dists - splits[j])))
      output[i, pull] <- 1
    }
  }
  return(output)
}

phyloRare <- function(comm.data, phylo, metric = 'PD', nrands = 1000,
                      proportion = FALSE) {
  mat <- as.matrix(comm.data)
  out.phylo <- rep(NA, nrow(mat))
  out.se <- out.phylo
  raremax <- min(rowSums(mat))
  for (i in 1:nrow(mat)) {
    site.taxa <- colnames(mat)[mat[i,] > 0]
    phylo.rand <- rep(NA, nrands)
    for (j in 1:nrands) {
      taxa <- sample(site.taxa, raremax)
      if(metric == 'PD') {
        phylo.rand[j] <- PD(phylo, taxa, proportion = proportion)
      } else {
        stop('Unknown phylo metric.')
      }
    }
    out.phylo[i] <- mean(phylo.rand)
    out.se[i] <- sd(phylo.rand)/sqrt(length(phylo.rand))
  }
  return(list(phylo = out.phylo, se = out.se))
}

## Jun's Functions
# PUT YOUR FUNCTIONS HERE JUN!