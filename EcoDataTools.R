## No copyright, no warranty
## Dominic John Bennett
## 25/05/2013
## Library of functions created for my Masters project (MRes Biodiversity Informatics
##  and Genomics , Imperial College London) of an ecological and phylogenetic slant

## Coordinate manipulation

DegToDec <- function(data) {
  # Take 'coordinate strings' in degrees (D M (S)) and convert to decimals.
  # Expects coordinate string such: ##'##'##'
  #   - where: # denotes any number of any length, ' denotes any non-numeric(s).
  # Function will stop if more than 3 numbers in coordinate string
  #
  # Args:
  #   data: a matrix or dataframe of coordinates (Latitude and Longitude).
  #
  # Returns:
  #   A dataframe of converted coordinates.
  
  # input checking
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("Input must be a matrix or dataframe, such:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  if(ncol(data)!= 2) {
    stop("Input must have two columns, such:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  two.numbers <- apply(data, 1:2, function(x) grepl("([^0-9][0-9]|[0-9][^0-9])", x))
  if(!all(two.numbers)) {
    print("Invalid data entries:")
    print(data[!two.numbers])
    stop("Input must be coordinate character strings with at least 2 numbers, e.g.:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  
  # vector inputs
  x <- data[ ,1] #latitude is x
  y <- data[ ,2] #longitude is y
  
  # find start and end of numbers in coordinate strings
  index.x1 <- gregexpr("[0-9\\-\\.]+", x)
  index.x2 <- gregexpr("([0-9][^0-9]|[0-9]$)", x)
  index.y1 <- gregexpr("[0-9\\-\\.]+", y)
  index.y2 <- gregexpr("([0-9][^0-9]|[0-9]$)", y)
  
  # vector outputs
  decimal.x <- rep(NA, nrow(data))
  decimal.y <- rep(NA, nrow(data))
  
  # loop through each row of coordinate string:
  #   - take the first number in string, add the second number/60
  #     and then the third/3600
  #    - stop if more than 3 numbers found.
  error.message1 <- "\nMore than 3 numbers in:\n data["
  error.message2 <- "\nTry checking for non-numeric characters in numbers e.g.:\n
  Here, there is an extra space between the 9 and 2 -- 51°28\'46.9 2\""
  for(i in 1:nrow(data)) {
    decimal.x[i] <- as.numeric(substr(x[i], index.x1[[i]][1], index.x2[[i]][1]))
    for(j in 2:length(index.x1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 1, ",", i,"] -- ",x[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.x <- as.numeric(substr(x[i], index.x1[[i]][j], index.x2[[i]][j]))
      decimal.x[i] <- decimal.x[i] + (next.x/division)
    }
    decimal.y[i] <- as.numeric(substr(y[i], index.y1[[i]][1], index.y2[[i]][1]))
    for(j in 2:length(index.y1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 2, ",", i,"] -- ",y[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.y <- as.numeric(substr(y[i], index.y1[[i]][j], index.y2[[i]][j]))
      decimal.y[i] <- decimal.y[i] + (next.y/division)
    }
  }
  
  # control for S and W
  decimal.x <- decimal.x * ifelse(grepl("(S|s)", x), -1, 1)
  decimal.y <- decimal.y * ifelse(grepl("(W|w)", y), -1, 1)
  
  # generate output
  decimals <- data.frame(Latitude = decimal.x, Longitude = decimal.y)
  return(decimals)
}

## Phylogenetic Community Analysis

## Libraries
library(picante)
library(caper)
library(lme4)

## Function
plotComm <- function(comm.data, phylo, groups = rep(1, nrow(comm.data))){
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
  variable.max <- (nrow(comm.data) * 10) #+ 150
  #variable.max <- ifelse(variable.max > 200, variable.max, 200) #min is 200
  plot(phylo, no.margin = T, show.tip.label = F, x.lim = c(0, variable.max))
  
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- comm.data/max(comm.data)
  
  # loop init
  ntips <- length(phylo$tip.label)
  spacing <- 10
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
    spacing <- spacing + 10
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