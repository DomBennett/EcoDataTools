## No copyright, no warranty
## Dominic John Bennett & Junying Lim
## Last update: 25/06/2013
## Functions for PD analysis BIG project 2013

## Libraries
library(picante) # Phylogenetic analytical tools
library(caper)
library(lme4)
library(stringr) #String manipulation tools
library(RCurl) ##GET and POST to URL queries functionality

## Dom's Functions:
extractEdges <- function(phylo, taxa, type = 1) {
  # Extract edges from a phylo object using 1 of 3 methods
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  taxa: vector of taxon names
  #  type:
  #     1 -- phylogeny consisting solely of the taxa, default
  #     2 -- edges from taxon tips to terminal node
  #     3 -- edges unique to taxa
  #
  # Return:
  #  vector of edges
  if (!type %in% c(1,2,3)) {
    stop("Type must be an integer: 1, 2 or 3.")
  }
  if (!is.vector(taxa) | !is.character(taxa)) {
    stop("Invalid or no taxa given.")
  }
  if (length(taxa) == length (phylo$tip.label)){
    return(phylo$edge)
  }
  if (type == 1 & length (taxa) == 1){
    stop("length(taxa) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  edges <- match (match (taxa, phylo$tip.label), phylo$edge[,2])
  end.nodes <- phylo$edge[edges, 1]
  term.node <- length (phylo$tip.label) + 1
  if (all(end.nodes %in% term.node)) {
    return(edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (phylo$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, phylo$edge[,2])
          end.node <- phylo$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, phylo$edge[,2])
        end.node <- phylo$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
      }
    }
    return (edges)
  }
}

commPD <- function(phylo, comm.data, type = 1, min.spp = 2,
                   taxon.names = colnames(comm.data)) {
  # Calculate Faith's PD for a community given a community phylogeny and a community
  #  matrix. Depends on extractEdges().
  #
  # Args:
  #  phylo: community phylogeny (ape class)
  #  comm.data: matrix of community data, species as cols sites as rows
  #  type: which branches to include for calculating PD, either 1, 2 or 3,
  #   default 2 -- see extractEdges()
  #  min.spp: minimum number of species to include for calculating PD, default 2
  #  taxon.names: if taxon.names are not the column names of comm.data specify
  #   here, else ignore.
  #
  # Return:
  # vector of PDs
  if (type == 1 & min.spp < 2) {
    stop("Cannot compute type 2 PD for fewer than 2 species")
  }
  calcPD <- function (row) {
    taxa <- taxon.names[as.logical(row)]
    edges <- extractEdges(phylo, taxa, type)
    return (sum(phylo$edge.length[edges]))
  }
  # add site names if none
  if (is.null(rownames(comm.data))) {
    rownames(comm.data) <- 1:nrow(comm.data)
  }
  # convert to incidence
  presences <- comm.data[comm.data > 0]
  comm.data[comm.data > 0] <- rep(1, length(presences))
  # drop sites with too few species
  row.sums <- rowSums(comm.data) > min.spp
  site.names <- rownames(comm.data)[row.sums]
  comm.data <- comm.data[row.sums, ]
  # calculate PD by site
  pds <- as.matrix(apply(comm.data, 1, calcPD))
  rownames(pds) <- site.names
  return (pds)
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
  edges <- extractEdges(phylo, phylo$tip.label[1], type = 2)
  phylo$edge.length <- phylo$edge.length/sum(phylo$edge.length[edges])
  # make all phylos the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(comm.data) < 20) {
    variable.max <- 1 + nrow(comm.data)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(comm.data)/10
    spacing.i <- 0.1 - 1/nrow(comm.data)
    # some stupid f***ing bizarre eqn! Why is it not linear!?!?!?!
    spacing.start <- 0.5 + spacing.i
  }
  plot(phylo, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max))
  
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- comm.data/max(comm.data)
  
  # loop init
  ntips <- length(phylo$tip.label)
  spacing <- spacing.start
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
    spacing <- spacing + spacing.i
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

genCommData <- function(phylo, focal, fact = 1, mean.incid, mean.abun = FALSE,
                        nsites = 1) {
  # Generate clustered/overdispersed data for testing community
  #  phylogenetic analyses. Return either incidence or abundance data.
  #
  # Args:
  #  phylo - phylo.object, on which the community data will be based
  #  focal - numeric index, indicating which tip to perform cluster/dispersion
  #  fact - scaling coefficient: larger the number the more pronounced the
  #   effect by a power law (hence 0 not allowed)
  #  mean.incid - the mean incidence of species in the community
  #  mean.abun - the mean abundance per site, if given output will be abundances
  #  nsites - number of sites
  #
  # Return:
  #  a matrix of community data
  invertVector <- function(dists) {
    # for reversing the probs for overdispersion
    u <- sort(unique(dists), TRUE)
    s <- sort(u)
    probs <- rep(NA, length(dists))
    for (i in 1:length(u)) {
      probs[u[i] == dists] <- s[i]
    }
    return (probs)
  }
  genAbuns <- function(row) {
    # for generating abundances for each row
    out.row <- rep(0, ntips)
    temp.probs <- probs
    temp.probs[row < 1] <- 0
    abundance <- abs(ceiling(rnorm(1, mean = mean.abun)))
    if (abundance == 0) {
      return (out.row)
    } else {
      abuns <- sample(1:ntips, abundance, prob = temp.probs, replace = TRUE)
      abuns <- table(abuns)
      out.row[as.numeric(names(abuns))] <- abuns
      return (out.row)
    }
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
    incidence <- abs(ceiling(rnorm(1, mean = mean.incid))) # avoid negative numbers
    output[i, ] <- ifelse(1:ntips %in% sample(1:ntips, incidence,
                                              prob = probs), 1, 0)
  }
  if (mean.abun != FALSE) {
    for (i in 1:nsites) {
      output[i, ] <- genAbuns(output[i, ])
    }
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

phyloRarefy <- function(comm.data, phylo, samp, metric = 'PD', nrands = 2000,
                        type = 1) {
  # Monte-Carlo method for rarefiying community data using common phylogenetic
  #  metrics. Depends on picante.
  #
  # Args:
  #  comm.data: community matrix, species as cols and sites as rows
  #  phylo: phylogeny (ape class)
  #  samp: sample size (minimum incidence of abundance)
  #  metric: either PD, PSV, PSE, PSR, PSC or PSD, default PD
  #  nrands: the number of iterations, default 2000
  #  types: if PD, type of PD to calculate, default 1
  #
  # Returns:
  #  matrix of rarefied PD by site and standard error
  # TODO(27/06/2013): unit test this and add test to wiki
  if (any(rowSums(comm.data) < samp)) {
    stop("Some sites have incidence/abundance less than samp")
  }
  # is it incidence data or abundance?
  if (any(comm.data > 1)) {
    abun <- TRUE
  } else {
    abun <- FALSE
  }
  # define internal functions for rapid vectorization
  # a lot of function definitions but avoids if statements for speed!!
  pluck <- function(row) {
    # this takes a community matrix and randomly distributes the occurences
    # of species using samp works for both incidence and abundance
    samples <- table(sample(1:length(row), samp, prob = row, replace = abun))
    sampled.row <- rep(0, length(row))
    sampled.row[as.numeric(names(samples))] <- samples
    return (sampled.row)
  } 
  PDrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(commPD(phylo, plucked, type, samp)))
  }
  PSVrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psv(plucked, phylo)[ ,1]))
  }
  PSErarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(pse(plucked, phylo)[ ,1]))
  }
  PSRrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psr(plucked, phylo)[ ,1]))
  }
  PSCrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psc(plucked, phylo)[ ,1]))
  }
  PSDrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psd(plucked, phylo)[ ,1]))
  }
  calcSE <- function (x) {
    return (sd(x)/sqrt(length(x)))
  }
  # now for the actual calculations ...
  mat <- as.matrix(comm.data)
  if (metric == 'PD') {
    mat.phylo <- sapply(1:nrands, PDrarefy)
  } else if (metric == 'PSV') {
    mat.phylo <- sapply(1:nrands, PSVrarefy)
  } else if (metric == 'PSE') {
    mat.phylo <- sapply(1:nrands, PSErarefy)
  } else if (metric == 'PSR') {
    mat.phylo <- sapply(1:nrands, PSRrarefy)
  } else if (metric == 'PSC') {
    mat.phylo <- sapply(1:nrands, PSCrarefy)
  } else if (metric == 'PSD') {
    mat.phylo <- sapply(1:nrands, PSDrarefy)
  } else {
    stop('Unknown metric. Use: PD, PSV, PSE, PSR, PSC or PSD')
  }
  out.phylo <- apply(mat.phylo, 1, mean)
  out.se <- apply(mat.phylo, 1, calcSE)
  out <- cbind(out.phylo, out.se)
  colnames(out) <- c(metric, paste0(metric, '.se'))
  return(out)
}

## Jun's Functions
phylotraitdist <- function(native, alien, phy, trait){
               #Print statements
	print(paste(dim(trait)[2], " traits found...", sep = ""))
	#Diagnostic checks on phylogeny
	taxon <- c(native, alien)
	if(sum(taxon %in% testphy$tip.label) != length(taxon)){
		print("Phylogeny does not contain all the taxa. Missing taxa will be excluded.")
	}
	tip <- phy$tip.label[!phy$tip.label %in% taxon] 
	trimmedphy <- drop.tip(phy, tip = tip) #Removing taxa that are not found
	#Diagnostic check on traits 
	if(sum(taxon %in% rownames(trait) != length(taxon))){
		print("Trait dataframe does not contain all the taxa. Missing taxa will be excluded.")
	}
	trait <- trait[rownames(trait) %in% taxon,] #Removing taxa that are not found
	alien <- intersect(intersect(alien, phy$tip.label), rownames(trait)) #Alien taxa found
	#Results vectors
	pnnd <- NULL
	traitdiff <- NULL
	#Distance matrices
	phydist <- cophenetic(trimmedphy) #Calculating phylogenetic distance matrix
	order <- match(taxon, colnames(phydist))
	phydist <- phydist[order,order]
	traitdist <- vector("list", length(colnames(trait)))
	for(j in 1:length(colnames(trait))){
		if(is.factor(trait[,j]) == TRUE){
			levels(trait[,j]) <- 1:length(levels(trait[,j])) #Convert levels of trait into numerical values
			temp <- decostand(dist(trait[,j], method = "manhattan", upper = TRUE, diag = TRUE), method = "pa")
			temp <- as.matrix(temp)
			rownames(temp) <- colnames(temp) <- rownames(trait)
			order <- match(taxon, colnames(temp))
			traitdist[[j]] <- temp[order,order]
		} else {
			temp <- dist(trait[,j], method = "manhattan", upper = TRUE, diag = TRUE)
			temp <- as.matrix(temp)
			rownames(temp) <- colnames(temp) <- rownames(trait)
			order <- match(taxon, colnames(temp))
			traitdist[[j]] <- temp[order,order]
		}
	} #Produces list of trait distance matrices
	#Find non-native index on cophenetic
	exclude <- colnames(phydist) %in% alien
	for(i in alien){
		temp <- phydist[!exclude,] #Remove all non-natives from matrix
		closest.nat <- which.min(temp[, colnames(phydist) %in% i]) #Find closest native by looking at minimum phylogenetic distance
		pnnd <- c(pnnd, temp[closest.nat, colnames(phydist) %in% i])
		for(j in 1:length(traitdist)){
			temp <- traitdist[[j]][!exclude,]
			traitdiff <- append(traitdiff, temp[closest.nat , colnames(traitdist[[j]]) %in% i])
		}
	}
	pnndresults <- data.frame(taxon = alien, pnnd = pnnd)
	traitresults <- matrix(data = traitdiff, nrow = length(alien), ncol = length(colnames(trait)), byrow = TRUE)
	colnames(traitresults) <- paste(colnames(trait), ".diff", sep = "")
	return(cbind(pnndresults, traitresults))
}

taxa2phylomatic <- function(taxa, output.dir, output.name, binom, storedtree){
  storedtree= "R20120829"
  informat = outformat = "newick"
  method = "phylomatic"
  url = 'http://phylodiversity.net/phylomatic/pmws'
  taxaformat = "slashpath"
  constraint.dir = file.path(output.dir, paste(output.name, "_tree.nwk", sep =""))
  unmatched.dir = file.path(output.dir, paste(output.name, "_unmatched.txt", sep = ""))
  # Phylomatic does not accept line breaks
  taxa <- str_replace_all(taxa, " ", "_")
  binom <- str_replace_all(binom, " ", "_")
  taxa <- paste(taxa, collapse = "\n")
  taxastring <- str_replace_all(str_replace_all(taxa, "/", "%2F"), "\n", "%0D")
  # Writing phylomatic query string
  phylomat.dat1 = paste("storedtree=", storedtree, "&", "informat=", informat, "&", "method=", method, "&", sep = "") 
  phylomat.dat2 = paste("taxaformat=", taxaformat, "&", "clean=true", "&", "outformat=", outformat, "&", sep = "")
  phylomat.dat3 = paste("taxa=", taxastring, sep = "")
  fullstr <- paste(phylomat.dat1, phylomat.dat2, phylomat.dat3, sep = "")
  # Parsing string to phylomatic using Curl
  print("Generating constraint tree using phylomatic...")
  temp <- postForm(uri = url, curl = getCurlHandle(), .opts = list(postfields = fullstr))
  # Cleaning up names (because tip labels may change in capitalization)
  for(i in 1:length(binom)){
    print(binom[i])
    temp <- gsub(pattern = binom[i], x = temp, replacement = binom[i], ignore.case = TRUE)
  }
  # Removing unmatched taxa
  matched <- gsub("\\[.+\\]", "", temp)
  # Export newick file
  print(paste("Writing constraint tree to ", constraint.dir))
  constraintnwk <-file(constraint.dir)
  writeLines(matched, constraintnwk)
  close(constraintnwk)
  # Exporting unmatched taxa to separte document
  m <- regexpr("\\[.+\\]", temp)
  unmatched <- regmatches(temp, m)
  print(unmatched)
  print(paste("Writing unmatched taxa to ", unmatched.dir))
  unmatchedfile <- file(unmatched.dir)
  writeLines(unmatched, unmatchedfile)
  close(unmatchedfile)
  #Defunct code. Does not handle large Request-URIs
  #urlquery <- paste(url, "?", phylomat.dat1, phylomat.dat2, phylomat.dat3, sep = "")
  #tt <- getURLContent(urlquery, curl = getCurlHandle())
  #phylomatic_nwk <- tt[[1]]
  #return(phylomatic_nwk)
  
  ## TO DO##
  #Argument for phylomatic tree to use
  #Argument to supply user trees
}

#newick and taxo2newick (Author: J. Lim)
#newick returns a list of taxa names into a polytomic newick tree
#taxo2newick returns a dataframe of taxa names into a newick tree,
#assuming monophyly in the penultimate taxonomic hierarchy (e.g., if data
#is on genera, then confamilial genera will be a monophyletic polytomy)
newick <- function(x){
  base <- x[1]
  if(length(x) > 1){
    for(i in 2:length(x)){
      base <- paste(base, x[i], sep = ", ")
    }
    return(paste("(", base, ")", sep = ""))
  } else {
    return(base)
  }
}

taxo2newick <-function(x, backbone = TRUE, backbonenewick, hightaxo, lowtaxo){
  taxo <- c("order", "class", "family")
  if(sum(match(tolower(names(x)), taxo, nomatch=0)) == 0){
    stop("Dataframe does not have taxonomic classifications")
  }else{
    hightaxoindex <- match(hightaxo, names(x))   #Find index of desired taxa
    hightaxolist <- as.vector(unique(x[, hightaxoindex]))   #Figure out what are the levels in the particular
    if(backbone == TRUE){
      base <- backbonenewick
    } else {
      base <- newick(hightaxolist)
    }
    for(i in 1:length(hightaxolist)){
      temp <- subset(x, x[,hightaxoindex] == hightaxolist[i]) #Subsets dataframe based on higher taxonomy to find the lower taxonomic levels
      lowtaxolist <- as.vector(unique(temp[,match(lowtaxo, names(temp))]))
      temp2 <- newick(lowtaxolist)
      base <- gsub(hightaxolist[i], temp2, base)
      print(base)
    }
  }
  return(base)
}
