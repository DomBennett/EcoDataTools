## Dominic Bennett
## 30/06/2013
## Generate random meta-data

## Functions
source("EcoDataTools.R")

## Dirs
dir <- file.path('cookbook', 'recipe1')

## Import Tree
big.phylo <- read.tree(file.path(dir, "big_phylo.txt"))

## Generate Data
# Parameters
nstudies <- 30 # number of studies to generate
nsites <- 40 # number of sites in each study
nspecies <- 20 # the size of the phylgeny in each study
nhtypes <- 2 # the number of habitat types to be compared

# Running
data <- list()
phylos <- list()
drop.n <- length(big.phylo$tip.label) - nspecies
facts <- split0(-2:2, nhtypes)
mean.incids <- split0(5:(nspecies/2), nhtypes)
for (i in 1:nstudies) {
  study.phylo <- drop.tip(big.phylo, sample(big.phylo$tip.label, drop.n))
  study.phylo$edge.length <- study.phylo$edge.length *
    sample(c(0.01, 0.1, 1, 10, 100), 1)
  study.abun <- mean.incids * sample(c(1, 2, 3, 20, 100), 1)
  study.focal <- sample(1:nspecies, 1)
  site.types <- sort(sample(1:nhtypes, nsites, replace = TRUE))
  site.nsites <- table(site.types)
  site.data <- genCommData(study.phylo, study.focal, facts[1], mean.incids[1],
                           study.abun[1], site.nsites[1])
  for (j in 2:nhtypes) {
    temp.site.data <- genCommData(study.phylo, study.focal, facts[j], mean.incids[j],
                                  study.abun[j], site.nsites[j])
    site.data <- rbind(site.data, temp.site.data)
  }
  rownames(site.data) <- site.types
  colnames(site.data) <- study.phylo$tip.label
  site.data <- as.matrix(na.omit(site.data))
  data <- c(data, list(site.data))
  phylos <- c(phylos, list(study.phylo))
}

## Sanity check
# plotting
study.samples <- sample(1:nstudies, 4)
par(mfrow = c(2,2))
for (i in 1:4) {
  plotComm(data[[study.samples[i]]], phylos[[study.samples[i]]],
           groups = rownames(data[[study.samples[i]]]))
}
# sums
#lapply(data, rowSums)
#lapply(data, function(x) sum(x > 0))
#lapply(phylos, function(x) sum(x$edge.length))