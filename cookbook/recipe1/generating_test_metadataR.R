## Dominic Bennett
## 30/06/2013
## Generate random meta-data

## Functions
source("EcoDataTools.R")

## Set seed
set.seed(4)

## Dirs
dir <- file.path('cookbook', 'recipe1')
outdir <- file.path('wiki', 'cookbook')

## Import Tree
big.phylo <- read.tree(file.path(dir, "big_phylo.txt"))

## Generate Data
# Parameters
nstudies <- 30 # number of studies to generate
nsites <- 40 # number of sites in each study
nspecies <- 20 # the size of the phylgeny in each study
nhstates <- 2 # the number of habitat states to be compared
facts <- split0(-2:2, nhstates) # the facts between habitat states
mean.incids <- sample(5:(nspecies/2), nhstates) # incidences between states

# Running
data <- list()
phylos <- list()
drop.n <- length(big.phylo$tip.label) - nspecies
for (i in 1:nstudies) {
  study.phylo <- drop.tip(big.phylo, sample(big.phylo$tip.label, drop.n))
  study.phylo$edge.length <- study.phylo$edge.length *
    sample(c(0.01, 0.1, 1, 10, 100), 1)
  study.abun <- mean.incids * sample(c(1, 2, 3, 20, 100), 1)
  study.focal <- sample(1:nspecies, 1)
  site.types <- sort(sample(1:nhstates, nsites, replace = TRUE))
  site.nsites <- table(site.types)
  site.data <- genCommData(study.phylo, study.focal, facts[1], mean.incids[1],
                           study.abun[1], site.nsites[1])
  for (j in 2:nhstates) {
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

# plotting
study.samples <- sample(1:nstudies, 4)
png("r1_random_data_sc.png", width = 720)
par(mfrow = c(2,2))
for (i in 1:4) {
  plotComm(data[[study.samples[i]]], phylos[[study.samples[i]]],
           groups = rownames(data[[study.samples[i]]]))
}
dev.off()
# sums
#lapply(data, rowSums)
#lapply(data, function(x) sum(x > 0))
#lapply(phylos, function(x) sum(x$edge.length))
save(list = c("data", "phylos"), file = file.path(dir, "test_metadata.RData"))