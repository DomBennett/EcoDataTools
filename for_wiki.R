## Dominic Bennett
## 21/07/2013
## Generate test data and images for wiki

## Dirs
source('EcoDataTools.R')

## set.seed
set.seed(4)

## extractEdges()
# showing the different types
plotEdges <- function (phylo, edges) {
  tip.cols <- ifelse(phylo$tip.label %in% taxa, "black", "grey")
  edge.lties <- ifelse(1:nrow(phylo$edge) %in% edges, 1, 3)
  plot.phylo(phylo, edge.lty = edge.lties, tip.color = tip.cols,
             show.tip.label = TRUE)
}
phylo <- read.tree(file.path("wiki","Phylocom_phylo.tre"))
phylo <- drop.tip(phylo, phylo$tip.label[1:12])
taxa <- phylo$tip.label[c(20,16:13,5:6)]
png(filename = file.path("wiki", "PD.png"), width = 1000, height = 800)
split.screen(c(2,2))
screen(1)
edges <- extractEdges(phylo, taxa, type = 1)
plotEdges(phylo, edges)
mtext("Type 1: phylogeny of the taxa", line = 2)
mtext(paste("PD:", sum(phylo$edge.length[edges]), "   min: 2"))
screen(2)
edges <- extractEdges(phylo, taxa, type = 2)
plotEdges(phylo, edges)
mtext("Type 2: branches from the taxa tips to the terminal node", line = 2)
mtext(paste("PD:", sum(phylo$edge.length[edges]), "   min: 1"))
screen(3)
edges <- extractEdges(phylo, taxa, type = 3)
plotEdges(phylo, edges)
mtext("Type 3: branches uniquely represented by the taxa", line = 2)
mtext(paste("PD:", sum(phylo$edge.length[edges]), "   min: 1"))
close.screen(al = TRUE)
dev.off()
# example usage
phylo <- read.tree(file.path("wiki","Phylocom_phylo.tre"))
taxa1 <- sample(phylo$tip.label, 16, prob = 1:length(phylo$tip.label))
taxa2 <- phylo$tip.label[!phylo$tip.label %in% taxa1]
edges1 <- extractEdges(phylo, taxa1, type = 3)
edges2 <- extractEdges(phylo, taxa2, type = 3)
tip.cols <- ifelse(phylo$tip.label %in% taxa1, "darkgreen", "darkred")
edges.cols <- rep("darkblue", nrow(phylo$edge))
edges.cols[1:nrow(phylo$edge) %in% edges1] <- rep("darkgreen", length(edges1))
edges.cols[1:nrow(phylo$edge) %in% edges2] <- rep("darkred", length(edges2))
png(filename = file.path("wiki", "PD_signal.png"), width = 1000, height = 800)
plot.phylo(phylo, edge.color = edges.cols, tip.color = tip.cols,
           show.tip.label = TRUE)
dev.off()

## Community plot images
phylo <- read.tree(file.path("wiki", "Phylocom_phylo.tre"))
png(filename = file.path("wiki","randCommData.png"))
plotComm(randCommData(phylo, 10, 5, pa = F), phylo)
dev.off()
png(filename = file.path("wiki", "evenCommData.png"))
plotComm(evenCommData(phylo, 10, 5), phylo)
dev.off()
png(filename = file.path("wiki", "genCommData_probs.png"), width = 960)
par(mar = c(5, 4, 4, 2), mfrow = c(1,2))
for (i in 1:2) {
  clust <- ifelse(i == 1, TRUE, FALSE)
  title <- ifelse(i ==1, '\'Clusteredness\'', '\'Overdispersed\'')
  focal.dists <- seq(1,10,1)
  fact <- c(0, 0.5,1,2,3)
  cols <- rainbow(length(fact), alpha = 0.9)
  plot(seq.int(0, 0.4, length.out = length(focal.dists)) ~ focal.dists,
       type = 'n', xlab = "Phylogenetic Distance from Focal",
       ylab = "P(Co-occurence)")
  for (j in 1:length(fact)) {
    temp.dists <- sort(focal.dists, decreasing = clust)
    probs <- temp.dists^fact[j]
    probs <- probs/sum(probs)
    lines(probs ~ focal.dists, col = cols[j], lwd = 2)
  }
  legend("topright", legend = fact, col = cols, pch = 19,
         title = title)
}
dev.off()
png(filename = file.path("wiki", "genCommData.png"))
type1 <- genCommData(phylo, focal = 16, fact = 1, mean.incid = 5, nsites = 10)
type2 <- genCommData(phylo, focal = 16, fact = 10, mean.incid = 5, nsites = 10)
type3 <- genCommData(phylo, focal = 16, fact = 0, mean.incid = 5, nsites = 10)
type4 <- genCommData(phylo, focal = 16, fact = -1, mean.incid = 5, nsites = 10)
type5 <- genCommData(phylo, focal = 16, fact = -10, mean.incid = 5, nsites = 10)
all <- rbind(type1, type2, type3, type4, type5)
types <- paste0('t', rep(1:5, each = 10))
plotComm(all, phylo, groups = types)
dev.off()
png(filename = file.path("wiki", "genCommData_abuns.png"))
type1 <- genCommData(phylo, focal = 16, fact = 1, mean.incid = 5, nsites = 10,
                     mean.abun = 20)
type2 <- genCommData(phylo, focal = 16, fact = -1, mean.incid = 5, nsites = 10,
                     mean.abun = 20)
all <- rbind(type1, type2)
types <- paste0('t', rep(1:5, each = 10))
plotComm(all, phylo, groups = types)
dev.off()

## Testing meanPhylo() 21/07/2013
# Generate distribution
nphylos <- 30
standard.dev <- 5
ori.phylo <- stree(32, 'balanced')
ori.phylo$edge.length <- rep(1, nrow(ori.phylo$edge))
phylo.dist <- list()
rand.phy.sizes <- ceiling(abs(rnorm(nphylos, sd = standard.dev)))
blength.modifier <- runif(nphylos, min = 0, max = 2)
for (i in 1:nphylos) {
  temp.phylo <- drop.tip(ori.phylo, sample(ori.phylo$tip.label,
                                            rand.phy.sizes[i]))
  temp.phylo$edge.length <- temp.phylo$edge.length * blength.modifier[i]
  phylo.dist <- c(phylo.dist, list(temp.phylo))
}
# Calculate mean tree
mean.phylo <- meanPhylo(phylo.dist)
topodist <- dist.topo(mean.phylo, unroot(ori.phylo))
bdist <- dist.topo(mean.phylo, unroot(ori.phylo), method = 'score')
png(filename = file.path("wiki", "meanPhylo_example.png"))
split.screen(c(1,2))
screen(1)
plot(mean.phylo, type = 'unrooted')
mtext("Mean phylogeny")
mtext(paste0("Topological Distance = ", topodist), side = 1)
mtext(paste0("Branch Len. Distance = ", signif(bdist, 2)), side = 1, line = 1)
screen(2)
plot(unroot(ori.phylo), type = 'unrooted')
mtext("Original phylogeny")
close.screen(all.screens = TRUE)
dev.off()

## Testing genNullDist() 22/07/2013
phylo <- stree(n = 32, "balanced", tip.label = as.character(1:32))
phylo$edge.length <- rep(1, length(phylo$edge)) # need to add lengths!
set.seed(4) # 2, 4, 55
scenario0 <- genCommData(phylo, focal = 26, fact = 0, mean.incid = 8,
                         mean.abun = 16, nsites = 40)
scenario1 <- rbind(genCommData(phylo, focal = 26, fact = 0, mean.incid = 12,
                               mean.abun = 24, nsites = 20),
                   genCommData(phylo, focal = 26, fact = 0, mean.incid = 2,
                               mean.abun = 4, nsites = 20))
scenario2 <- rbind(genCommData(phylo, focal = 26, fact = 0, mean.incid = 6,
                               mean.abun = 12, nsites = 20),
                   genCommData(phylo, focal = 26, fact = 2, mean.incid = 6,
                               mean.abun = 12, nsites = 20))
scenario3 <- rbind(genCommData(phylo, focal = 26, fact = -2, mean.incid = 6,
                               mean.abun = 12, nsites = 20),
                   genCommData(phylo, focal = 26, fact = 2, mean.incid = 6,
                               mean.abun = 24, nsites = 20))
scenario4 <- rbind(genCommData(phylo, focal = 26, fact = 2, mean.incid = 6,
                               mean.abun = 12, nsites = 20),
                   genCommData(phylo, focal = 6, fact = 2, mean.incid = 6,
                               mean.abun = 12, nsites = 20))
scenario5 <- rbind(genCommData(phylo, focal = 15, fact = 2, mean.incid = 4,
                               mean.abun = 8, nsites = 20) +
                     genCommData(phylo, focal = 17, fact = 2, mean.incid = 4,
                                 mean.abun = 8, nsites = 20),
                   genCommData(phylo, focal = 28, fact = 2, mean.incid = 4,
                               mean.abun = 8, nsites = 20) +
                     genCommData(phylo, focal = 4, fact = 2, mean.incid = 4,
                                 mean.abun = 8, nsites = 20))
htypes <- as.character(c(rep(1, 20), rep(2, 20)))
png(filename = file.path("wiki", "genNullDist_scenarios.png"))
par(mfrow = c(2,3))
plotComm(scenario0, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 0", 3, 1)
plotComm(scenario1, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 1", 3, 1)
plotComm(scenario2, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 2", 3, 1)
plotComm(scenario3, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 3", 3, 1)
plotComm(scenario4, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 4", 3, 1)
plotComm(scenario5, phylo, no.margin = FALSE, groups = htypes)
mtext("Scenario 5", 3, 1)
dev.off()