## Dominic Bennett
## 26/06/2013
## Generate test data and images for wiki

## Dirs
source('EcoDataTools.R')

## Calculating Faith's PD
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
mtext(paste("PD:", sum(phylo$edge.lengths[edges]), "   min: 2"))
screen(2)
edges <- extractEdges(phylo, taxa, type = 2)
plotEdges(phylo, edges)
mtext("Type 2: branches from the taxa tips to the terminal node", line = 2)
mtext(paste("PD:", sum(phylo$edge.lengths[edges]), "   min: 1"))
screen(3)
edges <- extractEdges(phylo, taxa, type = 3)
plotEdges(phylo, edges)
mtext("Type 3: branches uniquely represented by the taxa", line = 2)
mtext(paste("PD:", sum(phylo$edge.lengths[edges]), "   min: 1"))
close.screen(all.screens = TRUE)
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
type1 <- genCommData(phylo, focal = 16, nsites = 10, nspp = 5, fact = 1)
type2 <- genCommData(phylo, focal = 16, nsites = 10, nspp = 5, fact = 10)
type3 <- genCommData(phylo, focal = 16, nsites = 10, nspp = 5, fact = 0)
type4 <- genCommData(phylo, focal = 16, nsites = 10, nspp = 5, fact = -1)
type5 <- genCommData(phylo, focal = 16, nsites = 10, nspp = 5, fact = -10)
all <- rbind(type1, type2, type3, type4, type5)
types <- paste0('t', rep(1:5, each = 10))
plotComm(all, phylo, groups = types)
dev.off()