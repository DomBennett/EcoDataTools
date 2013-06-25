## Dominic Bennett
## 25/06/2013
## Generate test data and images for wiki

## Dirs + Libraries
source('functions.R')
setwd("wiki")

## Calculating Faith's PD
phylo <- read.tree("Phylocom_phylo.tre")
phylo <- drop.tip(phylo, phylo$tip.label[1:12])
taxa <- phylo$tip.label[c(20,16:13,5:6)]
png(filename = "PD.png", width = 1000, height = 800)
split.screen(c(2,2))
screen(1)
pd <- PD(phylo, taxa, type = 1, prop = FALSE, display = TRUE)
mtext("Type 1: phylogeny of the taxa", line = 2)
mtext(paste("PD:", pd, "   min: 2"))
screen(2)
pd <- PD(phylo, taxa, type = 2, prop = FALSE, display = TRUE)
mtext("Type 2: branches from the taxa tips to the terminal node", line = 2)
mtext(paste("PD:", pd, "   min: 1"))
screen(3)
pd <- PD(phylo, taxa, type = 3, prop = FALSE, display = TRUE)
mtext("Type 3: branches uniquely represented by the taxa", line = 2)
mtext(paste("PD:", pd, "   min: 1"))
close.screen(all.screens = TRUE)
dev.off()