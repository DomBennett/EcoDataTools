## Dominic Bennett
## 30/06/2013
## Comparing phylogenetic diversity between sites

## Functions
source('EcoDataTools.R')

## Set seed
set.seed(4)

## Dirs
dir <- file.path('cookbook', 'recipe1')
outdir <- file.path('wiki', 'cookbook')

## Load
load(file.path(dir, "test_metadata.RData"))

## rarefy and calculate PSE for each study
out <- list()
for (i in 1:length(data)) {
  temp.data <- data[[i]]
  temp.phylo <- phylos[[i]]
  incidence <- rowSums(ifelse(temp.data > 0, 1, 0))
  temp.data <- temp.data[incidence > 5, ]
  htypes <- rownames(temp.data)
  nhtypes <- length(unique(rownames(temp.data)))
  rownames(temp.data) <- NULL # phyloRarefy does not like non-unique names
  if (nhtypes > 1) {
    rarefied <- phyloRarefy(temp.data, temp.phylo, samp = 5, metric = 'PSE',
                            nrands = 1000)[,1]
    # TODO (30/06/2013): Does not work with PD?
    names(rarefied) <- htypes
    out <- c(out, list(rarefied))
  } else {
    next
  }
}

# plot
pse.values <- rle(unlist(out))$values
htypes <- factor(names(pse.values))
par(mar = c(5, 4, 4, 2), mfrow = c(1,1))
png("r1_pse_values.png", width = 720)
plot(pse.values ~ htypes, xlab = "Habitat Type", ylab = "PSE")
dev.off()

# report mean differences
