library(phytools)
library(geiger)
library(expm)
source("scale.tree.rates.R")
set.seed(35)
tree <- rcoal(20)
stretch.tree <- tree
stretch.tree$edge.length[20:38] <- stretch.tree$edge.length[20:38] * 10
plot(stretch.tree)
x <- sim.char(stretch.tree, par=matrix(c(-2,2,2,-2),2,2), nsim=1 , model="discrete")[,,1]
#x <- c(rep(1,10),rep(1:2,5))
#names(x) <- tree$tip.label
#shuffle
x <- sample(x)
plot(tree)

tree.scalar <- scaleTreeRates(tree = tree,
                              tip.states = x,
                              max.ratio = 2,
                              nbins = 10)

tree.scalar.stretch <- tree.scalar

tree.scalar.stretch$edge.length <- tree.scalar.stretch$edge.length * tree.scalar.stretch$scalar

plot(tree.scalar.stretch)

tree.scalar$scalar
colors <- c(viridis::viridis(n = sum(unique(tree.scalar$scalar) < 1),begin = 0,end = 0.2),
            "#000000",
            viridis::viridis(n = sum(unique(tree.scalar$scalar) > 1),begin = 0.8,end = 1))

names(colors) <- c(sort(unique(tree.scalar$scalar)))

plotScaledTree(tree = tree.scalar,
               colors = colors)

tree$edge
