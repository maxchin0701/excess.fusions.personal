library(phytools)
library(geiger)
library(expm)
library(RColorBrewer)
library(viridis)
source("scale.tree.rates.R")
set.seed(35)
tree <- rcoal(20)
stretch.tree <- tree
stretch.tree$edge.length[20:38] <- stretch.tree$edge.length[20:38] * 10
plot(stretch.tree)
x <- sim.char(stretch.tree, par=matrix(c(-2,2,2,-2),2,2), nsim=1 , model="discrete")[,,1]
plot(tree)

tree.scalar <- scaleTreeRates(tree = tree,
                                tip.states = x,
                                max.ratio = 1.5,
                                nbins = 30)



unique(tree.scalar$scalar)

plot.phyloscaled(tree.scalar)

# plot(tree.scalar, edge.color="red")
# 
# 
# 
# tree.scalar.stretch <- tree.scalar
# 
# tree.scalar.stretch$edge.length <- tree.scalar.stretch$edge.length * tree.scalar.stretch$scalar
# 
# plot(tree.scalar.stretch)
# 
# tree.scalar$scalar
# 
# slowclasses <- sum(unique(tree.scalar$scalar) < 1)
# fastclasses <- sum(unique(tree.scalar$scalar) > 1)
# brewer.pal(sum(c(slowclasses,1, fastclasses)), name="RdYlGn")
# 
# 
# cols <- as.numeric(as.factor(tree.scalar$scalar))
# shade <- brewer.pal(sum(c(slowclasses,1, fastclasses)), name="RdYlGn")
# plot.phylo(tree.scalar, edge.color = colors,edge.width = 1.75)
# edgelabels(tree.scalar$scalar)
# colors <- c(viridis(n = sum(unique(tree.scalar$scalar) < 1),begin = 1,end = 0.2, option="G"),
#             "#000000",
#             viridis(n = sum(unique(tree.scalar$scalar) > 1),begin = 1,end = 0.2,option = "F"))
# 
# names(colors) <- c(sort(unique(tree.scalar$scalar)))
# 
# plot.phyloscaled(tree = tree.scalar,
#                colors = colors)
# 
# tree$edge
