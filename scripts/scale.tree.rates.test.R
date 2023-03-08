library(phytools)
library(geiger)
source("scale.tree.rates.R")
set.seed(5)
tree <- rcoal(100)
stretch.tree <- tree
stretch.tree$edge.length[162:198] <- stretch.tree$edge.length[162:198] *10
plot(stretch.tree)
mat <- matrix(c(-2,2,0,0,0,0,0,0,0,0,
                2,-4,2,0,0,0,0,0,0,0,
                0,2,-4,2,0,0,0,0,0,0,
                0,0,2,-4,2,0,0,0,0,0,
                0,0,0,2,-4,2,0,0,0,0,
                0,0,0,0,2,-4,2,0,0,0,
                0,0,0,0,0,2,-4,2,0,0,
                0,0,0,0,0,0,2,-4,2,0,
                0,0,0,0,0,0,0,2,-4,2,
                0,0,0,0,0,0,0,0,2,-2),10,10,byrow = T)
x <- sim.char(stretch.tree, par=mat,nsim=1 , model="discrete")[,,1]
plot(tree)

tree.scalar <- scaleTreeRates(tree = tree,
                              tip.states = x,
                              max.ratios = c(1.5,1.5),
                              model = mat,
                              nbins = c(10,10),
                              max.transition = 1,
                              use.eigen=T,
                              use.expm=T,
                              print.rates=F)
tree.scalar$scalar
colors <- c(viridis::viridis(n = sum(unique(tree.scalar$scalar) < 1),begin = 0,end = 0.2),
            "#000000",
            viridis::viridis(n = sum(unique(tree.scalar$scalar) > 1),begin = 0.8,end = 1))

names(colors) <- c(sort(unique(tree.scalar$scalar)))


plotScaledTree(tree = tree.scalar,
               colors = colors)

tree.scalar$scalar

fit <- ace2(x=x, phy = tree.scalar,type="discrete",model=mat,use.eigen = T,use.expm = F)


tree.scalar.stretch <- tree.scalar

tree.scalar.stretch$edge.length <- tree.scalar.stretch$edge.length * 
  tree.scalar.stretch$scalar

plot(tree.scalar.stretch)

x

eigen(mat)
