library(phytools)
library(geiger)
library(expm)
library(RColorBrewer)
library(viridis)
library(ape)
source("funcitons.R")
# set.seed(35)
# tree <- rcoal(50)
# stretch.tree <- tree
# stretch.tree$edge.length[60:98] <- stretch.tree$edge.length[60:98] * 10
# plot(stretch.tree)
# x <- sim.char(stretch.tree, par=matrix(c(-2,2,2,-2),2,2), nsim=1 , model="discrete")[,,1]
set.seed(30)
tree <- rcoal(100)
plot(tree)
stretch.tree <- tree
stretch.tree$edge.length[108:198] <- stretch.tree$edge.length[108:198] * 10
plot(stretch.tree)
set.seed(1)
x <- sim.char(stretch.tree, par=Qmat, nsim=1 , model="discrete")[,,1]
x
names(x) <- tree$tip.label

plot(tree)
Qmat <- mat
diag(Qmat) <- -rowSums(Qmat)
mat <- matrix(c(rep(c(0,1,rep(0,8),2),9),0),nrow = 10,ncol = 10,byrow = T)
tree.scalar <- scaleTreeRates(tree = tree,
                                tip.states = x,
                                max.ratio = 1.5,
                                nbins = 10,
                              model=mat,
                              pi = "fitzjohn")

tree.scalar2 <- scaleTreeRates2(tree = tree,
                              tip.states = x,
                              max.ratio = 1.5,
                              nbins = 10)
unique(tree.scalar$scalar)

plot.phyloscaled(tree.scalar,palette="viridis",edge.width = 1.75,cex=1,show.tip.label=F)

plot.phyloscaled(tree.scalar2,palette = "viridis",edge.width = 1.75,cex=1)

plot(tree.scalar, edge.color="red")

plot(jitter(tree.scalar$scalar),jitter(tree.scalar2$scalar))

tree.scalar.stretch <- tree.scalar

tree.scalar.stretch$edge.length <- tree.scalar.stretch$edge.length * tree.scalar.stretch$scalar

plot(stretch.tree,show.tip.label = F)

plot(tree.scalar.stretch,show.tip.label = F)

tree$tip.label <- x
plot(tree)

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


yay <- 0
youre_fucked <- 0

for(i in 1:1000){
  print(paste0("ITERATION ",i))
  
  x <- sim.char(stretch.tree, par=Qmat, nsim=100, model="discrete",root=1)[,,1]
  x
  tree.scalar <- scaleTreeRates(tree = tree,
                                tip.states = x,
                                max.ratio = 1.5,
                                nbins = 10,
                                model=mat,
                                pi=c(1,rep(0,9)))
  
  
  if(mean(tree.scalar$scalar[1:107]) < mean(tree.scalar$scalar[108:198])){
    
    yay <- yay + 1
    
  } else {
    
    youre_fucked <- youre_fucked + 1
  }
  
  print(yay/(yay + youre_fucked))
}


tree.test
