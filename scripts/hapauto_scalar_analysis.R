#### PACKAGES ####
library(phytools)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)
tree <- read.tree("../data/trees/tree.nex")
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header=F))[-1,]
Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_hapauto.csv",
                           as.is=T,header=T))

#### PREPARE DATA ####
colnames(mat) <- 1:49
colnames(Qmat) <- 1:49
tip.states <- dat$hapauto - 1
names(tip.states) <- dat$tree.name

#### RUN SCALING ANALYSIS ####
scaled.tree <- scaleTreeRates(tree = tree,
                              tip.states = tip.states,
                              max.ratio = 2,
                              nbins=10,
                              model = mat,
                              fixedQ = Qmat,
                              pi="fitzjohn")

####P PLOT SCALED TREE AND DETERMINE CUTOFF ####
plot.phyloscaled(scaled.tree,cex=0.05,edge.width = 0.1)

plot(tree,cex=0.05,edge.width = 0.1)
edgelabels(cex=0.1,frame="none")

cut.tree <- extract.clade(tree,node=tree$edge[1120,2])
plot(cut.tree,cex=0.05,edge.width = 0.1)

#### SAVE SCALED TREE AND CUT TREE ####
write.nexus(cut.tree,file="../data/trees/cut.tree.nex")
save(scaled.tree,file="../outputs/scaled.tree.RData")

