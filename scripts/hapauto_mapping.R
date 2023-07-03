#### PACKAGES ####
library(phytools)
library(doSNOW)
library(viridis)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,3)]
tree <- read.tree("../data/trees/tree.nex")
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header = T))
Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_hapauto.csv",
                           as.is=T,header=T))

#### BUILD DATA MATRIX ####

#data.matrix
data.matrix <- matrix(0,nrow=nrow(dat),
                      ncol=max(dat$hapauto) - 1)

colnames(data.matrix) <- 1:49
rownames(data.matrix) <- dat$tree.name

#fill data matrix
for(i in 1:nrow(data.matrix)){
  data.matrix[i,as.numeric(dat$hapauto[i]) - 1] <- 1
}


#column and rownames for mat/Qmat
rownames(mat) <- 1:49
colnames(mat) <- 1:49
rownames(Qmat) <- 1:49
colnames(Qmat) <- 1:49

#### STOCHASTIC MAPPING ####
hists <- make.simmap2(tree = tree,
                     x = data.matrix,
                     model = mat,
                     nsim = 100,
                     Q = Qmat,
                     rejmax = 1000000,
                     rejint = 100000,
                     pi="fitzjohn",
                     monitor=T)

#### FIX STOCHASTIC MAPS ####
dat$sim.state <- dat$hapauto - 1
dat.for.fixing <- dat[,-2]
hists.fixed <- fix.simmap(hist,dat.for.fixing,mat)

#### SUMMARIZE STOCHASTIC MAPS ####
cols <- c(viridis(49))
names(cols) <- c(1:49)
plotSimmap(hist[[48]],col=cols,fsize = 0.05,lwd=1)
hists.summarized <- describe.simmap2(hists.fixed)

#### SAVE OUTPUTS ####
save(hists.fixed, file = "../outputs/hapauto_maps/hists.fixed.RData")
save(hists, file="../outputs/hapauto_maps/hists.RData")
save(hists.summarized, file = "../outputs/hapauto_maps/hists.summarized.RData")

