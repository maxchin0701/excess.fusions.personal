#### PACKAGES ####

library(phytools)
library(doSNOW)
library(viridis)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,4)]
tree <- read.tree("../data/trees/tree.nex")
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_SAF.csv",
                          as.is=T,header = T))
Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_SAF.csv",
                           as.is=T,header=T))

#### BUILD DATA MATRIX ####

#data.matrix
data.matrix <- matrix(0,nrow=nrow(dat),
                      ncol=max(dat$codedSCS) + 1)

colnames(data.matrix) <- 1:8
rownames(data.matrix) <- dat$tree.name

#fill data matrix
for(i in 1:nrow(data.matrix)){
  data.matrix[i,as.numeric(dat$codedSCS[i]) + 1] <- 1
}

#column and rownames for mat/Qmat
rownames(mat) <- 1:8
colnames(mat) <- 1:8
rownames(Qmat) <- 1:8
colnames(Qmat) <- 1:8

#### PERFORM STOCHASTIC MAPPING ####

hists <- make.simmap2(tree = tree,
                      x = data.matrix,
                      model = mat,
                      nsim = 100,
                      Q = Qmat,
                      rejmax = 1000000,
                      rejint = 100000,
                      pi=c(1,0,0,0,0,0,0,0),
                      monitor=T)

#### FIX SIMMAP ####

#loop through hists
for(i in 1:length(hists)){
  length <- 1/8
  states <- rep(length,8)
  names(states) <- c("1","2","3","4","5","6","7","8")
  hists[[i]]$maps[[1]] <- states
  hists[[i]]$mapped.edge[1,] <- rep(length,8)
} 

#### EXTRACT PROPORTION TIME IN NON ZERO SAF ####
hists.summarized <- describe.simmap2(hists)
prop.XY <- hists.summarized$times[,1]/hists.summarized$times[,9]
overalprop.xy <- sum(hists.summarized$times[,1])/sum(hists.summarized$times[,9])

#### SAVE OUTPUTS ####
save(hists.summarized, file = "../outputs/SAF_maps/hists.summarized.RData")
save(hists, file="../outputs/SAF_maps/hists.RData")



















