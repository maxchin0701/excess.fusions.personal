#### PACKAGES ####

library(phytools)
library(doSNOW)
library(viridis)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,4)]
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_SAF.csv",
                          as.is=T,header = T))
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#### LOOP THROUGH CLADES ####
for(i in 1:5){
  
  # load in subtree and Qmatrix
  split.tree <- read.tree(paste0("../data/trees/subtrees/tree_",clades[i],".nex"))
  Qmat <- as.matrix(read.csv(paste0("../data/transition_matrix/subtree_matrices/multiSAF/matrix_",clades[i],".csv"),
                            as.is=T,header = T))
  
  #Subset data
  split.data <- subset(dat, tree.name %in% split.tree$tip.label)
  
  #### BUILD DATA MATRIX ####
  #data.matrix
  data.matrix <- matrix(0,nrow=nrow(split.data),
                        ncol=8)
  
  colnames(data.matrix) <- 1:8
  rownames(data.matrix) <- split.data$tree.name
  
  #fill data matrix
  for(j in 1:nrow(data.matrix)){
    data.matrix[j,as.numeric(split.data$codedSCS[j]) + 1] <- 1
  }
  
  #change matri names
  rownames(mat) <- colnames(mat) <- rownames(Qmat) <- colnames(Qmat) <-1:8
  
  #### PERFORM STOCHASTIC MAPPING ####
  
  hists <- make.simmap2(tree = split.tree,
                        x = data.matrix,
                        model = mat,
                        nsim = 100,
                        Q = Qmat,
                        rejmax = 1000000,
                        rejint = 100000,
                        pi=c(1,0,0,0,0,0,0,0),
                        monitor=T)
  
  #### SUMMARIZE STOCHASTIC MAPS ####
  hists.summarized <- describe.simmap2(hists)
  
  #### SAVE OUTPUTS ####
  save(hists.summarized, file = paste0("../outputs/SAF_maps/subtrees/hists.",clades[i],".summarized.RData"))
  save(hists, file=paste0("../outputs/SAF_maps/subtrees/hists.",clades[i],".RData"))
}



