#### PACKAGES ####

library(phytools)
library(doSNOW)
library(viridis)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,3)]
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header = T))
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#### SET UP PARALLELIZATION ####
nClust <- 5

# Set up clusters, print intermediates to 
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

#### LOOP THROUGH CLADES ####
foreach(i=1:5,
        .verbose = T,
        .packages = c("phytools","maps","ape")) %dopar% {
  
  # load in subtree and Qmatrix
  split.tree <- read.tree(paste0("../data/trees/subtrees/tree_",clades[i],".nex"))
  Qmat <- as.matrix(read.csv(paste0("../data/transition_matrix/subtree_matrices/hapauto/matrix_",clades[i],".csv"),
                            as.is=T,header = T))
  
  #Subset data
  split.data <- subset(dat, tree.name %in% split.tree$tip.label)
  
  #### BUILD DATA MATRIX ####
  #data.matrix
  data.matrix <- matrix(0,nrow=nrow(split.data),
                        ncol=ncol(mat))
  
  colnames(data.matrix) <- 1:49
  rownames(data.matrix) <- split.data$tree.name
  
  #fill data matrix
  for(j in 1:nrow(data.matrix)){
    data.matrix[j,as.numeric(split.data$hapauto[j]) - 1] <- 1
  }
  
  #change matri names
  rownames(mat) <- colnames(mat) <- rownames(Qmat) <- colnames(Qmat) <-1:49
  
  #### PERFORM STOCHASTIC MAPPING ####
  
  hists <- make.simmap2(tree = split.tree,
                        x = data.matrix,
                        model = mat,
                        nsim = 100,
                        Q = Qmat,
                        rejmax = 1000000,
                        rejint = 100000,
                        pi="fitzjohn",
                        monitor=T)
  
  #### FIX STOCHASTIC MAPS ####
  split.data$sim.state <- split.data$hapauto - 1
  dat.for.fixing <- split.data[,-2]
  hists.fixed <- fix.simmap(hist,dat.for.fixing,mat)
  
  #### SUMMARIZE STOCHASTIC MAPS ####
  hists.summarized <- describe.simmap2(hists.fixed)
  
  #### SAVE OUTPUTS ####
  save(hists.summarized, file = paste0("../outputs/hapauto_maps/subtrees/hists.",clades[i],".summarized.RData"))
  save(hists, file=paste0("../outputs/hapauto_maps/subtrees/hists.",clades[i],".RData"))
  save(hists.fixed, file=paste0("../outputs/hapauto_maps/subtrees/hists.",clades[i],".fixed.RData"))
}

stopCluster(cl)
rm(cl)




s