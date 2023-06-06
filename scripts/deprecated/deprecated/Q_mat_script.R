# This script reads the transition matrix built in transition_matrix_script.R,
# along with the working mammalian karyotpye dataset with sim states and 
# mammalian tree. 

# Data is organized such that it is in the right format for the make.simmap 
# function of phytools. The sole purpose of running this simulation is to 
# obtain a Q matrix with reliable transition rate values, built by the function 
# prior to simulating trait history. The simulations themselves will not be 
# able to complete, thus the job can be cancelled once the Q matrix produced.

# All outputs of make.simmap will be outputted to a .csv file, which can then 
# be rearranged to build a Q matrix 

library(ape)
library(dplyr)
library(phytools)
library(devtools)
#

#Read tree and chromosome data
data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]
tree <- read.tree("../data/trees/4705sp_mean.nwk")

#Scale tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#Remove taxa with uncertain sim states
data <- subset(data,is.na(sim.state) == F)

#Convert species names to right format 
data$tree.name <- tolower(data$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
data$tree.name <- gsub(" ","_",data$tree.name)

#Trim phylogeny and chromosome data
data <- data[data$tree.name %in% tree$tip.label, ]
rownames(data) <- 1:nrow(data)
tree <- keep.tip(tree, data$tree.name)

#Establish range of haploid autosome numbers
chrmrng <- range(data$hapauto)

#Make data frame for sim.states
pruned.data <- as.data.frame(matrix(0, 
                                    length(tree$tip.label), 
                                    (chrmrng[2]-chrmrng[1] + 1) * 2))


row.names(pruned.data) <- tree$tip.label

#Fill pruned.data with sim.states
for(i in 1:nrow(pruned.data)){
  hit <- which(data$tree.name == row.names(pruned.data)[i])
  simcoding <- unique(data$sim.state[hit])
  weights <- 1/length(simcoding)
  pruned.data[i, simcoding] <- weights
}

#Remove unnecessary objects
rm(chrmrng,hit,i,simcoding,weights,data)

#Load model of chrom evolution (tmat)
tmat <- as.matrix(read.csv("../data/transition_matrix/matrix_mammalian.csv",
                           as.is=T,header = F))

#Make output file for simulation
sink(file = "../data/transition_matrix/output_Q.csv",split = T)

#Simulate (cancel once Q and pi have been produced)
hists <- make.simmap(tree,
                     x = as.matrix(pruned.data),
                     model = tmat,
                     pi = "fitzjohn",
                     Q = "mcmc",
                     nsim = 1)

#Save console output to file
sink()

?make.simmap
