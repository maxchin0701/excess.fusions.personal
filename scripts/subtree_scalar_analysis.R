#### PACKAGES ####

library(phytools)
library(geiger)
library(expm)
library(RColorBrewer)
library(viridis)
library(ape)
source("functions.R")

#### DATA ####
data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(3,5,9,11,12,15)]

clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#### DATA ORGANIZATION ####

#Remove taxa with uncertain sim states
data <- subset(data,is.na(sim.state) == F)

#Convert species names to right format 
data$tree.name <- tolower(data$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
data$tree.name <- gsub(" ","_",data$tree.name)

#Identify and remove multiple sim states for taxa (for now, troubleshoot in future)
data.dup <- data[duplicated(data$tree.name),]

#Remove duplicate elements
data.unique <- data[!duplicated(data$tree.name),]

data <- data.unique

rm(data.unique,data.dup)

#### LOOP THROUGH EACH SUBTREE ####

for(i in 1:5){
  
  print(paste0("subtree ",clades[i]))
  
  #### READ TREE AND QMAT ####
  tree <- read.tree(paste0("../data/trees/subtrees/tree_",clades[i],".nwk"))
  
  
  Qmat <- as.matrix(read.csv(paste0("../data/transition_matrix/subtree_matrices/matrix_",clades[i],".csv"),
                             as.is=T,header=T))
  
  #### BUILD TRANSITION MATRIX ####
  
  # #Set colnames to NULL
  # colnames(Qmat) <- NULL
  
  #Get rates
  asc <- Qmat[1,2]
  desc <- Qmat[2,1]
  SAF <- Qmat[2,(ncol(Qmat)/2 + 1)]
  Ro <- Qmat[(ncol(Qmat)/2 + 1),1]
  
  #Duplicate Qmat
  mat <- Qmat
  
  #Sub rates for parameter number
  mat[mat == asc] <- 1
  mat[mat == desc] <- 2
  mat[mat == SAF] <- 3
  mat[mat == Ro] <- 4
  
  #Set diagonals
  diag(mat) <- 0
  
  #Set colnames
  colnames(mat) <- 1:ncol(mat)
  colnames(Qmat) <- 1:ncol(Qmat)
  
  rm(Ro,desc,SAF,asc)
  
  #### TRIM DATASET AND BUILD STATES MATRIX ####
  
  #Trim phylogeny and chromosome data
  data.subset <- data[data$tree.name %in% tree$tip.label, ]
  
  #Establish range of haploid autosome numbers
  chrmrng <- range(data.subset$hapauto)
  
  #Determine sim.state
  for(j in 1:nrow(data.subset)){
    if(is.na(data.subset$hapauto[j]) == T ||
       data.subset$codedSCS[j] == "unclear"){
      data.subset$sim.state[j] <- NA
    } else if(data.subset$codedSCS[j] == "unfused"){
      data.subset$sim.state[j] <- data.subset$hapauto[j] - chrmrng[1] + 1
    } else if(data.subset$codedSCS[j] == "fused"){
      data.subset$sim.state[j] <- (data.subset$hapauto[j] - chrmrng[1] + 1) +
        (chrmrng[2] - chrmrng[1] + 1)
    }
  }
  
  rm(chrmrng)
  
  #### RUN SCALING ANALYSIS ####
  
  tip.states <- data.subset$sim.state
  names(tip.states) <- data.subset$tree.name
  
  scaled.tree <- scaleTreeRates(tree = tree,
                                tip.states = tip.states,
                                model = mat,
                                fixedQ = Qmat)
  
  #Get tree to right name
  tree.name <- paste0("scaled.tree.",clades[i])
  
  assign(tree.name,scaled.tree)
  
}

plot.phyloscaled(scaled.tree.carnivora,cex=0.1)
plot.phyloscaled(scaled.tree.artiodactyla,cex=0.1)
plot.phyloscaled(scaled.tree.yangochiroptera,cex=0.1)
plot.phyloscaled(scaled.tree.anomaluromorpha_castorimorpha_myomorpha,cex = 0.1)
plot.phyloscaled(scaled.tree.primatomorpha,cex = 0.1)


hists <- read.simmap(file="../data/simmap_out/simmap_out.nex", 
            format = "nexus",
            version = 1.5)

hists -> tree


