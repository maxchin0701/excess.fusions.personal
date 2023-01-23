#### PACKAGES ####

library(ape)
library(phytools)
library(igraph)
library(viridis)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]

clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

source("fix_simmap.R")

#### DATA ORGANIZATION ####

#Remove duplicate elements
data <- data[!duplicated(data$tree.name),]

#Convert species names to right format 
data$tree.name <- tolower(data$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
data$tree.name <- gsub(" ","_",data$tree.name)

#### LOOP THROUGH EACH SUBTREE ####

for(i in 1:5){
  
  #### LOAD IN TREE AND Q_MATRIX ####
  hists <- read.simmap(file=paste0("../data/subtree_simmap_out/subtree_",
                                   clades[i],
                                   "/simmap_out.nex"), 
                       format = "nexus",
                       version = 1.5)
  
  
  Qmat <- as.matrix(read.csv(paste0("../data/transition_matrix/subtree_matrices/matrix_",clades[i],".csv"),
                             as.is=T,header=T))
  
  #### BUILD TRANSITION MATRIX ####
  
  #Get rates
  asc <- Qmat[1,2]
  desc <- Qmat[2,1]
  SAF <- Qmat[2,(ncol(Qmat)/2 + 1)]
  Ro <- Qmat[(ncol(Qmat)/2 + 1),1]
  
  #Duplicate Qmat
  transition.matrix <- Qmat
  
  #Sub rates for parameter number
  transition.matrix[transition.matrix == asc] <- 1
  transition.matrix[transition.matrix == desc] <- 2
  transition.matrix[transition.matrix == SAF] <- 3
  transition.matrix[transition.matrix == Ro] <- 4
  
  #Set diagonals to 0
  diag(transition.matrix) <- 0
  
  #remove unnecesary objects
  rm(asc,desc,Ro,SAF,Qmat)
  
  
  #### TRIM DATASET ####
  
  #Trim phylogeny and chromosome data
  tip.dat <- data[data$tree.name %in% hists[[1]]$tip.label, ]
  
  #### FIX SIM.STATE ####
  
  #Establish range of haploid autosome numbers
  chrmrng <- range(tip.dat$hapauto)
  
  #Reassign new states
  for(j in 1:nrow(tip.dat)){
    if(is.na(tip.dat$hapauto[j]) == T ||
       tip.dat$codedSCS[j] == "unclear"){
      tip.dat$sim.state[j] <- NA
    } else if(tip.dat$codedSCS[j] == "unfused"){
      tip.dat$sim.state[j] <- tip.dat$hapauto[j] - chrmrng[1] + 1
    } else if(tip.dat$codedSCS[j] == "fused"){
      tip.dat$sim.state[j] <- (tip.dat$hapauto[j] - chrmrng[1] + 1) +
        (chrmrng[2] - chrmrng[1] + 1)
    }
  }
  
  #Trim down to necessary columns
  tip.dat <- tip.dat[,c(1,5)]
  
  #### RUN FIX.SIMMAP ####
  
  hists.fixed <- fix.simmap(hists = hists,
                            tips = tip.dat,
                            transition.matrix = transition.matrix)
  
  #### PLOT FIXED PHYLOS ####
  #Establish palette
  cols <- c(rainbow(ncol(transition.matrix)/2),
            viridis(ncol(transition.matrix)/2), 
            "black")
  names(cols) <- c(1:ncol(transition.matrix), "fail")
  
  #Plot each of 100 trees
  for(j in 1:100){
    
    cairo_pdf(paste0("../figures/subtree_simmap_fixed_plots/subtree_"
                     ,clades[i],
                     "/tree_",j,".pdf"))
    plotSimmap(hists.fixed[[j]],
               cols,
               fsize=sqrt(15/length(hists.fixed[[j]]$tip.label)),
               lwd = 1)
    dev.off()
    
  }
  
  #### SAVE FIXED SIMMAP ####
  
  write.simmap(hists.fixed,
               file=paste0("../data/subtree_simmap_out/subtree_",
                           clades[i],
                           "/simmap_fixed.nex"), 
               version = 1.5,
               format="nexus")
  
  
}


