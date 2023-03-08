library(phytools)
library(doSNOW)
library(viridis)

#### SIMMAP2 FUNCTIONS ####

source("make.Simmap2.R")

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]

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

rm(data.unique)


#### LOOP THROUGH EACH SUBTREE ####

for(i in 1:5){
  
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
  
  #Make data frame for sim.states
  pruned.data <- as.data.frame(matrix(0,
                                      length(tree$tip.label),
                                      (chrmrng[2]-chrmrng[1] + 1) * 2))
  
  row.names(pruned.data) <- tree$tip.label
  
  #Fill pruned.data with sim.states
  for(j in 1:nrow(pruned.data)){
    hit <- which(data.subset$tree.name == row.names(pruned.data)[j])
    simcoding <- unique(data.subset$sim.state[hit])
    weights <- 1/length(simcoding)
    pruned.data[j, simcoding] <- weights
  }
  
  #Remove unnecessary objects
  rm(chrmrng,hit,j,simcoding,weights)
  
  #Set colnames
  colnames(pruned.data) <- 1:ncol(pruned.data)
  colnames(mat) <- 1:ncol(pruned.data)
  colnames(Qmat) <- 1:ncol(pruned.data)
  
  
  #### SET UP PARALLELIZATION ####
  
  # Define number of clusters
  nClust <- 50
  
  # Set up clusters, print intermediates to 
  cl <- makeCluster(nClust, outfile = "")
  registerDoSNOW(cl)
  
  
  
  #### RUN SIMMAP ####
  sim.out <- foreach(j=1:100,
                     .verbose = T,
                     .packages = c("phytools","maps","ape")) %dopar% {
                       
                       result <- make.simmap2(tree = tree,
                                              x = as.matrix(pruned.data),
                                              model = mat,
                                              nsim = 1,
                                              rejmax = 5000000,
                                              rejint = 100000,
                                              Q = Qmat,
                                              pi="fitzjohn",
                                              monitor=F,
                                              parallel=c(j,100))
                     }
  
  #Close conection to cluster
  stopCluster(cl)
  rm(cl)
  
  #Set class of sim.out
  class(sim.out) <- c("multiSimmap","multiPhylo")

  #### PLOT MAPPED PHYLOGENYS AND SAVE ####
  
  #Establish palette
  cols <- c(rainbow(ncol(pruned.data)/2),viridis(ncol(pruned.data)/2), "black")
  names(cols) <- c(1:ncol(pruned.data), "fail")
  
  #Plot each of 100 trees
  for(j in 1:100){
    
    cairo_pdf(paste0("../figures/subtree_simmap_plots/subtree_",clades[i],"/tree_",j,".pdf"))
    plotSimmap(sim.out[[j]],
               cols,
               fsize=sqrt(15/length(sim.out[[j]]$tip.label)),
               lwd = 1)
    dev.off()
    
  }
  
  #### Save IMPORTANT OUTPUTS ####
  
  #Save count of transitions to object
  count <- describe.simmap(sim.out)$count
  
  colnames(count) <- gsub(",","_",colnames(count))
  
  #Save dwelling times for each tree
  for(j in 1:100){
    
    #Save times to object
    times <- describe.simmap(sim.out[[j]])$times
    
    #Save object
    write.csv(times, paste0("../data/subtree_simmap_out/subtree_",clades[i],"/dwelling_times/dwelling_times_",j,".csv"),quote=F,row.names=T)
  }
  
  #Summarize failed edge data
  failed.edges <- c()
  
  for(j in 1:100){
    sim.fail <- sim.out[[j]]$fail
    names(sim.fail) <- rep(j,length(sim.fail))
    failed.edges <- c(failed.edges,sim.fail)
    rm(sim.fail)
  }
  
  failed.edges <- factor(failed.edges)
  
  fail.df <- data.frame(summary(failed.edges))
  
  fail.df$edge <- rownames(fail.df)
  fail.df$n <- fail.df[,1]
  rownames(fail.df) <- NULL
  fail.df <- fail.df[,-1]
  
  #Save transitions
  write.csv(count,paste0("../data/subtree_simmap_out/subtree_",clades[i],"/transitions.csv"),quote=F,row.names=T)
  
  #Save failed tip data
  write.csv(fail.df,"../data/subtree_simmap_out/subtree_",clades[i],"/failed_edges.csv",quote=F,row.names=T)
  
  #Save simmap object
  write.simmap(sim.out,file="../data/subtree_simmap_out/subtree_",clades[i],"/simmap_out.nex",version = 1.5,format="nexus")
}



