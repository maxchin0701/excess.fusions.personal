library(phytools)
library(doSNOW)
library(viridis)

#### SIMMAP2 FUNCTIONS ####

source("make.Simmap2.R")

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]

tree <- read.tree("../data/trees/4705sp_mean.nwk")

mat <- as.matrix(read.csv("../data/transition_matrix/matrix_mammalian.csv",
                          as.is=T,header = F))

Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_mammalian.csv",
                           as.is=T,header=T))

#### DATA ORGANIZATION ####

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

#Scale tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))

# #Identify and remove multiple sim states for taxa (for now, troubleshoot in future)
# data.dup <- data[duplicated(data$tree.name),]
# 
# #Remove duplicate elements
# data.unique <- data[!duplicated(data$tree.name),]
# 
# #Use temp df to reorder tip data to match tree order
# data.alt <- data.unique
# 
# for(j in 1:length(tree$tip.label)){
#   index <- which(data.unique$tree.name == tree$tip.label[j])
#   data.alt[j,] <- data.unique[index,]
# }
# 
# data <- data.unique <- data.alt
# 
# rm(data.alt,data.unique)

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

#Set colnames
colnames(pruned.data) <- 1:98
colnames(mat) <- 1:98
colnames(Qmat) <- 1:98

#### SET UP PARALLELIZATION ####

# Define number of clusters
nClust <- 100

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
                              monitor=T,
                              parallel=c(j,100))
                   }

#Close cluster connection
stopCluster(cl)

#Set class of sim.out
class(sim.out) <- c("multiSimmap","multiPhylo")

#### PLOT MAPPED PHYLOGENYS AND SAVE ####

#Establish palette
cols <- c(rainbow(49),viridis(49), "black")
names(cols) <- c(1:98, "fail")

#Plot each of 50 trees
for(i in 1:100){
  
  cairo_pdf(paste0("../figures/simmap_plots/tree_",i,".pdf"))
  plotSimmap(sim.out[[i]],col=cols,fsize=0.03,lwd = 0.1)
  dev.off()
  
}

#### Save IMPORTANT OUTPUTS ####

#Save count of transitions to object
count <- describe.simmap(sim.out)$count

colnames(count) <- gsub(",","_",colnames(count))

#Save dwelling times for each tree
for(i in 1:100){
  
  #Save times to object
  times <- describe.simmap(sim.out[[i]])$times
  
  #Save object
  write.csv(times, paste0("../data/simmap_out/dwelling_times/dwelling_times_",i,".csv"),quote=F,row.names=T)
}

#Summarize failed edge data
failed.edges <- c()

for(i in 1:100){
  sim.fail <- sim.out[[i]]$fail
  names(sim.fail) <- rep(i,length(sim.fail))
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
write.csv(count,"../data/simmap_out/transitions.csv",quote=F,row.names=T)

#Save failed tip data
write.csv(fail.df,"../data/simmap_out/failed_edges.csv",quote=F,row.names=T)

#Save simmap object
write.simmap(sim.out,file="../data/simmap_out/simmap_out.nex",version = 1.5,format="nexus")






















