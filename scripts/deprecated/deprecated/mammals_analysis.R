# This script performs the empirical and theoretical analysis of sex-autosome 
# fusion in mammalian taxa

# Inputs are a working mammalian karyotype dataset with associated simulation
# states, a mammalian maximum likelihood phylogenetic tree, and a filled in 
# Q-matrix. 

# Simulations are performed using the sfreemap function of the sfreemap package

library(ape)
library(dplyr)
library(sfreemap)

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

#Identify and remove multiple sim states for taxa (for now, troubleshoot in future)
data.dup <- data[duplicated(data$tree.name),]

#Remove duplicate elements
data.unique <- data[!duplicated(data$tree.name),]


# #Convert to charcter
# data.unique$sim.state <- as.character((data.unique$sim.state))

#Convert data to vector form
data.vector <- c()

for(i in 1:nrow(data.unique)){
  data.vector <- c(data.vector,data.unique$sim.state[i])
}

#Set names
names(data.vector) <- data.unique$tree.name

#Reorder to match tip labels
data.vector <- data.vector[order(factor(names(data.vector),
                                        levels = tree$tip.label))]

#Create temp vector
tips <- c()

#Loop through and bind elements
for(i in 1:length(data.vector)){
  #Check for fused vs unfused
  if(data.vector[i] > 49){
    chrom <- paste0(as.numeric(data.vector[i] - 48),"n")
  } else {
    chrom <- paste0(data.vector[i]+1)
  }
  
  #Add species name
  tip <- paste0(names(data.vector[i])," (",chrom,")")
  
  #Bind to temporary vector
  tips <- c(tips,tip)
}

#Remove redundant variables
rm(tip,chrom)

# # #Establish range of haploid autosome numbers
# # chrmrng <- range(data$hapauto)
# # 
# # #Make data frame for sim.states
# # pruned.data <- as.data.frame(matrix(0, 
# #                                     length(tree$tip.label), 
# #                                     (chrmrng[2]-chrmrng[1] + 1) * 2))
# # 
# # 
# # row.names(pruned.data) <- tree$tip.label
# # 
# # #Fill pruned.data with sim.states
# # for(i in 1:nrow(pruned.data)){
# #   hit <- which(data$tree.name == row.names(pruned.data)[i])
# #   simcoding <- unique(data$sim.state[hit])
# #   weights <- 1/length(simcoding)
# #   pruned.data[i, simcoding] <- weights
# # }
# # 
# # #Remove unnecessary objects
# # rm(chrmrng,hit,i,simcoding,weights,data)
# 
# # 
# # tree.names <- data.frame(tree$tip.label)
# # 
# # #Keep only distinct rows
# # data <- data %>%
# #             distinct()
# # 
# # #Identify duplicate elements
# # data.dup <- data[duplicated(data$tree.name),]
# # 
# # #Remove duplicate elements
# # data.tree <- data[!duplicated(data$tree.name),]
# # 
# # #Bind species with duplicate elements to duplicated df
# # for(i in 1:nrow(data.tree)){
# #   if(data.tree$tree.name[i] %in% data.dup$tree.name){
# #     data.dup <- rbind(data.dup,
# #                       data.tree[i,])
# #   }
# # }
# # 
# # data.dup$tree.name <- as.factor(data.dup$tree.name)
# # 
# # #Use temp df to reorder tip data to match tree order
# # data.alt <- data.tree
# # 
# # for(j in 1:length(tree$tip.label)){
# #   index <- which(data.tree$tree.name == tree$tip.label[j])
# #   data.alt[j,] <- data.tree[index,]
# # }
# # 
# # data.tree <- data.alt
# # 
# # rm(data.alt)
# # 

#Plot tree

tree.alt <- tree

tree.alt$tip.label <- tips

plot(tree.alt, show.tip.label = T,
     #align.tip.label = T,
     cex = 0.02,
     no.margin = T,
     #label.offset = 0.01,
     edge.width = 0.00000000000001)

dist <- as.data.frame(cophenetic.phylo(tree.alt))

# tiplabels(tips, 
#           cex = 0.02,
#           align.tip.label=T,
#           #offset = 0.5,
#           frame = "none",)
# 
# ?plot.phylo
# #
# # #Make sim state matrix
# # sim.state <- matrix(0,nrow(data.tree),98)
# # rownames(sim.state) <- data.tree$tree.name
# # colnames(sim.state) <- 1:98

# # #Fill sim state matrix
# # for(k in 1:nrow(sim.state)){
# #   
# #   #Check if species has multiple states
# #   if(data.tree$tree.name[k] %in% data.dup$tree.name == F){
# #     
# #     #Fill with one if only single occurence
# #     sim.state[k,data.tree$sim.state[k]] <- 1
# #     
# #   } else {
# #     #Determine prob of each state
# #     prob <- 1/length(which(data.dup$tree.name==data.tree$tree.name[k]))
# #     
# #     #Subset only rows with right species name
# #     data.dup.subset <- subset(data.dup, tree.name == data.tree$tree.name[k])
# #     
# #     #Fill each sim.state with prob
# #     for(l in 1:nrow(data.dup.subset)){
# #       sim.state[k,data.dup.subset$sim.state[l]] <- prob
# #     }
# #   }
# # }

#Read in Q matrix
qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_mammalian.csv",
                           as.is=T,header = F))

#Set column names
colnames(qmat) <- 1:98

#Simulate
hist <- sfreemap(tree = tree,
                 tip_states = data.vector,
                 Q = as.matrix(qmat),
                 pi = "estimated")

summary(data.unique$sim.state)

tip_states <- build_states_matrix(tree$tip.label, data.test.vector, NULL)




