library(phytools)
library(ape)
library(expm)

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]

tree <- read.tree("../data/trees/4705sp_mean.nwk")

mat <- as.matrix(read.csv("../data/transition_matrix/matrix_mammalian.csv",
                           as.is=T,header=F))

Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_mammalian.csv",
                           as.is=T,header=T))


source("scale.tree.rates.R")

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


#Identify and remove multiple sim states for taxa (for now, troubleshoot in future)
data.dup <- data[duplicated(data$tree.name),]

#Remove duplicate elements
data.unique <- data[!duplicated(data$tree.name),]

#Use temp df to reorder tip data to match tree order
data.alt <- data.unique

for(j in 1:length(tree$tip.label)){
  index <- which(data.unique$tree.name == tree$tip.label[j])
  data.alt[j,] <- data.unique[index,]
}

data <- data.alt

rm(data.alt,data.unique,data.dup)

#### INITIAL FITMK #####
fit <- fitMkNew(tree,data$sim.state,model=mat,fixedQ = Qmat)

logLik.start <- fit$logLik



#### RUN SCALING ANALYSIS ####

tip.states <- data$sim.state
names(tip.states) <- data$tree.name

scaled.tree <- scaleTreeRates(tree = tree,
                              tip.states = tip.states,
                              max.ratio = 1.5,
                              nbins=10,
                              model = mat,
                              fixedQ = Qmat)


nb.tip <- length(tree$tip.label)
nb.node <- tree$Nnode

rate <- mat

nlminb(rep(0.1, length.out = 4), function(p) dev(p), 
       lower = rep(0, 4), upper = rep(1e+50, 4)) -> out

e1 <- tree$edge[, 1]
e2 <- tree$edge[, 2]

c(p, 0)[rate]

x <-data$sim.state
