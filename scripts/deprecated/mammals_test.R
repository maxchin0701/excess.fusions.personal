# make a small 10 tip tree (rcoal)
# make up some chrom numbers and fused unfused states
# get these in the same format including tip uncertainty (two dif chrom numb)

# This script creates a test dataset and tree, reading in a transition 
# matrix which describes the possible transitions. All analyses are tested
# here before being run on the working chromosomes dataset
library(ape)
library(chromePlus)
library(phytools)
library(sfreemap)

#Create test tree
tree.test <- rcoal(10,rooted=T)

#Scale tree to unit length
tree.test$edge.length <- tree.test$edge.length/max(branching.times(tree.test))

q <- matrix(0,ncol = 6,nrow = 6)

#Generate chromosome numbers
data.test.vector <- simChrom(tree.test, 
                             pars = c(1,1,1,1,0,0,0,0,1,1,25,1),
                             model="SAF",
                             limits = c(5,50),
                             verbose = T)

data.test.vector

#Build data
data.test <- data.frame(names(data.test.vector[[1]]),
                        data.test.vector[2],
                        data.test.vector[1])

#Remove redundant elements
rm(data.test.vector)

#Set colnames
colnames(data.test) <- c("species",
                         "chrom.num",
                         "state")

#Set col
data.matrix <- datatoMatrix(data.test,
                            c(2,50),
                            hyper = T)

#Make mkn model
model <- make.mkn(tree.test,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#Check argnames: 9506 different rates (98*97, expected)
length(argnames(model))

#Constrain model
model.con <- constrainMkn(data.matrix,
                          model,
                          hyper = T,
                          polyploidy = F,
                          constrain = list(drop.poly=T,
                                           drop.demi=T,
                                           nometa=T))

#Check argnames: 4 rates, model successfully constrained
argnames(model.con)

#Run mcmc
model.mcmc <- mcmc(model.con,
                   c(1,1,1,1),
                   prior=make.prior.exponential(r=0.5),
                   upper=c(100,100,100,100),
                   nsteps = 100,
                   1)

#Check for convergence (simple model, converges extremely fast)
plot(model.mcmc$i,model.mcmc$p,type="l")

#Plot to see distribution 
profiles.plot(model.mcmc["asc1"], col.line="red")
profiles.plot(model.mcmc["desc1"], col.line="red")


#Extract post burn portion
model.mcmc.postburn <- model.mcmc[25:100,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1),
            mean(model.mcmc.postburn$tran12),
            mean(model.mcmc.postburn$tran21))

names(params) <- data.test <- colnames(model.mcmc[,2:5])

#Convert chrom numbs + binary states to sim state
# sim.state.test <- c()
#
# for( i in 1:nrow(data.test)){
#   if(data.test$binary.state[i] == "0"){
#     sim.state.test <- c(sim.state.test,data.test$chrom.num[i])
#   } else {
#     sim.state.test <- c(sim.state.test,(data.test$chrom.num[i] + 5))
#   }
# }

#Name vector correctly
# names(sim.state.test) <- rownames(data.test)
# 
# #Create test dataset
# data.test <- data.frame(names(sim.state.test),
#                         sim.state.test)
# 
# #Assign column names
# colnames(data.test) <- c("tree.name","sim.state")

#Convert to charcter
# data.test$sim.state <- as.character((data.test$sim.state))

# #Convert data to vector form
# data.test.vector <- c()
# 
# for(i in 1:nrow(data.test)){
#   data.test.vector <- c(data.test.vector,data.test$sim.state[i])
# }
# 
# #Set names
# names(data.test.vector) <- data.test$tree.name
# 
# #Reorder to match tip labels
# data.test.vector <- data.test.vector[order(factor(names(data.test.vector),
#                                         levels = tree.test$tip.label))]
# 
# 
# #Create sim.state matrix
# pruned.data.test <- as.data.frame(matrix(0, 
#                                          length(tree.test$tip.label), 
#                                          10))
# 
# #Assign row names
# row.names(pruned.data.test) <- tree.test$tip.label
# 
# #Fill matrix
# for(i in 1:nrow(pruned.data.test)){
#   hit.test <- which(data.test$tree.name == row.names(pruned.data.test)[i])
#   simcoding.test <- unique(data.test$sim.state[hit.test])
#   weights.test <- 1/length(simcoding.test)
#   pruned.data.test[i, simcoding.test] <- weights.test
# }
# 
# rm(data.test,hit.test,simcoding.test,weights.test)
# 
# #Read matrix
# Tmat.test <- as.matrix(read.csv("../data/transition_matrix/transition.matrix.test.csv",
#                                 header=F))
# 
# #Convert NAs to 0
# Tmat.test[is.na(Tmat.test)] <- 0
# 
# #Set column names
# colnames(Tmat.test) <- 1:10
# 
# # ?fitMk
# # 
# model <- fitMk(tree = tree.test,
#                x = data.test.vector,
#                model = Tmat.test)
# 
# Qmat.alt <- as.Qmatrix(model)
# 
# #Simulate
# hists.test <- make.simmap(tree.test,
#                           x = as.matrix(pruned.data.test),
#                           model = Tmat.test,
#                           pi = c(0,0,0,0,1,0,0,0,0,0),
#                           nsim = 1)
# 
# #Save built Q matrix 
# Qmat.test <- hists.test$Q
# 
# plotSimmap(hists.test)
# 
# #Remove unnecessary objects
# rm(Tmat.test,pruned.data.test,hists.test,i)

# Qmat <- matrix(0,10,10)
# colnames(Qmat)<-rownames(Qmat) <- as.character(1:10)
# for(i in 1:10){#rows
#   for(j in 1:10){#cols
#     if(i == j+1){
#       Qmat[i,j] <- .5
#     }
#     if(j == i+1){
#       Qmat[i,j]<-.5
#     }
#   }
# }
# 
# diag(Qmat) <- -rowSums(Qmat)

#Simulate using sfreemap
hists.test2 <- sfreemap(tree = tree.test,
                        tip_states = data.test.vector,
                        Q = Qmat.test,
                        pi = "estimated",
                        model = "SYM")



describe.sfreemap(hists.test2)

?make.musse





