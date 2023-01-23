library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]
tree <- force.ultrametric(read.tree("../data/trees/4705sp_mean.nwk"))

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

data.unique <- data.alt

data <- data.unique[,-c(2,5)]

rm(data.alt,data.dup,data.unique,j,index)

#### MODEL ####

#Convert codedSCS column to right format
for(i in 1:nrow(data)){
  if(data$codedSCS[i] ==  "fused"){
    data$codedSCS[i] <- 0
  } else {
    data$codedSCS[i] <- 1
  }
}

data$codedSCS <- as.numeric(data$codedSCS)

#Convert to data matrix
data.matrix <- datatoMatrix(data,
                            c(2,50),
                            hyper = T)

#Make mkn model
model <- make.mkn(tree,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#Check argnames: 9506 different rates (98*97, expected)
argnames(model)

#Constrain model
model.con <- constrainMkn(data.matrix,
                          model,
                          hyper = T,
                          polyploidy = F,
                          verbose = T,
                          constrain = list(saf.model=T))

# Check args: 4 different rates in constrained model, output is as expected
argnames(model.con$`likelihood function`)

#Check parMat
parMat <- model.con$`parameter matrix`

#Run mcmc
model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                x.init=c(1,1,1,1),
                                prior=make.prior.exponential(r=2),
                                #upper=c(100,100,100,100),
                                nsteps = 100,
                                w=1)

plot(model.mcmc$i,model.mcmc$p,type = "l")

profiles.plot(model.mcmc["asc1"], col.line="red")

#Extract post burn portion
model.mcmc.postburn <- model.mcmc[25:100,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1),
            mean(model.mcmc.postburn$tranSAF),
            mean(model.mcmc.postburn$tranRo))

names(params) <- colnames(model.mcmc[,2:5])

#Sub into matrix
parMat[parMat == "asc1"] <- params[1]
parMat[parMat == "desc1"] <- params[2]
parMat[parMat == "tranSAF"] <- params[3]
parMat[parMat == "tranRo"] <- params[4]

#Convert to df
parMat <- as.data.frame(parMat)

#Mutate to right format
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

#Convert back to df
parMat <- as.data.frame(parMat)

#Set diagonals
diag(parMat) <- -rowSums(parMat)

#Save matrix
write.csv(parMat,
          paste0("../data/transition_matrix/Q_matrix_mammalian.csv"),
          row.names=F,col.names=F,quote=F)









