library(phytools)
library(dplyr)
library(ape)
library(devtools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../data/chromes/dat.csv",
                 as.is=T)[,c(1,4)]

#Create vector of clade names
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#### BUILD MODEL ####
parMat <- matrix("0",
                 nrow=max(dat$codedSCS) + 1,
                 ncol=max(dat$codedSCS) + 1)
colnames(parMat) <- 0:7
rownames(parMat) <- 0:7
mat <- parMat

for(i in 1:nrow(parMat)){
  if(i != nrow(parMat)){
    parMat[i,i+1] <- "SAF"
    mat[i,i+1] <- 1
  }
  if(i != 1){
    parMat[i,1] <- "Ro"
    mat[i,1] <- 2
  }
}

#build rate table
rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat), 
                                   3))
rate.table[, 1] <- rep(as.character(1:nrow(parMat)), 
                       each = ncol(parMat))
rate.table[, 2] <- rep(as.character(1:ncol(parMat)), nrow(parMat))
rate.table[, 3] <- as.character(c(t(parMat)))
rate.table <- rate.table[rate.table[, 1] != rate.table[, 
                                                       2], ]

#get formula for constraining
formulae <- c()
for (i in 1:nrow(rate.table)) {
  formulae[i] <- paste("q", rate.table[i, 1], rate.table[i, 
                                                         2], " ~ ", rate.table[i, 3], collapse = "", sep = "")
}

extras <- c("SAF", "Ro")

#### LOOP THROUGH SUBSET TREES, ESTIMATE RATES, BUILD MATRICES ####

#Loop
for(i in 1:5){
  
  # load in subtree
  split.tree <- read.tree(paste0("../data/trees/subtrees/tree_",clades[i],".nex"))
    
  #Subset data
  split.data <- subset(dat, tree.name %in% split.tree$tip.label)

  #### MODEL ####
  
  #### BUILD MODEL ####
  parMat <- matrix("0",
                   nrow=max(dat$codedSCS) + 1,
                   ncol=max(dat$codedSCS) + 1)
  colnames(parMat) <- 0:7
  rownames(parMat) <- 0:7

  for(j in 1:nrow(parMat)){
    if(j != nrow(parMat)){
      parMat[j,j+1] <- "SAF"
    }
    if(j != 1){
      parMat[j,1] <- "Ro"
    }
  }
  
  #make data matrix
  data.matrix <- matrix(0,nrow=nrow(split.data),
                        ncol=ncol(parMat))
  
  colnames(data.matrix) <- colnames(parMat)
  rownames(data.matrix) <- split.data$tree.name
  
  #fill data matrix
  for(j in 1:nrow(data.matrix)){
    data.matrix[j,as.numeric(split.data$codedSCS[j]) + 1] <- 1
  }
  
  #Make mkn model
  model <- make.mkn(split.tree,
                    data.matrix,
                    ncol(data.matrix),
                    strict=F,
                    control=list(method="ode"))

  #Check argnames: 9506 different rates (98*97, expected)
  model.con <- constrain(model, formulae = formulae, extra = extras)
  
  # Check args: 4 different rates in constrained model, output is as expected

  #Run mcmc
  model.mcmc <- diversitree::mcmc(lik=model.con,
                                  x.init=c(1,1),
                                  prior=make.prior.exponential(r=0.5),
                                  #upper=c(100,100,100,100),
                                  nsteps = 500,
                                  w=1)

  #profiles.plot(model.mcmc["asc1"], col.line="red")
  #### BUILD QMATRIX ####
  #Extract post burn portion
  model.mcmc.postburn <- model.mcmc[450:500,]

  #Get mean params
  params <- c(mean(model.mcmc.postburn$SAF),
              mean(model.mcmc.postburn$Ro))
  
  names(params) <- colnames(model.mcmc[,2:3])

  #Sub into matrix
  parMat[parMat == "SAF"] <- params[1]
  parMat[parMat == "Ro"] <- params[2]

  #Convert to df
  parMat <- as.data.frame(parMat)
  parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

  #Convert back to df
  parMat <- as.data.frame(parMat)

  #Set diagonals
  diag(parMat) <- -rowSums(parMat)

  #Save matrix
  write.csv(parMat,
            paste0("../data/transition_matrix/subtree_matrices/multiSAF/matrix_",clades[i],".csv"),
            row.names=F,quote=F)

}

