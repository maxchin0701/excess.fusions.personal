library(phytools)
library(dplyr)
library(ape)
library(devtools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../data/chromes/dat.csv",
                 as.is=T)[,c(1,3)]
tree <- read.tree("../data/trees/tree.nex")

#### CUT TREE ####

#Define edges of interest
edges <- c(181,448,796,1285,1718)

#Set up list
trees <- list()

#Loop through edges of interest
for(i in 1:length(edges)){

  trees[[i]] <- extract.clade(tree,
                             node=tree$edge[edges[i],2])

}

#Create vector of clade names
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#Remove redundant variables
rm(edges)


#### LOOP THROUGH SUBSET TREES, ESTIMATE RATES, BUILD MATRICES ####

#Loop
for(i in 1:5){
  
  #Subset data
  split.data <- subset(dat, tree.name %in% trees[[i]]$tip.label)

  # Subset tree to match data
  split.tree <- force.ultrametric(trees[[i]],
                                  method = "extend")

  #Scale tree to unit length
  split.tree$edge.length <- split.tree$edge.length/max(branching.times(split.tree))
  
  #### MODEL ####
  
  #Convert to data matrix
  data.matrix <- datatoMatrix(split.data,
                              c(2,50),
                              hyper = F)

  #Make mkn model
  model <- make.mkn(split.tree,
                    data.matrix,
                    ncol(data.matrix),
                    strict=F,
                    control=list(method="ode"))

  #Check argnames: 9506 different rates (98*97, expected)
  argnames(model)

  #Constrain model
  model.con <- constrainMkn(data.matrix,
                            model,
                            hyper = F,
                            polyploidy = F,
                            verbose = T,
                            constrain = list(drop.poly=T,
                                             drop.demi=T))

  # Check args: 4 different rates in constrained model, output is as expected
  argnames(model.con$`likelihood function`)

  #Check parMat
  parMat <- model.con$`parameter matrix`

  #Run mcmc
  model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                  x.init=c(1,1),
                                  prior=make.prior.exponential(r=10),
                                  #upper=c(100,100,100,100),
                                  nsteps = 500,
                                  w=1)

  #profiles.plot(model.mcmc["asc1"], col.line="red")
  #### BUILD QMATRIX ####
  #Extract post burn portion
  model.mcmc.postburn <- model.mcmc[450:500,]

  #Get mean params
  params <- c(mean(model.mcmc.postburn$asc1),
              mean(model.mcmc.postburn$desc1))
  
  names(params) <- colnames(model.mcmc[,2:3])

  #Sub into matrix
  parMat[parMat == "asc1"] <- params[1]
  parMat[parMat == "desc1"] <- params[2]

  #Convert to df
  parMat <- as.data.frame(parMat)
  parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

  #Convert back to df
  parMat <- as.data.frame(parMat)

  #Set diagonals
  diag(parMat) <- -rowSums(parMat)

  #Save matrix
  write.csv(parMat,
            paste0("../data/transition_matrix/subtree_matrices/hapauto/matrix_",clades[i],".csv"),
            row.names=F,quote=F)

  #Save tree newick
  write.tree(split.tree,
             paste0("../data/trees/subtrees/tree_",clades[i],".nex"))
  
}

