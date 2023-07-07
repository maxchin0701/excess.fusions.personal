library(phytools)
library(dplyr)
library(ape)
library(devtools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../data/chromes/dat.csv",
                 as.is=T)
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
rm(edges, i)


#### LOOP THROUGH SUBSET TREES, ESTIMATE RATES, BUILD MATRICES ####

#Loop
for(i in 1:5){
  
  # #Subset data
  # split.data <- subset(data, tree.name %in% trees[[i]]$tip.label)[,-c(1,3,6)]
  # 
  #Subset tree to match data
  split.tree <- force.ultrametric(trees[[i]],
                                  method = "extend")

  #Scale tree to unit length
  split.tree$edge.length <- split.tree$edge.length/max(branching.times(split.tree))

  #Save tree PDF
  cairo_pdf(paste0("../figures/subtrees/",clades[i],"_tree.pdf"))

  plot(split.tree, show.tip.label = T,
       cex = sqrt(15/length(split.tree$tip.label)),
       no.margin = T,
       edge.width = 1)

  dev.off()

  #Convert codedSCS column to right format
  for(j in 1:nrow(split.data)){
    if(split.data$codedSCS[j] ==  "fused"){
      split.data$codedSCS[j] <- 0
    } else {
      split.data$codedSCS[j] <- 1
    }
  }

  #Change to numeric
  split.data$codedSCS <- as.numeric(split.data$codedSCS)

  #Convert to data matrix
  data.matrix <- datatoMatrix(split.data,
                              c(min(split.data$hapauto),
                                max(split.data$hapauto)),
                              hyper = T)

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

  #Plot log probabilities and save
  cairo_pdf(paste0("../data/transition_matrix/subtree_matrices/mcmc_plots/",clades[i],".pdf"))
  plot(model.mcmc$i,model.mcmc$p,type = "l",
       main = paste0("MCMC logP ", clades[i]))
  dev.off()

  #profiles.plot(model.mcmc["asc1"], col.line="red")

  #Extract post burn portion
  model.mcmc.postburn <- model.mcmc[50:100,]

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
            paste0("../data/transition_matrix/subtree_matrices/matrix_",clades[i],".csv"),
            row.names=F,quote=F)

  #Save tree newick
  write.tree(split.tree,
             paste0("../data/trees/subtrees/tree_",clades[i],".nwk"))
  
}

