library(phytools)
library(dplyr)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(3,5,9,11,12,15)]
tree <- read.tree("../data/trees/4705sp_mean.nwk")

#### DATA ORGANIZATION ####

# plot(tree, show.tip.label = T,
#      cex = 0.008,
#      no.margin = T,
#      edge.width = 0.00000000000001)

#dist.test.complete <- distTips(tree,tips = c("otomys_irroratus","chrotomys_silaceus"),
#method="patristic")

#Remove taxa with uncertain sim states
data <- subset(data,is.na(sim.state) == F)

#Convert genus/species names to right format 
data$tree.name <- tolower(data$tree.name)
data$Genus <- tolower(data$Genus)

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

data <- data.unique

rm(data.alt)

rm(data.dup,index,j,data.unique)

#### SPLIT TIP LABELS
#### TALLY SPECIES/GENUS ####

#Tally by factor level
tally.genus <- as.data.frame(table(as.factor(data$Genus)))

#Check tally
sum(tally.genus$Freq)

#Assign colnames
colnames(tally.genus) <- c("genus","freq")

#### ASSIGN CUTOFF AND IDENTIFY GENERA WHICH QUALIFY ####

#Assign cutoff
cutoff <- 5

#Identify genera
genera <- as.vector(tally.genus$genus[which(tally.genus$freq>=cutoff)])

# #Remove unnecesary objects
rm(tally.genus)

#### LOOP THROUGH QUALIFYING GENERA, SUBSET TREES, ESTIMATE RATES, BUILD MATRICES ####

#Loop
for(i in genera){

  #Subset data
  genera.data <- subset(data, Genus=="microtus")[,-c(1,3,6)]
  
  #Subset tree to match data
  genera.tree <- force.ultrametric(keep.tip(tree, genera.data$tree.name))
  
  #Scale tree to unit length
  genera.tree$edge.length <- genera.tree$edge.length/max(branching.times(genera.tree))
  
  #Save tree
  cairo_pdf(paste0("../figures/genera_trees/",i,"_tree.pdf"))
  
  plot(genera.tree, show.tip.label = T,
       cex = 1,
       no.margin = T,
       edge.width = 1)
  
  dev.off()
  
  if(length(unique(genera.data$hapauto))==1){
    print(paste0("genera ",i, " does not exhibit variation in chromosome numbers"))
    next
  }
  
  #Convert codedSCS column to right format
  for(j in 1:nrow(genera.data)){
    if(genera.data$codedSCS[j] ==  "fused"){
      genera.data$codedSCS[j] <- 0
    } else {
      genera.data$codedSCS[j] <- 1
    }
  }
  
  #Change to numeric
  genera.data$codedSCS <- as.numeric(genera.data$codedSCS)
  
  #Convert to data matrix
  data.matrix <- datatoMatrix(genera.data,
                              c(min(genera.data$hapauto),
                                max(genera.data$hapauto)),
                              hyper = T)
  
  #Make mkn model
  model <- make.mkn(genera.tree,
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
  model.mcmc <- mcmc(model.con$`likelihood function`,
                     c(1,1,1,1),
                     prior=make.prior.exponential(r=0.5),
                     #upper=c(100,100,100,100),
                     nsteps = 1000,
                     1)
  
  #Plot log probabilities and save
  cairo_pdf(paste0("../data/transition_matrix/genera_matrices/mcmc_plots/",i,".pdf"))
  plot(model.mcmc$i,model.mcmc$p,type = "l",
       main = paste0("MCMC logP ", i))
  dev.off()
  
  #profiles.plot(model.mcmc["asc1"], col.line="red")
  
  #Extract post burn portion
  model.mcmc.postburn <- model.mcmc[500:1000,]
  
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
  
  #Save matrix
  write.csv(parMat,
            paste0("../data/transition_matrix/genera_matrices/matrix_",i,".csv"),
            row.names=F,quote=F)
  
}

hist(rexp(100,0.5))









