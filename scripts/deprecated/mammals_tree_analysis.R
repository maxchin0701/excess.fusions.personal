library(ape)
library(dplyr)
library(plyr)
library(phytools)
library(coda)
library(diversitree)
library(chromePlus)
library(lattice)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]
tree <- read.tree("../data/trees/4705sp_mean.nwk")

#### DATA ORGANIZATION ####

# plot(tree, show.tip.label = T,
#      cex = 0.008,
#      no.margin = T,
#      edge.width = 0.00000000000001)

#Scale tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#dist.test.complete <- distTips(tree,tips = c("otomys_irroratus","chrotomys_silaceus"),
                               #method="patristic")

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

#Use temp df to reorder tip data to match tree order
data.alt <- data.unique
 
for(j in 1:length(tree$tip.label)){
 index <- which(data.unique$tree.name == tree$tip.label[j])
 data.alt[j,] <- data.unique[index,]
}
 
data.unique <- data.alt
 
rm(data.alt)

#### MAKE PHYLOGENY WITH TREE NAMES + HAPAUTO ####

#Convert data to vector form
data.vector <- c()

for(i in 1:nrow(data.unique)){
  data.vector <- c(data.vector,data.unique$sim.state[i])
}

#Set names
names(data.vector) <- data.unique$tree.name

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

#Create temporary alternate tree
tree.alt <- tree

#Rename tip labels with tips vector
tree.alt$tip.label <- tips

#Assign tree.alt back to tree, remove redundant tree.alt (keep original tree for later)
tree.orig <- tree

tree <- tree.alt

rm(tree.alt)

#Plot tree and save as PDF
cairo_pdf("../figures/chrom_phylogeny.pdf")

plot(tree, show.tip.label = T,
     cex = 0.02,
     no.margin = T,
     edge.width = 0.00000000000001)

dev.off()

#### ANALYZE CHROM-DIST RATIO OF ADJACENT TIPS ####

#Get distance between all tips
dist <- as.data.frame(cophenetic.phylo(tree))

#Extract distances of interest(tips which are adjacent)
chrom.dist <- data.frame(NULL)

#Loop through dist rows to build chrom.dist 
for(i in 1:(nrow(dist)-1)){

  #Extract distance of interest
  dist.adj <- dist[i,i+1]
  
  #Calculate difference in chromosome numbers
  if(data.unique$hapauto[i+1] >= data.unique$hapauto[i]){
    chrom.diff <- data.unique$hapauto[i+1] - data.unique$hapauto[i]
  } else {
    chrom.diff <- data.unique$hapauto[i] - data.unique$hapauto[i+1]
  }
  
  #divide by distance
  chrom.dist.ratio <- chrom.diff/dist.adj
  
  #Bind to df
  chrom.dist <- rbind(chrom.dist,
                      c(data.unique$tree.name[i],
                        data.unique$tree.name[i+1],
                        data.unique$hapauto[i],
                        data.unique$hapauto[i+1],
                        dist.adj,
                        chrom.diff,
                        chrom.dist.ratio)
                      )
}

#Remove unnecessary objects
rm(chrom.diff,chrom.dist.ratio,dist.adj,index,name1,name2)

#Name columns
colnames(chrom.dist) <- c("species1",
                          "species2",
                          "species1.hapauto",
                          "species2.hapauto",
                          "tree.dist",
                          "chrom.diff",
                          "chrom.dist.ratio")

#Convert to right format
chrom.dist$species1.hapauto <- as.numeric(chrom.dist$species1.hapauto)
chrom.dist$species2.hapauto <- as.numeric(chrom.dist$species2.hapauto)
chrom.dist$tree.dist <- as.numeric(chrom.dist$tree.dist)
chrom.dist$chrom.diff <- as.numeric(chrom.dist$chrom.diff)
chrom.dist$chrom.dist.ratio <- as.numeric(chrom.dist$chrom.dist.ratio)

#Invert rows to make comparisons easier
chrom.dist <- chrom.dist[order(nrow(chrom.dist):1),]
row.names(chrom.dist) <- NULL

#Scatterplot of tree distance vs chr change
cairo_pdf(filename = "../figures/chr_dist_scatter.pdf")

plot(x = chrom.dist$chrom.diff,
     y = chrom.dist$tree.dist,
     xlab = "Chromosome number change",
     ylab = "Pairwise distance (scaled branch length)",
     main = "Relationship between pairwise distance and Chr. number change in adjacent taxa",
     cex.main = 0.75,
     cex.lab=0.75)

dev.off()

#Histogram of chr.change-distance ratios
cairo_pdf(filename = "../figures/chr_dist_hist.pdf")

hist(chrom.dist$chrom.dist.ratio,
     main = "Distribution of chromosome change - distance ratios among adjacent taxa",
              xlab = "Chromosome change divided by pairwise distance",
              cex.main = 0.8,
              cex.lab = 0.8)

dev.off()

#### IDENTIFY OUTLIERS #####

#Check HPD intervals of chrom.dist.ratio
upper <- HPDinterval(as.mcmc(chrom.dist$chrom.dist.ratio))[2]

#Identify taxa with chrom.dist.ratio > upper HPD
to.rm <- which(chrom.dist$chrom.dist.ratio > upper)

to.rm <- chrom.dist[to.rm,]

#Make vector of names to remove
names.rm <- c("dipodomys_heermanni","sicista_strandi","sigmodon_hispidus",
              "akodon_serrensis","akodon_orophilus","akodon_montensis",
              "eligmodontia_puerulus","taterillus_pygargus","meriones_crassus",
              "rattus_rattus","chrotomys_silaceus","otomys_irroratus",
              "stenocephalemys_albipes","mastomys_erythroleucus","eliurus_minor",
              "eothenomys_proditor","ellobius_talpinus","microtus_oregoni",
              "microtus_canicaudus","microtus_montanus","microtus_dogramacii",
              "octodon_degus","blarina_hylophaga","sorex_araneus")



#### MODEL ####

#Create dataframe with removed names
data.model <- data.unique[,-c(2,5)]

#Match tree to df
tree.model <- force.ultrametric(keep.tip(tree.orig, data.model$tree.name))

#rm(tree.orig)

#Convert codedSCS column to right format
for(i in 1:nrow(data.model)){
  if(data.model$codedSCS[i] ==  "fused"){
    data.model$codedSCS[i] <- 0
  } else {
    data.model$codedSCS[i] <- 1
  }
}

data.model$codedSCS <- as.numeric(data.model$codedSCS)

#Convert to data matrix
data.matrix <- datatoMatrix(data.model,
                            c(2,50),
                            hyper = T)

#Make mkn model
model <- make.mkn(tree.model,
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
                   prior=make.prior.exponential(r=5),
                   #upper=c(100,100,100,100),
                   nsteps = 100,
                   1)

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

#### GET TIP RATES AND PLOT ONTO PHYLOGENY #####

#Convert tip states to named vector
tip.states <- as.matrix(data.unique[,-c(2,4,5)])[,2]
names(tip.states) <- data.unique$tree.name

#Convert tip.states vector to numeric
tip.states <- as.numeric(tip.states)

#Run ancestral character estimation
fit <- ace(tip.states,tree.orig,type = "continuous")

#Extract state probabilities for each node
node.probs <- as.data.frame(fit$ace)

#Plot tree
plot(tree.orig,show.tip.label = F,edge.width = .1)
tiplabels(data.unique$hapauto, frame = "none",cex=.1)
nodelabels(node.probs$`fit$ace`,frame = "none",cex=0.1)

#Get tip-node distances
node.dist <- dist.nodes(tree.orig)

node.dist <- node.dist[1:950,951:1898]

#Make tip rates df
tip.rates <- data.frame(NULL)

#Loop through all tip names
for(i in 1:length(tree.orig$tip.label)){
  
  # #Extract distance between tip and all other tips/nods
  # node.dist <- as.data.frame(distNodes(tree,node=tree$tip.label[i]))
  # 
  # #Remove rows which have distance of 0(no separation == same node)
  # node.dist <- node.dist[-which(node.dist$time == 0),]
  
  #Extract closest node to tip
  adj.node <- as.numeric(colnames(node.dist)[which(node.dist[i,] == 
                                                  min(node.dist[i,]))])
  
  #Extract distance between closest node and tip
  adj.node.dist <- node.dist[i,paste0("",adj.node,"")]
  
  #Extract simulated state for that node
  adj.state <- node.probs[paste0("",adj.node,""),1]
  
  #Extract difference between tip state and node state
  state.diff <- abs(adj.state - data.unique$hapauto[i])
  
  #Bind to df
  tip.rates <- rbind(tip.rates,
                     c(data.unique$tree.name[i],
                       as.numeric(adj.node),
                       as.numeric(data.unique$hapauto[i]),
                       as.numeric(adj.state),
                       as.numeric(state.diff),
                       as.numeric(adj.node.dist),
                       as.numeric(state.diff/adj.node.dist)))
}

#Set column names
colnames(tip.rates) <- c("species",
                         "closest.node",
                         "tip.state",
                         "node.state",
                         "chrom.change",
                         "pair.dist",
                         "tip.rate")

#Get columns into right format
tip.rates$closest.node <- as.numeric(tip.rates$closest.node)
tip.rates$tip.state <- as.numeric(tip.rates$tip.state)
tip.rates$node.state <- as.numeric(tip.rates$node.state)
tip.rates$chrom.change <- as.numeric(tip.rates$chrom.change)
tip.rates$pair.dist <- as.numeric(tip.rates$pair.dist)
tip.rates$tip.rate <- as.numeric(tip.rates$tip.rate)

#Create vector of tip rates
tip.rates.vector <- as.vector(tip.rates$tip.rate)

names(tip.rates.vector) <- tip.rates$species

#Plot tree with bars
cairo_pdf("../figures/phylogeny_tip_rates.pdf")

plotTree.wBars(tree.orig,
               x = tip.rates.vector,
               show.tip.label = F,
               method = "plotTree",
               width=0.1,
               col = "darkslategray4",
               no.margin=T,
               border = NA,
               lwd=0.001)
tiplabels(tree$tip.label, frame = "none",cex=.02,offset = 0.01)
edgelabels(cex = 0.03, width = 0.001,frame="none")

dev.off()

getParent(tree.orig,node = 4)





#### BUILD SIMPLIFIED MODEL ####

#Remove unnecessary columns
data.simp <- data.unique[,-c(2,4,5)]

#Make data.matrix
data.matrix.simp <- datatoMatrix(data.simp,
                                  c(2,50),
                                  hyper = F)

#Build model
model.simp <- make.mkn(force.ultrametric(tree.orig),
                        data.matrix.simp,
                        ncol(data.matrix.simp),
                        strict=F,
                        control=list(method="ode"))

#Check argnames
argnames(model.simp)

#Constrain
model.simp.con <- constrainMkn(data.matrix.simp,
                               model.simp,
                               hyper = F,
                               polyploidy = F,
                               verbose = T,
                               constrain = list(drop.poly=T,
                                                drop.demi=T))

#Check argnames
argnames(model.simp.con$`likelihood function`)

#Check parmat
parmat.simp <- model.simp.con$`parameter matrix`

#Run MCMC
model.mcmc.simp <- mcmc(model.simp.con$`likelihood function`,
                        c(1,1),
                        prior=make.prior.exponential(r=5),
                        #upper=c(100,100,100,100),
                        nsteps = 100,
                        1)

#Plot to check for convergence
plot(model.mcmc.simp$i,model.mcmc$p,type = "l")

profiles.plot(model.mcmc.simp["asc1"], col.line="red")

#Extract post burn portion
simp.mcmc.postburn <- model.mcmc.simp[25:100,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1))

names(params) <- colnames(model.mcmc[,2:3])
















   