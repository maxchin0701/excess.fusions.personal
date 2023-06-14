#### LOAD PACKAGES ####
library(ape)
library(phytools)
library(viridis) 

#### LOAD IN TEST STOCHASTIC MAP ####
load("../outputs/test.hist.RData")
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,3)]

#### SUBSET TREE TO JUST ARTIODACTYLS ####

#subset map
hist.subset <- extract.clade.simmap(tree = hists.fixed.alt[[1]],
                                    node=hists.fixed.alt[[1]]$edge[448,2])

dat <- dat[which(dat$tree.name %in% hist.subset$tip.label),]

#plot
cols <- c(viridis(49))
names(cols) <- c(1:49)
plotSimmap(hists.fixed.alt[[1]],col=cols)

#### GET NODE STATES ####

#vector to store
node.states <- c()

for(i in 143:283){
  
  #get edge which descends from node
  desc.edge <- which(hist.subset$edge[,1] == i)[1]
  
  #get initial state of descendent edge
  node.state <- names(hist.subset$maps[[desc.edge]][1])
  
  #check that it matches with last state of leading edge
  if(i != 143){
    
    lead.edge <- which(hist.subset$edge[,2] == i)
    
    lead.state <- names(hist.subset$maps[[lead.edge]][length(hist.subset$maps[[lead.edge]])])
    
    if(lead.state != node.state){
      print("warning: node state does not match last state of leading edge")
    }
  }
  node.states[i-142] <- node.state
}

#### PLOT ####
plotSimmap(hist.subset,col=cols,fsize=0.1,lwd=0.1,pts=T)
plotSimmap(hist.subset,col=cols,fsize=0.0000001,lwd=0.1)
nodelabels(node.states,frame="none",cex=0.6)
tiplabels(dat$hapauto,frame="none",cex=0.6)




