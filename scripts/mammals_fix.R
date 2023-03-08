#### PACKAGES ####

library(ape)
library(phytools)
library(igraph)


#### LOAD IN DATA ####

hists <- read.simmap(file="../data/simmap_out/simmap_out.nex", 
                     format = "nexus",
                     version = 1.5)

transition.matrix <- as.matrix(read.csv("../data/transition_matrix/matrix_mammalian.csv",
                                        as.is = T,
                                        header = F))

tip.dat <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                    as.is=T)[,c(5,15)]

source("fix_simmap.R")

#### PRUNE TIP DATA ####

#Remove duplicate elements
tip.dat <- tip.dat[!duplicated(tip.dat$tree.name),]

#Convert species names to right format 
tip.dat$tree.name <- tolower(tip.dat$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
tip.dat$tree.name <- gsub(" ","_",tip.dat$tree.name)

#Prune to tree
tip.dat <- tip.dat[tip.dat$tree.name %in% hists[[1]]$tip.label, ]


#### RUN FIX.SIMMAP ####

hists.fixed <- fix.simmap(hists = hists,
                          tips = tip.dat,
                          transition.matrix = transition.matrix)


#### PLOT ####

#Establish palette
cols <- c(rainbow(49),viridis(49), "black")
names(cols) <- c(1:98, "fail")

#Plot each of 100 fixed trees
for(i in 1:100){
  
  cairo_pdf(paste0("../figures/fixed_simmap_plots/tree_",i,".pdf"))
  plotSimmap(hists.fixed[[i]],col=cols,fsize=0.03,lwd = 0.1)
  dev.off()
  
}

#### SAVE SIMMAP OBJECT ####

write.simmap(hists.fixed,
             file="../data/simmap_out/simmap_fixed.nex",
             version = 1.5,
             format="nexus")






