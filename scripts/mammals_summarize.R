#### PACKAGES ####

library(ape)
library(phytools)
library(doSNOW)


#### LOAD DATA ####

hists <- read.simmap(file=paste0("../data/simmap_out/simmap_fixed.nex"), 
                     format = "nexus",
                     version = 1.5)

source("describe.simmap2.R")

#### SUMMARIZE TRANSITION COUNTS ####

#Get summary using describe.simmap
hists.summary <- describe.simmap2(hists)

#Extract transition counts
count <- hists.summary$count

#Change colnames so saving as .csv doesnt mess up
colnames(count) <- gsub(",","_",colnames(count))

#### SAVE TRANSITION TIMES FOR EACH TREE ####

for(i in 1:100){
  
  #Save times to object
  times <- describe.simmap2(hists[[i]])$times
  
  #Save object
  write.csv(times,
            paste0("../data/simmap_out/dwelling_times_fixed/dwelling_times_",i,".csv"),
            quote=F,
            row.names=T)
}

#### SAVE TRANSITION COUNTS ####

write.csv(count,
          paste0("../data/simmap_out/transitions_fixed.csv"),
          quote=F,
          row.names=T)



     