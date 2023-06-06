#### PACKAGES ####

library(phytools)
library(evobiR)
library(coda)
library(ggplot2)

#### LOAD IN DATA ####

dat <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                as.is=T)[,c(5,9,11,12,15)]

clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

hists <- read.simmap(file="../data/simmap_out/simmap_out.nex", 
                     format = "nexus",
                     version = 1.5)

source("pfsa2.R")
source("describe.simmap2.R")

#### DATA ORGANIZATION ####

#Remove duplicate elements
dat <- dat[!duplicated(dat$tree.name),]

#Convert species names to right format 
dat$tree.name <- tolower(dat$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
dat$tree.name <- gsub(" ","_",dat$tree.name)

#Prune to tree
dat <- dat[dat$tree.name %in% hists[[1]]$tip.label, ]

dat.all <- dat

#### FIX SCS ####

#Reset row names
row.names(dat.all) <- NULL

#Fix XO
dat.all$sex.chromosome.system <- gsub("XO females and XY males",
                                  "XY",
                                  dat.all$sex.chromosome.system)

dat.all$sex.chromosome.system <- gsub("XO",
                                  "XY",
                                  dat.all$sex.chromosome.system)

#Fix gazella X/A/Y
dat.all$sex.chromosome.system[926] <- "XXY"

#Fix Muntiacus X/A/Y
dat.all$sex.chromosome.system <- gsub("X/A, X/A/Y",
                                  "XYY",
                                  dat.all$sex.chromosome.system)

#Fix tragelauphus X/Y/A
dat.all$sex.chromosome.system <- gsub("X/X, X/Y/A",
                                  "XXY",
                                  dat.all$sex.chromosome.system)

#Fix nanger rejuvenation
dat.all$sex.chromosome.system <- gsub("rujuvenation of autosomes to sex chromosomes",
                                  "XXYY",
                                  dat.all$sex.chromosome.system)

#Fix lasiopodomys
dat.all$sex.chromosome.system[946] <- "XY"

#Fix tachyglossus
dat.all$sex.chromosome.system[948] <- "XXXXXYYYY"

#Fix X1Y1Y2
dat.all$sex.chromosome.system <- gsub("X1Y1Y2",
                                  "XYY",
                                  dat.all$sex.chromosome.system)

#Fix X1X2Y
dat.all$sex.chromosome.system <- gsub("X1X2Y1Y2",
                                  "XXYY",
                                  dat.all$sex.chromosome.system)

#Fix X1X2Y
dat.all$sex.chromosome.system <- gsub("X1X2Y",
                                  "XXY",
                                  dat.all$sex.chromosome.system)

#Factor
dat.all$sex.chromosome.system <- as.factor(dat.all$sex.chromosome.system)


#### LOOP THROUGH EACH SUBTREE ####

#Dataframe for results
prop.subtrees <- as.data.frame(matrix(c(clades,rep(NA,30)),
                                      nrow=5,
                                      ncol=7,
                                      byrow = F))

colnames(prop.subtrees) <- c("tree",
                             "obs.mean",
                             "null.mean",
                             "obs.lower",
                             "obs.upper",
                             "null.lower",
                             "null.upper")

for(i in 1:5){
  
  #### LOAD IN HISTS ####
  hists <- read.simmap(file=paste0("../data/subtree_simmap_out/subtree_",
                                   clades[i],
                                   "/simmap_fixed.nex"), 
                       format = "nexus",
                       version = 1.5)
  
  Qmat <- as.matrix(read.csv(paste0("../data/transition_matrix/subtree_matrices/matrix_",clades[i],".csv"),
                             as.is=T,header=T))
  
  
  #### TRIM DATASET ####
  
  #Trim phylogeny and chromosome data
  dat <- dat.all[dat.all$tree.name %in% hists[[1]]$tip.label, ]
  
  #### FIX SIM.STATE ####
  #Establish range of haploid autosome numbers
  chrmrng <- range(dat$hapauto)
  
  #Reassign new states
  for(j in 1:nrow(dat)){
    if(is.na(dat$hapauto[j]) == T ||
       dat$codedSCS[j] == "unclear"){
      dat$sim.state[j] <- NA
    } else if(dat$codedSCS[j] == "unfused"){
      dat$sim.state[j] <- dat$hapauto[j] - chrmrng[1] + 1
    } else if(dat$codedSCS[j] == "fused"){
      dat$sim.state[j] <- (dat$hapauto[j] - chrmrng[1] + 1) +
        (chrmrng[2] - chrmrng[1] + 1)
    }
  }
  
  #### DESCRIBE SIMMAP AND GET TR COUNT####
  
  #Summarize
  hists.summarized <- describe.simmap2(hists)
  
  #Get transitions
  trans <-  hists.summarized$count
  
  #### EXTRACT COLUMNS ASSOCIATED WITH TRANSITIONS OF INTEREST ####
  
  #Get column names
  names <- colnames(trans[,2:ncol(trans)])
  
  #Split by commas
  split.names <- strsplit(names,",")
  
  #Create vectors of rows
  SAF.cols <- c()
  AAF.cols <- c()
  
  #Establish cutoff between fused and unfused
  cutoff <- nrow(Qmat)/2
  
  #Loop through
  for(j in 1:length(split.names)){
    if(as.numeric(split.names[[j]][1]) == 
       as.numeric(split.names[[j]][2]) - (cutoff - 1) &&
       as.numeric(split.names[[j]][1]) <= cutoff &&
       as.numeric(split.names[[j]][2]) >= (cutoff + 1)){
      SAF.cols <- c(SAF.cols,j+1) 
    } else if(as.numeric(split.names[[j]][1]) == as.numeric(split.names[[j]][2])+1 &&
              as.numeric(split.names[[j]][1]) != (cutoff + 1)){
      AAF.cols <- c(AAF.cols,j+1)
    }
  }
  
  #Check colnames
  colnames(trans[,SAF.cols])
  colnames(trans[,AAF.cols])
  
  #### SUM SAF AND AA COUNTS ####
  
  #Get rowSums
  SAF.counts <- rowSums(trans[,SAF.cols])
  AAF.counts <- rowSums(trans[,AAF.cols])
  
  #Divide
  obspropSAF <- SAF.counts/(AAF.counts + SAF.counts)
  
  #### EXPERIMENTAL SAF ####
  
  expSA <- c()
  
  for(j in 1:100){
    
    #print iterations
    print(paste0("subtree ",clades[i],",tree ",j))
    
    #read in times
    times <- describe.simmap2(hists[[j]])$times[2,]
    
    #Vector to store individual proportions
    expSa.prop <- c()
    
    for(k in 1:(length(times)-1)){
      
      #Get state
      state <- as.numeric(names(times)[k])
      
      #Check if state is fused or unfused
      codedSCS <- ifelse(test = state > cutoff, 
                         yes = "fused", 
                         no = "unfused")
      
      if(codedSCS == "unfused"){
        
        #Convert state to Diploid autosome
        Da <- (state + min(dat$hapauto) - 1) * 2
        
        #Calculate proportion SAF
        expSa.prop[k] <- Pfsa2(Da = Da,
                                 scs = "XY",
                                 mud = 0.5) * times[k]
        
      } else {
        
        #Convert state to diploid autosome
        Da <- (state - (cutoff + 1) + min(dat$hapauto)) * 2
        
        #Check associated states
        scs <- as.character(unique(dat$sex.chromosome.system[
          which(dat$sim.state == state)]))
        
        if(length(scs) > 1){
          
          #Get total number of species in state
          state.total <- length(which(dat$sim.state == state))
          
          #Set proportion equal to 0
          prop <- 0
          
          for(l in scs){
            
            #Get numner of species with state and scs
            state.scs <- length(which(dat$sim.state == state &
                                        dat$sex.chromosome.system == l))
            
            #Calculate proportion and add to sum for state
            prop <- prop + as.numeric(Pfsa2(Da = Da,
                                            scs = l,
                                            mud = 0.5) * 
                                        times[k] * 
                                        state.scs/state.total)
            
          }
          
          #Bind to expSa.prop
          expSa.prop[k] <- prop
          
        } else if(length(scs) == 1){
          
          #calculate prop if single scs
          expSa.prop[k] <- Pfsa2(Da = Da,
                                   scs = scs,
                                   mud = 0.5) * times[k]
          
        } else {
          
          #Create variable for range
          state.range <- c()
          
          #Get range around state
          if(state < (cutoff + 6)){
            
            #Set lower bound
            state.range[1] <- cutoff + 1
            
            #Set upper bound
            state.range[2] <- state + 5
            
          } else if(state > (cutoff * 2) - 5){
            
            #Set lower bound
            state.range[1] <- state-5
            
            #Set upper bound
            state.range[2] <- cutoff * 2
            
          } else {
            
            #Set lower bound
            state.range[1] <- state - 5
            
            #Set upper bound
            state.range[2] <- state + 5
            
          }
          
          if(length(which(dat$sim.state %in% state.range[1]:state.range[2])) != 0){
            
            #Extract scs for data with range
            scs  <- as.character(unique(dat$sex.chromosome.system[which(
              dat$sim.state %in% state.range[1]:state.range[2])]))
            
            prop <- 0
            
            states.total <- length(which(dat$sim.state %in% state.range[1]:state.range[2]))
            
            for(l in scs){
              
              #Get number of fused for sex chromosome system
              states.scs <- length(which(dat$sim.state %in%  
                                           state.range[1]:state.range[2] & 
                                           dat$sex.chromosome.system == l))
              
              #Calculate proportion and add to sum
              prop <- prop + as.numeric(Pfsa2(Da = Da,
                                              scs = l,
                                              mud = 0.5) * 
                                          times[k] * 
                                          states.scs/states.total)
              
            }
            
            #Bind to expSa.prop
            expSa.prop[k] <- prop
            
          } else {
            
            #Calculate proportion and bind
            expSa.prop[k] <- Pfsa2(Da = Da,
                                     scs = "XY",
                                     mud = 0.5) * times[k]
            
          }
          
          # #Get total number of fused
          # total <- length(which(dat$codedSCS == "fused"))
          
          #Variable to store proportion
          
        }
      }
      
      #Set names
      names(expSa.prop[[k]]) <- state
      
    }
    
    #Sum proportions for each state
    expSA[j] <- sum(as.numeric(expSa.prop))
  }
  
  #### EXTRACT MEAN AND INTERVALS ####
  
  #Calculate and save
  prop.subtrees [i,2] <- mean(obspropSAF)
  prop.subtrees [i,3] <- mean(expSA)
  prop.subtrees [i,4] <- HPDinterval(as.mcmc(obspropSAF))[1]
  prop.subtrees [i,5] <- HPDinterval(as.mcmc(obspropSAF))[2]
  prop.subtrees [i,6] <- HPDinterval(as.mcmc(expSA))[1]
  prop.subtrees [i,7] <- HPDinterval(as.mcmc(expSA))[2]
  
  #### PLOT ####
  
  #Make plot
  overlap_plot <- ggplot()+
    geom_density(mapping = aes(obspropSAF,
                               fill="Observed"),
                 alpha=0.5,
                 bw=0.005)+
    geom_density(mapping=aes(expSA,
                             fill="Null"),
                 alpha=0.5,
                 bw=0.005)+
    ggtitle(paste0("Overlap in distribution of null vs. observed SAF"))+
    labs(subtitle = paste0(clades[i]))+
    scale_y_continuous("Density")+
    scale_fill_viridis_d()+
    xlab("Proportion sex-autosome fusion")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background= element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(face = "bold",
                                    size = 15),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold",
                                    size = 17,
                                    hjust=0.5))
  
  #Display
  plot(overlap_plot)
  
  #Save
  ggsave(overlap_plot,
         filename = paste0("../figures/subtree_overlap/overlap_",
                           clades[i],".pdf"),
         width = 7,
         height = 7,
         units = "in")
  
  
}

##### SAVE STATS ####

write.csv(prop.subtrees,
          paste0("../data/subtree_simmap_out/prop_subtrees.csv"),
          quote=F,
          row.names=T)



