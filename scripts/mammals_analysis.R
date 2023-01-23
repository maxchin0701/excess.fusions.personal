#### PACKAGES ####

library(phytools)
library(evobiR)
library(coda)
library(ggplot2)

#### LOAD IN DATA ####
trans <- read.csv("../data/simmap_out/transitions_fixed.csv",
                  as.is = T,check.names = F)

dat <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                as.is=T)[,c(5,9,11,12,15)]

hists <- read.simmap(file="../data/simmap_out/simmap_out.nex", 
                     format = "nexus",
                     version = 1.5)

source("pfsa2.R")

#### PRUNE TIP DATA ####

#Remove duplicate elements
dat <- dat[!duplicated(dat$tree.name),]

#Convert species names to right format 
dat$tree.name <- tolower(dat$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
dat$tree.name <- gsub(" ","_",dat$tree.name)

#Prune to tree
dat <- dat[dat$tree.name %in% hists[[1]]$tip.label, ]

#### FIX SCS ####

#Reset row names
row.names(dat) <- NULL

#Fix XO
dat$sex.chromosome.system <- gsub("XO females and XY males",
                                  "XY",
                                  dat$sex.chromosome.system)

dat$sex.chromosome.system <- gsub("XO",
                                  "XY",
                                  dat$sex.chromosome.system)

#Fix gazella X/A/Y
dat$sex.chromosome.system[926] <- "XXY"

#Fix Muntiacus X/A/Y
dat$sex.chromosome.system <- gsub("X/A, X/A/Y",
                                  "XYY",
                                  dat$sex.chromosome.system)

#Fix tragelauphus X/Y/A
dat$sex.chromosome.system <- gsub("X/X, X/Y/A",
                                  "XXY",
                                  dat$sex.chromosome.system)

#Fix nanger rejuvenation
dat$sex.chromosome.system <- gsub("rujuvenation of autosomes to sex chromosomes",
                                  "XXYY",
                                  dat$sex.chromosome.system)

#Fix lasiopodomys
dat$sex.chromosome.system[946] <- "XY"

#Fix tachyglossus
dat$sex.chromosome.system[948] <- "XXXXXYYYY"

#Fix X1Y1Y2
dat$sex.chromosome.system <- gsub("X1Y1Y2",
                                  "XYY",
                                  dat$sex.chromosome.system)

#Fix X1X2Y
dat$sex.chromosome.system <- gsub("X1X2Y1Y2",
                                  "XXYY",
                                  dat$sex.chromosome.system)

#Fix X1X2Y
dat$sex.chromosome.system <- gsub("X1X2Y",
                                  "XXY",
                                  dat$sex.chromosome.system)

#Factor
dat$sex.chromosome.system <- as.factor(dat$sex.chromosome.system)

#### EXTRACT COLUMNS ASSOCIATED WITH TRANSITIONS OF INTEREST ####

#Get column names
names <- colnames(trans[,3:ncol(trans)])

#Split by commas
split.names <- strsplit(names,"_")

#Create vectors of rows
SAF.cols <- c()
AAF.cols <- c()

#Loop through
for(i in 1:length(split.names)){
  if(as.numeric(split.names[[i]][1]) == as.numeric(split.names[[i]][2])-48 &&
     as.numeric(split.names[[i]][1]) <= 49 &&
     as.numeric(split.names[[i]][2]) >= 50){
    SAF.cols <- c(SAF.cols,i+2)
  } else if(as.numeric(split.names[[i]][1]) == as.numeric(split.names[[i]][2])+1 &&
            as.numeric(split.names[[i]][1]) != 50){
    AAF.cols <- c(AAF.cols,i+2)
  }
}

colnames(trans[SAF.cols])

#### SUM SAF AND AA COUNTS ####

#Get rowSums
SAF.counts <- rowSums(trans[,SAF.cols])
AAF.counts <- rowSums(trans[,AAF.cols])

#Divide
obspropSAF <- SAF.counts/(AAF.counts + SAF.counts)

#### EXPERIMENTAL SAF ####

expSA <- c()

for(i in 1:100){
  
  #Read in times
  times <- read.csv(paste0("../data/simmap_out/dwelling_times_fixed/dwelling_times_",
                           i,".csv"),
                    as.is = T,check.names = F)[2,]
  
  #Vector to store individual proportions
  expSa.prop <- c()
  
  for(j in 2:(ncol(times)-1)){
    
    #Get state
    state <- as.numeric(colnames(times)[j])
    
    #Check if state is fused or unfused
    codedSCS <- ifelse(test = state > 49, 
                  yes = "fused", 
                  no = "unfused")
    
    if(codedSCS == "unfused"){
      
      #Convert state to Diploid autosome
      Da <- (state + 1) * 2
      
      #Calculate proportion SAF
      expSa.prop[j-1] <- Pfsa2(Da = Da,
                          scs = "XY",
                          mud = 0.5) * times[j]
      
    } else {
      
      #Convert state to diploid autosome
      Da <- (state - 48) * 2
      
      #Check associated states
      scs <- as.character(unique(dat$sex.chromosome.system[
        which(dat$sim.state == state)]))
      
      if(length(scs) > 1){
        
        #Get total number of species in state
        state.total <- length(which(dat$sim.state == state))
        
        #Set proportion equal to 0
        prop <- 0
        
        for(k in scs){

          #Get numner of species with state and scs
          state.scs <- length(which(dat$sim.state == state &
                                      dat$sex.chromosome.system == k))
          
          #Calculate proportion and add to sum for state
          prop <- prop + as.numeric(Pfsa2(Da = Da,
                                           scs = k,
                                           mud = 0.5) * 
                                      times[j] * 
                                      state.scs/state.total)
          
        }
        
        #Bind to expSa.prop
        expSa.prop[j-1] <- prop
        
      } else if(length(scs) == 1){
        
        #calculate prop if single scs
        expSa.prop[j-1] <- Pfsa2(Da = Da,
                                 scs = scs,
                                 mud = 0.5) * times[j]
        
      } else {
        
        #Create variable for range
        state.range <- c()
        
        #Get range around state
        if(state < 55){
          
          #Set lower bound
          state.range[1] <- 50
          
          #Set upper bound
          state.range[2] <- state + 5
          
        } else if(state > 93){
          
          #Set lower bound
          state.range[1] <- state-5
          
          #Set upper bound
          state.range[2] <- 98
          
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
          
          for(k in scs){

            #Get number of fused for sex chromosome system
            states.scs <- length(which(dat$sim.state %in%  
                                        state.range[1]:state.range[2] & 
                                        dat$sex.chromosome.system == k))
            
            #Calculate proportion and add to sum
            prop <- prop + as.numeric(Pfsa2(Da = Da,
                                            scs = k,
                                            mud = 0.5) * 
                                        times[j] * 
                                        states.scs/states.total)
            
          }
          
          #Bind to expSa.prop
          expSa.prop[j-1] <- prop
          
        } else {
          
          #Calculate proportion and bind
          expSa.prop[j-1] <- Pfsa2(Da = Da,
                                   scs = "XY",
                                   mud = 0.5) * times[j]
          
        }
        
        # #Get total number of fused
        # total <- length(which(dat$codedSCS == "fused"))
        
        #Variable to store proportion
        
      }
    }
    
    #Set names
    names(expSa.prop[[j-1]]) <- state
    
  }
  
  #Sum proportions for each state
  expSA[i] <- sum(as.numeric(expSa.prop))
}

#Extract mean and intervals
mean(obspropSAF)
mean(expSA)
HPDinterval(as.mcmc(obspropSAF))
HPDinterval(as.mcmc(expSA))

#### PLOT ####

#Plot
overlap_plot <- ggplot()+
                  geom_density(mapping = aes(obspropSAF,
                                             fill="Observed"),
                               alpha=0.5)+
                  geom_density(mapping=aes(expSA,
                                           fill="Expected"),
                               alpha=0.5)+
                  ggtitle(paste0("Overlap in distribution of expected vs. observed SAF"))+
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

plot(overlap_plot)



