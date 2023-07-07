#### PACKAGES ####

library(phytools)
library(evobiR)
library(coda)
source("functions.R")

clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#### LOOP THROUGH
for(i in 1:5){
  load(file=paste0("../outputs/hapauto_maps/subtrees/hists.",clades[i],".summarized.RData"))
  hapauto.summarized <- hists.summarized
  load(file=paste0("../outputs/SAF_maps/subtrees/hists.",clades[i],".summarized.RData"))
  SAF.summarized <- hists.summarized
  
  rm(hists.summarized)
  
  #### EXTRACT OBSERVED TRANSITIONS ####
  
  #get column names
  SAF.names <- colnames(SAF.summarized$count[,2:ncol(SAF.summarized$count)])
  hapauto.names <- colnames(hapauto.summarized$count[,2:ncol(hapauto.summarized$count)])
  
  #Split by commas
  SAF.split.names <- strsplit(SAF.names,",")
  hapauto.split.names <- strsplit(hapauto.names,",")
  
  #Create vectors of rows
  SAF.cols <- c()
  hapauto.cols <- c()
  
  #Loop through
  for(j in 1:length(hapauto.split.names)){
    if(j <= length(SAF.split.names) &&
       as.numeric(SAF.split.names[[j]][1]) == as.numeric(SAF.split.names[[j]][2])-1){
      SAF.cols <- c(SAF.cols,j+1)
    } 
    
    if(as.numeric(hapauto.split.names[[j]][1]) == as.numeric(hapauto.split.names[[j]][2])+1){
      hapauto.cols <- c(hapauto.cols,j+1)
    }
  }
  
  #### SUM SAF AND AA COUNTS ####
  
  #Get rowSums
  SAF.counts <- rowSums(SAF.summarized$count[,SAF.cols])
  total.counts <- rowSums(hapauto.summarized$count[,hapauto.cols])
  
  #Divide
  obspropSAF <- SAF.counts/(total.counts)
  
  #### NULL SAF ####
  
  #null proportions vector
  expSA <- c()
  
  for(j in 1:100){
    
    #print iterations
    print(paste0("Tree ",j))
    
    #get current times
    times <- hapauto.summarized$times[j,1:(ncol(hapauto.summarized$times) - 1)]/
      hapauto.summarized$times[j,ncol(hapauto.summarized$times)]
    
    #vector for current map proportions
    expSA.cur <- c()
    
    #loop through states
    for(k in 1:length(times)){
      
      #convert state to diploid autosome
      Da <- (k + 1) * 2
      
      #Calculate proportion SAF
      expSA.cur[k] <- Pfsa2(Da = Da,
                               scs = "XY",
                               mud = 0.5) * times[k]
      
      names(expSA.cur)[k] <- Da
      
    }
    
    expSA[j] <- sum(expSA.cur)
  }
  
  #### SUMMARIZE AND SAVE ####
  #Dataframe for raw proportions
  raw.props <- as.data.frame(matrix(NA,
                                    nrow=200,
                                    ncol=2,
                                    byrow = F))
  
  colnames(raw.props) <- c("proportion","category")
  
  raw.props[1:100,] <- cbind(obspropSAF,
                             rep("Observed",100))
  
  raw.props[101:200,] <- cbind(expSA,
                             rep("Null",100))
  
  
  #Dataframe for results
  hpd.intervals <- as.data.frame(matrix(NA,
                                nrow=4,
                                ncol=3,
                                byrow = F))
  
  colnames(hpd.intervals) <- c("x",
                       "y",
                       "category")
  hpd.intervals$category <- c(rep("Observed",2),
                      rep("Null",2))
  hpd.intervals$y <- c(rep(-25,2),
                       rep(-50,2))
  
  hpd.intervals$x <- c(HPDinterval(as.mcmc(obspropSAF)),
                       HPDinterval(as.mcmc(expSA)))
  
  write.csv(hpd.intervals,
            paste0("../outputs/subtrees/HPD_",clades[i],"_intervals.csv"),
            quote=F,
            row.names=T)
  
  write.csv(raw.props,
            paste0("../outputs/subtrees/proportions_",clades[i],"_raw.csv"),
            quote=F,
            row.names=T)
}





