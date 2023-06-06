#This script adds a column to the working mammalian chromosome dataset which
#codifies the simulation states associated with taxa

#Read chromosome data
data <- read.csv("../data/chromes/mammal_chroms_working.csv",
                 as.is=T)#[,c(6,12,13)]

#Remove first row
data <- data[,-1]

#Convert codedSCS to factor
data$codedSCS <- as.factor(data$codedSCS)

#Determine sim.state
for(i in 1:nrow(data)){
  if(is.na(data$hapauto[i]) == T ||
     data$codedSCS[i] == "unclear"){
    data$sim.state[i] <- NA
  } else if(data$codedSCS[i] == "unfused"){
    data$sim.state[i] <- data$hapauto[i] - 1 
  } else if(data$codedSCS[i] == "fused"){
    data$sim.state[i] <- data$hapauto[i] + 48
  }
}

#Save data
write.csv(data,"../data/chromes/mammal_chroms_sim_state.csv",row.names = F)


