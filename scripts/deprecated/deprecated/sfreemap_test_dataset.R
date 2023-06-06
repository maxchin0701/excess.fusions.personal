# Functions from the package sfreemap are tested here on example datasets
# provided

# Focus was placed on confirming correct function given the provision of a 
# pre-determined Q matrix and testing functionality given Q matrix states
# which are not represented in tip states

library(phytools)
library(sfreemap)
library(parallel)

library(devtools)

#Load tips and trees
tips <- sfreemap.corals.tips

trees <- sfreemap.corals.trees

#Perform analysis
map <- sfreemap(trees[[1]],
                tips,
                prior="equal",
                parallel=F)

#Output: As intended

#Save built Q matrix
Qmat.ex <- map$Q 

#Rerun analysis with provided Q matrix
map <- sfreemap(trees,
                tips,
                Qmat.ex,
                parallel=F)

#Output: As intended

#Add third state (All of these are originally "colonial")
tips[c(10,20,30,40,50,60)] <- "test"

#Run analysis
map <- sfreemap(trees,
                tips,
                parallel=F)

map.3 <- sfreemap::sfreemap(trees,
                tips,
                parallel=F)

#Output: as intended

Qmat.ex.alt <- map[[1]]$Q

#Set third state back to original
tips[c(10,20,30,40,50,60)] <- "colonial"

#Run analysis with "mismatched" matrix and tip states
map.2 <- sfreemap(trees,
                tips,
                Qmat.ex.alt,
                parallel=F,
                poss.states=c("colonial","solitary","test"))

#Error: tb[, , i] %*% liks[des1, ] : non-conformable arguments

#Plot
plot_distribution_tree(map_posterior_distribution(map[[1]],
                                                  map,
                                                  parallel=F),
                       state = "solitary",
                       tip_states = tips)

plot_distribution_tree(map_posterior_distribution(map.3[[1]],
                                                  map.3,
                                                  parallel=F),
                       state = "solitary",
                       tip_states = tips)


#### SAVE OUTPUTS ####

#Output Q matrix
out.Q <- as.data.frame(map$Q)

#Output priors
out.pi <- as.data.frame(map$prior)
out.pi <- as.data.frame(cbind(rownames(out.pi),out.pi[,1]))

colnames(out.pi) <- c("state","priors")
rownames(out.pi) <- NULL

#Output logL
out.logL <- as.data.frame(map$logL)
colnames(out.logL) <- "logL"

#Output mapped dwelling times
out.times <- as.data.frame(map$mapped.edge)
out.times <- as.data.frame(cbind(rownames(out.times),out.times))

colnames(out.times)[1] <- "edge"
rownames(out.times) <- NULL

#Output mapped transitions
out.trans <- as.data.frame(map$mapped.edge.lmt)
out.trans <- as.data.frame(cbind(rownames(out.trans),out.trans))

colnames(out.trans)[1] <- "edge"
rownames(out.trans) <- NULL



#Save all
write.csv(out.Q,"Q",row.names = F,quote=F)
write.csv(out.pi,"priors",row.names = F,quote=F)
write.csv(out.logL,"logL",row.names = F,quote=F)
write.csv(out.times,"times",row.names = F,quote=F)
write.csv(out.trans,"trans",row.names = F,quote=F)





