library(phytools)
library(abc)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]
data.sim <- read.csv("../data/chromes/simChrom_data.csv")
tree <- read.tree("../data/trees/4705sp_mean.nwk")

#### DATA ORGANIZATION ####

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

data <- data.unique <- data.alt

rm(data.alt,data.unique)

#### CALCULATE EMPIRICAL SUMMARY STATS ####

#Calculate margin of error
margin <- qt(0.975,df=914)*sd(data$hapauto)/sqrt(914)

#Fill vector with stats
stats.emp <- c(mean(data$hapauto),
               sd(data$hapauto)/mean(data$hapauto),
               mean(data$hapauto) - margin,
               mean(data$hapauto) + margin)

#Name vector (should be same as in data.sim)
colnames(data.sim)

names(stats.emp) <- c("hapauto.mean",
                      "hapauto.cv",
                      "hapauto.95..lower",
                      "hapauto.95..upper")

#### RUN ABC ####

#Split into params and sumstats
params.sim <- data.sim[,c(1,2,4)]

stats.sim <- data.sim[,5:8]

#Run ABC
abc.out <- abc(target = stats.emp,
               param = params.sim,
               sumstat = stats.sim,
               tol=0.1,
               method = "rejection")

summary(abc.out)

params
