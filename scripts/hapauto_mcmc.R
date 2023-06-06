#### LOAD PACKAGES ####
library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../data/chromes/dat.csv",
                 as.is=T)[,c(1,3)]
tree <- force.ultrametric(read.tree("../data/trees/tree.nex"))

#### MODEL ####

#Convert to data matrix
data.matrix <- datatoMatrix(dat,
                            c(2,50),
                            hyper = F)

#Make mkn model
model <- make.mkn(tree,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#constrain model
model.con <- constrainMkn(data.matrix,
                          model,
                          hyper = F,
                          polyploidy = F,
                          verbose = T,
                          constrain = list(drop.poly=T,
                                           drop.demi=T))

#run MCMC
model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                x.init=c(1,1),
                                prior=make.prior.exponential(r=0.5),
                                #upper=c(100,100,100,100),
                                nsteps = 500,
                                w=1)

#### BUILD QMATRIX ####

#Extract post burn portion
model.mcmc.postburn <- model.mcmc[450:500,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1))

names(params) <- colnames(model.mcmc[,2:3])

#Sub into matrix
parMat <- model.con$`parameter matrix`
mat <- parMat
parMat[parMat == "asc1"] <- params[1]
parMat[parMat == "desc1"] <- params[2]
mat[mat == "asc1"] <- 1
mat[mat == "desc1"] <- 2

#convert to numeric
parMat <- as.data.frame(parMat)
mat <- as.data.frame(mat)
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)
mat <- sapply(mat[,1:ncol(mat)],as.numeric)

#calculate diagonal
diag(parMat) <- -rowSums(parMat)

#### SAVE QMATRIX ####
write.csv(parMat,
          paste0("../data/transition_matrix/Q_matrix_hapauto.csv"),
          row.names=F,quote=F)
write.csv(mat,
          paste0("../data/transition_matrix/transition_matrix_hapauto.csv"),
          row.names=F,quote=F)
