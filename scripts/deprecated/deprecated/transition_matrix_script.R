# This file reads an empty transition matrix built from mammalian karyoptype
# data and builds/outputs a filled matrix from which a Q matrix can be built 
# (see Q_mat_script)

#Load matrix
matrix <- read.csv("../data/transition_matrix/matrix_mammalian_empty.csv",header=F)

#matrix first element = 0
matrix[1,1] <- 0

#Rownames and colnames
rownames(matrix) <- matrix[,1]

colnames(matrix) <- matrix[1,]

#Loop through rows and columns to fill
for(i in 2:nrow(matrix)){
  for(j in 2:ncol(matrix)){
    
    #AA Fisions
    if(matrix[1,j]==matrix[i,1]+1){
      matrix[i,j] <- 1
    
      #AA Fusions  
    } else if(matrix[1,j]+1==matrix[i,1]){
      matrix[i,j] <- 2
      
      #SA Fusions
    } else if(matrix[1,j]==matrix[i,1]+49){
      matrix[i,j] <- 3
      
      #Neo XY -> XY
    } else if(matrix[1,j]==matrix[i,1]-50){
      matrix[i,j] <- 4
    } else {
      matrix[i,j] <- 0
    }
  } 
}

#Remove first column and row
matrix <- matrix[-1,-1]

#Save file
write.table(matrix,"../data/transition_matrix/matrix_mammalian.csv",sep=',',col.names = F,row.names = F)






