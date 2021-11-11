#freq <- c(374,602,170,64,18,255,139,71,4,23,42,55,2,6,17,53)
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0)
#freq <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1)

NI <- ifelse(floor(sqrt(length(freq)))
             <ceiling(sqrt(length(freq))),
             stop(),sqrt(length(freq)))
row <- gl(NI,NI,length=NI^2)
col <- gl(NI,1,length=NI^2)
u <- c(1:NI)
sample <- data.frame(freq,row,col)

### SI ###
array_cs <- array (0, dim =c(NI^2, (NI-1)))
index_i <- index_j <- 0

index_i <- index_i + 1
for (i in 1:NI^2) {
  index_j <- index_j + 1
  for (j in 1:(NI-1)) {
    
    
    print(i)
    #print(index_i)
    print("--------")
    print(j)
    #print(index_j)
    print("--------")
    
    if (i==index_i) {
      array_cs[i, j] <- array_cs[i, j] + 1
    } 
    if (j==index_j) {
      array_cs[i, j] <- array_cs[i, j] + 1
    } 
  }
}
### SU ###