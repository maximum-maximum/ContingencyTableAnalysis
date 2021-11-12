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
array_si <- array(0, dim=c(NI^2, (NI-1)))

k <- 1
for (i in 1:NI) {
  for (j in 1:(NI-1)) {
    if (i <= (NI-1)) {
      array_si[k, i] <- array_si[k, i] + 1  
    }
    if (j <= (NI-1)) {
      array_si[k, j] <- array_si[k, j] + 1
    }
    k <- k + 1
  }
}


### SU ###
array_su <- array(0, dim=c(NI^2, (NI-1)))

l <- 1
for (i in 1:NI) {
  for (j in 1:(NI-1)) {
    if (i <= (NI-1)) {
      array_su[l, i] <- array_su[l, i] + 1  
    }
    if (j <= (NI-1)) {
      array_su[l, j] <- array_su[l, j] + 1
    }
    l <- l + 1
  }
}

array_su2 <- c(1:NI %x% 1:NI)
