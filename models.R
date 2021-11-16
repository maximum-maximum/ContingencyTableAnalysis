library(Rsolnp)

##### data samples #####
## Yamamoto-Tomizawa2010 Table1
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0) 

## Yamamoto-Tomizawa2010 Table2
freq2 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1) 

## Tominaga1979 (Tahata2016, Tahata-Sudo-Arimoto)
freq3 <- c(374, 602, 170, 64, 18, 255, 139, 71, 4, 23, 42, 55, 2, 6, 17, 53) 



model = function(freq) {
  NI <- ifelse(floor(sqrt(length(freq)))
              <ceiling(sqrt(length(freq))),
              stop(),sqrt(length(freq)))
  row <- gl(NI,NI,length=NI^2)
  col <- gl(NI,1,length=NI^2)
  u <- c(1:NI)
  sample <- data.frame(freq, row, col)


  
  ##### define design matrices #####
  array1 <- array(0, dim=c(NI^2, (NI-1)))
  k <- 1
  for (i in 1:NI) {
    for (j in 1:(NI-1)) {
      if (i <= (NI-1)) {
        array1[k, i] <- array1[k, i] + 1
      }
      if (j <= (NI-1)) {
        array1[k, j] <- array1[k, j] + 1
      }
      k <- k + 1
    }
  }
  
  array2 <- c(1:NI %x% 1:NI)
  
  array2star <- array2
  for (i in 1:NI) {
    array2star[i+NI*(i-1)] <- 0
  }
  
  array3 <- array(0, dim=c(NI,NI,NI-1))
  for (k in 1:(NI-1)) {
    for (i in 1:NI) {
      for (j in 1:NI) {
        array3[i, j, k] <- i^k - j^k
      }
    }
  }
  f <- list()
  for (k in 1:(NI-1)) {
    f[[k]] <- c(aperm(array3[,,k]))
  }
  
  array4 <- array(0, dim=c(NI^2, NI))
  for (i in 1:NI) {
    for (j in i:NI) {
      array4[j+NI*(j-1), i] <- 1
      break
    }
  }
  
  
  
  ##### bind matrices #####
  ### SI 
  array_si <- array1
  
  
  ### SU 
  array_su <- cbind(array1, array2)


  ### LSQIk
  array_lsqi1 <- cbind(array1, f[[1]], array4)
  array_lsqi2 <- cbind(array1, f[[2]], array4)
  array_lsqi3 <- cbind(array1, f[[3]], array4)
  
  
  ### LSQUk
  array_lsqu1 <- cbind(array1, f[[1]], array4, array2star)
  array_lsqu2 <- cbind(array1, f[[2]], array4, array2star)
  array_lsqu3 <- cbind(array1, f[[3]], array4, array2star)


  ### SQU
  array_squ <- cbind(array1, array4, array2star)


  ### S
  array_s <- array(0,dim=c(1,NI,NI))
  s <- c()
  for(i in 1:NI){
    for(j in 1:NI){
      if(i==j){
        array_s[1,i,j]=1
        s <-cbind(s,c(array_s[1,,]))
      }else if(i<j){
        array_s[1,i,j]=array_s[1,j,i]=1
        s <- cbind(s,c(array_s[1,,]))
      }
      array_s <- array(0,dim=c(1,NI,NI))
    }
  }
  

  
  ##### result #####
  m <- list()
  m <- append(m, list(SI = glm(freq~array_si, family=poisson, data=sample)))
  m <- append(m, list(SU = glm(freq~array_su, family=poisson, data=sample)))
  m <- append(m, list(S = glm(freq~s, family=poisson, data=sample)))
  m <- append(m, list(LSQI1 = glm(freq~array_lsqi1, family=poisson, data=sample)))
  m <- append(m, list(LSQI2 = glm(freq~array_lsqi2, family=poisson, data=sample)))
  m <- append(m, list(LSQI3 = glm(freq~array_lsqi3, family=poisson, data=sample)))
  m <- append(m, list(LSQU1 = glm(freq~array_lsqu1, family=poisson, data=sample)))
  m <- append(m, list(LSQU2 = glm(freq~array_lsqu2, family=poisson, data=sample)))
  m <- append(m, list(LSQU3 = glm(freq~array_lsqu3, family=poisson, data=sample)))
  m <- append(m, list(SQU = glm(freq~array_squ, family=poisson, data=sample)))
  
  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[1]]+theta_lsquk+array_psi, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[2]]+theta_lsquk+array_psi, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[3]]+theta_lsquk+array_psi, family=poisson, data=sample)))
  
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[1]]+theta, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[2]]+theta, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[3]]+theta, family=poisson, data=sample)))

  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[1]]+theta_lsquk+psi, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[1]]+theta, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[2]]+theta_lsquk+psi, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[2]]+theta, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[3]]+theta_lsquk+psi, family=poisson, data=sample)))
  # m <- append(m, list(LSQUk_ver2 = glm(freq~cs+f2[[3]]+theta, family=poisson, data=sample)))
  
  #m <- append(m, list(LSQIk = glm(freq~array_si+f[[1]]+psi, family=poisson, data=sample)))
  #m <- append(m, list(LSQIk = glm(freq~array_si+f[[2]]+psi, family=poisson, data=sample)))
  #m <- append(m, list(LSQIk = glm(freq~array_si+f[[3]]+psi, family=poisson, data=sample)))
  
  #m <- append (m, list (SItheta = glm (freq~array_si+theta, family=poisson, data=sample)))
  
  for (i in 1:length(m)) {
    print(m[[i]]$deviance)
  }
  return (m)
}