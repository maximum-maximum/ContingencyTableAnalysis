library(Rsolnp)

freq <- c(374,602,170,64,18,255,139,71,4,23,42,55,2,6,17,53) # Tahata2016
freq2 <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0) # Yamamoto, Tomizawa 2010 Table1
freq3 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1) # Yamamoto, Tomizawa 2010 Table2


model = function(freq) {
  NI <- ifelse(floor(sqrt(length(freq)))
              <ceiling(sqrt(length(freq))),
              stop(),sqrt(length(freq)))
  row <- gl(NI,NI,length=NI^2)
  col <- gl(NI,1,length=NI^2)
  u <- c(1:NI)
  sample <- data.frame(freq,row,col)


  
  ##### SI #####
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
  
  

  ##### SU #####
  theta <- c(1:NI %x% 1:NI)
  array_su <- cbind(array_si, theta)


  
  ##### LSQUk ver.1(別表現) #####
  ## Defines xl (l = 1, . . . , r −1), 
  ## and f[[l]] is the vector corresponding to xl
  array_f <- array(0, dim=c(NI,NI,NI-1))
  for (k in 1:(NI-1)) {
    for (i in 1:NI) {
      for (j in 1:NI) {
        array_f[i, j, k] <- i^k - j^k
      }
    }
  }
  f <- list()
  for (k in 1:(NI-1)) {
    f[[k]] <- c(aperm(array_f[,,k]))
  }

  theta_lsquk <- theta
  psi <- c()
  k <- 1
  for (i in 1:NI){
    for (j in 1:NI){
      if (i == j) {
        theta_lsquk[k] <- 0
        psi[k] <- 1
      } else {
        psi[k] <- 0
      }
      k <- k + 1
    }
  }

  
  
  ##### LSQUk ver.2 #####
  array_cs <- array(0, dim=c(NI,NI))
  for (i in 1:NI) {
    for (j in 1:NI) {
      if (i < j) {
        array_cs[i, j] <- 1
      }
    }
  }
  cs <- c(aperm(array_cs))
  
  ## Defines xl (l = 1, . . . , r −1), 
  ## and f[[l]] is the vector corresponding to xl
  array_f2 <- array(0, dim=c(NI,NI,NI-1))
  for (k in 1:(NI-1)) {
    for (i in 1:NI) {
      for (j in 1:NI) {
        array_f2[i, j, k] <- u[j]^k - u[i]^k
      }
    }
  }
  f2 <- list()
  for (k in 1:(NI-1)) {
    f2[[k]] <- c(aperm(array_f[,,k]))
  }
  
  
  
  ### Answer ###
  m <- list ()
  m <- append(m, list(SI = glm(freq~array_si, family=poisson, data=sample)))
  m <- append(m, list(SU = glm(freq~array_su, family=poisson, data=sample)))
  
  m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[1]]+theta_lsquk+psi, family=poisson, data=sample)))
  m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[2]]+theta_lsquk+psi, family=poisson, data=sample)))
  m <- append(m, list(LSQUk_ver1 = glm(freq~array_si+f[[3]]+theta_lsquk+psi, family=poisson, data=sample)))
  
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
  
  return (m)
}

model(freq)