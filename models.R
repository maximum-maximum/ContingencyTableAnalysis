library(Rsolnp)

##### Data samples #####
## Yamamoto-Tomizawa2010 Table1
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0) 

## Yamamoto-Tomizawa2010 Table2
freq2 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1) 

## Tominaga1979 (Tahata2016)
freq3 <- c(374, 602, 170, 64, 18, 255, 139, 71, 4, 23, 42, 55, 2, 6, 17, 53) 

## Smith2006 (Tahata-Sudo-Arimoto)
freq4 <- c(98, 150, 135, 53, 37, 131, 133, 43, 9, 16, 33, 15, 4, 1, 4, 21)



model = function(freq) {
  NI <- ifelse(floor(sqrt(length(freq)))
              <ceiling(sqrt(length(freq))),
              stop(),sqrt(length(freq)))
  row <- gl(NI, NI, length=NI^2)
  col <- gl(NI, 1, length=NI^2)
  u <- c(1:NI)
  sample <- data.frame(freq, row, col)


  
  ##### Define design matrices #####
  array1 <- array(0, dim=c(NI^2, (NI-1)))
  k <- 1
  for (i in 1:NI) {
    for (j in 1:NI) {
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
  
  
  
  ##### Bind matrices #####
  ### SI 
  array_si <- array1
  
  
  ### SU 
  array_su <- cbind(array1, array2)


  ### SQU
  array_squ <- cbind(array1, array4, array2star)
  
  
  ### S
  array_s <- array(0, dim=c(1,NI,NI))
  s <- c()
  for (i in 1:NI) {
    for (j in 1:NI) {
      if (i == j) {
        array_s[1, i, j] <- 1
        s <- cbind(s, c(array_s[1,,]))
      } else if (i < j) {
        array_s[1, i, j] <- array_s[1, j, i] <- 1
        s <- cbind(s, c(array_s[1,,]))
      }
      array_s <- array(0, dim=c(1,NI,NI))
    }
  }
  
  
  ### LSQIk
  for (i in 1:(NI-1)) assign(paste('array_lsqi', i, sep=''), cbind(array1, f[[i]], array4)) 
  
  
  ### LSQUk
  for (i in 1:(NI-1)) assign(paste('array_lsqu', i, sep=''), cbind(array1, f[[i]], array4, array2star)) 
  
  
  
  ##### Show results #####
  m <- list()
  m <- append(m, list(SI = glm(freq~array_si, family=poisson, data=sample)))
  m <- append(m, list(SU = glm(freq~array_su, family=poisson, data=sample)))
  m <- append(m, list(SQU = glm(freq~array_squ, family=poisson, data=sample)))
  m <- append(m, list(S = glm(freq~s, family=poisson, data=sample)))

  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste('freq~array_lsqi', i, sep='')), family=poisson, data=sample)))
    names(m)[length(m)] <- paste('LSQI', i, sep='')  
  }

  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste('freq~array_lsqu', i, sep='')), family=poisson, data=sample)))
    names(m)[length(m)] <- paste('LSQU', i, sep='')  
  }
  
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
  
  df <- G2 <- c()
  for (i in 1:length(m)) {
    df <- append(df, m[[i]]$df.residual)
    G2 <- append(G2, m[[i]]$deviance)
    # print(paste(names(m)[i], m[[i]]$deviance), quote=F)
  }
  result <- data.frame(model=names(m), df=df, G2=G2)
  return (result)
}