library(Rsolnp)
#freq <- c(374,602,170,64,18,255,139,71,4,23,42,55,2,6,17,53)
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0)
freq3 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1)

# model = function(freq) {
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


  ##### LSQUk #####
  ## Defines x0, and cs is the vector corresponding to x0.
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
  array_f <- array(0, dim=c(NI,NI,NI-1))
  for (k in 1:(NI-1)) {
    for (i in 1:NI) {
      for (j in 1:NI) {
        array_f[i, j, k] <- u[j]^k - u[i]^k
      }
    }
  }
  f <- list()
  for (k in 1:(NI-1)) {
    f[[k]] <- c(aperm(array_f[,,k]))
    f
  }
  
  ## In a manner similar to Lawal (2001, 2004), 
  ## the following code defines XS, and s is the matrix corresponding to XS
  array_s <- array(0, dim=c(1,NI,NI))
  s <- c()
  for (i in 1:NI) {
    for (j in 1:NI) {
      if (i == j) {
        array_s[1, i, j] <- 1
        s <- cbind(s, c(array_s[1,,]))
      } else if (i < j){
        array_s[1, i, j] <- array_s[1, j, i] <- 1
        s <- cbind(s, c(array_s[1,,]))
      }
      array_s <- array(0, dim=c(1,NI,NI))
    }
  }


  ### Answer ###
  m <- list ()
  m <- append(m, list (SI = glm (freq~array_si, family = poisson, data = sample )))
  m <- append(m, list (SU = glm (freq~array_su, family = poisson, data = sample )))
  #m <- append (m, list (SItheta = glm (freq~array_si+theta, family = poisson, data = sample )))
  LSk <- ELSk <- list()
  str <- str2 <- c()
  
  ans_si <- ans_su <- c()
# }