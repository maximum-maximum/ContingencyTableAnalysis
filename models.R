##### delete existing objects #####
rm(list = ls())



##### data samples #####
## Yamamoto-Tomizawa2010 Table1
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0) 

## Yamamoto-Tomizawa2010 Table2
freq2 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1) 

## Tominaga1979 (Tahata2016)
freq3 <- c(374, 602, 170, 64, 18, 255, 139, 71, 4, 23, 42, 55, 2, 6, 17, 53) 

## Smith2006, cross-classification of GSS (Tahata-Sudo-Arimoto)
freq4 <- c(98, 150, 135, 53, 37, 131, 133, 43, 9, 16, 33, 15, 4, 1, 4, 21)



model <- function(freq) {
  NI <- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), stop(), sqrt(length(freq)))
  row <- gl(NI, NI, length=NI^2)
  col <- gl(NI, 1, length=NI^2)
  sample <- data.frame(freq, row, col)


  
  ##### define design matrices #####
  array1 <- array(0, dim=c(NI^2, (NI-1)))
  k <- 1
  for (i in 1:NI) {
    for (j in 1:NI) {
      if (i <= (NI-1)) array1[k, i] <- array1[k, i] + 1
      if (j <= (NI-1)) array1[k, j] <- array1[k, j] + 1
      k <- k + 1
    }
  }
  
  
  array2 <- c(1:NI %x% 1:NI)
  
  
  array2star <- array2
  for (i in 1:NI) array2star[i+NI*(i-1)] <- 0
  
  
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


  ### SQI
  array_sqi <- cbind(array1, array4)
  
  
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
  
  
  ### LSIk
  ff <- c()
  for (i in 1:(NI-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsi', i), cbind(array1, ff))  
  }
  
  
  ### LSUk
  ff <- c()
  for (i in 1:(NI-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsu', i), cbind(array1, ff, array2)) 
  }
  
  
  ### LSQIk
  ff <- c()
  for (i in 1:(NI-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsqi', i), cbind(array1, ff, array4))  
  }
  
  
  ### LSQUk
  ff <- c()
  for (i in 1:(NI-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsqu', i), cbind(array1, ff, array4, array2star)) 
  }
  
  
  ### MEk
  library(Rsolnp)
  ConstFunc_MEk = function(p){
    const <- c()
    sum <- sum(p)
    const <- append(const,sum-1)
    l <- ifelse(length(solution_MEk)>=(NI-1),length(solution_GGMk) + 1,length(solution_MEk) + 1)
    for(k in 1:l) const <- append(const,sum((f[[k]])*p))
    return(const)
  }
  eq.valueMEk <- list()
  for(i in 1:(NI-1)) eq.valueMEk <- append(eq.valueMEk,list(rep(0,i+1)))
  zero <- function(x) return (x[x>0])
  ObjFunc = function(p) return(-sum(freq*log(p)))
  saturated_model = function(freq) return(-sum(zero(freq)*log(zero(freq/sum(freq)))))
  eq.LB <- rep(0,length(freq))
  p0 <- rep(1/length(freq), length(freq))
  solution_MEk <- list()
  for(i in 1:(NI-1)){
    solution_MEk <- append(solution_MEk, list(solnp(p0, fun=ObjFunc, eqfun=ConstFunc_MEk, eqB=eq.valueMEk[[i]], LB=eq.LB)))
  }
  ans_MEk <- c()
  for(i in 1:(NI-1)){
    ans_MEk <- append(ans_MEk, -2*(saturated_model(freq) - solution_MEk[[i]]$value[length(solution_MEk[[i]]$value)]))
  }
  
  
  
  ##### show results #####
  m <- list()
  m <- append(m, list(SI=glm(freq~array_si, family=poisson, data=sample)))
  m <- append(m, list(SU=glm(freq~array_su, family=poisson, data=sample)))
  m <- append(m, list(SQI=glm(freq~array_sqi, family=poisson, data=sample)))
  m <- append(m, list(SQU=glm(freq~array_squ, family=poisson, data=sample)))
  m <- append(m, list(S=glm(freq~s, family=poisson, data=sample)))

  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsi', i)), family=poisson, data=sample)))
    names(m)[length(m)] <- paste0('LSI', i)
  }
  
  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsu', i)), family=poisson, data=sample)))
    names(m)[length(m)] <- paste0('LSU', i)  
  }
  
  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsqi', i)), family=poisson, data=sample)))
    names(m)[length(m)] <- paste0('LSQI', i)
  }

  for (i in 1:(NI-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsqu', i)), family=poisson, data=sample)))
    names(m)[length(m)] <- paste0('LSQU', i)  
  }
  
  df <- G2 <- AIC <- Pvalue <- code <- c()
  signif.code <- ''
  for (i in m) {
    df <- append(df, i$df.residual)
    G2 <- append(G2, round(i$deviance, digits=3))
    AIC <- append(AIC, round(i$aic, digits=3))
    Pvalue <- append(Pvalue, round(1-pchisq(i$deviance, i$df.residual), digits=4))
    p <- round(1-pchisq(i$deviance, i$df.residual), digits=4)
    if (p < 0.05) signif.code <- paste0(signif.code, "*") 
    if (p < 0.01) signif.code <- paste0(signif.code, "*") 
    if (p < 0.001) signif.code <- paste0(signif.code, "*") 
    code <- append(code, signif.code)
  }
  result <- data.frame(model=names(m), df=df, G2=G2, AIC=AIC, Pvalue=Pvalue, code=code)
  for (i in 1:(NI-1)) result <- rbind(result, list(paste0('ME',i), i, round(ans_MEk[i], digits=3), '', '', ''))
  return (result)
}