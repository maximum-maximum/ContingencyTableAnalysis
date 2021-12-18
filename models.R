##### delete existing objects #####
rm(list = ls(all.names = TRUE))



##### data samples #####
## Yamamoto-Tomizawa2010 Table1
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0) 

## Yamamoto-Tomizawa2010 Table2
freq2 <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1) 

## Tominaga1979 (Tahata2016)
freq3 <- c(374, 602, 170, 64, 18, 255, 139, 71, 4, 23, 42, 55, 2, 6, 17, 53) 

## Smith2006, cross-classification of GSS (Tahata-Sudo-Arimoto)
freq4 <- c(98, 150, 135, 53, 37, 131, 133, 43, 9, 16, 33, 15, 4, 1, 4, 21)



##### initailize global objects #####
modelResults <- list()
inputData <- c()
r <- 0



model <- function(freq) {
  r <<- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), stop(), sqrt(length(freq)))
  inputData <<- freq

  
  ##### define design matrices #####
  array1 <- array(0, dim=c(r^2, (r-1)))
  k <- 1
  for (i in 1:r) {
    for (j in 1:r) {
      if (i <= (r-1)) array1[k, i] <- array1[k, i] + 1
      if (j <= (r-1)) array1[k, j] <- array1[k, j] + 1
      k <- k + 1
    }
  }
  
  
  array2 <- c(1:r %x% 1:r)
  
  
  array2star <- array2
  for (i in 1:r) array2star[i+r*(i-1)] <- 0
  
  
  array3 <- array(0, dim=c(r,r,r-1))
  for (k in 1:(r-1)) {
    for (i in 1:r) {
      for (j in 1:r) {
        array3[i, j, k] <- i^k - j^k
      }
    }
  }
  f <- list()
  for (k in 1:(r-1)) {
    f[[k]] <- c(aperm(array3[,,k]))
  }
  
  
  array4 <- array(0, dim=c(r^2, r))
  for (i in 1:r) {
    for (j in i:r) {
      array4[j+r*(j-1), i] <- 1
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
  array_s <- array(0, dim=c(1,r,r))
  s <- c()
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
        array_s[1, i, j] <- 1
        s <- cbind(s, c(array_s[1,,]))
      } else if (i < j) {
        array_s[1, i, j] <- array_s[1, j, i] <- 1
        s <- cbind(s, c(array_s[1,,]))
      }
      array_s <- array(0, dim=c(1,r,r))
    }
  }
  
  
  ### LSIk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsi', i), cbind(array1, ff))  
  }
  
  
  ### LSUk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsu', i), cbind(array1, ff, array2)) 
  }
  
  
  ### LSQIk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsqi', i), cbind(array1, ff, array4))  
  }
  
  
  ### LSQUk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0('array_lsqu', i), cbind(array1, ff, array4, array2star)) 
  }
  
  
  ### MEk
  library(Rsolnp)
  constraintFunc = function(p){
    const <- c()
    const <- append(const, sum(p)-1)
    l <- ifelse(length(solnpList)>=(r-1),length(solnpList_GGMk) + 1,length(solnpList) + 1)
    for(k in 1:l) const <- append(const,sum(f[[k]]*p))
    return(const)
  }
  equationValue <- list()
  for(i in 1:(r-1)) equationValue <- append(equationValue, list(rep(0,i+1)))
  removeZero <- function(freq) return (freq[freq>0])
  objectFunc <- function(p) return(-sum(freq*log(p)))
  fullModel = function(freq) return(-sum(removeZero(freq)*log(removeZero(freq/sum(freq)))))
  paramLowerBound <- rep(0, length(freq))
  p0 <- rep(1/length(freq), length(freq))
  solnpList <- list()
  for(i in 1:(r-1)){
    solnp <- solnp(p0, fun=objectFunc, eqfun=constraintFunc, eqB=equationValue[[i]], LB=paramLowerBound)
    solnpList <- append(solnpList, list(solnp))
  }
  MEkG2 <- MEkAIC <- c()
  for(i in 1:(r-1)){
    MEkG2 <- append(MEkG2, -2*(fullModel(freq) - solnpList[[i]]$value[length(solnpList[[i]]$value)]))
    MEkAIC <- append(MEkAIC, -2*(fullModel(freq) - solnpList[[i]]$value[length(solnpList[[i]]$value)]))
  }
  
  
  #objectFunc = function(p) return(-sum(freq*log(p)))
  # constraintFunc = function(p){
  #   MEkConstraintsList <- c()
  #   MEkConstraintsList <- append(MEkConstraintsList, sum(p)-1)
  #   l <- r-1
  #   constraintsList <- list()
  #   for (k in 1:l) {
  #     tmp <- c()
  #     for (i in 1:r) {
  #       for (j in 1:r) {
  #         tmp <- c(tmp, i^k-j^k)
  #       }
  #     }
  #     constraintsList[[k]] <- tmp
  #     MEkConstraintsList <- append(MEkConstraintsList, sum(tmp*p))
  #   }
  #   return(MEkConstraintsList)
  # }
  
  
  ##### show results #####
  m <- list()
  m <- append(m, list(SI=glm(freq~array_si, family=poisson, data=list(freq))))
  m <- append(m, list(SU=glm(freq~array_su, family=poisson, data=list(freq))))
  m <- append(m, list(SQI=glm(freq~array_sqi, family=poisson, data=list(freq))))
  m <- append(m, list(SQU=glm(freq~array_squ, family=poisson, data=list(freq))))
  m <- append(m, list(S=glm(freq~s, family=poisson, data=list(freq))))

  for (i in 1:(r-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsi', i)), family=poisson, data=list(freq))))
    names(m)[length(m)] <- paste0('LSI', i)
  }
  
  for (i in 1:(r-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsu', i)), family=poisson, data=list(freq))))
    names(m)[length(m)] <- paste0('LSU', i)  
  }
  
  for (i in 1:(r-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsqi', i)), family=poisson, data=list(freq))))
    names(m)[length(m)] <- paste0('LSQI', i)
  }

  for (i in 1:(r-1)) {
    m <- append(m, list(glm(as.formula(paste0('freq~array_lsqu', i)), family=poisson, data=list(freq))))
    names(m)[length(m)] <- paste0('LSQU', i)  
  }
  
  df <- G2 <- AIC <- Pvalue <- code <- c()
  for (i in m) {
    p <- round(1 - pchisq(i$deviance, i$df.residual), digits=4)
    signif.code <- ''
    for (alpha in c(0.05, 0.01, 0.001)) {
      if (p < alpha) signif.code <- paste0(signif.code, "*")   
    }
    
    df <- append(df, i$df.residual)
    G2 <- append(G2, round(i$deviance, digits=3))
    AIC <- append(AIC, round(i$aic, digits=3))
    Pvalue <- append(Pvalue, p)
    code <- append(code, signif.code)
  }
  result <- data.frame(model=names(m), df=df, G2=G2, AIC=AIC, Pvalue=Pvalue, code=code)
  for (i in 1:(r-1)) {
    result <- rbind(result, list(paste0('ME',i), i, round(MEkG2[i], digits=3), '', '', ''))
  }  
  modelResults <<- m
  cat("\n")
  print(result)
  cat("-----\n")
  cat("Signif. codes:  0  '***'  0.001  '**'  0.01  '*'  0.05  ''\n")
}

detail <- function(model) {
  selectedModelResult <- modelResults[[model]]
  fittingValue <- round(fitted(selectedModelResult), 3)
  resultMatrix <- t(matrix(paste0(inputData,' (',fittingValue,')'),r,r))
  
  print(summary(selectedModelResult))
  cat('Data:\n')
  print(resultMatrix)
  cat('(The parenthesized values are the MLEs of expected frequencies under the selected model)')
}