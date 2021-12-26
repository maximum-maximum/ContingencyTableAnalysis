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
globalAnalysResults <- list()
inputData <- c()
r <- 0



model <- function(freq, sort=FALSE) {
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
  arraySI <- array1
  
  
  ### SU 
  arraySU <- cbind(array1, array2)


  ### SQI
  arraySQI <- cbind(array1, array4)
  
  
  ### SQU
  arraySQU <- cbind(array1, array4, array2star)
  
  
  ### S
  s <- array(0, dim=c(1,r,r))
  arrayS <- c()
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
        s[1, i, j] <- 1
        arrayS <- cbind(arrayS, c(s[1,,]))
      } else if (i < j) {
        s[1, i, j] <- s[1, j, i] <- 1
        arrayS <- cbind(arrayS, c(s[1,,]))
      }
      s <- array(0, dim=c(1,r,r))
    }
  }
  
  
  ### LSIk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0("arrayLSI", i), cbind(array1, ff))
  }
  
  
  ### LSUk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0("arrayLSU", i), cbind(array1, ff, array2))
  }
  
  
  ### LSQIk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0("arrayLSQI", i), cbind(array1, ff, array4))
  }
  
  
  ### LSQUk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0("arrayLSQU", i), cbind(array1, ff, array4, array2star))
  }


  ### LSk
  ff <- c()
  for (i in 1:(r-1)) {
    ff <- cbind(ff, f[[i]])
    assign(paste0("arrayLS", i), cbind(arrayS, ff))
  }
  

  ##### analyze with each model #####
  analysResults <- list()
  models <- c("SI", "SU", "SQI", "SQU", "S", "LSI", "LSU", "LSQI", "LSQU", "LS")
  for (modelIndex in 1:length(models)) {
    if (modelIndex <= which(models == "S")) {
      formula <- as.formula(paste0("freq~array", models[modelIndex]))
      glm <- glm(formula, family=poisson, data=list(freq))
      analysResults <- append(analysResults, list(glm))
      names(analysResults)[length(analysResults)] <- models[modelIndex]
    } else {
      for (i in 1:(r-1)) {
        formula <- as.formula(paste0("freq~array", models[modelIndex], i))
        glm <- glm(formula, family=poisson, data=list(freq))
        analysResults <- append(analysResults, list(glm))
        names(analysResults)[length(analysResults)] <- paste0(models[modelIndex], i)
      }
    }
  }
  
  
  ### MEk
  library(Rsolnp)
  #source("solnp.R")
  constraintFunc <- function(p){
    MEkConstraints <- c()
    MEkConstraints <- append(MEkConstraints, sum(p)-1)
    
    l <- length(solnpList) + 1
    for (k in 1:l) {
      momentConstraints <- c()
      for (i in 1:r) {
        for (j in 1:r) {
          momentConstraints <- c(momentConstraints, i^k-j^k)
        }
      }
      MEkConstraints <- append(MEkConstraints, sum(momentConstraints*p))
    }
    return(MEkConstraints)
  }
  equationValue <- list()
  for (i in 1:(r-1)) equationValue <- append(equationValue, list(rep(0, i+1)))
  removeZero <- function(freq) return(freq[freq > 0])
  objectFunc <- function(p) return(-sum(freq*log(p)))
  fullModel <- function(freq) return(-sum(removeZero(freq)*log(removeZero(freq/sum(freq)))))
  paramLowerBound <- rep(0, length(freq))
  p0 <- rep(1/length(freq), length(freq))
  solnpList <- list()
  for (i in 1:(r-1)) {
    solnp <- solnp(p0, fun=objectFunc, eqfun=constraintFunc, eqB=equationValue[[i]], LB=paramLowerBound)
    solnpList <- append(solnpList, list(solnp))
  }
  
  constMolecule <- 0
  constDenominator <- 0
  for (i in 1:sum(freq)) constMolecule <- constMolecule + log(i)
  for (i in removeZero(freq)) {
    for (j in 1:i) {
      constDenominator <- constDenominator + log(j)
    }
  }
  
  for (i in (r-1):1) {
    G2 <- -2*(fullModel(freq) - solnpList[[i]]$value[length(solnpList[[i]]$value)])
    maxLogLikeli <- constMolecule - constDenominator + (-solnpList[[i]]$value[length(solnpList[[i]]$value)])
    paramSize <- r^2 - 1 - i
    AIC <- -2*maxLogLikeli + 2*paramSize
    
    fittingValue <- round(sum(freq)*(solnpList[[i]]$pars), 3)
    resultMatrix <- t(matrix(paste0(freq," (",fittingValue,")"), r, r))
    
    analysResults <- append(analysResults, list(list(deviance=G2, df.residual=i, aic=AIC, result=resultMatrix)))
    names(analysResults)[length(analysResults)] <- paste0("ME", i)
    
    # fullMaxLogLikeli <- constMolecule - constDenominator + (-fullModel(freq))
    # AIC2 <- -2*fullMaxLogLikeli + G2 +2*paramSize
    # print(AIC)
    # print(AIC2)
    # cat("\n")
  }
  
  
  
  ##### show results #####
  anothNames <- dfs <- G2s <- AICs <- pValues <- codes <- c()
  anothNameTargetModels <- paste0(c("LSI", "LSU", "LSQI", "LSQU", "LS", "ME"), r-1)
  anothNameModels <- paste0("(", c("I", "U", "QI", "QU", "QS", "MH"), ")")
  for (model in analysResults) {
    modelName <- names(analysResults[match(list(model), analysResults)])
    anothName <- ""
    for (i in 1:length(anothNameTargetModels)) {
      if (modelName == anothNameTargetModels[i]) anothName <- anothNameModels[i]
    }
    p <- round(1 - pchisq(model$deviance, model$df.residual), 4)
    code <- ""
    if (0.05 < p && p < 0.1) code <- "."
    else {
      for (alpha in c(0.05, 0.01, 0.001)) {
        if (p < alpha) code <- paste0(code, "*")
      }
    }
    
    anothNames <- append(anothNames, anothName)
    dfs <- append(dfs, model$df.residual)
    G2s <- append(G2s, round(model$deviance, 3))
    AICs <- append(AICs, round(model$aic, 3))
    pValues <- append(pValues, p)
    codes <- append(codes, code)
  }
  resultForDisplay <- data.frame(model=names(analysResults), anothName=anothNames, df=dfs, G2=G2s, AIC=AICs, pValue=pValues, code=codes)
  names(resultForDisplay)[6] <- "Pr(>G2)"
  names(resultForDisplay)[2] <- names(resultForDisplay)[7] <- ""
  
  globalAnalysResults <<- analysResults
  cat("\n")
  
  if (is.character(sort)) print(resultForDisplay[order(resultForDisplay[[sort]]),])
  else print(resultForDisplay)
  
  cat("---\n")
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
}



detail <- function(model) {
  selectedModelResult <- globalAnalysResults[[model]]
  modelName <- names(globalAnalysResults[match(list(selectedModelResult), globalAnalysResults)])
  cat("Model:", modelName, "\n")
  if (length(selectedModelResult) == 30) {
    fittingValue <- round(fitted(selectedModelResult), 3)
    resultMatrix <- t(matrix(paste0(inputData," (",fittingValue,")"), r, r))
    
    print(summary(selectedModelResult))
    cat("Data:\n")
    print(resultMatrix)
  } else print(selectedModelResult)
  cat("(The parenthesized values are the MLEs of expected frequencies under the selected model)")
}