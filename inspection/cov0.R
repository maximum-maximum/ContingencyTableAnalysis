library(Rsolnp)

tom = c(64,57,57,72,36,21,94,94,105,141,97,71,58,54,65,77,54,54,46,40,60,94,78,71)
tom2 <- c(122,30,20,472,226,51,66,704,306,115,96,1072,130,59,38,501,50,31,15,249)

test <- function(freq, row, col){
  print(matrix(freq, row, col, T))
  
  cov0constraintFunc <- function(p){
    cov0Constraints <- c()
    cov0Constraints <- append(cov0Constraints, sum(p)-1)
    
    pMatrix <- matrix(p, row, col, T)
    
    multi <- c()
    for (i in 1:row) {
      for (j in 1:col) {
        multi <- append(multi, rowSums(pMatrix)[i] * colSums(pMatrix)[j])
      }
    }
    cov0Constraints <- append(cov0Constraints, sum(c(1:row %x% 1:col) * (p-multi)))
    return(cov0Constraints)
  }

  removeZero <- function(freq) return(freq[freq > 0])
  objectFunc <- function(p) return(-sum(freq*log(p)))
  fullModel <- function(freq) return(-sum(removeZero(freq)*log(removeZero(freq/sum(freq)))))
  paramLowerBound <- rep(0, length(freq))
  p0 <- rep(1/length(freq), length(freq))

  solnp <- solnp(p0, fun=objectFunc, eqfun=cov0constraintFunc, eqB=c(0,0), LB=paramLowerBound)

  moleculeOfConst <- denominatorOfConst <- 0
  for (i in 1:sum(freq)) moleculeOfConst <- moleculeOfConst + log(i)
  for (i in removeZero(freq)) {
    for (j in 1:i) {
      denominatorOfConst <- denominatorOfConst + log(j)
    }
  }
  G2 <- 2*((-fullModel(freq)) - (-solnp$value[length(solnp$value)]))
  print(G2)
}