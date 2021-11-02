#  EXAMPLE 6.  Test Marginal Homogeneity in a 3x3 Table
# 
#  Data Source:  Table 10.16, Agresti, 445:2002.
#
#  y <- Y ~ MP(nu, p |strata=1, fixed="all");  i.e. Y ~ multinomial.
#
#  Specifically, 
#  y <- Y ~ multinomial(160, p=(p[1,1],...,p[3,3]))
#
#  GOAL:    Test H0:  p[1,+] = p[+,1], p[2,+] = p[+,2], p[3,+] = p[+,3] vs.  H1: not H0.

d <- scan(what=list(Siskel="",Ebert="",count=0))
Pro Pro     64
Pro Mixed    9
Pro Con     10
Mixed Pro   11
Mixed Mixed 13
Mixed Con    8
Con Pro     13
Con Mixed    8
Con Con     24


d <- data.frame(d)

h.fct <- function(p){
    p.Siskel <- M.fct(d$Siskel)%*%p
    p.Ebert  <- M.fct(d$Ebert)%*%p
    as.matrix(c(p.Siskel[-3] - p.Ebert[-3]))
}
a1 <- mph.fit(y=d$count, h.fct=h.fct)
mph.summary(a1,T)

#  Suppose that we wish to report on the observed and fitted marginal probabilities
#  

L.fct <- function(p) {
    p.Siskel <- M.fct(d$Siskel)%*%p
    p.Ebert <- M.fct(d$Ebert)%*%p
    L <- as.matrix(c(p.Siskel,p.Ebert))
    rownames(L) <- c(paste(sep="", "P(Siskel=",levels(d$Siskel),")"), 
                     paste(sep="", "P(Ebert=", levels(d$Ebert),")"))
    L
}
a2 <- mph.fit(y=d$count,h.fct=h.fct,L.fct=L.fct,X=diag(6))
mph.summary(a2,T)

#  M.fct(factor)%*%p  gives the marginal probabilities corresponding to the 
#  levels of 'factor'.   The  marginal probabilities are ordered by the levels of 'factor'.
#
#  Alternatively,  in this rectangular table setting, we can find the marginal 
#  probabilities using the apply(...) function.  In this case, the marginal 
#  probabilities are ordered as they are entered in the data set.
#

h.fct <- function(p) {
    p <- matrix(p,3,3,byrow=T)
    p.Siskel <- apply(p,1,sum)
    p.Ebert <- apply(p,2,sum)
    as.matrix(c(p.Siskel[-3] - p.Ebert[-3]))
}

L.fct <- function(p) {
    p <- matrix(p,3,3,byrow=T)
    p.Siskel <- apply(p,1,sum)
    p.Ebert <- apply(p,2,sum)
    L <- as.matrix(c(p.Siskel,p.Ebert))
    rownames(L) <- c("P(Siskel=Pro)","P(Siskel=Mixed)","P(Siskel=Con)",
                     "P(Ebert=Pro)", "P(Ebert=Mixed)", "P(Ebert=Con)")
    L
}
b <- mph.fit(y=d$count,h.fct=h.fct,L.fct=L.fct,X=diag(6))

