#  EXAMPLE 4.  Test of Independence in a 2x2 Table.  
#
#  y=(y[1,1],y[1,2],y[2,1],y[2,2]) = (25,18,13,21) <- Y ~ MP(nu, p|strata=1, fixed="all");  i.e. Y ~ multinomial
# 
#  In other symbols,
#  y =(y[1,1],y[1,2],y[2,1],y[2,2]) <- Y ~ multinomial(77, p=(p[1,1],p[1,2],p[2,1],p[2,2]))
#
#  GOAL:    Test H0: p[1,1]*p[2,2]/p[1,2]/p[2,1] = 1    vs. H1: not H0.

d <- scan(what=list(A="",B="",count=0))

#1 1 25
#1 2 18
#2 1 13
#2 2 21

d <- data.frame(d)

a1 <- mph.fit(y=d$count, h.fct=function(p){log(p[1]*p[4]/p[2]/p[3])})

# Alternative specifications of independence....
a2 <- mph.fit(y=d$count, h.fct=function(p){p <- matrix(p,2,2,byrow=T); log(p[1,1]*p[2,2]/p[1,2]/p[2,1])})
a3 <- mph.fit(y=d$count, h.fct=function(p){p[1]*p[4]/p[2]/p[3] - 1})
a4 <- mph.fit(y=d$count, h.fct=function(p){p[1]/(p[1]+p[2]) - p[3]/(p[3]+p[4])})
a5 <- mph.fit(y=d$count, L.fct="logm", X=model.matrix(~A+B,data=d))


# Suppose we wished to output observed and fitted values of log OR, OR, and P(B=1|A=1)-P(B=1|A=2)...

L.fct <- function(p) {
    L <- as.matrix(c(
               log(p[1]*p[4]/p[2]/p[3]),
               p[1]*p[4]/p[2]/p[3],
               p[1]/(p[1]+p[2]) - p[3]/(p[3]+p[4])
         ))
    rownames(L) <- c("log OR","OR","P(B=1|A=1) - P(B=1|A=2)")
    L
}

a6 <- mph.fit(y=d$count,h.fct=function(p){log(p[1]*p[4]/p[2]/p[3])}, L.fct=L.fct,X=diag(3))

#  Unrestricted Model...
b <- mph.fit(y=d$count, L.fct=L.fct, X=diag(3))

mph.summary(a6,T)
mph.summary(b,T)
