# EXAMPLE 1.  Test whether a binomial probability equals 0.5
#
#  y=(15, 22) <-  Y ~ MP(nu,p|strata=1,fixed="all");  i.e. Y ~ multinomial
#  
#  In other symbols,
#  
#     y=(15, 22) <-  Y=(Y[1], Y[2])  ~ multinomial(37, p=(p[1],p[2])).
#  
#  GOAL:    Test H0:  p[1] = 0.5   vs.  H1:  not H0.

a1 <- mph.fit(y=c(15, 22), constraint=function(p){p[1]-0.5})

#Alternative specifications...
a2 <- mph.fit(y=c(15, 22), constraint=function(p){p[1]-p[2]})
a3 <- mph.fit(y=c(15, 22), constraint=function(p){log(p[1]/p[2])})
a4 <- mph.fit(y=c(15, 22), constraint=function(m){m[1]-m[2]},h.mean=T)
a5 <- mph.fit(y=c(15, 22), link=function(p){p}, X=matrix(1,2,1))
a6 <- mph.fit(y=c(15, 22), link="logm",X=matrix(1,2,1))

#  
#  Alternatively, assume that 
#    
#      y=(15, 22) <- Y ~ MP(nu, p|strata=1,fixed="none"); i.e. Y ~ indep Poisson.
#  
#  In other symbols,
#      y=(15, 22) <- Y = (Y[1],Y[2]),  where  Y[i] indep ~  Poisson(nu p[i]), i=1,2.
#
#  GOAL:     Test H0:  p[1] = 0.5  vs.  H1: not H0.  

b1 <- mph.fit(y=c(15, 22), constraint=function(p){p[1]-0.5},fixed.strata="none")

mph.summary(a1,T)
mph.summary(b1,T)
