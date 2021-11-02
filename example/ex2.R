#  EXAMPLE 2. Test whether a multinomial probability vector is uniform. 
#             Test whether a multinomial probability vector equals a specific value.  
#
#   y <- Y=(Y[1],...,Y[6]) ~ MP(nu, p|strata=1, fixed="all");  i.e. Y ~ multinomial.
#
#   In other symbols,
#
#   y <- Y ~ multinomial(15, p=(p[1],...,p[6]))
#
#   GOAL:   Test H0: p[1]=p[2]=...=p[6]   vs. H1:  not H0.
#

y <- rmultinom(1,15,rep(1,6)) 
a1 <- mph.fit(y,L.fct = function(p){p}, X=matrix(1,6,1), y.eps=0.1)

#Alternative specification...
a2 <- mph.fit(y,h.fct = function(p){as.matrix(p[-6]-1/6)},y.eps=0.1)

mph.summary(a1,T); mph.summary(a2,T)

#Test whether p = (1, 2, 3, 1, 2, 3)/12...
#
p0 <- c(1,2,3,1,2,3)/12
b <-  mph.fit(y,h.fct= function(p){as.matrix(p[-6]-p0[-6])},y.eps=0.1)
mph.summary(b,T)
