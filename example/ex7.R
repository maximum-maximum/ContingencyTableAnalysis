# EXAMPLE 7.   Log-linear Model for 2x2x2 Table 
#
# Data Source:  Table 8.16, Agresti 347:2002
# 
#  y <- Y ~ MP(nu, p| strata=1, fixed="all");  i.e. Y ~ multinomial.
#
#  Specifically, 
#
#  y <- Y ~ multinomial(621, p).
#
# The counts in y are cross-classification counts for variables 
# G=Gender,  I=Information Opinion,   H= Health Opinion.  
#
# GOAL:  Fit the loglinear models [GI, GH, IH] and [G, IH].
#  

d <- scan(what=list(G ="",I="",H="",count=0))
Male   Support Support  76
Male   Support Oppose  160
Male   Oppose  Support   6
Male   Oppose  Oppose   25
Female Support Support 114
Female Support Oppose  181
Female Oppose  Support  11
Female Oppose  Oppose   48

d <- data.frame(d)

#Fit loglinear model [GI, GH, IH]...

a1 <- mph.fit(y=d$count, link="logm", X=model.matrix(~G+I+H + G:I + G:H + I:H, data=d))

#Fit loglinear model [G, IH]...

a2 <- mph.fit(y=d$count, link="logm", X=model.matrix(~G+I+H + I:H,data=d)) 


#Different Sampling Distribution Assumptions:
#
#Alternatively,  assume 
#   y <- Y ~ MP(nu, p| strata=1, fixed="none");    that is,  Y ~ indep Poisson.
#   
#   In other symbols, 
#   y <- Y,  where Y[i] indep ~  Poisson(m[i] = nu p[i]).
#                     Here,  nu is the unknown expected sample size.

b2 <- mph.fit(y=d$count, link="logm", X=model.matrix(~G+I+H + I:H, data=d), fixed="none")

#Alternatively,  assume 
#   y <- Y ~ MP(nu, p| strata=Gender, fixed="all").   That is, Y ~ prod multinomial.
#
#   In other symbols,   
#   y <- Y = (Y[1,1,1],Y[1,1,2],...,Y[2,2,2]),  where  (Y[i,1,1],...,Y[i,2,2]) indep ~ multinomial(n[i], p[i,,]).
#             Here,  p[i,j,k] = P(I=j,H=k|G=i) and    n[1]=267 and n[2]=354 are the a priori fixed 
#             samples sizes for males and females.  

c2 <- mph.fit(y=d$count, link="logm", X=model.matrix(~G+I+H + I:H, data=d), strata=d$G)

#Alternatively, assume
#   y <- Y ~ MP(nu, p| strata=Gender, fixed="none").  That is, Y ~ prod Poisson.
#
#   In other symbols, 
#   y <- Y = (Y[1,1,1],Y[1,1,2],...,Y[2,2,2]),  where Y[i,j,k] indep ~ Poisson(m[i,j,k] = nu[i] p[i,j,k]).
#             Here,  p[i,j,k] = P(I=j, H=k| G=i) and  nu[1] and nu[2] are the unknown expected 
#             sample sizes for males and for females.  

d2 <- mph.fit(y=d$count, link="logm", X=model.matrix(~G+I+H + I:H, data=d), strata=d$G, fixed="none")

  
cbind(a2$m,b2$m,c2$m,d2$m, sqrt(diag(a2$covm)),sqrt(diag(b2$covm)),sqrt(diag(c2$covm)),sqrt(diag(d2$covm)))
cbind(a2$p,b2$p,c2$p,d2$p, sqrt(diag(a2$covp)),sqrt(diag(b2$covp)),sqrt(diag(c2$covp)),sqrt(diag(d2$covp)))
