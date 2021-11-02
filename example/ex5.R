#  EXAMPLE 5.  Test of Independence in a 4x4 Table.  (Using Log-Linear Model.) 
#
#  Data Source:  Table 2.8, Agresti, 57:2002. 
#
#  y <- Y ~ MP(nu, p|strata=1, fixed="all");  i.e. Y ~ multinomial
# 
#  In other symbols,
#  y <- Y ~ multinomial(96, p=(p[1,1],p[1,2],p[2,1],p[2,2]))
#
#  GOAL:    Test H0: p[1,1]*p[2,2]/p[1,2]/p[2,1] = 1    vs. H1: not H0.
#

d <- scan(what=list(Income="",JobSatisf="",count=0))
<15 VD 1
<15 LD 3
<15 MS 10
<15 VS 6
15-25 VD 2
15-25 LD 3
15-25 MS 10
15-25 VS 7
25-40 VD 1
25-40 LD 6
25-40 MS 14
25-40 VS 12
>40 VD 0
>40 LD 1
>40 MS 9
>40 VS 11

d <- data.frame(d)

a <- mph.fit(y=d$count, link="logp", X=model.matrix(~Income+JobSatisf,data=d))
mph.summary(a)

#Alternatively,
b <- mph.fit(y=d$count, link="logm", X=model.matrix(~Income+JobSatisf,data=d))
mph.summary(b)
