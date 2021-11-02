# EXAMPLE 8.  Fit Linear-by-Linear Log-Linear Model
#
# Data Source:  Table 8.15,  Agresti, 345:2002
#
# y <- Y ~ MP(nu, p|strata=1, fixed="all");  i.e. Y ~ multinomial
#
# Specifically,
# y <- Y ~ multinomial(1425, p)
#
# GOAL:  Assess the fit of the linear-by-linear log-linear model.
#

d <- scan(what=list(Schooling="",Abortion = "", count=0))
<HS  Disapprove  209
<HS  Middle      101
<HS  Approve     237
HS   Disapprove  151
HS   Middle      126
HS   Approve     426
>HS  Disapprove   16
>HS  Middle       21
>HS  Approve     138

Schooling.score <- -1*(d$Schooling=="<HS") + 0*(d$Schooling=="HS") + 1*(d$Schooling==">HS")
Abortion.score  <- -1*(d$Abortion=="Disapprove") + 0*(d$Abortion=="Middle") +  1*(d$Abortion=="Approve")

d <- data.frame(d,Schooling.score,Abortion.score)

a <- mph.fit(y=d$count, link="logm",X=model.matrix(~Schooling + Abortion + Schooling.score:Abortion.score,data=d))
mph.summary(a,T)

