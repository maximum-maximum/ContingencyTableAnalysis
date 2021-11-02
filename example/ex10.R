#  EXAMPLE 10.   Cumulative Logit Model
#
#  Data Source:  Table 7.19, Agresti, 306:2002.
#
#  y <- Y ~ MP(nu, p| strata=Therapy*Gender, fixed="all");  i.e. Y ~ prod multinomial
#
#  Here,  y[i,j,k] is the cross-classification count corresponding to 
#  Therapy=i, Gender=j, Response=k.
#
#  The table probabilities are defined as p[i,j,k] = P(Response=k|Therapy=i,Gender=j).
#
#  Goal:  Fit the cumulative logit proportional odds model that includes 
#         the main effect of Therapy and Gender.
#

d <- scan(what=list(Therapy="",Gender="",Response="",count=0))
Sequential Male Progressive 28
Sequential Male NoChange 45
Sequential Male Partial 29
Sequential Male Complete 26
Sequential Female Progressive 4
Sequential Female NoChange 12
Sequential Female Partial 5
Sequential Female Complete 2
Alternating Male Progressive 41
Alternating Male NoChange 44
Alternating Male Partial 20
Alternating Male Complete 20
Alternating Female Progressive 12
Alternating Female NoChange 7
Alternating Female Partial 3
Alternating Female Complete 1

d <- data.frame(d)
strata <- paste(sep="",d$Therapy,".",d$Gender)
d <- data.frame(d,strata)

d3 <- subset(d,Response!="Complete")
levels(d3$Response) <- c(NA,"NoChange","Partial","Progressive")


L.fct <- function(p) {
   p <- matrix(p,4,4,byrow=T)
   clogit <- c()
   for (s in 1:4) {
     clogit <- c(clogit, 
         log(sum(p[s,1])  /sum(p[s,2:4])),
         log(sum(p[s,1:2])/sum(p[s,3:4])),
         log(sum(p[s,1:3])/sum(p[s,4]))
     )
   }
   L <- as.matrix(clogit)
   rownames(L) <- c(paste(sep="","log odds(R < ",2:4,"|",d3$strata,")"))

   L                 
}

a <- mph.fit(d$count,link=L.fct,X=model.matrix(~ -1 + Response + Therapy + Gender,data=d3) ,strata=strata)


#  Fit the related non-proportional odds cumulative logit model
b <- mph.fit(d$count,link=L.fct,
             X=model.matrix(~ Response+Response*Therapy + Response*Gender-1-Therapy-Gender,data=d3),strata=strata)

mph.summary(a,T)
mph.summary(b,T)
