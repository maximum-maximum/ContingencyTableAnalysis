#  EXAMPLE 9.  Marginal Standardization of a Contingency Table  
#
#  Data Source:  Table 8.15, Agresti 345:2002.   
#
#  GOAL:   
#     For a two-way table,  find the standardized values of y, say y*, 
#     that  satisfy  (i)  y* has the same odds ratios as y and 
#                    (ii) y* has row and column totals equal to 100.  
#  
#  Note:  This is equivalent to the problem of finding the fitted 
#         values for the following model...
#             x <- Y ~ multinomial(n, p=(p[1,1],...,p[3,3]))
#                  p[1,+] = p[2,+] = p[3,+] = p[+,1] = p[+,2] = p[+,3] = 1/3
#                  p[1,1]*p[2,2]/p[2,1]/p[1,2] = or[1,1]
#                  p[1,2]*p[2,3]/p[2,2]/p[1,3] = or[1,2]
#                  p[2,1]*p[3,2]/p[3,1]/p[2,2] = or[2,1]
#                  p[2,2]*p[3,3]/p[3,2]/p[2,3] = or[2,2],
#                  where  or[i,j] = y[i,j]*y[i+1,j+1]/y[i+1,j]/y[i,j+1] are 
#                  the observed (y) odds ratios.
#         If m is the vector of fitted values, then y* = m*300/sum(m)
#         are the standardized values of y.  
#         Here  x can be any vector of 9 counts.
#         Choosing x so that the sum is 300 leads to sum(m) = 300, so that 
#         y* = m in this case.  
#      
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

d <- data.frame(d)

h.fct <- function(p) {
   p.Schooling <- M.fct(d$Schooling)%*%p
   p.Abortion  <- M.fct(d$Abortion )%*%p
   p <- matrix(p,3,3,byrow=T)
   as.matrix(c(
                p.Schooling[-3]-1/3,p.Abortion[-3]-1/3,
                p[1,1]*p[2,2]/p[2,1]/p[1,2] - 209*126/151/101,
                p[1,2]*p[2,3]/p[2,2]/p[1,3] - 101*426/126/237,
                p[2,1]*p[3,2]/p[3,1]/p[2,2] - 151*21/16/126,
                p[2,2]*p[3,3]/p[3,2]/p[2,3] - 126*138/21/426
             ))
}

b <- mph.fit(y=d$count,h.fct=h.fct)  
ystar <-   b$m*300/sum(b$m) 
matrix(round(ystar,1),3,3,byrow=T)
 

x <- c(rep(33,8),36)
b <- mph.fit(y=x,h.fct=h.fct) 
ystar <- b$m
matrix(round(ystar,1),3,3,byrow=T)
