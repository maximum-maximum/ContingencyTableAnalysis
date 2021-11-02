#  EXAMPLE 3.  Test whether a multinomial probability vector satisfies a particular constraint.
#
#   Data Source:  Agresti 25:2002.
#
#   y=(30,63,63) <- Y ~ MP(nu, p| strata=1, fixed="all");  i.e. Y ~ multinomial.
#
#   In other symbols,
#
#   y=(30,63,63)  <- Y ~ multinomial(156, p=(p[1],p[2],p[3]))
#
#   GOAL:    Test H0: p[1]+p[2] = p[1]/(p[1]+p[2])   vs. H1: not H0.
#

y <- c(30, 63, 63)
h.fct <- function(p) {
    (p[1] + p[2]) - p[1]/(p[1]+p[2])
}
a <- mph.fit(y, h.fct=h.fct)
mph.summary(a,T)

