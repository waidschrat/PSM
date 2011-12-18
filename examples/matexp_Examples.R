rm(list=ls())

# R CMD SHLIB matexp.f - TEST
# dyn.load("~/PSM/Rpackages/Development/matexp/src/matexp.so")



detach(package:matexp)
detach(package:mexp)

hostname <- "gridterm"

library(matexp,lib.loc=paste("~/PSM/Rpackages/",hostname,sep=""))
library(mexp,  lib.loc=paste("~/PSM/Rpackages/",hostname,sep=""))


M <- 3
(x <- matrix( 1:(M**2) , ncol=M))
d <- 0
(x <- matrix( c(2,0,1,2-d) , ncol=2))


matexp(x,order=8)
mexp(  x,order=8)


# Identity test
matexp(x,order=8) %*% matexp(-x,order=8)
mexp(x,order=8) %*% mexp(-x,order=8)


# Perf Test
T1 <- system.time(
for(i in 1:10000) {
  matexp(x,order=8) })
T2 <- system.time(
for(i in 1:10000) {
  mexp(x,order=8) })

# user cpu, system cpu, elapsed, subproc1, subproc2 times
T1
T2
T1[1]/T2[1]


# Test from mexp
test1 <- t(matrix(c(
                    4, 2, 0,
                    1, 4, 1,
                    1, 1, 4), 3, 3))
matexp(test1)
## Results on Power Mac G3 under Mac OS 10.2.8
##                    [,1]               [,2]               [,3]
## [1,] 147.86662244637000 183.76513864636857  71.79703239999643
## [2,] 127.78108552318250 183.76513864636877  91.88256932318409
## [3,] 127.78108552318204 163.67960172318047 111.96810624637124
## -- these agree with ward (1977, p608)
##
## A naive alternative to mexp, using spectral decomposition:

mexp2 <- function(matrix){
  z <- eigen(matrix,sym=FALSE)
  Re(z$vectors %*% diag(exp(z$values)) %*%
     solve(z$vectors))
}
try(
    mexp2(test1)
    ) ## now gives an error from solve !
##
## older result was
##                   [,1]                [,2]               [,3]
##[1,] 147.86662244637003  88.500223574029647 103.39983337000028
##[2,] 127.78108552318220 117.345806155250600  90.70416537273444
##[3,] 127.78108552318226  90.384173332156763 117.66579819582827
## -- hopelessly inaccurate in all but the first column.
##
##

## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
test2 <- t(matrix(c(
                    29.87942128909879, .7815750847907159, -2.289519314033932,
                    .7815750847907159, 25.72656945571064,  8.680737820540137,
                    -2.289519314033932, 8.680737820540137,  34.39400925519054),
                  3, 3))
matexp(test2)
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692357 -18231880972009844 -30475770808580828
##[2,] -18231880972009852  60605228702227024 101291842930256144
##[3,] -30475770808580840 101291842930256144 169294411240859072
## -- which agrees with Ward (1977) to 13 significant figures
mexp2(test2)
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692405 -18231880972009100 -30475770808580196
##[2,] -18231880972009160  60605228702221760 101291842930249376
##[3,] -30475770808580244 101291842930249200 169294411240850880
## -- in this case a very similar degree of accuracy.
##


## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
test3 <- t(matrix(c(
                    -131, 19, 18,
                    -390, 56, 54,
                    -387, 57, 52), 3, 3))
matexp(test3)
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608.
mexp2(test3)
##                   [,1]               [,2]                [,3]
##[1,] -1.509644158796182 0.3678794391103086 0.13533528117547022
##[2,] -5.632570799902948 1.4715177585023838 0.40600584352641989
##[3,] -4.934938326098410 1.1036383173309319 0.54134112676302582
## -- in this case, a similar level of agreement with Ward (1977).
##

## ----------------------------
## Test case 4 from Ward (1977)
## ----------------------------
test4 <-
  structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
            .Dim = c(10, 10))
attributes(matexp(test4))
max(abs(matexp(test4) - mexp2(test4)))
##[1] 8.746826694186494e-08
## -- here mexp2 is accurate only to 7 d.p., whereas mexp
##    is correct to at least 14 d.p.
##
## Note that these results are achieved with the default
## settings order=8, method="Pade" -- accuracy could
## presumably be improved still further by some tuning
## of these settings.
