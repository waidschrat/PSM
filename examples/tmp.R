
XX = 0:200

x <- 14
n <- 27000
p <- x/n

prbs <- dbinom(XX, n, p)

plot(XX,prbs,type="s")

