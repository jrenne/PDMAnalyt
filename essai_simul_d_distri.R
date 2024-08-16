
Pi <- .03
Dy <- .02
zeta <- exp(-Pi-Dy)
q <- .06
beta <- .1
d_bar <- .8
s_star <- d_bar*(zeta*(1+q) - 1 - beta)
sigma <- .02

chi <- .95

n <- 10000

r <- q * zeta * d_bar
r_1 <- r
d <- d_bar
d_1 <- d_bar
d_2 <- d_bar

all_d <- d
all_r <- r

for(i in 1:n){
  r <- q * zeta * (d_1 - chi * zeta * d_2) + zeta * chi * r_1
  d <- zeta * d_1 - s_star - beta * d_1 - sigma * rnorm(1) + r
  
  all_d <- c(all_d,d)
  all_r <- c(all_r,r)
  
  d_1 <- d
  d_2 <- d_1
  r_1 <- r
}

plot(all_d,type="l")
plot(all_r,type="l")

print(c(mean(all_d),sd(all_d)))

plot(density(all_d))

