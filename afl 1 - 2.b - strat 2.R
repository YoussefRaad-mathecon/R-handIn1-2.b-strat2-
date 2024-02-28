set.seed(2024)
rUS <- 0.03
rJ <- 0
sigmaX <- c(0.1, 0.02)
sigmaJ <- c(0, 0.25)
Y0 <- 1/100
X0 <- Y0
S0 <- 30000
K <- S0
T <- 2

muJ <- 0
muX <- 0

Nhedge <- 252*8*T
Nrep <- 10000

d1 <- function(S, K, rJ, sigmaX, sigmaJ, timetomat) {
  (log(S/K) + (rJ - t(sigmaX) %*% sigmaJ + norm(sigmaJ, type = "2")^2 / 2) * timetomat) / 
    (sqrt(timetomat) * norm(sigmaJ, type = "2"))
}
d2 <- function(S, K, rJ, sigmaX, sigmaJ, timetomat) {
  (log(S/K) + (rJ - t(sigmaX) %*% sigmaJ - norm(sigmaJ, type = "2")^2 / 2) * timetomat) / 
    (sqrt(timetomat) * norm(sigmaJ, type = "2"))
}



gdelta <- function(S, Y0, X0, K, rUS, rJ, sigmaX, sigmaJ, timetomat) {
  d1 <- d1(S, K, rJ, sigmaX, sigmaJ, timetomat)
  g <- Y0 * exp((rJ - t(sigmaX) %*% sigmaJ - rUS) * timetomat) * (pnorm(d1) - 1)
  deltaQP <- g/X0
  return(deltaQP)
}



FQP <- function(S, Y0, K, rUS, rJ, sigmaX, sigmaJ, timetomat) {
  d1 <- d1(S, K, rJ, sigmaX, sigmaJ, timetomat)
  d2 <- d2(S, K, rJ, sigmaX, sigmaJ, timetomat)
  priceQP <- Y0 * exp(-rUS * timetomat) * (K * pnorm(-d2) - exp((rJ - t(sigmaX) %*% sigmaJ) * timetomat)* S * pnorm(-d1))
}

St <- rep(S0, length = Nrep)
Xt <- rep(X0, length = Nrep)
dt <- T/Nhedge


initial_outlay <- c(FQP(S0, Y0, K, rUS, rJ, sigmaX, sigmaJ, timetomat = T))
Vpf <- rep(initial_outlay, length = Nrep)
a <- gdelta(St, Y0, Xt, K, rUS, rJ, sigmaX, sigmaJ, timetomat = T)
b <- (-a)*St
c <- Vpf - a*St*Xt - b*Xt

for (i in 2:Nhedge) {
  W <- cbind(rnorm(Nrep), rnorm(Nrep))
  St <- St * exp((muJ + rJ - t(sigmaX) %*% sigmaJ - 0.5 * norm(sigmaJ, type = "2")^2) * dt + sqrt(dt)*c(W%*%sigmaJ))
  Xt <- Xt * exp((muX + rUS - rJ - 0.5 * norm(sigmaX, type = "2")^2) * dt + sqrt(dt)*c(W%*%sigmaX))
  Vpf <- a * St * Xt + b * exp(dt * rUS)*Xt + c * exp(dt*rUS)
  a <- gdelta(St, Y0, Xt, K, rUS, rJ, sigmaX, sigmaJ, timetomat = (T - (i - 1) * dt))
  b <- (-a)*St
  c <- Vpf - a*St*Xt - b*Xt
}


W <- cbind(rnorm(Nrep), rnorm(Nrep))
St <- St * exp((muJ + rJ - t(sigmaX) %*% sigmaJ - 0.5 * norm(sigmaJ, type = "2")^2) * dt + sqrt(dt)*c(W%*%sigmaJ))
Xt <- Xt * exp((muX + rUS - rJ - 0.5 * norm(sigmaX, type = "2")^2) * dt + sqrt(dt)*c(W%*%sigmaX))

Stmax <- 0:max(St)

plot(St, Vpf, col=c("red", "blue", "green", "orange"), 
     xlab="S(T)", ylab="Value of hedge portfolio", 
     xlim=c(0, max(St)), main = "Value of hedge portfolio vs. S(T)", 
     xaxt='n') 
text(100, 200, paste("Hedge point count =", Nhedge), adj=-1)
lines(Stmax, Y0 * pmax(K - Stmax, 0), type="l", lwd=3)


ticks <- seq(0, max(St), by=5000) 

axis(1, at=ticks, labels=format(ticks, big.mark=",", scientific=FALSE))


