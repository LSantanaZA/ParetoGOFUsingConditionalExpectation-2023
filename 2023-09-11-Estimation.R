
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ PARAMETER ESTIMATION ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ppareto <- function(x,beta,sigma=1){
  out = 1-(x/sigma)^-beta
  return(out)
}

qpareto <- function(p,beta,sigma=1){
  out = sigma*(1-p)^(-1/beta)
  return(out)
}

dpareto <- function(x,beta,sigma=1){
  out = beta*sigma^beta/(x^beta+1)
  return(out)
}

rpareto <- function(n,beta,sigma=1){
  X = sigma*(1-runif(n))^(-beta^-1)
  return(X)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ MLE sigma known ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMLE1 <- function(X){
  n   = length(X)
  out = list(betaH=n/sum(log(X)))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ MLE sigma unknown ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMLE2 <- function(X){
  n      = length(X)
  sigmaH = min(X)
  betaH  = n/sum(log(X/sigmaH))
  out    = list(betaH=betaH,sigmaH=sigmaH)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ MME sigma known ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMME1 <- function(X){
  out = list(betaH=mean(X)/(mean(X)-1))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ MME sigma unknown ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMME2 <- function(X){
  betaH  = (n*mean(X)-min(X))/(n*(mean(X)-min(X)))
  sigmaH = mean(X)*(betaH-1)/betaH
  out    = list(betaH=betaH,sigmaH=sigmaH)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#--------------------------------------------------------------------------------------------
#Minimum Distance Estimator (LS)
ParetoMDE.LS <- function(x) {
  n <- length(x)
  MD.LS <- function(param) {
    beta <- param[1]
    ans1 <- (1 / n) * sum((ppareto(x, beta, 1) - (1:n) / (n + 1)) ^ 2)
    return(ans1)
  }
  b <- mean(x)/(mean(x)-1) #MME estimator as "first guess"
  ans2 <- optim(c(b), MD.LS, method = "BFGS")
  ans3 <- list(betaH=ans2$par[1],convergence=ans2$convergence)
  return(ans3)
}
#--------------------------------------------------------------------------------------------
#Minimum Distance Estimator (CvM)
ParetoMDE.CvM <- function(x) {
  n <- length(x)
  MD.CvM <- function(param) {
    beta <- param[1]
    ans1 <-   1/(12*n) + sum((ppareto(x,beta,1)-(2*(1:n)-1)/(2*n))^2)
    return(ans1)
  }
  b <- mean(x)/(mean(x)-1) #MME estimator as "first guess"
  ans2 <- optim(c(b), MD.CvM, method = "BFGS")
  ans3 <- list(betaH=ans2$par[1],convergence=ans2$convergence)
  return(ans3)
}
#--------------------------------------------------------------------------------------------
#Minimum Distance Estimator (Phi=KL)
ParetoMDE.Phi.kl <- function(x) {
  n <- length(x)
  MD.Phi.kl <- function(param) {
    beta <- param[1]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dpareto(x, beta, 1) / fhat(x)
    mean(phi.kl(fdg))
  }
  b <- mean(x)/(mean(x)-1) #MME estimator as "first guess"
  ans2 <- optim(c(b), MD.Phi.kl, method = "Nelder-Mead")
  ans3 <- list(betaH = ans2$par[1], convergence=ans2$convergence)
  return(ans3)
}
#--------------------------------------------------------------------------------------------
#Minimum Distance Estimator (Phi=ChiSq)
ParetoMDE.Phi.chisq <- function(x) {
  n <- length(x)
  MD.Phi.chisq <- function(param) {
    beta <- param[1]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dpareto(x, beta, 1) / fhat(x)
    gdf <- fdg ^ (-1)
    mean(gdf * phi.chisq(fdg))
  }
  b <- mean(x)/(mean(x)-1) #MME estimator as "first guess"
  ans2 <- optim(c(b), MD.Phi.chisq, method = "Nelder-Mead")
  ans3 <- list(betaH = ans2$par[1], convergence=ans2$convergence)
  return(ans3)
}
#--------------------------------------------------------------------------------------------
#Minimum Distance Estimator (Phi=Total Variation)
ParetoMDE.Phi.tv <- function(x) {
  n <- length(x)
  MD.Phi.tv <- function(param) {
    beta <- param[1]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dpareto(x, beta, 1) / fhat(x)
    gdf <- fdg ^ (-1)
    mean(gdf * phi.tv(fdg))
  }
  b <- mean(x)/(mean(x)-1) #MME estimator as "first guess"
  ans2 <- optim(c(b), MD.Phi.tv, method = "Nelder-Mead")
  ans3 <- list(betaH = ans2$par[1], convergence=ans2$convergence)
  return(ans3)
}


#--------------------------------------------------------------------------------------------
# if(FALSE){
#   n     = 5000
#   beta  = 3
#   sigma = 2
#   X     = sort(rpareto(n,beta,sigma))
#   
#   betaH = ParetoMLE1(X/sigma)
#   Y     = (X/sigma)^betaH
#   betaH
#   ParetoMLE1(Y)
#   
#   est = ParetoMLE2(X)
#   Y   = (X/est$sigmaH)^est$betaH
#   est
#   ParetoMLE2(Y)
#   
#   betaH = ParetoMME1(X/sigma)
#   Y     = (X/sigma)^betaH
#   betaH
#   ParetoMME1(Y)
#   
#   est = ParetoMME2(X)
#   Y   = (X/est$sigmaH)^est$betaH
#   est
#   ParetoMLE2(Y)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
