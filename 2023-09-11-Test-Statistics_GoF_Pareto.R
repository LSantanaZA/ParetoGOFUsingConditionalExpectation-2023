# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New classes of tests for the Pareto distribution based on a conditional expectation characterisation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# T. Nombebe, J.S. Allison, L. Santana and I.J.H. Visagie
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Abstract
# We propose and investigate new classes of goodness-of-fit tests for the Pareto Type I distribution based on a characterisation involving a conditional expectation. A Monte Carlo power study is included in order to assess the finite sample performance of the newly developed tests for four different estimation methods. The results from the simulation study shows that the newly proposed tests are competitive in terms of power performance when compared to some existing tests. It also shows that the majority of tests produce their highest powers when the unknown shape parameter is estimated by the method of moments. A practical example, where we consider the annual salaries of English Premier League football players for two consecutive seasons, is also included to illustrate the use of the newly proposed test


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~# Thobeka 1 (L_{n,s,a}^{[1]}) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First test statistic for the third paper
Tns1 <- function(X,beta,a,s){
  X <- sort(X)
  n <- length(X)
  cumul <- 0
  for(j in 1:n){
    for(k in 1:n){
      term1 <- (X[j]^s) * (X[k]^s)/(n^2) 
      term2 <- (1 - min(X[j],X[k])^(1-a) )/(a-1)
      term3 <- 2*(X[j]^(s-beta-a +1) - 1)/(s - beta - a + 1)
      term4 <- 1/(2*s-2*beta - a + 1)   
      cumul <- cumul + term1* ( term2  - term3 - term4)
    }
  }
  return(cumul)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~# Thobeka 2 (L_{n,s,a}^{[2]}) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Second test statistic for the third paper
Tns2 <- function(X,beta,a,s){
  X <- sort(X)
  n <- length(X)
  cumul <- 0
  for(j in 1:n){
    for(k in 1:n){
      term1 <- (X[j]^s) * (X[k]^s)/(n^2) 
      term2 <- (1 - min(X[j],X[k])^(1-a) )/(a-1)
      term3 <- 2*(min(X[j],X[k])^(s-a +1) - 1)/(n*(s -  a + 1))
      term4 <- (min(X[j],X[k])^(2*s-a+1) - 1)/(n^2*(2*s - a + 1))   
      cumul <- cumul + term1* ( term2  - term3 + term4)
    }
  }
  return(cumul)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~# Thobeka 3 (T_{n,s,a}^{[1]}) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Third test statistic for the third paper

Tns3 <- function(X,beta,a,s){
  X <- sort(X)
  n <- length(X)
  term1 <- 0
  for(j in 1:n){
    for(k in 1:n){
      term1.1 <- (X[j]^s) * (X[k]^s)/(n^2*(a-1))
      term1.2 <- (1 - min(X[j],X[k])^(1-a) )
      term1 <- term1 + term1.1*term1.2 
    }
  }
  
  term2 <- sum((2*beta/(n*(beta-s)*(s-beta-a+1)))*(X^s*(X^(s-beta-a+1) -1)))
  term3 <- beta^2/((beta-s)^2*(2*s-2*beta-a+1))
  ans <- term1 - term2 - term3
  return(ans)
} 

# #Fourth test for the third paper
# 
# Tns4 <- function(X,beta,a,ss){ #NOTE: h should be calculated inside the function 
#   X <- sort(X)
#   n <- length(X)
#   h <- 1.06*sd(X)*n^(-1/5) # normal scaled rule of thumb (Silverman))
#   #h = 1.06*S(X)*n^(-0.2)
#   tmp_A <- tmp_B <- tmp_C <- 0
#   for(j in 1:n){
#     for(k in 1:n){
#       tmp_A <- tmp_A +  ((X[j]^ss) * (X[k]^ss) * (1  - min(X[j],X[k])^(1-a) ))/(n^2*(a-1))
#       tmp_B <- tmp_B +  ((X[j]^ss)*( (min(X[j],X[k] + h))^(ss+2-a)   - (max(X[k] - h,1))^(ss-a+2)) )/(h*2*n^2*(beta-ss)*(ss-a+2))
#       tmp_C <- tmp_C +  ( (min(X[j] + h, X[k] + h))^(2*ss+3-a) - (max(X[j] - h, X[k] - h,1))^(2*ss+3-a)  )/((h^2)*4*n^2*((beta-ss)^2)*(2*ss+3-a))
#     }
#   }
#   out <- tmp_A - tmp_B + tmp_C
#   return(out)
# }
# 
# #Fifth test for the third paper
# 
# Tns5 <- function(X,beta,a,ss){ 
#   X    <- sort(X)
#   n    <- length(X)
#   fh   <- approxfun(density(X))
#   cuml <- 0
#   for(kk in 1:n){
#     cuml <- cuml + (X[kk]^(-a))*
#       (
#         mean((X^ss)*as.numeric(X > X[kk])) - (fh(X[kk])*(X[kk]^(ss+1)))/(beta-ss)
#       )^2
#   }
#   return(cuml/n)
# }
# 
# #Sixth test for the third paper
# 
# Tns6 <- function(X,beta,a,ss){ 
#   X    <- sort(X)
#   n    <- length(X)
#   fh   <- approxfun(density(X))
#   cuml <- 0
#   for(kk in 1:n){
#     cuml <- cuml + (X[kk]^(-a))*
#       (
#         beta*mean((X^ss)*as.numeric(X > X[kk])) - mean(X^ss)*fh(X[kk])*X[kk]^(ss+1)
#       )^2
#   }
#   return(cuml/n)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~ (18) Lethani et al (L_{n,m,a}) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sn2 <- function(X,beta,m=3,a=2){ #beta is not actually used, but included to conform to the argument list of all other statistics
  n = length(X)
  X = sort(X)
  sumj = 0
  for (j in 1:n){
    sumk = 0 
    for (k in 1:n){
      vj = ((n-j+1)^m-(n-j)^m)/n^m
      vk = ((n-k+1)^m-(n-k)^m)/n^m
      termA = exp(-((X[j]^(1/m)-X[k]^(1/m))^2)/(4*a))
      termB = 2*n*vj*exp(-(X[j]-X[k]^(1/m))^2/(4*a))
      termC = n^2*vj*vk*exp(-(X[j]-X[k])^2/(4*a))
      sumk = sumk + termA-termB+termC
    }
    sumj = sumj + sumk
  }
  out <- sqrt(pi/a)*sumj/n^2
  return(out)
}


#############################################################################################



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (1) Kolmogorov-Smirnov ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KSn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  FH  = ppareto(X,betaH,sigmaH)
  KSp = max((1:n)/n-FH)
  KSm = max(FH-(0:(n-1))/n)
  out = max(KSp,KSm)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (2) Cram?r-von Mises ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CVn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  FH  = ppareto(X,betaH,sigmaH)
  out = sum((FH-(2*(1:n)-1)/(2*n))^2)+1/(12*n)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (3) Anderson-Darling ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ADn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  z   = ppareto(X,betaH,sigmaH)
  T1  = 2*(1:n)-1
  T2  = log(z) + log(1-z[n:1])
  out = -n-mean(T1*T2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~ (4) Modified Anderson-Darling ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  z   = ppareto(X,betaH,sigmaH)
  S1  = sum(z)
  T1  = 2-(2*(1:n)-1)/n
  T2  = log(1-z)
  S2  = sum(T1*T2)
  out = n/2-2*S1-S2
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (5) Zhang's A ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZAn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  j   = 1:n
  FH  = ppareto(X,betaH,sigmaH)
  out = -sum(log(FH)/(n-j+0.5)+log(1-FH)/(j-0.5))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (6) Zhang's B ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ZBn <- function(X,betaH,sigmaH=1){
#   n   = length(X)
#   j   = 1:n
#   FH  = ppareto(X,betaH,sigmaH)
#   arg = (1/FH-1)/((n-0.5)/(j-0.75)-1)
#   out = sum((log(arg))^2)
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (7) Zhang's C ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ZCn <- function(X,betaH,sigmaH=1){
#   n   = length(X)
#   j   = 1:n
#   FH  = ppareto(X,betaH,sigmaH)
#   T1  = n*(j-0.5)/(n-j+0.5)^2*log((j-0.5)/(n*FH))
#   T2  = n/(n-j+0.5)*log((n-j+0.5)/(n*(1-FH)))
#   out = 2*sum(T1+T2)
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (8) Kullback-Leibler ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# H <- function(X,m){
#   n    = length(X)
#   j    = 1:n
#   ind1 = pmin(j+m,n)
#   ind2 = pmax(j-m,1)
#   out  = mean(log(n/(2*m)*(X[ind1]-X[ind2])))
#   return(out)
# }
# 
# KLn <- function(X,betaH,m,sigmaH=1){
#   n   = length(X)
#   out = -H(X,m)-log(betaH)-betaH*log(sigmaH)+(betaH+1)*mean(log(X))
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (9) Hellinger distance ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# HDn <- function(X,m,betaH,sigmaH=1){
#   n    = length(X)
#   j    = 1:n
#   ind1 = pmin(j+m,n)
#   ind2 = pmax(j-m,1)
#   num  = (n/(2*m)*(X[ind1]-X[ind2]))^(-0.5)-sqrt(dpareto(X,betaH,sigmaH))
#   den  = (n/(2*m)*(X[ind1]-X[ind2]))^(-1)
#   out  = sum(num^2/den)/(2*n)
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Phi-divergence ~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fhHat <- function(X){
#   n   = length(X)
#   s   = sd(X)
#   h   = 1.06*s*n^(-1/5)
#   out = rep(0,n)
#   for(j in 1:n){
#     out[j] = mean(dnorm((X[j]-X)/h))/h
#   }
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (10) Kullback-Leibler ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DKn <- function(X,betaH,sigmaH=1){
#   fH  = fhHat(X)
#   out = mean(log(fH/dpareto(X,betaH,sigmaH)))
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (11) Meintanis (2009a) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# MEnOLD <- function(X,a,betaH,sigmaH=1){
#   n  = length(X)
#   FH = ppareto(X,betaH,sigmaH)
#   S1 = matrix(0,n,n)
#   for (j in 1:n){
#     for (k in 1:n){
#       S1[j,k] = 2*a/((FH[j]-FH[k])^2+a^2)
#     }
#   }
#   S2 = rep(0,n)
#   for (j in 1:n){
#     S2[j] = atan(FH[j]/a)+atan((1-FH[j])/a)
#   }
#   out = 2*n*(2*atan(1/a)-a*log(1+1/a^2))+sum(S1)/n-4*sum(S2)
#   return(out)
# }
# 
# MEn <- function(X,a,betaH,sigmaH=1){
#   n  = length(X)
#   FH = ppareto(X,betaH,sigmaH)
#   S1 = rep(0,n)
#   k  = 1:n
#   for (j in 1:n){
#     S1[j] = sum(2*a/((FH[j]-FH)^2+a^2))
#   }
#   S2  = atan(FH[k]/a)+atan((1-FH[k])/a)
#   out = 2*n*(2*atan(1/a)-a*log(1+1/a^2))+sum(S1)/n-4*sum(S2)
#   return(out)
# }


Gn2 <- function(X,betaH,a){  #Rename it to G2 and get the order of arguments the same as the rest
  sigmaH = 1
  n      = length(X)
  FH     = ppareto(X,betaH,sigmaH)
  S1     = rep(0,n)
  k      = 1:n
  for (j in 1:n){
    S1[j] = sum(2*a/((FH[j]-FH)^2+a^2))
  }
  S2  = atan(FH[k]/a)+atan((1-FH[k])/a)
  out = 2*n*(2*atan(1/a)-a*log(1+1/a^2))+sum(S1)/n-4*sum(S2)
  return(out)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (12) Meintanis (2009b) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Ia0 <- function(X,a){
#   out = 1/(a+log(X))
# }
# 
# Ia1 <- function(X,a){
#   out = (1-a-log(X))/(a+log(X))^2
# }
# 
# Ia2 <- function(X,a){
#   out = (2-2*a+a^2+2*(a-1)*log(X)+(log(X))^2)/(a+log(X))^3
# }
# 
# M2nOld <- function(X,a,betaH,sigmaH=1){
#   n  = length(X)
#   X  = X/sigmaH
#   S1 = matrix(0,n,n)
#   for (j in 1:n){
#     for(k in 1:n){
#       S1[j,k] = Ia0(X[j]*X[k],a)
#     }
#   }
#   S2 = matrix(0,n,n)
#   for (j in 1:n){
#     for(k in 1:n){
#       S2[j,k] = Ia2(X[j]*X[k],a)
#     }
#   }
#   S3 = matrix(0,n,n)
#   for (j in 1:n){
#     for(k in 1:n){
#       S3[j,k] = Ia1(X[j]*X[k],a)
#     }
#   }
#   S4 = rep(0,n)
#   for (j in 1:n){
#     S4[j] = Ia0(X[j],a)
#   }
#   S5 = rep(0,n)
#   for (j in 1:n){
#     S5[j] = Ia1(X[j],a)
#   }
#   T1  = ((betaH+1)^2*sum(S1)+sum(S2)+2*(betaH+1)*sum(S3))/n
#   T2  = (n*betaH*Ia0(1,a)-2*(betaH+1)*sum(S4)-2*sum(S5))*betaH
#   out = T1+T2
# }
# 
# M2n <- function(X,a,betaH,sigmaH=1){
#   n  = length(X)
#   X  = X/sigmaH
#   S1 = rep(0,n)
#   for (j in 1:n){
#     S1[j] = sum(Ia0(X[j]*X,a))
#   }
#   S2 = rep(0,n)
#   for (j in 1:n){
#     S2[j] = sum(Ia2(X[j]*X,a))
#   }
#   S3 = rep(0,n)
#   for (j in 1:n){
#     S3[j] = sum(Ia1(X[j]*X,a))
#   }
#   S4  = Ia0(X,a)
#   S5  = Ia1(X,a)
#   T1  = ((betaH+1)^2*sum(S1)+sum(S2)+2*(betaH+1)*sum(S3))/n
#   T2  = (n*betaH*Ia0(1,a)-2*(betaH+1)*sum(S4)-2*sum(S5))*betaH
#   out = T1+T2
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~ (13) Obradovic et al. (2015) ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# OJnMarius <- function(X){
#   n    = length(X)
#   N    = choose(n,2)
#   tmp0 = combn(X,2)
#   tmp1 = sum(rank(c(X,apply(cbind(tmp0[1,]/tmp0[2,],tmp0[2,]/tmp0[1,]),1,max)))[1:n])
#   tmp3 = 1/(2*n*N)*(2*tmp1-(n+1)*(N+n))
#   ans  = abs(tmp3)
#   return(ans)
# }
# 
# MnOld <- function(X,t){
#   n = length(X)
#   S = matrix(0,n,n)
#   for (j in 1:(n-1)){
#     for (k in (j+1):n){
#       S[j,k] = (max(X[j]/X[k],X[k]/X[j])<=t)
#     }
#   }
#   out = sum(S)/choose(n,2)
#   return(out)
# }
# 
# Mn <- function(X,t){
#   n = length(X)
#   S = rep(0,n)
#   for (j in 1:(n-1)){
#     k    = (j+1):n
#     S[j] = sum(pmax(X[j]/X[k],X[k]/X[j])<=t)
#   }
#   out = sum(S)/choose(n,2)
#   return(out)
# }
# 
# MnVec <- function(X){
#   n = length(X)
#   S = matrix(0,n,n)
#   for (j in 1:(n-1)){
#     k = (j+1):n
#     for (l in 1:n){
#       S[j,l] = sum(pmax(X[j]/X[k],X[k]/X[j])<=X[l])
#     }
#   }
#   out = colSums(S)/choose(n,2)
#   return(out)
# }
# 
# OJn <- function(X){
#   n   = length(X)
#   Mn  = MnVec(X)
#   out = mean(Mn-(1:n)/n)
#   return(abs(out))
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ (14) Allison et al. (1) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Delta2Old <- function(n){
#   out = rep(0,n)
#   S   = matrix(0,n,n)
#   for (l in 1:n){
#     for (j in 1:n){
#       for (k in 1:n){
#         S[j,k] = (min(j,k)<=l)
#       }
#     }
#     out[l] = sum(S)
#   }
#   return(out)
# }
# 
# Delta2 <- function(n){
#   out = rep(0,n)
#   S   = rep(0,n)
#   for (l in 1:n){
#     for (j in 1:n){
#       S[j] = sum((pmin(j,1:n)<=l))
#     }
#     out[l] = sum(S)
#   }
#   return(out)
# }
# 
# Delta3 <- function(n){
#   out = rep(0,n)
#   S   = matrix(0,n,n)
#   for (l in 1:n){
#     for (j in 1:n){
#       for (k in 1:n){
#         S[j,k] = sum(pmin(j,k,1:n)<=l)
#       }
#     }
#     out[l] = sum(S)
#   }
#   return(out)
# }
# 
# A1nOld <- function(X,m){
#   n = length(X)
#   if(m == 2){ Delta = Delta2(n) }
#   if(m == 3){ Delta = Delta3(n) }
#   S1 = rep(0,n)
#   S  = rep(0,n)
#   for (l in 1:n){
#     for (j in 1:n){
#       S1[j] = (X[j]^(1/m)<X[l])
#     }
#     S[l] = mean(S1)
#   }
#   out = mean(S-Delta/n^m)
#   return(out)
# }
# 
# A1n <- function(X,m){
#   n = length(X)
#   if(m == 2){ Delta = Delta2(n) }
#   if(m == 3){ Delta = Delta3(n) }
#   S = rep(0,n)
#   for (l in 1:n){
#     S[l] = mean(X^(1/m)<X[l])
#   }
#   out = mean(S-Delta/n^m)
#   return(abs(out))
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ (15) Allison et al. (2) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# A2n <- function(X,m){
#   n = length(X)
#   if(m == 2){ Delta = Delta2(n) }
#   if(m == 3){ Delta = Delta3(n) }
#   S = rep(0,n)
#   for (l in 1:n){
#     S[l] = mean(X^(1/m)<X[l])
#   }
#   out = mean((S-Delta/n^m)^2)
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~ (16) Obradovic (1) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Gn <- function(x,X){
#   gx <- function(x,Xj,Xk,Xl){
#     md  = median(c(Xj,Xk,Xl))
#     mn  = min(c(Xj,Xk,Xl))
#     out = (md/mn<=x)
#     return(out)
#   }
#   n  = length(X)
#   S1 = matrix(0,n,n)
#   for (j in 1:(n-2)){
#     for (k in (j+1):(n-1)){
#       S1[j,k] = gx(x,X[j],X[k],X[n])*(n-k)
#     }
#   }
#   S1  = sum(S1)
#   S2a = matrix(0,n,n)
#   for (k in 1:(n-1)){
#     for (l in (k+1):n){
#       S2a[k,l] = gx(x,X[k],X[l],X[l])
#     }
#   }
#   S2a = sum(S2a)
#   S2b = matrix(0,n,n)
#   for (l in 1:(n-1)){
#     for (k in (l+1):n){
#       S2b[k,l] = gx(x,X[k],X[l],X[l])
#     }
#   }
#   S2b = sum(S2b)
#   S3  = rep(0,n)
#   for (l in 1:n){
#     S3[l] = gx(x,X[l],X[l],X[l])
#   }
#   S3  = sum(S3)
#   out = (6*S1+3*(S2a+S2b)+S3)/n^3
#   return(out)
# }
# 
# Hn <- function(n){
#   out = rep(0,n)
#   S   = rep(0,n)
#   for (l in 1:n){
#     for (j in 1:n){
#       S[j] = sum((pmin(j,1:n)<=l))
#     }
#     out[l] = sum(S)/n^2
#   }
#   return(out)
# }
# 
# O1n <- function(X){
#   n = length(X)
#   G = rep(0,n)
#   for(j in 1:n){
#     G[j] = Gn(X[j],X)
#   }
#   out = mean(G-Hn(n))
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~ (17) Obradovic (2) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Jn <- function(x,X){
#   jx <- function(x,Xj,Xk,Xl){
#     mx  = max(c(Xj,Xk,Xl))
#     md  = median(c(Xj,Xk,Xl))
#     out = (mx/md<=x)
#     return(out)
#   }
#   n  = length(X)
#   S1 = matrix(0,n,n)
#   for (l in 3:n){
#     for (k in 2:(l-1)){
#       S1[k,l] = (X[l]/X[k]<x)*(k-1)
#     }
#   }
#   S1  = sum(S1)
#   S2a = matrix(0,n,n)
#   for (k in 1:(n-1)){
#     for (l in (k+1):n){
#       S2a[k,l] = jx(x,X[k],X[l],X[l])
#     }
#   }
#   S2a = sum(S2a)
#   S2b = matrix(0,n,n)
#   for (l in 1:(n-1)){
#     for (k in (l+1):n){
#       S2b[k,l] = jx(x,X[k],X[l],X[l])
#     }
#   }
#   S2b = sum(S2b)
#   S3  = rep(0,n)
#   for (l in 1:n){
#     S3[l] = jx(x,X[l],X[l],X[l])
#   }
#   S3  = sum(S3)
#   out = (6*S1+3*(S2a+S2b)+S3)/n^3
#   return(out)
# }
# 
# Kn <- function(x,X){
#   kx <- function(x,Xj,Xk,Xl){
#     md  = median(c(Xj,Xk,Xl))
#     mn  = (min(c(Xj,Xk,Xl)))
#     out = ((md/mn)^2<=x)
#     return(out)
#   }
#   n  = length(X)
#   S1 = matrix(0,n,n)
#   for (j in 1:(n-2)){
#     for (k in (j+1):(n-1)){
#       S1[j,k] = kx(x,X[j],X[k],X[n])*(n-k)
#     }
#   }
#   S1  = sum(S1)
#   S2a = matrix(0,n,n)
#   for (k in 1:(n-1)){
#     for (l in (k+1):n){
#       S2a[k,l] = kx(x,X[k],X[l],X[l])
#     }
#   }
#   S2a = sum(S2a)
#   S2b = matrix(0,n,n)
#   for (l in 1:(n-1)){
#     for (k in (l+1):n){
#       S2b[k,l] = kx(x,X[k],X[l],X[l])
#     }
#   }
#   S2b = sum(S2b)
#   S3  = rep(0,n)
#   for (l in 1:n){
#     S3[l] = kx(x,X[l],X[l],X[l])
#   }
#   S3  = sum(S3)
#   out = (6*S1+3*(S2a+S2b)+S3)/n^3
#   return(out)
# }
# 
# O2n <- function(X){
#   n   = length(X)
#   out = rep(0,n)
#   for(j in 1:n){
#     out[j] = Jn(X[j],X)-Kn(X[j],X)
#   }
#   out = mean(out)
#   return(out)
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
