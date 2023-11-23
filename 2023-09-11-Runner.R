# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SIMULATION PROGRAM FOR ~~~~~~~~~~~~~~~
#
# New classes of tests for the Pareto distribution based on a 
# conditional expectation characterisation
#
# T. Nombebe, J.S. Allison, L. Santana and I.J.H. Visagie
#
# 2023
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial settings ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
cat("\014")
set.seed(1234)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(stargazer)
source("2023-09-11-Estimation.R")
source("2023-09-11-Test-Statistics_GoF_Pareto.R")
source("2023-09-11-DistsForSim.R")

starttime <- format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss")

#Calculate the power via WarpSpeed
powersWS <- function(TS,alpha){
  TS_X = TS[,1]
  TS_S = TS[,2]
  cv   = quantile(TS_S,1-alpha)
  out  = mean(TS_X>cv)
  return(out)
}

#evaluate strings as functions
my.eval  <- function(str) eval(str2expression(str))
my.eval.vec  <- function(str) eval(str2expression(str))
my.eval.vec  <- Vectorize(my.eval.vec)

#Calculate the power (plain)
powersPlain <- function(TS,alpha,cv){
  out  = mean(TS>cv)
  return(out)
}

#Count how many tests there are, unpack the names and parameter settings.
unpack.test <- function(parms){
  namer <- NULL
  parm0 <- NULL
  ntests <- 0
  for(ii in 1:length(parms)){
    if(is.list(parms[[ii]])){
      tmp <- prod(unlist(lapply(parms[[ii]],length)))
      namer[(ntests+1):(ntests+tmp)] <- names(parms[ii])
      ntests <- ntests + tmp
      parm0 <- c(parm0,as.list(as.data.frame(t(expand.grid(parms[[ii]])))))
    }else if(is.null(parms[[ii]])){ 
      ntests <- ntests + 1
      namer[ntests] <- names(parms[ii])
      parm0 <- c(parm0,list(NULL))
    } else{
      warning("Parameter settings must be a list or a NULL.")
    }
  }
  names(parm0) <- NULL
  parm0text    <- lapply(parm0,function(x){ if(is.null(x)) return("") else return(paste0(": (", paste0(x,collapse=", "),")"))})
  return(list(ntests=ntests,namer=namer,parm=parm0,parmtext=parm0text))
}

#Estimate the parameters using any of the methods. Needs edits!
EstParm <- function(X,meth) switch(meth,
                                   MME        =ParetoMME1(         X),
                                   MLE        =ParetoMLE1(         X),
                                   MDELS      =ParetoMDE.LS(       X),
                                   MDECvM     =ParetoMDE.CvM(      X),
                                   MDEPhikl   =ParetoMDE.Phi.kl(   X),
                                   MDEPhichisq=ParetoMDE.Phi.chisq(X),
                                   MDEPhitv   =ParetoMDE.Phi.tv(   X)
                                   )$betaH

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample size
n        = 30 #20 #30
# distributions used
Range    = c(2,4,5,6,7,8,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40) #1:16 #2:36 #1:4 use pareto, gamma, weibull, lognormal, lfr
# estimatino method used
meth     = "MME" #"MLE" #"MME" #"MLE" or "MDELS" or "MDECvM" or "MDEPhikl" or "MDEPhichisq" or "MDEPhitv" or something else (see the function "EstParm")
# number of simulations for power
MC       = 1e2#1e4
# number of simulations for critical value calculation
MCC      = 5e2#5e4
# significance value for tests
alpha    = 0.05
# Names of the test statistic functions:
tst.names= c(
   "KSn" ,
   "CVn" ,
   "ADn" ,
   "MAn" ,
   "ZAn" ,
  # "ZBn" ,
  # "ZCn" ,
  # "DKn" ,
  # "KLn" ,
  "Tns1",
  "Tns2",
  "Tns3",
  "Sn2",
  #"Tns4",
  #"Tns5",
  #"Tns6",
  "Gn2"
)
# Additional Parameter settings for each test (MUST be a list or a NULL):
parms    = list( 
  NULL,                                              #KSn
  NULL,                                              #CVn
  NULL,                                              #ADn
  NULL,                                              #MAn
  NULL,                                              #ZAn
  # NULL,                                            #ZBn
  # NULL,                                            #ZCn
  # NULL,                                            #DKn
  # list(c(1,10)),                                   #KLn
  list(sort(unique(c(3))),sort(unique(c(0.1)))),     #Tns1
  list(sort(unique(c(3))),sort(unique(c(0.1)))),     #Tns2
  list(sort(unique(c(3))),sort(unique(c(0.1)))),     #Tns3
  list(sort(unique(c(3))),sort(unique(c(2)))),       #Sn2  (m and a)
  #list(sort(unique(c(3))),sort(unique(c(0.25)))),   #Tns5
  #list(sort(unique(c(3))),sort(unique(c(0.25))))    #Tns6
  list(sort(unique(c(2))))                           #Gn2  (a)
  )   
names(parms) = tst.names
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
timestarted  <- format(Sys.time(),"%Y-%m-%d-%HH%MM%SS")
filenamecsv  <- paste0("Pareto_",timestarted,"_n=",n,"_MC=",MC,"_meth=",meth,"_dists=",paste0(Range,collapse=","),".csv")
filenamerdat <- paste0("Pareto_",timestarted,"_n=",n,"_MC=",MC,"_meth=",meth,"_dists=",paste0(Range,collapse=","),".Rdata")
# filenamecsv  <- paste0("Pareto_",timestarted,"_n=",n,"_MC=",MC,"_meth=",meth,"_dists=",min(Range),"-",max(Range),".csv")
# filenamerdat <- paste0("Pareto_",timestarted,"_n=",n,"_MC=",MC,"_meth=",meth,"_dists=",min(Range),"-",max(Range),".Rdata")
fileCV       <- paste0("MLE_CV_n=",n,"_MCC=",MCC,"_alpha=",alpha,"_tsts=",paste0(tst.names,collapse=","),".txt")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unpack   = unpack.test(parms)
ntests   = unpack$ntests   # number of tests
newname  = unpack$namer    # names of tests (repetition for parameter combinations) : length = ntests
newparm  = unpack$parm     # parameter combinations: length = ntests: Vector of parameters
newparmtx= unpack$parmtext # parameter combinations: length = ntests: TEXT VERSION FOR PRINTING

Results  = matrix(nrow=length(Range),ncol=ntests+1)  # (0, number of distributions, number of tests + 1)
dimnames(Results) = list(paste("dist=",Range),c("",paste0(newname,newparmtx)))

nm    = numeric(length(Range))
tme1  = proc.time()

count = 0

pb    = txtProgressBar(min=0,max=length(Range),style=3)

DistNumcntr <- 1
for (DistNum in Range){#START FOR LOOP (DistNum in Range)
  p            = numeric(ntests) 
  if(meth != "MLE"){ #START IF MME vs MLE vs OTHERS
    TS           = array(dim=c(MC,2,ntests))
    dimnames(TS) = list(1:MC,c("X","Xstar"),paste0(newname,newparmtx))  
    
    for (j in 1:MC){#START FOR LOOP (j in 1:MC)
      simdata         = SimFromDist(DistNum,n)
      X               = sort(simdata$X)
      nm[DistNumcntr] = simdata$nm
      betaH           = EstParm(X,meth) 
      Xstar           = sort(rpareto(n,betaH))
      betaHS          = EstParm(Xstar,meth)
      
      for(k in 1:ntests){
        if(is.null(newparm[k][[1]])){
          TS[j,,k] = c(my.eval(newname[k])(    X,betaH),                                        
                       my.eval(newname[k])(Xstar,betaHS))                     #Calculate test statistics that have no additional parameters
        }else if(length(newparm[k][[1]])==1){
          TS[j,,k] = c(my.eval(newname[k])(    X,betaH ,newparm[[k]][1]),
                       my.eval(newname[k])(Xstar,betaHS,newparm[[k]][1]))     #Calculate test statistics that have ONE additional parameter
        }else if(length(newparm[k][[1]])==2){
          TS[j,,k] = c(my.eval(newname[k])(    X,betaH ,newparm[[k]][1],newparm[[k]][2]),
                       my.eval(newname[k])(Xstar,betaHS,newparm[[k]][1],newparm[[k]][2]))  #Calculate test statistics that have TWO additional parameters
        }  
      }  
    }#END FOR LOOP (j in 1:MC)
    for(k in 1:ntests)
      p[k]  = powersWS(TS[,,k],alpha)
    
  }else{ # ELSE IF MLE ESTIMATION
    
    
    # ----------
    #CHECK IF THE CRITICAL VALUE CSV FILE EXISTS. IF YES, THEN SKIP THE CV CALCULATION. IF NO, CALCULATE THE CV.
    suppressWarnings(
      cvs <- try(dget(fileCV),silent=TRUE)
    )
    if(inherits(cvs, "try-error")){ #START IF CRITICAL VALUE CSV FILE EXISTS OR NOT
      TS.cv           = array(dim=c(MCC,ntests))
      dimnames(TS.cv) = list(1:MCC,paste0(newname,newparmtx))  
      cvs             = numeric(ntests)
      for (j in 1:MCC){#START FOR LOOP (j in 1:MC)
        simdata         = SimFromDist(4,n)
        X               = sort(simdata$X)
        betaH           = EstParm(X,meth) 
        Y               = X^betaH 
        for(k in 1:ntests){
          if(is.null(newparm[k][[1]])){
            TS.cv[j,k] = my.eval(newname[k])(Y,1)                                  #Calculate test statistics that have no additional parameters
          }else if(length(newparm[k][[1]])==1){ 
            TS.cv[j,k] = my.eval(newname[k])(Y,1 ,newparm[[k]][1])                 #Calculate test statistics that have ONE additional parameter
          }else if(length(newparm[k][[1]])==2){
            TS.cv[j,k] = my.eval(newname[k])(Y,1 ,newparm[[k]][1],newparm[[k]][2]) #Calculate test statistics that have TWO additional parameters
          }  
        }  
      }#END FOR LOOP (j in 1:MC)
      
      for(k in 1:ntests){
        cvs[k]  = quantile(TS.cv[,k],1-alpha)
      }  
      
      names(cvs) <- paste0(newname,newparmtx)
      dput(cvs,file=fileCV) #CREATE THE CSV FILE CONTAINING THE CRITICAL VALUES
      
    } #END IF CRITICAL VALUE CSV FILE EXISTS OR NOT
    # ----------
    
    
    
    TS           = array(dim=c(MC,ntests))
    dimnames(TS) = list(1:MC,paste0(newname,newparmtx))  
    for (j in 1:MC){#START FOR LOOP (j in 1:MC)
      simdata         = SimFromDist(DistNum,n)
      X               = sort(simdata$X)
      nm[DistNumcntr] = simdata$nm
      betaH           = EstParm(X,meth) 
      Y               = X^betaH
      for(k in 1:ntests){
        if(is.null(newparm[k][[1]])){
          TS[j,k] = my.eval(newname[k])(Y,1)                                  #Calculate test statistics that have no additional parameters
        }else if(length(newparm[k][[1]])==1){ 
          TS[j,k] = my.eval(newname[k])(Y,1 ,newparm[[k]][1])                 #Calculate test statistics that have ONE additional parameter
        }else if(length(newparm[k][[1]])==2){
          TS[j,k] = my.eval(newname[k])(Y,1 ,newparm[[k]][1],newparm[[k]][2]) #Calculate test statistics that have TWO additional parameters
        }  
      }  
    }#END FOR LOOP (j in 1:MC)
    
    for(k in 1:ntests){
      p[k]  = powersPlain(TS[,k],alpha,cvs[k])
    }  
    
    
  }  #END IF MME vs MLE vs OTHERS
  
  Results[DistNumcntr,] = c(DistNum,p)
  count = count+1
  # ~~~~~~~~~~~~~~~ SAVE SAVE SAVE SAVE SAVE ~~~~~~~~~~~~~~~~~
  write.csv(Results,file=filenamecsv)
  save.image(filenamerdat)
  # ~~~~~~~~~~~~~~~ SAVE SAVE SAVE SAVE SAVE ~~~~~~~~~~~~~~~~~
  DistNumcntr <- DistNumcntr + 1
  setTxtProgressBar(pb,count)
}#END FOR LOOP (DistNum in Range)

Results

tme2 = proc.time()
tme  = tme2[3]-tme1[3]
################################################################


