# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRACTICAL DATA APPLICATION EXAMPLE FOR:
#
# New classes of tests for the Pareto distribution based on a conditional expectation characterisation
#
# T. Nombebe, J.S. Allison, L. Santana and I.J.H. Visagie
#
# 2023
#
# SEE PAGE 14
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


graphics.off()
cat("\014")
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("2023-09-11-Estimation.R")
source("2023-09-11-Test-Statistics_GoF_Pareto.R")
source("2023-09-11-DistsForSim.R")
set.seed(11803371)
B       = 1e5
ntests  = 10
alpha   = 0.1

#  EPL 2021 to 2022 DataSet
X1 <- EPL21to22 <- c( 26800000,20800000,19500000,18200000,17680000,16900000,15600000,15600000,15080000,15080000,14144000,13000000,13000000,11440000,10400000,10400000,10400000,10400000,10400000,10310000)
#EPL 2022 to 2023 DataSet
X2 <- EPL22to23 <- c(20800000,19500000,19500000,18200000,18200000,17680000,16900000,15600000,15600000,15470000,15340000,15080000,13780000,13000000,13000000,13000000,12480000,11700000,11440000,11440000,10400000,10400000,10400000,10400000,10400000,10400000,10400000,10400000)
X1 <- sort(X1)/1e7
X2 <- sort(X2)/1e7
X <- data.frame(Earnings=c(X1,X2),Date=c(rep("2021-2022",length(X1)),rep("2022-2023",length(X2))))
n <- table(X$Date)

pdf(file = paste0("TN3-PracData-Boxplots-",format(Sys.time(), "%Y%m%d-%Hh%M"),".pdf"),
    width=11, height=9)
boxplot(X$Earnings~factor(X$Date),main="",ylab="Earnings (GBP)",xlab="")
dev.off()


pdf(file = paste0("TN3-PracData-Violinplots-",format(Sys.time(), "%Y%m%d-%Hh%M"),".pdf"),
    width=11, height=9)
library("ggplot2")
ggplot(X,aes(x=factor(Date),y=Earnings)) +
  geom_hline(yintercept=1,lty=2)+
  geom_boxplot(width=0.35,aes(fill=factor(Date)))+
  geom_violin(width=0.5,aes(fill=factor(Date),alpha=0.25))+
  xlab("")+
  ylab("(Rescaled) Earnings (GBP)")+
  ggtitle("")+
  labs(fill = "Users By labs")+
  guides(fill="none",alpha="none")+
  theme_classic()+
  theme_update(axis.title=element_text(size=26),
        axis.text=element_text(size=22))
dev.off()

betahat.ML <- betahat.MM <- numeric(2)
betahat.ML[1] <- 1/mean(log(X1))
betahat.ML[2] <- 1/mean(log(X2))
betahat.MM[1] <- mean(X1)/(mean(X1) - 1)
betahat.MM[2] <- mean(X2)/(mean(X2) - 1)

pPareto <- function(x,betahat){  1 - x^(-betahat)}
qPareto <- function(p,betahat){  (1 - p)^(-1/betahat)}

MLE.line <- c("grey0",3,1) #col, lwd, lty
MME.line <- c("grey40",3,2) #col, lwd, lty



pdf(file = paste0("TN3-PracData-EDFplots-",format(Sys.time(), "%Y%m%d-%Hh%M"),".pdf"),
    width=11, height=9)
par(mfcol=c(1,2))
par(mar=c(4.2,4.6,0.1,0.5))
plot(ecdf(X1),main="",xlab="Scaled EPL earnings, 2021-2022 (GBP)",ylab=expression(P(X < x)),cex.lab=1.4)
xseq <- seq(1, max(X1)*1.05,length=1e3)
lines(xseq,pPareto(xseq,betahat.ML[1]),
      col=MLE.line[1],
      lwd=as.numeric(MLE.line[2]),
      lty=as.numeric(MLE.line[3]))
lines(xseq,pPareto(xseq,betahat.MM[1]),
      col=MME.line[1],
      lwd=as.numeric(MME.line[2]),
      lty=as.numeric(MME.line[3]))
s1 <- bquote(MLE~(hat(beta)[MLE]==~.(round(betahat.ML[1],2))) )
s2 <- bquote(MME~(hat(beta)[MME]==~.(round(betahat.MM[1],2))) )
legend("bottomright", legend=c("EDF",
                           s1, # paste0("MLE",round(betahat.ML[1],2)),
                           s2),  # paste0("MME",round(betahat.MM[1],2))
       cex=1.1,
       box.col="white", 
       inset=c(0.05,0.05),
       col = c("black",MLE.line[1],MME.line[1]), 
       lwd = c(1,as.numeric(MLE.line[2]),as.numeric(MME.line[2])),
       lty = c(1,as.numeric(MLE.line[3]),as.numeric(MME.line[3])),
       pch=c(16,NA,NA))

plot(ecdf(X2),main="",xlab="Scaled EPL earnings, 2022-2023 (GBP)",ylab=expression(P(X < x)),cex.lab=1.4)
xseq <- seq(1, max(X2)*1.05,length=1e3)
lines(xseq,pPareto(xseq,betahat.ML[2]),
      col=MLE.line[1],
      lwd=as.numeric(MLE.line[2]),
      lty=as.numeric(MLE.line[3]))
lines(xseq,pPareto(xseq,betahat.MM[2]),
      col=MME.line[1],
      lwd=as.numeric(MME.line[2]),
      lty=as.numeric(MME.line[3]))

s3 <- bquote(MLE~(hat(beta)[MLE]==~.(round(betahat.ML[2],2))) )
s4 <- bquote(MME~(hat(beta)[MME]==~.(round(betahat.MM[2],2))) )
legend("bottomright", legend=c("EDF",
                           s3, # paste0("MLE",round(betahat.ML[1],2)),
                           s4 ), # paste0("MME",round(betahat.MM[1],2))
       cex=1.1,
       box.col="white" , 
       inset=c(0.05,0.05),
       col = c("black",MLE.line[1],MME.line[1]), 
       lwd = c(1,as.numeric(MLE.line[2]),as.numeric(MME.line[2])),
       lty = c(1,as.numeric(MLE.line[3]),as.numeric(MME.line[3])),
       pch=c(16,NA,NA))
dev.off()


#-----------------------------------------------------









(TS_KS1_   = c(KSn(X1,betahat.ML[1]),KSn(X1,betahat.MM[1])))
(TS_KS2_  = c(KSn(X2,betahat.ML[2]),KSn(X2,betahat.MM[2])))

(TS_CV1_   = c(CVn(X1,betahat.ML[1]),CVn(X1,betahat.MM[1])))
(TS_CV2_   = c(CVn(X2,betahat.ML[2]),CVn(X2,betahat.MM[2])))

(TS_AD1_   = c(ADn(X1,betahat.ML[1]),ADn(X1,betahat.MM[1])))
(TS_AD2_   = c(ADn(X2,betahat.ML[2]),ADn(X2,betahat.MM[2])))

(TS_MA1_   = c(MAn(X1,betahat.ML[1]),MAn(X1,betahat.MM[1])))
(TS_MA2_   = c(MAn(X2,betahat.ML[2]),MAn(X2,betahat.MM[2])))

(TS_ZA1_   = c(ZAn(X1,betahat.ML[1]),ZAn(X1,betahat.MM[1])))
(TS_ZA2_   = c(ZAn(X2,betahat.ML[2]),ZAn(X2,betahat.MM[2])))

(TS_G21_   = c(Gn2(X1,betahat.ML[1],a=2),Gn2(X1,betahat.MM[1],a=2)))
(TS_G22_   = c(Gn2(X2,betahat.ML[2],a=2),Gn2(X2,betahat.MM[2],a=2)))

(TS_NA1_   = c(Sn2(X1,betahat.ML[1],m=3,a=2),Sn2(X1,betahat.MM[1],m=3,a=2)))
(TS_NA2_   = c(Sn2(X2,betahat.ML[2],m=3,a=2),Sn2(X2,betahat.MM[2],m=3,a=2)))

(TS_L11_   = c(Tns1(X1,betahat.ML[1],s=0.1,a=3),Tns1(X1,betahat.MM[1],s=0.1,a=3)))
(TS_L12_   = c(Tns1(X2,betahat.ML[2],s=0.1,a=3),Tns1(X2,betahat.MM[2],s=0.1,a=3)))

(TS_L21_   = c(Tns2(X1,betahat.ML[1],s=0.1,a=3),Tns2(X1,betahat.MM[1],s=0.1,a=3)))
(TS_L22_   = c(Tns2(X2,betahat.ML[2],s=0.1,a=3),Tns2(X2,betahat.MM[2],s=0.1,a=3)))

(TS_T11_   = c(Tns3(X1,betahat.ML[1],s=0.1,a=3),Tns3(X1,betahat.MM[1],s=0.1,a=3)))
(TS_T12_   = c(Tns3(X2,betahat.ML[2],s=0.1,a=3),Tns3(X2,betahat.MM[2],s=0.1,a=3)))

pvalcalc <- function(MAT,TS)
  c(mean(MAT[,1] > TS[1]),mean(MAT[,2] > TS[2]))

p1  = p2 = matrix(numeric(ntests*2),ncol=2) 
TS_KS1   = TS_KS2 =
  TS_CV1 = TS_CV2 =
  TS_AD1 = TS_AD2 =
  TS_MA1 = TS_MA2 =
  TS_ZA1 = TS_ZA2 =
  TS_G21 = TS_G22 =
  TS_NA1 = TS_NA2 =
  TS_L11 = TS_L12 =
  TS_L21 = TS_L22 =
  TS_T11 = TS_T12 = matrix(0,B,2)

betaHML <- numeric(2)
betaHMM <- numeric(2)

for (j in 1:B){#START FOR LOOP (j in 1:B)
  XML1        = sort(rpareto(n[1],betahat.ML[1]))
  XML2        = sort(rpareto(n[2],betahat.ML[2]))
  
  XMM1        = sort(rpareto(n[1],betahat.MM[1]))
  XMM2        = sort(rpareto(n[2],betahat.MM[2]))
  
  betaHML[1]  = ParetoMLE1(XML1)$betaH
  betaHML[2]  = ParetoMLE1(XML2)$betaH 
  
  betaHMM[1]  = ParetoMME1(XMM1)$betaH
  betaHMM[2]  = ParetoMME1(XMM2)$betaH
  
  TS_KS1[j,]  = c(KSn( XML1, betaHML[1]),          KSn( XMM1, betaHMM[1]))
  TS_CV1[j,]  = c(CVn( XML1, betaHML[1]),          CVn( XMM1, betaHMM[1]))
  TS_AD1[j,]  = c(ADn( XML1, betaHML[1]),          ADn( XMM1, betaHMM[1]))
  TS_MA1[j,]  = c(MAn( XML1, betaHML[1]),          MAn( XMM1, betaHMM[1]))
  TS_ZA1[j,]  = c(ZAn( XML1, betaHML[1]),          ZAn( XMM1, betaHMM[1]))
  TS_G21[j,]  = c(Gn2( XML1, betaHML[1],a=2),      Gn2( XMM1, betaHMM[1],a=2))
  TS_NA1[j,]  = c(Sn2( XML1, betaHML[1],m=3,a=2),  Sn2( XMM1, betaHMM[1],m=3,a=2))
  TS_L11[j,]  = c(Tns1(XML1, betaHML[1],s=0.1,a=3),Tns1(XMM1, betaHMM[1],s=0.1,a=3))
  TS_L21[j,]  = c(Tns2(XML1, betaHML[1],s=0.1,a=3),Tns2(XMM1, betaHMM[1],s=0.1,a=3))
  TS_T11[j,]  = c(Tns3(XML1, betaHML[1],s=0.1,a=3),Tns3(XMM1, betaHMM[1],s=0.1,a=3))
  
  TS_KS2[j,]  = c(KSn( XML2, betaHML[2]),          KSn( XMM2, betaHMM[2]))
  TS_CV2[j,]  = c(CVn( XML2, betaHML[2]),          CVn( XMM2, betaHMM[2]))
  TS_AD2[j,]  = c(ADn( XML2, betaHML[2]),          ADn( XMM2, betaHMM[2]))
  TS_MA2[j,]  = c(MAn( XML2, betaHML[2]),          MAn( XMM2, betaHMM[2]))
  TS_ZA2[j,]  = c(ZAn( XML2, betaHML[2]),          ZAn( XMM2, betaHMM[2]))
  TS_G22[j,]  = c(Gn2( XML2, betaHML[2],a=2),      Gn2( XMM2, betaHMM[2],a=2))
  TS_NA2[j,]  = c(Sn2( XML2, betaHML[2],m=3,a=2),  Sn2( XMM2, betaHMM[2],m=3,a=2))
  TS_L12[j,]  = c(Tns1(XML2, betaHML[2],s=0.1,a=3),Tns1(XMM2, betaHMM[2],s=0.1,a=3))
  TS_L22[j,]  = c(Tns2(XML2, betaHML[2],s=0.1,a=3),Tns2(XMM2, betaHMM[2],s=0.1,a=3))
  TS_T12[j,]  = c(Tns3(XML2, betaHML[2],s=0.1,a=3),Tns3(XMM2, betaHMM[2],s=0.1,a=3))
  
}#END FOR LOOP (j in 1:B)

TS_names <- c(
  "KS",
  "CV",
  "AD",
  "MA",
  "ZA",
  "G2",
  "NA",
  "L1",
  "L2",
  "T1")

dimnames(p1) <- list(TS_names,c("MLE '21-'22","MME '21-'22"))
dimnames(p2) <- list(TS_names,c("MLE '22-'23","MME '22-'23"))

p1[1,]  = pvalcalc(TS_KS1,TS_KS1_)
p1[2,]  = pvalcalc(TS_CV1,TS_CV1_)
p1[3,]  = pvalcalc(TS_AD1,TS_AD1_)
p1[4,]  = pvalcalc(TS_MA1,TS_MA1_)
p1[5,]  = pvalcalc(TS_ZA1,TS_ZA1_)
p1[6,]  = pvalcalc(TS_G21,TS_G21_)
p1[7,]  = pvalcalc(TS_NA1,TS_NA1_)
p1[8,]  = pvalcalc(TS_L11,TS_L11_)
p1[9,]  = pvalcalc(TS_L21,TS_L21_)
p1[10,] = pvalcalc(TS_T11,TS_T11_)

p2[1,]  = pvalcalc(TS_KS2,TS_KS2_)
p2[2,]  = pvalcalc(TS_CV2,TS_CV2_)
p2[3,]  = pvalcalc(TS_AD2,TS_AD2_)
p2[4,]  = pvalcalc(TS_MA2,TS_MA2_)
p2[5,]  = pvalcalc(TS_ZA2,TS_ZA2_)
p2[6,]  = pvalcalc(TS_G22,TS_G22_)
p2[7,]  = pvalcalc(TS_NA2,TS_NA2_)
p2[8,]  = pvalcalc(TS_L12,TS_L12_)
p2[9,]  = pvalcalc(TS_L22,TS_L22_)
p2[10,] = pvalcalc(TS_T12,TS_T12_)
cbind(p1,p2)

rp1 <- round(p1*100,0)
rp2 <- round(p2*100,0)
alph<- alpha*100
rgb <- c(0.749,0.749,0.749)
txt <- "
\\begin{table}[!htbp!] 
  \\centering 
  \\caption{Practical data example: $p$-values for testing Paretoness of the EPL earnings above 10 million GBP for the 2021-2022 and 2022-2023 seasons (values less than $\\alpha=0.1$ are highlighted).\\label{tbl:pracdata}} 
      \\begin{tabular}{l|cc|cc} 
\\cline{2-5} 
\\multicolumn{1}{c}{}& \\multicolumn{2}{c|}{2021--2022} & \\multicolumn{2}{c}{2022--2023} \\\\ 
\\cline{2-5} 
\\multicolumn{1}{c}{}&   MLE      &   MME      &   MLE     &    MME \\\\
\\hline  "
for(ii in 1:nrow(rp1)){
  txt <- paste0(txt,"$",TS_names[ii],"$ & ", 
                ifelse(rp1[ii,1]<alph,paste0(" \\cellcolor[rgb]{",rgb[1],",",rgb[2],",",rgb[3],"}",rp1[ii,1]),rp1[ii,1]) ," & ", 
                ifelse(rp1[ii,2]<alph,paste0(" \\cellcolor[rgb]{",rgb[1],",",rgb[2],",",rgb[3],"}",rp1[ii,2]),rp1[ii,2]) ," & ", 
                ifelse(rp2[ii,1]<alph,paste0(" \\cellcolor[rgb]{",rgb[1],",",rgb[2],",",rgb[3],"}",rp2[ii,1]),rp2[ii,1]) ," & ", 
                ifelse(rp2[ii,2]<alph,paste0(" \\cellcolor[rgb]{",rgb[1],",",rgb[2],",",rgb[3],"}",rp2[ii,2]),rp2[ii,2]) , "\\\\ \n")
}  
txt <- paste0(txt,"\\hline
\\end{tabular}
\\end{table} ")
cat(txt)




