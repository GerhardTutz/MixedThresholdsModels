


library("pracma")
library(splines)
library(tidyr)
library(spatstat) 
library("foreign") 
library("ordinal")


install.packages("Matrix")
install.packages("lme4")


library("lme4")
library("Matrix")



str(sleepstudy)



##### with thresholds


#str(sleepstudy)

dats<-sleepstudy

library(tidyr)

dataswide <- spread(dats, Days, Reaction)

datm<- t(dataswide[2:11])


#lin <- "lin"  
lin <- "log"
#lin <- "log1"
#lin <- "logit"

respf<- "NV"
#respf<- "logistic"
#respf<- "Gumbel"
#respf<- "Gompertz"

indicator<- "C"
#indicator<- "D"


### always run (1) and (2) after selection of lin and respfct   

## run: Difficultyfunctions.R
## run: Responsefunctions.R
source("Difficultyfunctions")
source("Responsefunctions")



####(1) define function for response function

if(respf=="NV"){respfctd <- respfctdNV
respfctc<-respfctcNV
respfctder<-respfctderNV}
if(respf=="logistic"){respfctd <- respfctdlogit
respfctc<-respfctclogit
respfctder<-respfctderlogit}
if(respf=="Gumbel"){respfctd <- respfctdGumbel
respfctc<-respfctcGumbel
respfctder<-respfctderGumbel}
if(respf=="Gompertz"){respfctd <- respfctdGompertz
respfctc<-respfctcGompertz
respfctder<-respfctderGompertz}

####(2) define function for difficulty function

if(lin =="lin") {diffunct<-diffunctlin
derdiffunct<-derdiffunctlin}
if(lin =="log") {diffunct<-diffunctlog
derdiffunct<-derdiffunctlog}
if(lin =="log1") {diffunct<-diffunctlog1
derdiffunct<-derdiffunctlog1}



### fits with basis program 

## run: ProgramsBSplinesPureResp.R
## run: ProgramsBSplinesPenalties.R
source("ProgramsBSplinesPureResp")
source("ProgramsBSplinesPenalties")


I<- dim(datm)[1]
indicatorsingle  <- "C"  #### discrete or continuous
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)

cov0<-matrix(0,1,18)


### with log y

dum<-c(-110,20,0)
par<-c(rep(dum,I),1)

ThrS0 <-ThreshbasisRespCov(datm,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=par, 
                           startind = "yes",der="der",cov=cov0,hessiander="FALSE")




#### Plots for  estimates

parmatrest <- ThrS0$parmatrix[,1:2] ###
dat<-datm 

min <- min(dat)-0.5
max <-  max(dat)+0.5
max <-500

pcdum <-c ('1','2','3','4','5','6','7','8','9')
plotitems<-c(3,5,9)

theta <-1# specify theta

ylims <- c(0,0.025)
y <- seq(min,max,(max-min)/45)


#### density

prob <- 0*y  ### only dummy
delta<-0*y  ## dummy

par(cex.axis=1.0,cex.main=1.2,cex=1.2,lwd=0.7) 

i<-1
for (l in 1:length(y))prob[l]<- {respfctd(theta - diffunct(parmatrest[i,1],parmatrest[i,2],
                                                           y[l]))*derdiffunct(parmatrest[i,1],parmatrest[i,2],y[l])
}

plot(y,prob,cex.lab=1.0, main ="Densities",xlab ="y", type="b",
     ylab="", pch =pcdum[1],ylim=ylims )

for (i in plotitems){
  for (l in 1:length(y))prob[l]<- {respfctd(theta - diffunct(parmatrest[i,1],parmatrest[i,2],y[l]))*
      derdiffunct(parmatrest[i,1],parmatrest[i,2],y[l])
  } 
  sum(prob)*(y[2]-y[1])
  lines(y,prob,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2, type="b",
        pch =pcdum[i])
}







