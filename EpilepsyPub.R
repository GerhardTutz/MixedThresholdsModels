

#### fits for epilepsy data from masss

library("MASS")
install.packages('TMB', type = 'source')
library("glmmTMB")
library("bbmle") ## for AICtab
library("ggplot2")
library("sjlabelled")



library(tidyr)
library(spatstat)

#epil2
epil   ## from mass
summary(epil)


#### subject 49 with response 102
epil$y[193]
epil$subject[193]
epil$subject[193:196]

##set without subject 49
epilnew<-epil [c(1:192,197:236),]

## preparing data set
resp<-epilnew[c('subject','y','period')]
resp<-as.data.frame(resp)
dataresp <- spread(resp, key="period", value="y")
dataresp<-t(dataresp)
dataresp<-dataresp[2:5,]
dim(epilnew)

dum<-seq(1,dim(epil)[1]-4,4)
covdum<-epilnew[dum,]
cov<-covdum[c('trt','base','age')]
covnum<-as_numeric(cov)

covnum<-t(as.matrix(covnum))
covnum[1,]<-covnum[1,]-1   ## trt o 1


# data
dataresp  ## period x subjects
covnum    ## cov x subjects

#dim(covnum)



########################
#### repfct and respfct have to be chosen

## run: Difficultyfunctions.R
## run: PResponsefunctions.R
source("Difficultyfunctions")
source("Responsefunctions")


### selection of difficulty function (lin) and response function (respfct)

## choose lin, that is, difficulty function
lin <- "log1"

## choose respf

#respf<- "NV"
#respf<- "logistic"   
#respf<- "Gompertz"
respf<- "Gumbel"

# choose continuous or discrete
#indicator<- "C"
indicator<- "D0"

### always run difffct and respfct after selection of lin and respfct:

if(lin =="lin") {diffunct<-diffunctlin
derdiffunct<-derdiffunctlin}
if(lin =="log") {diffunct<-diffunctlog
derdiffunct<-derdiffunctlog}
if(lin =="log1") {diffunct<-diffunctlog1
derdiffunct<-derdiffunctlog1}
if(lin =="logit") {diffunct<-diffunctlogit
derdiffunct<-derdiffunctlogit}

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




I<- dim(dataresp)[1]
indicatorsingle  <- "D0"  #### discrete or continuous
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)



####################################
## 

## run: ProgramsBSplinesPureResp.R
## run: ProgramsBSplinesPenalties.R
source("ProgramsBSplinesPureResp")
source("ProgramsBSplinesPenalties")



#### with global program

ThrcovGlob <-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=pargl, 
                                      startind = "no",der="nonder",cov=covnum,hessiander="FALSE", maxit=1000)

ThrcovGlob$parmatrix 
ThrcovGlob$Loglik

##repeat
start<-ThrcovGlob$par
maxit<-1000
ThrcovGlob<-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=start, 
                                     startind = "yes",der="nonder",cov=covnum,hessiander="FALSE", maxit=maxit)
ThrcovGlob$parmatrix 
ThrcovGlob$Loglik
ThrcovGlob$AIC



### without treatment

ThrcovGlobtreat<-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=start, 
                                          startind = "no",der="nonder",cov=covnum[2:3,],hessiander="FALSE", maxit=maxit)
ThrcovGlobtreat$parmatrix 
ThrcovGlobtreat$Loglik
ThrcovGlobtreat$AIC



##repeat
start<-ThrcovGlobtreat$par
maxit<-1000
ThrcovGlobtreat<-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=start, 
                                          startind = "yes",der="nonder",cov=covnum[2:3,],hessiander="FALSE", maxit=maxit)
ThrcovGlobtreat$parmatrix 
ThrcovGlobtreat$Loglik
ThrcovGlobtreat$AIC


### without base or age

covnow<-covnum[c(1,3),]  ## base
covnow<-covnum[c(1,2),]  ## age

ThrcovGlobbase<-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=start, 
                                         startind = "no",der="nonder",cov=covnow,hessiander="FALSE", maxit=maxit)
ThrcovGlobbase$parmatrix 
ThrcovGlobbase$Loglik
ThrcovGlobbase$AIC



##repeat
start<-ThrcovGlobbase$par
maxit<-1000
ThrcovGlobbase<-ThreshbasisRespCovGlobal(dataresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=start, 
                                         startind = "yes",der="nonder",cov=covnow,hessiander="FALSE", maxit=maxit)
ThrcovGlobbase$parmatrix 
ThrcovGlobbase$Loglik
ThrcovGlobbase$AIC

bas<- -604.678
val<- -605.169
dif<-2*(bas-val)
1-pchisq(dif, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE)



### fit with glmmTMB
(fm <- glmmTMB(y ~ trt+base + age +  (1|subject),
               data=epilnew, family=nbinom2))

(fm <- glmmTMB(y ~ trt+base + age +  (1|subject),
               data=epilnew, family=poisson()))



################################################################
################################################################


