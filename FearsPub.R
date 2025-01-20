
##### fit fears data without and with covariates


library("pracma")
library(splines)
library(tidyr)
library(spatstat) 
library("foreign") 

###### Fears data   

#### load fears200
datfears<-readRDS("fears200")
head(datfears)

datfearsresp<-t(datfears[,1:5]-1)  ### -1
datfearscov<-t(datfears[,6:10])



#### repfct and respfct have to be chosen

### selection of difficulty function (lin) and response function (respfct)

## choose lin, that is, difficulty function

#lin <- "log1"
#lin <- "lin"
lin<- "logit"  


## choose respf

#respf<- "NV"
respf<- "logistic"   


# choose continuous or discrete
#indicator<- "C"
indicator<- "D"



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



########################################
###  if(lin =="logit")  difficulty function has to defined anew
### for continuous symmetric, for discrete assymmetric


## thresholds and adaptation:
min<- 0  
max <- 6
c<- 0.5

minnew<-min-c #*width
maxnew<-max+c  #*width

minnewcat<-min-2*c #*width
maxnewcat<-max+0.01 

if(indicator=="C"){
  
  diffunctlogit <-function(intdiff,slopediff,y){
    #minnew<--0.5
    #maxnew<-6.5
    r<- intdiff+slopediff*log((y-minnew)/(maxnew-y))
    return(r)  }
  derdiffunctlogit <-function(intdiff,slopediff,y){
    #minnew<--0.5
    #maxnew<-6.5
    r<- slopediff*(maxnew-minnew)/((y-minnew)*(maxnew-y))
    return(r)  }
}

if(indicator=="D"){
  diffunctlogit <-function(intdiff,slopediff,y){
    #minnew<--2
    #maxnew<-6.01 
    r<- intdiff+slopediff*log((y-minnewcat)/(maxnewcat-y))
    return(r)   }
  derdiffunctlogit <-function(intdiff,slopediff,y){
    #minnew<--2
    #maxnew<-6.01 
    r<- slopediff*(maxnewcat-minnewcat)/((y-minnewcat)*(maxnewcat-y))
    return(r)  }
}  


if(lin =="logit") {diffunct<-diffunctlogit
derdiffunct<-derdiffunctlogit}


########
I<- dim(datfearsresp)[1]
indicatorsingle  <- "D0"  #### discrete or continuous
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)
indicator<-indicatorvec





## run: Difficultyfunctions.R
## run: Responsefunctions.R
source("Difficultyfunctions")
source("Responsefunctions")
## run: ProgramsBSplinesPureResp.R
## run: ProgramsBSplinesPenalties.R
source("ProgramsBSplinesPureResp")
source("ProgramsBSplinesPenalties")


## without covariates  

cov0<-matrix(0,1,dim(datfearsresp)[2])  ## empty covariates
Thrw1<- ThreshbasisRespCov(datfearsresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start, 
                           startind = "no",der="der",cov0,hessiander="FALSE")
Thrw1$parmatrix


####  with covariates

datfearscov[1,]<- (datfearscov[1,]-mean(datfearscov[1,]))/std(datfearscov[1,]) #age standardized
cov<-datfearscov[1:4,]

ThrCnM <-ThreshbasisRespCov(datfearsresp,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start, 
                            startind = "no",der="der",cov,hessiander="TRUE")

round(ThrCnM$parmatrix, digits=3)  
round(ThrCnM$itemmatrixtval, digits=3)
