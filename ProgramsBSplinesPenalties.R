

############### functions with covariate penalties

ThreshbasisRespCovPen <- function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der,
                                  cov,hessiander, penvector,pnorm,pentype,cepsilon,maxit){
  
  #### first run diffunct with indicator lin, not basis, and respf
  #### with covariates in fixed, now +x*gamma!
  #### hessiander :true or false for computation with derivative
  
  #### datb I x P matrix responses
  ####  basis: splines for splines, lin for linear,
  #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
  #### cov            covariates (cov x person matrix)
  #### covind          if "yes"  covariates included
  #### ord:           order of B-Splines
  #### numknotsstart: initial number of basis functions
  #### lambdpenin:    penalizes violation of increase 
  #### lambdc:        penalizes variation of coefficients 0: varying coefficients
  #### pencov         penalizes covariates
  #### penvector      vector of penalty parameters for variables
  #### pnprm          obsolete, only in derivative (=2)
  #### pentype: ridge or grlasso or grlassosel
  #### cepsilon: only for grlassosel, approximation (=0 no approximation)
  #### maxit          number of iterations 
  
  #### derivative (der)  only for numknots=2 (linear)
  
  datb<- dat
  
  #Memlist <- list("PenM1" = I)
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  itemsums <- rowSums(datb)/P
  numcov<-dim(cov)[1]
  
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  if(basis =="splines"){
    dif<- maxresp-minresp
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    ### determine number of knots and basis
    sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
    
    for (i in 1:I){for (p in 1:P){
      bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    #fact <-1
    #parp   <- seq(1,numknots,1)
    #par <- parp-itemsums[1]
    #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
    #par <- c(par,1)
    
    fact <-1
    parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
    par <- parp
    for (l in 2:I) par <- c(par,parp)
    par <- c(par,1)
  }### end splines
  #length(par)
  
  
  
  if(basis !="splines"){
    
    
    #for (i in 1:numcov)indicator<-c(indicator,"C")  #### !!!
    
    numknots <- 2+numcov
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,diffunct(0,1,datb[i,p]),-cov[,p])}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,derdiffunct(0,1,datb[i,p]),0*cov[,p])}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      #basdummy[i,p,]<- c(1,datdum)
      basdummy[i,p,]<- c(1,diffunct(0,1,datdum),-cov[,p])  ### here covariates
    }}
    
    ### extension covariates
    #basderiv[,,3]
    #basdummy[,,2]
    par1 <- -itemsums  ### 
    par <- c(par1[1],1,0*cov[,p]) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],1,0*cov[,p])
    par <- c(par,1) ### stdmixt
  } # end not splines
  
  
  
  
  
  if(startind == "yes") par <- start 
  
  
  if(der!="der") {fittrw <- optim(par, LoglikIntCov, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                                  I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                                  lower = -Inf, upper = Inf, hessian = FALSE,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,
                                  pentype=pentype,cepsilon=cepsilon,control = list(maxit=maxit)) }
  
  if(der=="der") {fittrw <- optim(par, LoglikIntCov, gr = derLoglikBsplIPen,lambdpenin ,lambdc ,numknots=numknots,
                                  dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                  lower = -Inf, upper = Inf, hessian = hessiander,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon,
                                  control = list(maxit=maxit)) }
  
  
  #LoglikIntExtBsplIN(par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
  #                     indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  ll <- length(par)
  itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,I)) ##I x numknots matrix
  stdmixt <- fittrw$par[ll]
  
  
  itemmatrixstd<-  0
  itemmatrixtval<-0
  
  if(hessiander=="TRUE"){
    covm<-diag(solve(fittrw$hessian))
    stdev<-sqrt(covm)
    itemmatrixstd<-  t(matrix(stdev[1:ll-1],numknots,I))
    itemmatrixtval<-itemmatrixest/itemmatrixstd
  }
  
  #loglikorig<--fittrw$value
  
  derivfin <-0  ## dummy
  #derivative here
  
  derivfin <-derLoglikBsplIPen(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
                               indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,
                               numcov=numcov,
                               penvector=penvector,pnorm=pnorm)
  
  ## logorig
  loglikorig <- -LoglikIntExtBsplIN(fittrw$par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  ### positive stddev
  stdmixt<-abs(stdmixt) 
  fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
  ####
  AIC <-  -2*loglikorig +2*length(par)
  BIC <-  -2*loglikorig +2*log(length(par))
  
  #########
  newList <- list("parmatrix" = itemmatrixest,"stdmixt"=stdmixt,"par"=fittrw$par, "numbasis"=numknots, 
                  "Loglikwithpen"=-fittrw$value, "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                  "AIC"=AIC, "BIC"=BIC,"derivative"=derivfin,"itemmatrixstd"=itemmatrixstd,
                  "itemmatrixtval"=itemmatrixtval)
  return(newList)
  
}
########################################################



LoglikIntCov <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                        indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
                        numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon){
  
  ### item specific
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  
  # numcov penalty on covariates
  # penvector new!!
  # pnorm: chosen norm: pnorm=2 is Euclidean-squared
  # cepsilon: modification of norm,  =0 no modification
  #print(par)
  
  P <- dim(dat)[2]
  
  numpar <- length(par) 
  sdnormal<- par[numpar]
  #paroriginal <- par[1:numpar-1]
  
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = sdnormal,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots)
    
    #if (int<(10)^{-300}) int<-(10)^{-300} #log((10)^{-300}) 
    ind[p]<-log(int)
    #int<-abs(int)  # 
    
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
    
  } #  end p
  
  dimBspl<-length(par)-1
  parBspl <- par[1:dimBspl]
  itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
  
  
  ####penalty increase
  #lambdpenin <- 0
  if (lambdpenin >0){
    
    upp<- numknots-1
    peninc<-0
    for (i in 1:I){
      for (l in 1:upp ){dif  <-itemmatrix[i,l+1]- itemmatrix[i,l]
      peninc <-peninc + dif^2 *indicatt(-dif , 0)}}
    
    #pen <- dif%*%dif
    sum <- sum -  lambdpenin*peninc ### + weil negativ , jetzt - ????
    
  } #end if
  ############
  
  #### penalty common difficulties
  #lambdc <- 3
  if (lambdc >0){
    upp1<- I-1
    upp2<- numknots-1
    penc<-0
    for (i in 1:upp1){
      for (l in 1:upp2 ){dif  <-(itemmatrix[i+1,l+1]- itemmatrix[i+1,l])-(itemmatrix[i,l+1]- itemmatrix[i,l])
      penc <- penc +dif*dif}}
    
    ##### funktioniert: with matrix penc = pencc jetzt auskommentiert
    
    #upp1<- I-1
    #upp2<- numknots-1
    #counth <- 0
    #for (i in 1:upp1){
    #  for (l in 1:upp2 ){
    #  counth <- counth +1
    #  penM <- matrix(0,1,length(par))
    #  penM[(i-1)*numknots+l+1]<-1
    #  penM[(i-1)*numknots+l]<--1
    #  penM[(i)*numknots+l+1]<--1
    #  penM[(i)*numknots+l]<- 1
    #  if (counth <= 1) pennew <- penM
    # if(counth > 1) pennew <- rbind(pennew,penM)
    # }}
    #pencc <- t(pennew%*%par)%*% (pennew%*%par) 
    
    ####ende auskommentiert
    
    if (numknots == 2) {penc<-0
    for (i in 2:I) penc<- penc + (itemmatrix[i,2]-itemmatrix[i-1,2])^2
    
    ## alternative pennew=penc
    M<- matrix(0,I-1,I)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,I-2))
    M[I-1,]<-c(rep(0,I-2),d1)
    I2 <- I-2
    for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    pennew <- t(itemmatrix[,2])%*%t(M)%*%M%*%itemmatrix[,2]
    ###
    }
    
    
    sum <- sum -  lambdc*penc ### 
    #sum <- sum -  lambdc*pennew 
  } #end if
  
  
  #### penalty covariates
  pencov<-sum(penvector)
  red<-numknots-numcov+1
  if (pencov >0){
    penm<-as.matrix(itemmatrix[,red:numknots])
    #print("numknots")
    #print(numknots)
    
    #print(penm)   #################################################!
    #pemmdif<-penm[2,]-penm[1,]  ## difference matrix
    #for (i in 3:I) pemmdif<-rbind(pemmdif,penm[i,]-penm[i-1,])
    
    ### new generate pairs
    cpairs<-I*(I-1)/2
    M<-matrix(0,cpairs,I)
    count<-0
    for (i in 1:(I-1)){
      for (j in (i+1):I) {count<-count+1
      M[count,i]<-1
      M[count,j]<--1}} 
    diffm<-M%*%penm
    #diffv<-as.vector(diffm)
    #############
    
    
    if(pentype=="ridge"){
      #normvect<-Norm(pemmdif[,1], p = pnorm)
      #for (v in 2:numcov){normvect<-c(normvect,Norm(pemmdif[,v], p = pnorm))}   
      
      #normvect<-(Norm(diffm[,1], p = pnorm))^pnorm
      # for (v in 2:numcov){normvect<-c(normvect,(Norm(diffm[,v], p = pnorm))^pnorm)}
      normvect<-t(diffm[,1])%*%diffm[,1]
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,t(diffm[,v])%*%diffm[,v])}   
      
      sum <- sum -  penvector%*%normvect
    } ## end ridge
    
    if(pentype=="grlasso"){
      
      c<-cepsilon   ## with c, modified norm!
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}
      
      
      sum <- sum -  penvector%*%normvect
    } ## end lasso
    
    if(pentype=="grlassosel"){
      
      c<-cepsilon   ## with c, modified norm!
      
      # first term
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}
      
      sum <- sum -  penvector%*%normvect#*(cpairs^(1/2))
      
      
      # second term
      #normvect<-(t(penm[,1])%*%penm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(penm[,v])%*%penm[,v])^(1/2))}
      
      
      #print("v")
      #print(v)
      
      normvect<-NormShiftc((t(penm[,1])%*%penm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(penm[,v])%*%penm[,v]),c))
      #print("v")
      #print(v,numcov)
      }
      
      sum <- sum -  penvector%*%normvect#*(numcov^(1/2))
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)+(t(penm[1,])%*%penm[1,])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2)+(t(penm[v,])%*%penm[v,])^(1/2))}   
      
      
    } ## end grlassosel
    
  } ## end pen cov
  
  
  
  sum<- -sum
  return(sum)}
###############

derLoglikBsplIPen <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                             indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,numcov=numcov,
                             penvector=penvector,pnorm=pnorm){
  ### derivative, only continuous
  ### only ridge
  #parmatr <- matrix(par,I,2)
  # derivfin <-derLoglikBsplI(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
  #indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  
  slope <- 1 ### slope fixed  weight on theta!
  P <- dim(dat)[2]
  
  lpar <- length(par)
  der <- matrix(0,lpar,1)
  
  parmatr <- t(matrix(par[1:lpar-1],numknots,I))
  stdmixt <- par[lpar]
  
  # basis functions prodfctIntExtBsplINVec 
  
  for (p in 1:P){
    datitem <- dat[,p]
    
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots) 
    #### modified for global
    #Ilimit<-I
    #if(lpar <= I+1)Ilimit<-1
    ######
      
      for (l in 1:I){
      for (j in 1:numknots){ 
        comp <- l
        pers<-p
        #count <- (l-1)*numknots+j
        count <- j
        #der[l,1] <-der[l,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
        #                    indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp)/int } 
        der[count,1] <-der[count,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,
                                                  numknots=numknots,I=I,par=par, slope=slope,indicator=indicator,order = 10, bas=bas,
                                                  basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp, j=j)/int } # end j 
    } #end l
  } # end p
  
  
  
  
  
  
  #### now numeric !!!!!
  deltan<- 0.00005
  vect1<- par
  vect1[lpar] <- vect1[lpar]+deltan
  #loglikorig  <- LoglikN(par,dat,I,indicator, lin=lin)
  #loglikorignow <- LoglikN(vect1,dat,I,indicator, lin=lin)
  
  loglikorig<-LoglikIntExtBsplIN(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                 indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  loglikorignow<-LoglikIntExtBsplIN(vect1,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  derivnow <- (loglikorignow-loglikorig)/deltan
  der[lpar,1]<- -derivnow
  
  
  ##### end numeric
  
  
  
  #lambdc <- 3
  
  if (lambdpenin >0){
    numknots
    M<- matrix(0,numknots-1,numknots)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,numknots-2))
    numknots1 <- numknots-1
    if(numknots >2)for(l in 2:numknots1) M[l,]<- c(rep(0,l-1),d1,rep(0,numknots-l-1))
    dum <- diag(1,I)
    Mtot <- kronecker(dum,M)
    Mtot<- cbind(Mtot,rep(0,dim(Mtot)[1]))
    pen1 <- 2*t(Mtot)%*%(Mtot%*%par*indicattmatrix(-Mtot%*%par,0))
    
    der <- der - lambdpenin*pen1
  }
  ## end penin
  #dim(Mtot)
  
  if (lambdc >0){
    
    if (numknots > 2){
      upp1<- I-1
      upp2<- numknots-1
      counth <- 0
      for (i in 1:upp1){
        for (l in 1:upp2 ){
          counth <- counth +1
          penM <- matrix(0,1,length(par))
          penM[(i-1)*numknots+l+1]<-1
          penM[(i-1)*numknots+l]<--1
          penM[(i)*numknots+l+1]<--1
          penM[(i)*numknots+l]<- 1
          if (counth <= 1) pennew <- penM
          if(counth > 1) pennew <- rbind(pennew,penM)
        }}
      pencc <- t(pennew%*%par)%*% (pennew%*%par) 
      
      deriv1 <- 2*t(pennew)%*%pennew%*%par
      der <- der -  lambdc*deriv1
    }
    
    if (numknots == 2) {
      M<- matrix(0,I-1,I)
      d1 <- c(-1,1)
      M[1,]<- c(d1,rep(0,I-2))
      M[I-1,]<-c(rep(0,I-2),d1)
      I2 <- I-2
      for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
      
      deriv <- 2*t(M)%*%M%*%parmatr[,2]
      for (l in  1:I){der[l*2,1]<-der[l*2,1] -lambdc*deriv[l]}
      
    }} #end if
  
  pencov<-sum(penvector)
  
  
  ####################### new penalty covariates
  
  if (pencov >0){
    red<-numknots-numcov+1
    penm<-as.matrix(parmatr[,red:numknots])
    
    #if (numknots == 2) {
    #M<- matrix(0,I-1,I)
    #d1 <- c(-1,1)
    #M[1,]<- c(d1,rep(0,I-2))
    # M[I-1,]<-c(rep(0,I-2),d1)
    # I2 <- I-2
    # for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    
    ### new generate pairs
    cpairs<-I*(I-1)/2
    M<-matrix(0,cpairs,I)
    count<-0
    for (i in 1:(I-1)){
      for (j in (i+1):I) {count<-count+1
      M[count,i]<-1
      M[count,j]<--1}} 
    #diffm<-M%*%penm 
    
    #}
    
    
    #fact1<-(t(M%*% as.matrix(penm[,1]))%*%(M%*% as.matrix(penm[,1])))^(-1/2)/2
    #deriv1 <- penvector[1]*fact1[1]*(2*t(M)%*%M%*%penm[,1])
    deriv1 <- penvector[1]*(2*t(M)%*%M%*%penm[,1])
    
    seq<-c(2:numcov)
    if (numcov >1)for (v in seq){          ####### 24.2
      # fact1<-(t(M%*% as.matrix(penm[,v]))%*%(M%*% as.matrix(penm[,v])))^(-1/2)/2
      # derivn <- penvector[v]*fact1[1]*(2*t(M)%*%M%*%penm[,v])
      derivn <- penvector[v]*(2*t(M)%*%M%*%penm[,v])
      deriv1<-cbind(deriv1,derivn)   
    } ## end if numcov
    
    derivm<-cbind(matrix(0,I,2),deriv1)
    derivv<- as.matrix(c(t(derivm)))
    
    
    #der <-  der -  as.matrix(c(rep(0,((numknots-numcov)*I)),derivv,0))
    der<- der-as.matrix(c(derivv,0))
    
  }  ### end pencov
  
  
  
  
  der <- -der
  return(der)}
###############



##################################################

prodfctIntExtBsplINderiv1Pen <-function(theta,datitem=datitem,numknots=numknots,I=I,par=par,slope=slope,
                                        indicator=indicator,order = order, bas,basderiv,basdummy,
                                        pers=pers,comp=comp, j=j,numcov=numcov,penvector=penvector,pnorm=pnorm){
  #### with std mixt only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  #####  j basis function
  
  prod <-1    
  set <- seq(1,I,1)
  
  
  dimBspl<-length(par)-1
  parBspl <- par[1:dimBspl]
  itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
  
  
  nonset <- set[-comp] 
  
  for (i in nonset){
    iteml <- itemmatrix[i,]
    item=i
    if(indicator[i]=="C"){prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,
                                                           basdummy=basdummy,item=item,pers=pers)} 
    if(indicator[i]!="C"){prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,
                                                                basdummy=basdummy,item=item,pers=pers,indicator[i])}
    #if(indicator[i]=="D1"){prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,
    #                                                             basdummy=basdummy,item=item,pers=pers,indicator[i])}
  }
  
  ## der f
  iteml <-itemmatrix[comp,]
  deltaf <-  bas[comp,pers,]%*%iteml
  prodderiv <- basderiv[comp,pers,]%*%iteml
  
  
  #if(indicator[comp]=="C"){prod <- prod*(basderiv[comp,pers,j]*respfctd(theta-deltaf)+
  #                                         (theta-deltaf)*respfctd(theta-deltaf)*bas[comp,pers,j]*prodderiv)}
  if(indicator[comp]=="C"){prod <- prod*(basderiv[comp,pers,j]*respfctd(theta-deltaf)-
                                           respfctder(theta-deltaf)*bas[comp,pers,j]*prodderiv)}
  
  if(indicator[comp]=="D1"){ 
    if(datitem[comp] <= 1){
      val <-  bas[comp,pers,]%*%iteml
      prod <- prod*respfctd(slope*theta - val)*bas[comp,pers,j]}
    if(datitem[comp] > 1){
      val <-  bas[comp,pers,]%*%iteml
      valprev <-  basdummy[comp,pers,]%*%iteml
      prod <- prod*(-respfctd(slope*theta - valprev)*basdummy[comp,pers,j]+
                      respfctd(slope*theta - val)*bas[comp,pers,j])}
  }## end D1 
  
  if(indicator[comp]=="D0"){ 
    if(datitem[comp] <= 0){
      val <-  bas[comp,pers,]%*%iteml
      prod <- prod*respfctd(slope*theta - val)*bas[comp,pers,j]}
    if(datitem[comp] > 0){
      val <-  bas[comp,pers,]%*%iteml
      valprev <-  basdummy[comp,pers,]%*%iteml
      prod <- prod*(-respfctd(slope*theta - valprev)*basdummy[comp,pers,j]+
                      respfctd(slope*theta - val)*bas[comp,pers,j])}
  }## end D1 
  
  
  return(prod)}

#########################



#####################
NormShiftc <-function(x,c){
  ## used argument x^Tx
  y  <-sqrt(abs(x)+c)-sqrt(c)
  return(y)}

c<-.1 
xpl<-seq(-3,3,.01)
ypl<-NormShiftc(xpl[1],c)
for(i in 2:length(xpl))ypl<-c(ypl,NormShiftc(xpl[i],c))
plot(xpl,ypl)


###########################################################
###### for loss estimation 

LoglikCovLoss <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                         indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
                         numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon){
  
  ### item specific
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  
  # numcov penalty on covariates
  # penvector new!!
  # pnorm: chosen norm: pnorm=2 is Euclidean-squared
  # cepsilon: modification of norm,  =0 no modification
  
  
  P <- dim(dat)[2]
  
  numpar <- length(par) 
  sdnormal<- par[numpar]
  #paroriginal <- par[1:numpar-1]
  
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = sdnormal,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots)
    ind[p]<-log(int)
    #int<-abs(int)  # 
    
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
    
  } #  end p
  
  
  return(sum)}
################################################

PosteriorEstimatesN  <-function(grid,dat,I,basis=basis,indicatorvec=indicatorvec,ord=ord,cov= cov,
                                numknotsstart=numknotsstart,parmatrest=parmatrest,stdest=stdest){ 
  
  
  ####computes centered posterior person parameters with covariates
  ######   parmatrest: estimated parameters
  
  
  numgrid <- length(grid)
  dens <- matrix(0,numgrid,1)
  P<-dim(dat)[2]
  esttheta <- matrix(0,P,1)
  
  
  ### general than splines
  datb<- dat
  
  #Memlist <- list("PenM1" = I)
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  itemsums <- rowSums(datb)/P
  numcov<-dim(cov)[1]
  
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  if(basis =="splines"){
    dif<- maxresp-minresp
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    ### determine number of knots and basis
    sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
    
    for (i in 1:I){for (p in 1:P){
      bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    #fact <-1
    #parp   <- seq(1,numknots,1)
    #par <- parp-itemsums[1]
    #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
    #par <- c(par,1)
    
    fact <-1
    parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
    par <- parp
    for (l in 2:I) par <- c(par,parp)
    par <- c(par,1)
  }### end splines
  #length(par)
  
  
  
  if(basis !="splines"){
    
    
    #for (i in 1:numcov)indicator<-c(indicator,"C")  #### !!!
    
    numknots <- 2+numcov
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,diffunct(0,1,datb[i,p]),-cov[,p])}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,derdiffunct(0,1,datb[i,p]),0*cov[,p])}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      #basdummy[i,p,]<- c(1,datdum)
      basdummy[i,p,]<- c(1,diffunct(0,1,datdum),-cov[,p])  ### here covariates
    }}
    
    ### extension covariates
    #basderiv[,,3]
    #basdummy[,,2]
    par1 <- -itemsums  ### 
    par <- c(par1[1],1,0*cov[,p]) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],1,0*cov[,p])
    par <- c(par,1) ### stdmixt
  } # end not splines
  
  
  
  
  
  
  ### end general
  
  for (p  in 1:P){obs<- datb[,p]
  for (l in 1:numgrid){
    theta <- grid[l]
    #dens[l,1]<-prodfct(theta,dat=obs,I=I,parmatr=parmatrest,slope=1,indicator=indicator, lin =lin)*dnorm(theta,0,stdest)
    dens[l,1]<-prodfctIntExtBsplINVec(theta,dat=obs,I=I,par=parmatrest,slope=1,
                                      indicator=indicatorvec,bas,basderiv,basdummy,
                                      pers=p,numknots=numknots)*dnorm(theta,0,stdest)
  }
  num<-  which.max(dens)
  esttheta[p]<- grid[num]
  }
  
  #### centering
  esttheta<-esttheta-mean(esttheta)
  
  return(esttheta)
}
###### end fct  
##########################################################


ThreshbasisRespCovPenItSpec <- function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der,
                                  cov,hessiander, penvector,penall,pnorm,pentype,cepsilon,maxit){
  
  #### first run diffunct with indicator lin, not basis, and respf
  #### with covariates in fixed, now +x*gamma!
  #### hessiander :true or false for computation with derivative
  
  #### penall =  "all" also slopes and intercepts penalized
  #              penvector longer
  
  #### dat  I x P matrix responses
  ####  basis: splines for splines, lin for linear,
  #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
  #### cov            covariates as array var x pers x item  !!!!!!!!!  array
   
  #### ord:           order of B-Splines
  #### numknotsstart: initial number of basis functions
  #### lambdpenin:    penalizes violation of increase 
  #### lambdc:        penalizes variation of coefficients 0: varying coefficients
  #### pencov         penalizes covariates
  #### penvector      vector of penalty parameters for variables
  #### pnprm          obsolete, only in derivative (=2)
  #### pentype: ridge or grlasso or grlassosel
  #### cepsilon: only for grlassosel, approximation (=0 no approximation)
  #### maxit          number of iterations 
  
  #### derivative (der)  only for numknots=2 (linear)
  
  datb<- dat
  
  #Memlist <- list("PenM1" = I)
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  itemsums <- rowSums(datb)/P
  numcov<-dim(cov)[1]
  #dim(covwine)[1]
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  if(basis =="splines"){
    dif<- maxresp-minresp
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    ### determine number of knots and basis
    sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
    
    for (i in 1:I){for (p in 1:P){
      bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    #fact <-1
    #parp   <- seq(1,numknots,1)
    #par <- parp-itemsums[1]
    #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
    #par <- c(par,1)
    
    fact <-1
    parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
    par <- parp
    for (l in 2:I) par <- c(par,parp)
    par <- c(par,1)
  }### end splines
  #length(par)
  
  
  
  if(basis !="splines"){
    
    
    #for (i in 1:numcov)indicator<-c(indicator,"C")  #### !!!
    
    numknots <- 2+numcov
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,diffunct(0,1,datb[i,p]),-cov[,p,i])}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,derdiffunct(0,1,datb[i,p]),0*cov[,p,i])}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      #basdummy[i,p,]<- c(1,datdum)
      basdummy[i,p,]<- c(1,diffunct(0,1,datdum),-cov[,p,i])  ### here covariates
    }}
    
    ### extension covariates
    #basderiv[,,3]
    #basdummy[,,2]
    par1 <- -itemsums  ### 
    par <- c(par1[1],1,0*cov[,p,1]) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],1,0*cov[,p,i])
    par <- c(-10,1,0*cov[,p,1]) ### 1:slope difficulty!!!
    for(i in 2:I) par <- c(par,-10,1,0*cov[,p,i])
    par <- c(par,2) ### stdmixt
  } # end not splines
  
  
  
  
  
  if(startind == "yes") par <- start 
  
  if(penall!="all"){
  if(der!="der") {fittrw <- optim(par, LoglikIntCov, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                                  I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                                  lower = -Inf, upper = Inf, hessian = FALSE,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,
                                  pentype=pentype,cepsilon=cepsilon,control = list(maxit=maxit)) }
  
  if(der=="der") {fittrw <- optim(par, LoglikIntCov, gr = derLoglikBsplIPen,lambdpenin ,lambdc ,numknots=numknots,
                                  dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                  lower = -Inf, upper = Inf, hessian = hessiander,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon,
                                  control = list(maxit=maxit)) }
  
  }
  if(penall=="all"){
    if(der!="der") {fittrw <- optim(par, LoglikIntCovAll, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                                    I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf, hessian = FALSE,
                                    numcov=numcov,penvector=penvector,pnorm=pnorm,
                                    pentype=pentype,cepsilon=cepsilon,control = list(maxit=maxit)) }
    
    if(der=="der") {fittrw <- optim(par, LoglikIntCovAll, gr = derLoglikBsplIPen,lambdpenin ,lambdc ,numknots=numknots,
                                    dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                    lower = -Inf, upper = Inf, hessian = hessiander,
                                    numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon,
                                    control = list(maxit=maxit)) }
    
  }
  
  
  
  
  
  ll <- length(par)
  itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,I)) ##I x numknots matrix
  stdmixt <- fittrw$par[ll]
  
  
  itemmatrixstd<-  0
  itemmatrixtval<-0
  
  if(hessiander=="TRUE"){
    covm<-diag(solve(fittrw$hessian))
    stdev<-sqrt(covm)
    itemmatrixstd<-  t(matrix(stdev[1:ll-1],numknots,I))
    itemmatrixtval<-itemmatrixest/itemmatrixstd
  }
  
  #loglikorig<--fittrw$value
  
  derivfin <-0  ## dummy
  #derivative here
  
  
  ### derivative omitted
  #derivfin <-derLoglikBsplIPen(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
  #                             indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,
  #                             numcov=numcov,
  #                             penvector=penvector,pnorm=pnorm)
  
  ## logorig
  loglikorig <- -LoglikIntExtBsplIN(fittrw$par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  #-LoglikIntExtBsplIN(par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
  #                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  
  ### positive stddev
  stdmixt<-abs(stdmixt) 
  fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
  ####
  AIC <-  -2*loglikorig +2*length(par)
  BIC <-0
  #BIC <-  -2*loglikorig +2*log(length(par))
  
  #########
  newList <- list("parmatrix" = itemmatrixest,"stdmixt"=stdmixt,"par"=fittrw$par, "numbasis"=numknots, 
                  "Loglikwithpen"=-fittrw$value, "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                  "AIC"=AIC, "BIC"=BIC,"derivative"=derivfin,"itemmatrixstd"=itemmatrixstd,
                  "itemmatrixtval"=itemmatrixtval)
  return(newList)
  
}
########################################################


LoglikIntCovAll <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                        indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
                        numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon){
  
  ### item specific all parameters penalized also slopes and intercepts
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  
  # numcov 
  # penvector new!!
  # pnorm: chosen norm: pnorm=2 is Euclidean-squared
  # cepsilon: modification of norm,  =0 no modification
  #print(par)
  
  P <- dim(dat)[2]
  
  numpar <- length(par) 
  sdnormal<- par[numpar]
  #paroriginal <- par[1:numpar-1]
  
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = sdnormal,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots)
    #if (int<(10)^{-300}) int<-(10)^{-300} #log((10)^{-300}) 
    ind[p]<-log(int)
    #int<-abs(int)  # 
    
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
    
  } #  end p
  
  dimBspl<-length(par)-1
  parBspl <- par[1:dimBspl]
  itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
  
  
  ####penalty increase
  #lambdpenin <- 0
  if (lambdpenin >0){
    
    upp<- numknots-1
    peninc<-0
    for (i in 1:I){
      for (l in 1:upp ){dif  <-itemmatrix[i,l+1]- itemmatrix[i,l]
      peninc <-peninc + dif^2 *indicatt(-dif , 0)}}
    
    #pen <- dif%*%dif
    sum <- sum -  lambdpenin*peninc ### + weil negativ , jetzt - ????
    
  } #end if
  ############
  
  #### penalty common difficulties
  #lambdc <- 3
  if (lambdc >0){
    upp1<- I-1
    upp2<- numknots-1
    penc<-0
    for (i in 1:upp1){
      for (l in 1:upp2 ){dif  <-(itemmatrix[i+1,l+1]- itemmatrix[i+1,l])-(itemmatrix[i,l+1]- itemmatrix[i,l])
      penc <- penc +dif*dif}}
    
    ##### funktioniert: with matrix penc = pencc jetzt auskommentiert
    
    #upp1<- I-1
    #upp2<- numknots-1
    #counth <- 0
    #for (i in 1:upp1){
    #  for (l in 1:upp2 ){
    #  counth <- counth +1
    #  penM <- matrix(0,1,length(par))
    #  penM[(i-1)*numknots+l+1]<-1
    #  penM[(i-1)*numknots+l]<--1
    #  penM[(i)*numknots+l+1]<--1
    #  penM[(i)*numknots+l]<- 1
    #  if (counth <= 1) pennew <- penM
    # if(counth > 1) pennew <- rbind(pennew,penM)
    # }}
    #pencc <- t(pennew%*%par)%*% (pennew%*%par) 
    
    ####ende auskommentiert
    
    if (numknots == 2) {penc<-0
    for (i in 2:I) penc<- penc + (itemmatrix[i,2]-itemmatrix[i-1,2])^2
    
    ## alternative pennew=penc
    M<- matrix(0,I-1,I)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,I-2))
    M[I-1,]<-c(rep(0,I-2),d1)
    I2 <- I-2
    for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    pennew <- t(itemmatrix[,2])%*%t(M)%*%M%*%itemmatrix[,2]
    ###
    }
    
    
    sum <- sum -  lambdc*penc ### 
    #sum <- sum -  lambdc*pennew 
  } #end if
  
  
  #### penalty covariates
  pencov<-sum(penvector)
  red<-numknots-numcov+1
  if (pencov >0){
    penm<-as.matrix(itemmatrix[,red:numknots])
    penm<-as.matrix(itemmatrix) ########  !!
    #print("numknots")
    #print(numknots)
    
    #print(penm)   #################################################!
    #pemmdif<-penm[2,]-penm[1,]  ## difference matrix
    #for (i in 3:I) pemmdif<-rbind(pemmdif,penm[i,]-penm[i-1,])
    
    ### new generate pairs
    cpairs<-I*(I-1)/2
    M<-matrix(0,cpairs,I)
    count<-0
    for (i in 1:(I-1)){
      for (j in (i+1):I) {count<-count+1
      M[count,i]<-1
      M[count,j]<--1}} 
    diffm<-M%*%penm
    #diffv<-as.vector(diffm)
    #############
    
    
    if(pentype=="ridge"){
      #normvect<-Norm(pemmdif[,1], p = pnorm)
      #for (v in 2:numcov){normvect<-c(normvect,Norm(pemmdif[,v], p = pnorm))}   
      
      #normvect<-(Norm(diffm[,1], p = pnorm))^pnorm
      # for (v in 2:numcov){normvect<-c(normvect,(Norm(diffm[,v], p = pnorm))^pnorm)}
      normvect<-t(diffm[,1])%*%diffm[,1]
      for (v in 2:numknots){normvect<-c(normvect,t(diffm[,v])%*%diffm[,v])}  ####! 
      
      sum <- sum -  penvector%*%normvect
    } ## end ridge
    
    if(pentype=="grlasso"){
      
      c<-cepsilon   ## with c, modified norm!
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      for (v in 2:numknots){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}  ####!
      
      
      sum <- sum -  penvector%*%normvect
    } ## end lasso
    
    if(pentype=="grlassosel"){
      
      c<-cepsilon   ## with c, modified norm!
      
      # first term
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      for (v in 2:numknots){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}  ####!
      
      sum <- sum -  penvector%*%normvect#*(cpairs^(1/2))
      
      
      # second term
      #normvect<-(t(penm[,1])%*%penm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(penm[,v])%*%penm[,v])^(1/2))}
      
      
      #print("v")
      #print(v)
      
      normvect<-NormShiftc((t(penm[,1])%*%penm[,1]),c)
      for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(penm[,v])%*%penm[,v]),c))
      #print("v")
      #print(v,numcov)
      }
      
      sum <- sum -  penvector%*%normvect#*(numcov^(1/2))
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)+(t(penm[1,])%*%penm[1,])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2)+(t(penm[v,])%*%penm[v,])^(1/2))}   
      
      
    } ##
    
  } ## end pen cov
  
  
  
  sum<- -sum
  return(sum)}
###############


ThreshbasisRespCovItemsonly <- function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der,
                                  cov,hessiander, penvector,pnorm,pentype,cepsilon,maxit){
  
  #### only item specific covariates, only non-splines modified
  
  #### first run diffunct with indicator lin, not basis, and respf
  #### with covariates in fixed, now +x*gamma!
  #### hessiander :true or false for computation with derivative
  
  #### datb I x P matrix responses
  ####  basis: splines for splines, lin for linear,
  #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
  #### covitem            covariates (variables x item)
  #### covind          if "yes"  covariates included
  #### ord:           order of B-Splines
  #### numknotsstart: initial number of basis functions
  #### lambdpenin:    penalizes violation of increase 
  #### lambdc:        penalizes variation of coefficients 0: varying coefficients
  #### pencov         penalizes covariates
  #### penvector      vector of penalty parameters for variables
  #### pnprm          obsolete, only in derivative (=2)
  #### pentype: ridge or grlasso or grlassosel
  #### cepsilon: only for grlassosel, approximation (=0 no approximation)
  #### maxit          number of iterations 
  
  #### derivative (der)  only for numknots=2 (linear)
  
  
  datb<- dat
  
  #Memlist <- list("PenM1" = I)
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  itemsums <- rowSums(datb)/P
  numcov<-dim(cov)[1]
  
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  if(basis =="splines"){
    dif<- maxresp-minresp
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    ### determine number of knots and basis
    sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
    
    for (i in 1:I){for (p in 1:P){
      bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    #fact <-1
    #parp   <- seq(1,numknots,1)
    #par <- parp-itemsums[1]
    #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
    #par <- c(par,1)
    
    fact <-1
    parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
    par <- parp
    for (l in 2:I) par <- c(par,parp)
    par <- c(par,1)
  }### end splines
  #length(par)
  
  
  
  if(basis !="splines"){
    
    
    #for (i in 1:numcov)indicator<-c(indicator,"C")  #### !!!
    
    numknots <- 2+numcov
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,diffunct(0,1,datb[i,p]),-cov[,i])}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,derdiffunct(0,1,datb[i,p]),0*cov[,i])}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      #basdummy[i,p,]<- c(1,datdum)
      basdummy[i,p,]<- c(1,diffunct(0,1,datdum),-cov[,i])  ### here covariates
    }}
    
    ### extension covariates
    #basderiv[,,3]
    #basdummy[,,2]
    par1 <- -itemsums  ### 
    par <- c(par1[1],1,0*cov[,i]) ### 1:slope difficulty
    #for(i in 2:I) par <- c(par,par1[i],1,0*cov[,p])
    par <- c(par,1) ### stdmixt
  } # end not splines
  
  
  
  
  
  if(startind == "yes") par <- start 
  
  
  if(der!="der") {fittrw <- optim(par, LoglikIntCov, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                                  I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                                  lower = -Inf, upper = Inf, hessian = FALSE,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,
                                  pentype=pentype,cepsilon=cepsilon,control = list(maxit=maxit)) }
  
  #if(der=="der") {fittrw <- optim(par, LoglikIntCov, gr = derLoglikBsplIPen,lambdpenin ,lambdc ,numknots=numknots,
  #                                dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
  #                                lower = -Inf, upper = Inf, hessian = hessiander,
  #                                numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon,
  #                                control = list(maxit=maxit))  }
  ### mod
  if(der=="der") {fittrw <- optim(par, LoglikIntCov, gr = derLoglikBsplIPen,lambdpenin ,lambdc ,numknots=numknots,
                                  dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                  lower = -Inf, upper = Inf, hessian = hessiander,
                                  numcov=numcov,penvector=penvector,pnorm=pnorm,
                                  control = list(maxit=maxit))  }
  
  
  #LoglikIntCov(par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
  #                        indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
  #                        numcov=numcov,penvector=penvector,pnorm=2,pentype=pentype,cepsilon=cepsilon) 
  ll <- length(par)
  itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,1)) ##I x numknots matrix
  stdmixt <- fittrw$par[ll]
  
  
  itemmatrixstd<-  0
  itemmatrixtval<-0
  stdstdmixt<-0
  
  if(hessiander=="TRUE"){
    covm<-diag(solve(fittrw$hessian))
    stdev<-sqrt(covm)
    itemmatrixstd<-  t(matrix(stdev[1:ll-1],numknots,1))
    itemmatrixtval<-itemmatrixest/itemmatrixstd
    stdstdmixt<-stdev[ll]
    }
  
  #loglikorig<--fittrw$value
  
  derivfin <-0  ## dummy
  #derivative here
  
  derivfin <-derLoglikBsplIPen(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
                               indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,
                               numcov=numcov,
                               penvector=penvector,pnorm=pnorm)
  
  ## logorig
  loglikorig <- -LoglikIntExtBsplIN(fittrw$par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  ### positive stddev
  #stdmixt<-abs(stdmixt) 
  #fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
  ####
  AIC <-  -2*loglikorig +2*length(par)
  BIC <- 0# -2*loglikorig +2*log(length(par))
  
  #########
  newList <- list("parmatrix" = itemmatrixest,"itemmatrixstd"=itemmatrixstd,"stdmixt"=stdmixt,
                  "stdstdmixt"=stdstdmixt,"par"=fittrw$par, "numbasis"=numknots, 
                  "Loglikwithpen"=-fittrw$value, "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                  "AIC"=AIC, "BIC"=BIC,"derivative"=derivfin,
                  "itemmatrixtval"=itemmatrixtval)
  return(newList)
  
}
############################


###### Programs global. not item specific estimates for all parameters  



ThreshbasisRespCovGlobaln <- 
  function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der,
           cov,hessiander){   
    
    #### covariates not item-specific!  
    
    ### only nonder
    
    #### first run diffunct with indicator lin, not basis, and respf
    #### with covariates in fixed, now +x*gamma!
    #### cov Var x persons, only person specific covariates 
    
    #### hessiander :true or false for computation with derivative
    #### datb I x P matrix responses
    ####  basis: splines for splines, lin for linear,
    #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
    #### ord:           order of B-Splines
    #### numknotsstart: initial number of basis functions
    #### lambdpenin:    penalizes violation of increase 
    #### lambdc:        penalizes variation of coefficients 0: varying coefficients
    
    
    datb<- dat
    
    #Memlist <- list("PenM1" = I)
    
    I <- dim(datb)[1]  
    P <- dim(datb)[2]
    itemsums <- rowSums(datb)/P
    numcov<-dim(cov)[1]
    #dim(covwine)[1]
    
    minresp <- min(datb)
    maxresp <- max(datb)
    
    if(basis =="splines"){
      dif<- maxresp-minresp
      knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
      
      ### determine number of knots and basis
      sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
      numknots <- dim(sdim)[2]
      
      bas <- array(0, dim=c(I,P,numknots))
      basderiv <- array(0, dim=c(I,P,numknots))
      basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
      
      for (i in 1:I){for (p in 1:P){
        bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
      
      for (i in 1:I){for (p in 1:P){
        basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
      
      for (i in 1:I){for (p in 1:P){
        datdum <- datb[i,p]-1
        basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
      
      #fact <-1
      #parp   <- seq(1,numknots,1)
      #par <- parp-itemsums[1]
      #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
      #par <- c(par,1)
      
      fact <-1
      parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
      par <- parp
      for (l in 2:I) par <- c(par,parp)
      par <- c(par,1)
    }### end splines
    #length(par)
    
    
    
    if(basis !="splines"){
      
      
      #for (i in 1:numcov)indicator<-c(indicator,"C")  #### !!!
      
      numknots <- 2+numcov
      bas <- array(0, dim=c(I,P,numknots))
      basderiv <- array(0, dim=c(I,P,numknots))
      basdummy <- array(0, dim=c(I,P,numknots)) 
      
      
      for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,diffunct(0,1,datb[i,p]),-cov[,p,i])}
      for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,derdiffunct(0,1,datb[i,p]),0*cov[,p,i])}
      
      for (i in 1:I){for (p in 1:P){
        datdum <- datb[i,p]-1
        #basdummy[i,p,]<- c(1,datdum)
        basdummy[i,p,]<- c(1,diffunct(0,1,datdum),-cov[,p,i])  ### here covariates
      }}
      
      ### extension covariates
      #basderiv[,,3]
      #basdummy[,,2]
      par1 <- -itemsums  ### 
      par <- c(par1[1],1,0*cov[,p,1]) ### 1:slope difficulty
      #for(i in 2:I) par <- c(par,par1[i],1,0*cov[,p,i])
      #par <- c(-10,1,0*cov[,p,1]) ### 1:slope difficulty!!!
      #for(i in 2:I) par <- c(par,-10,1,0*cov[,p,i])
      par <- c(par,2) ### stdmixt
    } # end not splines
    
    #par <-c(-1,sl,0,0,0,2)
    
    pentype<-"ridge"  ## dummies
    cepsilon<-.05     ##dummies
    
    if(startind == "yes") par <- start 
    
    
    if(der!="der") {fittrw <- optim(par, LoglikIntCovGlobal, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                                    I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf, hessian = FALSE,
                                    numcov=numcov,penvector=penvector,pnorm=pnorm,
                                    pentype=pentype,cepsilon=cepsilon,control = list(maxit=maxit)) }
    
    if(der=="der") {fittrw <- optim(par, LoglikIntCovGlobal, gr = derLoglikIntCovGlobal,lambdpenin ,lambdc ,numknots=numknots,
                                    dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                    lower = -Inf, upper = Inf, hessian = hessiander,
                                    numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon,
                                    control = list(maxit=maxit)) }
    
    
    ll <- length(par)
    #itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,I)) ##I x numknots matrix
    itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,1))
    stdmixt <- fittrw$par[ll]
    
    
    itemmatrixstd<-  0
    itemmatrixtval<-0
    
    if(hessiander=="TRUE"){
      covm<-diag(solve(fittrw$hessian))
      stdev<-sqrt(covm)
      itemmatrixstd<-  stdev[1:(ll-1)]
      itemmatrixtval<-itemmatrixest/itemmatrixstd
      stdstdmixt<-stdev[ll]
    }
    
    #loglikorig<--fittrw$value
    
    derivfin <-0  ## dummy
    #derivative here
    
    
    ### derivative omitted
    #derivfin <-derLoglikBsplIPen(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
    #                             indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,
    #                             numcov=numcov,
    #                             penvector=penvector,pnorm=pnorm)
    
    ## logorig
    #loglikorig <- -LoglikIntExtBsplIN(fittrw$par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
    #                                  indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
    #-LoglikIntExtBsplIN(par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
    #                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
    loglikorig <- -fittrw$value 
    ### positive stddev
    stdmixt<-abs(stdmixt) 
    fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
    ####
    AIC <-  -2*loglikorig +2*length(par)
    BIC <-0
    #BIC <-  -2*loglikorig +2*log(length(par))
    
    #########
    newList <- list("parmatrix" = itemmatrixest,"stdmixt"=stdmixt,"stdstdmixt"=stdstdmixt,
                    "par"=fittrw$par, "numbasis"=numknots, 
                    "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                    "AIC"=AIC, "BIC"=BIC,"derivative"=derivfin,"itemmatrixstd"=itemmatrixstd,
                    "itemmatrixtval"=itemmatrixtval)
    return(newList)
    
  }
########################################################

#LoglikIntCov <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
#                        indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
#                       numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon){

LoglikIntCovGlobal <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                              indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy, 
                              numcov=numcov,penvector=penvector,pnorm=pnorm,pentype=pentype,cepsilon=cepsilon){  
  
  ### not item specific parameters
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  
  # numcov penalty on covariates
  # penvector new!!
  # pnorm: chosen norm: pnorm=2 is Euclidean-squared
  # cepsilon: modification of norm,  =0 no modification
  #print(par)
  
  P <- dim(dat)[2]
  
  numpar <- length(par) 
  sdnormal<- par[numpar]
  #paroriginal <- par[1:numpar-1]
  
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVecGlobal, mu = 0, sd = sdnormal,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots)
    
    #if (int<(10)^{-300}) int<-(10)^{-300} #log((10)^{-300}) 
    ind[p]<-log(int)
    #int<-abs(int)  # 
    
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
    
  } #  end p
  
  dimBspl<-length(par)-1
  parBspl <- par[1:dimBspl]
  itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
  
  
  ####penalty increase
  #lambdpenin <- 0
  if (lambdpenin >0){
    
    upp<- numknots-1
    peninc<-0
    for (i in 1:I){
      for (l in 1:upp ){dif  <-itemmatrix[i,l+1]- itemmatrix[i,l]
      peninc <-peninc + dif^2 *indicatt(-dif , 0)}}
    
    #pen <- dif%*%dif
    sum <- sum -  lambdpenin*peninc ### + weil negativ , jetzt - ????
    
  } #end if
  ############
  
  #### penalty common difficulties
  #lambdc <- 3
  if (lambdc >0){
    upp1<- I-1
    upp2<- numknots-1
    penc<-0
    for (i in 1:upp1){
      for (l in 1:upp2 ){dif  <-(itemmatrix[i+1,l+1]- itemmatrix[i+1,l])-(itemmatrix[i,l+1]- itemmatrix[i,l])
      penc <- penc +dif*dif}}
    
    ##### funktioniert: with matrix penc = pencc jetzt auskommentiert
    
    #upp1<- I-1
    #upp2<- numknots-1
    #counth <- 0
    #for (i in 1:upp1){
    #  for (l in 1:upp2 ){
    #  counth <- counth +1
    #  penM <- matrix(0,1,length(par))
    #  penM[(i-1)*numknots+l+1]<-1
    #  penM[(i-1)*numknots+l]<--1
    #  penM[(i)*numknots+l+1]<--1
    #  penM[(i)*numknots+l]<- 1
    #  if (counth <= 1) pennew <- penM
    # if(counth > 1) pennew <- rbind(pennew,penM)
    # }}
    #pencc <- t(pennew%*%par)%*% (pennew%*%par) 
    
    ####ende auskommentiert
    
    if (numknots == 2) {penc<-0
    for (i in 2:I) penc<- penc + (itemmatrix[i,2]-itemmatrix[i-1,2])^2
    
    ## alternative pennew=penc
    M<- matrix(0,I-1,I)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,I-2))
    M[I-1,]<-c(rep(0,I-2),d1)
    I2 <- I-2
    for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    pennew <- t(itemmatrix[,2])%*%t(M)%*%M%*%itemmatrix[,2]
    ###
    }
    
    
    sum <- sum -  lambdc*penc ### 
    #sum <- sum -  lambdc*pennew 
  } #end if
  
  
  #### penalty covariates
  pencov<-sum(penvector)
  red<-numknots-numcov+1
  if (pencov >0){
    penm<-as.matrix(itemmatrix[,red:numknots])
    #print("numknots")
    #print(numknots)
    
    #print(penm)   #################################################!
    #pemmdif<-penm[2,]-penm[1,]  ## difference matrix
    #for (i in 3:I) pemmdif<-rbind(pemmdif,penm[i,]-penm[i-1,])
    
    ### new generate pairs
    cpairs<-I*(I-1)/2
    M<-matrix(0,cpairs,I)
    count<-0
    for (i in 1:(I-1)){
      for (j in (i+1):I) {count<-count+1
      M[count,i]<-1
      M[count,j]<--1}} 
    diffm<-M%*%penm
    #diffv<-as.vector(diffm)
    #############
    
    
    if(pentype=="ridge"){
      #normvect<-Norm(pemmdif[,1], p = pnorm)
      #for (v in 2:numcov){normvect<-c(normvect,Norm(pemmdif[,v], p = pnorm))}   
      
      #normvect<-(Norm(diffm[,1], p = pnorm))^pnorm
      # for (v in 2:numcov){normvect<-c(normvect,(Norm(diffm[,v], p = pnorm))^pnorm)}
      normvect<-t(diffm[,1])%*%diffm[,1]
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,t(diffm[,v])%*%diffm[,v])}   
      
      sum <- sum -  penvector%*%normvect
    } ## end ridge
    
    if(pentype=="grlasso"){
      
      c<-cepsilon   ## with c, modified norm!
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}
      
      
      sum <- sum -  penvector%*%normvect
    } ## end lasso
    
    if(pentype=="grlassosel"){
      
      c<-cepsilon   ## with c, modified norm!
      
      # first term
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2))}   
      
      normvect<-NormShiftc((t(diffm[,1])%*%diffm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(diffm[,v])%*%diffm[,v]),c))}
      
      sum <- sum -  penvector%*%normvect#*(cpairs^(1/2))
      
      
      # second term
      #normvect<-(t(penm[,1])%*%penm[,1])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(penm[,v])%*%penm[,v])^(1/2))}
      
      
      #print("v")
      #print(v)
      
      normvect<-NormShiftc((t(penm[,1])%*%penm[,1]),c)
      if(numcov > 1)for (v in 2:numcov){normvect<-c(normvect,NormShiftc((t(penm[,v])%*%penm[,v]),c))
      #print("v")
      #print(v,numcov)
      }
      
      sum <- sum -  penvector%*%normvect#*(numcov^(1/2))
      
      #normvect<-(t(diffm[,1])%*%diffm[,1])^(1/2)+(t(penm[1,])%*%penm[1,])^(1/2)
      #for (v in 2:numcov){normvect<-c(normvect,(t(diffm[,v])%*%diffm[,v])^(1/2)+(t(penm[v,])%*%penm[v,])^(1/2))}   
      
      
    } ## end grlassosel
    
  } ## end pen cov
  
  
  
  sum<- -sum
  return(sum)}
###############

###########################
#prodfctIntExtBsplINVec  <-function(theta,datitem=datitem,I=I,par=par,slope=slope,
#                                   indicator=indicator, bas,basderiv,basdummy,pers=pers,numknots=numknots){
#
prodfctIntExtBsplINVecGlobal  <-function(theta,datitem=datitem,I=I,par=par,slope=slope,
                                         indicator=indicator, bas,basderiv,basdummy,pers=pers,numknots=numknots){  
  
  #### with std mixt only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  prod <-1    
  
  dimBspl<-length(par)-1
  parBspl <- par[1:dimBspl]
  #itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
  itemmatrix<-  parBspl
  for (i in 1:I){
    #iteml <- itemmatrix[i,]
    iteml <- parBspl
    item=i
    if(indicator[i]=="C"){
      prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                       item=item,pers=pers)} 
    if(indicator[i]!="C"){
      prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                            item=item,pers=pers,indicator[i])}
    
    #if(indicator[i]=="D0"){
    #  prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
    #                                        item=item,pers=pers,indicator[i])}
    #if(indicator[i]=="D1"){
    #  prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
    #                                        item=item,pers=pers,indicator[i])} ##idetisch mit D0
  }
  
  return(prod)}  
#########################

# same as 
#derLoglikBsplIPen <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
#                             indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,numcov=numcov,
#                             penvector=penvector,pnorm=pnorm){
## only pentype = "ridge", cepsilon = 0.05 added


derLoglikIntCovGlobal <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                 indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,numcov=numcov,
                                 penvector=penvector,pnorm=pnorm,pentype = "ridge", cepsilon = 0.05){
  
  ### derivative, only continuous
  ### only ridge
  #parmatr <- matrix(par,I,2)
  # derivfin <-derLoglikBsplI(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
  #indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  
  slope <- 1 ### slope fixed  weight on theta!
  P <- dim(dat)[2]
  
  lpar <- length(par)
  der <- matrix(0,lpar,1)
  
  parmatr <- t(matrix(par[1:lpar-1],numknots,I))
  stdmixt <- par[lpar]
  
  # basis functions prodfctIntExtBsplINVec 
  
  for (p in 1:P){
    datitem <- dat[,p]
    
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots) 
    #### modified for global
    #Ilimit<-I
    #if(lpar <= I+1)Ilimit<-1
    ######
    
    for (l in 1:I){
      for (j in 1:numknots){ 
        comp <- l
        pers<-p
        #count <- (l-1)*numknots+j
        count <- j
        #der[l,1] <-der[l,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
        #                    indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp)/int } 
        der[count,1] <-der[count,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,
                                                  numknots=numknots,I=I,par=par, slope=slope,indicator=indicator,order = 10, bas=bas,
                                                  basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp, j=j)/int } # end j 
    } #end l
  } # end p
  
  
  
  
  
  
  #### now numeric !!!!!
  deltan<- 0.00005
  vect1<- par
  vect1[lpar] <- vect1[lpar]+deltan
  #loglikorig  <- LoglikN(par,dat,I,indicator, lin=lin)
  #loglikorignow <- LoglikN(vect1,dat,I,indicator, lin=lin)
  
  loglikorig<-LoglikIntExtBsplIN(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                 indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  loglikorignow<-LoglikIntExtBsplIN(vect1,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  derivnow <- (loglikorignow-loglikorig)/deltan
  der[lpar,1]<- -derivnow
  
  
  ##### end numeric
  
  
  
  #lambdc <- 3
  
  if (lambdpenin >0){
    numknots
    M<- matrix(0,numknots-1,numknots)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,numknots-2))
    numknots1 <- numknots-1
    if(numknots >2)for(l in 2:numknots1) M[l,]<- c(rep(0,l-1),d1,rep(0,numknots-l-1))
    dum <- diag(1,I)
    Mtot <- kronecker(dum,M)
    Mtot<- cbind(Mtot,rep(0,dim(Mtot)[1]))
    pen1 <- 2*t(Mtot)%*%(Mtot%*%par*indicattmatrix(-Mtot%*%par,0))
    
    der <- der - lambdpenin*pen1
  }
  ## end penin
  #dim(Mtot)
  
  if (lambdc >0){
    
    if (numknots > 2){
      upp1<- I-1
      upp2<- numknots-1
      counth <- 0
      for (i in 1:upp1){
        for (l in 1:upp2 ){
          counth <- counth +1
          penM <- matrix(0,1,length(par))
          penM[(i-1)*numknots+l+1]<-1
          penM[(i-1)*numknots+l]<--1
          penM[(i)*numknots+l+1]<--1
          penM[(i)*numknots+l]<- 1
          if (counth <= 1) pennew <- penM
          if(counth > 1) pennew <- rbind(pennew,penM)
        }}
      pencc <- t(pennew%*%par)%*% (pennew%*%par) 
      
      deriv1 <- 2*t(pennew)%*%pennew%*%par
      der <- der -  lambdc*deriv1
    }
    
    if (numknots == 2) {
      M<- matrix(0,I-1,I)
      d1 <- c(-1,1)
      M[1,]<- c(d1,rep(0,I-2))
      M[I-1,]<-c(rep(0,I-2),d1)
      I2 <- I-2
      for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
      
      deriv <- 2*t(M)%*%M%*%parmatr[,2]
      for (l in  1:I){der[l*2,1]<-der[l*2,1] -lambdc*deriv[l]}
      
    }} #end if
  
  pencov<-sum(penvector)
  
  
  ####################### new penalty covariates
  
  if (pencov >0){
    red<-numknots-numcov+1
    penm<-as.matrix(parmatr[,red:numknots])
    
    #if (numknots == 2) {
    #M<- matrix(0,I-1,I)
    #d1 <- c(-1,1)
    #M[1,]<- c(d1,rep(0,I-2))
    # M[I-1,]<-c(rep(0,I-2),d1)
    # I2 <- I-2
    # for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    
    ### new generate pairs
    cpairs<-I*(I-1)/2
    M<-matrix(0,cpairs,I)
    count<-0
    for (i in 1:(I-1)){
      for (j in (i+1):I) {count<-count+1
      M[count,i]<-1
      M[count,j]<--1}} 
    #diffm<-M%*%penm 
    
    #}
    
    
    #fact1<-(t(M%*% as.matrix(penm[,1]))%*%(M%*% as.matrix(penm[,1])))^(-1/2)/2
    #deriv1 <- penvector[1]*fact1[1]*(2*t(M)%*%M%*%penm[,1])
    deriv1 <- penvector[1]*(2*t(M)%*%M%*%penm[,1])
    
    seq<-c(2:numcov)
    if (numcov >1)for (v in seq){          ####### 24.2
      # fact1<-(t(M%*% as.matrix(penm[,v]))%*%(M%*% as.matrix(penm[,v])))^(-1/2)/2
      # derivn <- penvector[v]*fact1[1]*(2*t(M)%*%M%*%penm[,v])
      derivn <- penvector[v]*(2*t(M)%*%M%*%penm[,v])
      deriv1<-cbind(deriv1,derivn)   
    } ## end if numcov
    
    derivm<-cbind(matrix(0,I,2),deriv1)
    derivv<- as.matrix(c(t(derivm)))
    
    
    #der <-  der -  as.matrix(c(rep(0,((numknots-numcov)*I)),derivv,0))
    der<- der-as.matrix(c(derivv,0))
    
  }  ### end pencov
  
  
  
  
  der <- -der
  return(der)}
###############











### alt
NormShifted <-function(x,eps){
  y<- 2*x**2
  
  if(abs(x)>eps)y<-sqrt(abs(x)-eps)+ 2*eps**2
  return(y)}

eps<-.05 
xpl<-seq(-3,3,.01)
ypl<-NormShifted(xpl[1],eps)
for(i in 2:length(xpl))ypl<-c(ypl,NormShifted(xpl[i],eps))
plot(xpl,ypl) 


