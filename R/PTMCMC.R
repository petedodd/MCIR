## new as yet undocumented one with adaptive metropolis, and PTMCMC

## CONTENTS:
## MHstepper
## MHadapt
## (tr)
## (getEmpSig)
## PTMCMC
## PTMCMCcontinue


## metropolis hastings
MHstepper <- function(x0,sigprop,Tchain,LL){
    DIM <- length(x0)
    X <- matrix(NA,nrow=DIM,ncol=Tchain)
    LLchn <- rep(NA,Tchain)
    Vbef <- LL(x0)
    LLchn[1] <- Vbef
    X[,1] <- x0
    acc <- 0
    if(is.null(dim(sigprop))){
        cat('making SIGMA...\n')
        sigmat <- diag(sigprop,DIM)
    } else {
        sigmat <- sigprop
    }
    for( i in 2:(Tchain)){
        if(!i%%ceiling(Tchain/10)) cat(round(1e2*i/Tchain),' % complete \n')
        x <- X[,i-1] + mvtnorm::rmvnorm(n=1,mean=rep(0,DIM),sigma=sigmat) #matrix(sigprop * rnorm(DIM),nrow=DIM,ncol=1)
        Vnow <- LL(x)
        test <- ( log(runif(1)) <= (Vnow - Vbef) )
        if(test){
            X[,i] <- x
            Vbef <- Vnow
            acc <- acc + 1
        } else {
            X[,i] <- X[,i-1]
        }
        LLchn[i] <- Vbef
    }    
    ## end
    acc <- round(1e2*acc/Tchain)
    cat('acc = ',acc,'\n')
    list(X = X,LLchn = LLchn,acc=acc)
}

## MH with heuristic adaptation
MHadapt <- function(x0,sigprop,Tchain,LL,lots=10,nadapt=1e3,targetacc=.4){
    nsmall <- ceiling(nadapt/lots)
    s1 <- sigprop
    sz <- az <- rep(NA,lots)
    XZ <- LZ <- list()
    for(i in 1:lots){
        go <- MHstepper(x0=x0,sigprop=s1,Tchain=nsmall,LL=LL)
        az[i] <- go$acc
        XZ[[i]] <- go$X
        LZ[[i]] <- go$LLchn
        cat('i = ',i,', acc = ',az[i],' %\n')
        sz[i] <- s1
        s1 <- s1*((go$acc*1e-2)/targetacc)
    }
    XZ <- do.call('cbind',XZ)
    list(sigprop=s1,acz=az,sz=sz,XZ=XZ,LZ=LZ)
}

## testmc <- MHadapt(x0=xx0,.5,1e3,LL=function(x)LAISLL(matrix(x,ncol=1)))

## trace function
tr <- function(x)sum(diag(x))           #trace

## covariance calculation, including shrinkage estimators
getEmpSig <- function(X,shrink='none',verbose=FALSE){
    DIM <- nrow(X)
    n <- ncol(X)
    if(any(is.na(X))){
        sna <- which.min(apply(X,2,function(x)any(!is.na(x)))) #start of NAs
    } else {sna <- ncol(X)+1}
    sna <- sna-1
    XS <- X[,1:sna]                     #X w/o the NAs
    mn <- rowMeans(XS)
    XS <- XS - matrix(mn,ncol = sna,nrow=nrow(XS),byrow=FALSE) #centre
    S <- XS %*% t(XS) / sna                                   #empiric
    if(shrink=='none') return(S)
    F <- diag(tr(S)/DIM,DIM)            #shrink target
    if(shrink=='OAS'){                  #
        ## OAS shrinkage -- from scikit-learn
        rho <- .1
        alpha <- mean(S^2)
        mu <- tr(S)/DIM
        num <- alpha + mu^2
        den <- (sna + 1) * (alpha - mu^2/DIM)
        rho <- 1
        if(den>0) rho <- min(num/den,1)
        SS <- (1-rho) * S + rho * diag(mu,DIM)
        if(verbose) cat('rho = ',rho,'\n')
    }
    return(SS)
}


## adaptive-metropolis
AMMHstepper <- function(x0,sigprop,Tchain,LL,adaptstart=1e2,lag=500,bet=.1,
                        cool=10,colp=1){
    DIM <- length(x0)
    cool <- (cool/Tchain)        #length-scale for cooling, input cool=no e-folds 
    X <- matrix(NA,nrow=DIM,ncol=Tchain)
    LLchn <- rep(NA,Tchain)
    Vbef <- LL(x0)
    LLchn[1] <- Vbef
    scale <- scale0 <- abs(Vbef/10)                       #scaling for cooling
    if(cool==0) scale <- scale0 <- 1                      #turn off cooling
    Vbef <- Vbef/scale
    X[,1] <- x0
    acc <- 0
    if(is.null(dim(sigprop))){
        cat('making SIGMA...\n')
        sigmat <- diag(sigprop,DIM)
    } else {
        sigmat <- sigprop
    }
    for( i in 2:(Tchain)){
        if(!i%%ceiling(Tchain/10)){
            cat(round(1e2*i/Tchain),' % complete \n')
            cat('scale = ',scale,'\n')
            cat('LL = ',LLchn[i-1],'\n')
        }
        x <- X[,i-1] + mvtnorm::rmvnorm(n=1,mean=rep(0,DIM),sigma=sigmat) #matrix(sigprop * rnorm(DIM),nrow=DIM,ncol=1)
        Vnow <- LL(x)/scale
        test <- ( log(runif(1)) <= (Vnow - Vbef) )
        if(test){
            X[,i] <- x
            Vbef <- Vnow
            acc <- acc + 1
        } else {
            X[,i] <- X[,i-1]
        }
        LLchn[i] <- Vbef * scale
        scale <- scale0*exp(-i*cool*colp) +(1-exp(-i*cool))^colp
        if(i > adaptstart){
            if(runif(1)<bet){
                sigmat <- diag(1e-2/DIM,DIM) #safe one
            } else {
                sigmat <- 2.38^2*getEmpSig(X[,max(1,i-lag):i],shrink='OAS')/DIM
                if(tr(sigmat)<1e-6) sigmat <- diag(1e-2/DIM,DIM) #safe one
            }
            ## print(sigmat)
        }
    }    
    ## end
    acc <- round(1e2*acc/Tchain)
    cat('acc = ',acc,'\n')
    cat('final scale = ',scale,'\n')
    list(X = X,LLchn = LLchn,acc=acc,sigmat=sigmat)
}

## testmc <- AMMHstepper(x0=sola$par,.2,Tchain=3e3,LL=function(x)LAISLL(matrix(x,ncol=1)),adaptstart = 5e3,cool=10,colp=1)

## parallel tempered MCMC with adaptive metropolis
PTMCMC <- function(x0,sigprop,ntemps=10,nwalkers=10,niter=1e3,swapn=1,tempfac=2,
                   adaptstart = Inf,lag=250,bet=.1, #adaptation
                   LL){
    nchains <- ntemps*nwalkers
    if(is.null(dim(x0))){
        DIM <- length(x0)
        x0 <- matrix(rep(x0,nchains) + rnorm(n=DIM*nchains,sd=mean(sigprop)),
                     nrow=DIM,ncol=nchains)
    } else { DIM <- nrow(x0)}
    X <- array(NA,dim=c(DIM,nchains,niter))
    LLchn <- array(NA,dim=c(nchains,niter))
    stemps <- tempfac^((0:(ntemps-1))/2)
    cat('temperature ladder = ',stemps,'\n')
    temps <- rep(stemps,each = nwalkers)
    X[,,1] <- x0
    Vbef <- apply(x0,2,LL)
    Vbef <- Vbef / temps
    LLchn[,1] <- Vbef
    acc <- rep(0,ntemps)
    acc2 <- 0
    if(is.null(dim(sigprop))){
        if(length(sigprop)==1){
            cat('growing siprop...\n')
            sigprop <- rep(sigprop,ntemps)
            sigprop <- sigprop * stemps^1
        }
        cat('making SIGMA...\n')
        sigmat <- array(NA,dim=c(DIM,DIM,ntemps))
        for(i in 1:ntemps)
            sigmat[,,i] <- diag(sigprop[i],DIM)
    } else {
        if(!all(dim(sigprop)==c(DIM,DIM,ntemps))){
            stop('If sigprop array, dims must be npars,npars,ntemps!')
        }
        sigmat <- sigprop
        cat('using SIGMA supplied...\n')
    }
    x <- array(NA,dim=c(DIM,nchains))
    pflag <- TRUE                       #for displaying adaptation
    for( i in 2:niter){
        if(!i%%ceiling(niter/10)) cat(round(1e2*i/niter),' % complete \n')
        ## normal proposing
        for(j in 1:ntemps)
            x[,(1+(j-1)*nwalkers):(j*nwalkers)] <- t( mvtnorm::rmvnorm(n=nwalkers,mean=rep(0,DIM),sigma=sigmat[,,j]) ) 
        x <- X[,,i-1] + x
        Vnow <- apply(x,2,LL)
        Vnow <- Vnow / temps
        test <- ( log(runif(nchains)) <= (Vnow - Vbef) ) #accept
        X[,,i] <- X[,,i-1]              #previous
        if(any(is.na(test))) print(test)
        X[,test,i] <- x[,test]          #set to new
        for(j in 1:ntemps)
            acc[j] <- acc[j] + sum(test[1:nwalkers + (j-1)*nwalkers])/nwalkers
        Vbef[test] <- Vnow[test]
        ## print(mean(test))
        ## swapping
        if(!i%%swapn){
            swap <- sample(nchains,2)        #propose swap
            Vnow <- Vbef
            ## cat(i,' :------\n')
            ## print(swap)
            ## print(Vnow[swap])
            Vnow[rev(swap)] <- Vnow[swap] * temps[swap]/ temps[rev(swap)]  #real LL
            test <- ( log(runif(1)) <= sum(Vnow[swap] - Vbef[swap]) ) #accept
            ## print(Vnow[swap])
            ## print(Vbef[swap])
            if(test){
                ## print(X[1,swap,i])
                X[,rev(swap),i] <- X[,swap,i]
                ## print(X[1,swap,i])
                ## xtmp <- X[,swap[1],i]
                ## X[,swap[1],i] <- X[,swap[2],i]
                ## X[,swap[2],i] <- xtmp
                Vbef[swap] <- Vnow[swap]
                acc2 <- acc2 + 1
            }
        }
        LLchn[,i] <- Vbef
        ## adaptation
        if(i > adaptstart){
            if(pflag){
                cat('Starting adaptive phase...\n')
                pflag <- FALSE
            }
            for(j in 1:ntemps){
                if(runif(1)<bet){
                    sigmat[,,j] <- diag(1e-2/DIM,DIM) #safe one
                } else {
                    xtmp <- X[,(1+(j-1)*nwalkers):(j*nwalkers), max(1,i-lag):i]
                    xtmp <- matrix(xtmp,nrow=DIM,ncol=prod(dim(xtmp)[2:3]))
                    sigmat[,,j] <- 2.38^2*getEmpSig(xtmp,shrink='OAS')/DIM
                    if(tr(sigmat[,,j])<1e-7) sigmat[,,j] <- diag(1e-2/DIM,DIM) #safe one
                }
            }
            ## print(sigmat)
        }
    }    
    ## end
    acc <- round(1e2*acc/niter)
    acc2 <- round(1e2*acc2*swapn/niter)
    print(acc)
    ## cat('acc = ',acc,'\n')
    cat('acc2 = ',acc2,'\n')
    list(X = X,LLchn = LLchn,acc=acc,acc2=acc2,sigmat=sigmat,
         ntemps=ntemps,nwalkers=nwalkers,niter=niter,swapn=swapn,tempfac=tempfac,
         lag=lag,bet=bet,adaptstart=adaptstart,LL=LL)
}


## for carrying on PTMCMC
PTMCMCcontinue <- function(previous,niter=NULL){
    if(is.null(niter)) niter <- previous$niter
    xinit <- previous$X[,,dim(previous$X)[3]]
    ##  do new run
    newrun <- PTMCMC(x0=xinit,sigprop=previous$sigmat,
                     ntemps=previous$ntemps,nwalkers=previous$nwalkers,
                     niter=niter,
                     swapn=previous$swapn,tempfac = previous$tempfac,
                     adaptstart = previous$adaptstart,lag=previous$lag,
                     bet=previous$bet,                 
                     LL=previous$LL)
    ## join
    cat('joining to previous output...\n')
    ## chain updates
    dz1 <- dim(previous$X)
    dz2 <- dim(newrun$X)
    Xboth <- array(NA,dim=c(dz1[1:2],dz1[3]+dz2[3]))
    LLchnb <- array(NA,dim=dim(Xboth)[2:3])
    Xboth[,,1:dz1[3]] <- previous$X
    Xboth[,,(dz1[3]+1):dim(Xboth)[3]] <- newrun$X
    LLchnb[,1:dz1[3]] <- previous$LLchn
    LLchnb[,(dz1[3]+1):dim(Xboth)[3]] <- newrun$LLchn
    newrun$X <- Xboth
    newrun$LLchn <- LLchnb
    ## other updates
    newrun$acc <- (newrun$acc*dz1[3] + previous$acc*dz2[3])/dim(Xboth)[3]
    newrun$acc2  <- (newrun$acc2*dz1[3] + previous$acc2*dz2[3])/dim(Xboth)[3]
    newrun$niter <- newrun$niter + dz1[3]
    ## return
    return(newrun)
}



## ## tests...

## NW <- 20
## NT <- 1
## NI <- .25*1e3
## cat('evaluations = ',NW*NT*NI,'\n')

## ## todo
## xx0 <- t(randomLHS(n=NW*NT,k=nrow(rngs))-.5)*1e-1
## ## x0 <- rep(0,nrow(rngs))

## testmc <- PTMCMC(x0=xx0,sigprop=1e-1/nrow(rngs),
##                  ntemps=NT,nwalkers=NW,niter=NI,swapn=10,tempfac = 2,
##                  adaptstart = 1e2/NW,lag=250,bet=.1,                 
##                  LL=function(x)LAISLL(matrix(x,ncol=1)))


## mcl <- list()
## for(i in 1:NW)
##     mcl[[i]] <- mcmc(t(testmc$X[,i,]))

## xyplot(as.mcmc.list(mcl),layout=c(2,7))


## burnin <- 5e2
## lowtemp <- testmc$X[,,-c(1:burnin)]
## lowtemp <- matrix(c(lowtemp[,1:NW,]),nrow=nrow(rngs),ncol=NW*dim(lowtemp)[3])
## lowtemp <- t(lowtemp)
## if(nrow(lowtemp)>1e3){
##     pairs(lowtemp[sample(nrow(lowtemp),1e3),],cex=.1)
## } else{ pairs(lowtemp,cex=.1) }
