## todo:
## - correct progress thing, which is wrong
## - collect more data for INS post-processing
## - INS postprocessing
enest <- function(Nlive,niter,DIM,expand=1.5,toteval,rejtol=1e-3,LL){
    initlive <- lhs::randomLHS(n=Nlive,k=DIM)    #initial points
    initlogl <- apply(initlive,1,LL)
    ord <- order(initlogl)
    liveall <- matrix(nrow=Nlive+niter,ncol=DIM)
    LLall <- rep(NA,Nlive+niter)
    liveall[1:Nlive,] <- initlive[ord,]
    LLall[1:Nlive] <- initlogl[ord]
    maxit <- ceiling(1/rejtol)
    neval <- Nlive
    kmax <- 1
    accrec <- rep(NA,niter)
    INSL <- INSX <- INSV <- INSM <- list()
    cat('initial live points calculated...\n')
    cat('...beginning iterations...\n')
    for(i in 1:niter){
        if(!i%%(ceiling(niter/10)))
            cat('iteration ',i,' / ',niter,'...\n')
        k <- 0
        accept <- FALSE
        LIVE <- liveall[(i):(Nlive+i-1),]
        mn <- colMeans(LIVE)
        V <- chol(cov(LIVE))
        INSM[[i]] <- mn
        INSV[[i]] <- V
        INSX[[i]] <- INSL[[i]] <- list()
        while(k<maxit & !accept){
            xnw <- mn + rnorm(DIM) %*% V * expand
            if(any(xnw>1)|any(xnw<0)){
                bad <- xnw>1 | xnw<0
                xnw[bad] <- runif(sum(bad))
            } #stop('out of U!')
            LLnew <- LL(xnw)
            INSL[[i]][[k]] <- LLnew
            INSX[[i]][[k]] <- xnw
            neval <- neval + 1
            k <- k + 1
            kmax <- max(kmax,k)
            if(LLnew > LLall[i]) break;
        }
        if(k>maxit) break;
        if(neval>toteval) break;
        LLall[i+Nlive] <- LLnew
        liveall[i+Nlive,] <- xnw
        accrec[i] <- k
    }
    ord <- order(LLall)
    list(X=liveall[ord,],LLall=LLall[ord],kmax=kmax,neval=neval,
         Nlive=Nlive,accrec=accrec,
         INSM=INSM,INSV=INSV,INSX=INSX,INSL=INSL)
}


postprocess <- function(inpt,INS=FALSE){
    Nlive <- inpt$Nlive
    if(!INS){
        niter <- sum(!is.na(inpt$LLall))
        P <- inpt$X[1:niter,]
        Lz <- inpt$LLall[1:niter]
        Xz <- exp(-(1:niter)/Nlive)
        wz <- -diff(Xz)
        wz <- c(1-Xz[1],wz)
        LLmax <- max(Lz)
        wz <- exp(Lz-LLmax) * wz
        Z <- sum(wz) + mean(exp(Lz-LLmax)) * rev(Xz)[1]
        wz <- wz / Z
        Z <- Z * exp(LLmax)
    } else {                            #INS version
        
    }
    list(Z=Z,w=wz,X=P,LLmax=LLmax)
}


## logit <- function(x)log(x/(1-x))
## ilogit <- function(x)1/(1+exp(-x))


## LLtest <- function(x){
##     y <- logit(x)
##     ans <- -sum(y^2/2) -sum(log(x*(1-x)))
##     ## ans <- -sum(y^2/2 + y - 2*log(1+exp(-y))) #SEEMS WRONG
##     ans
## }


## LLtestv <- function(x){
##     sapply(x,LLtest)
## }

## LLtestve <- function(x){
##     exp(sapply(x,LLtest))
## }

## integrate(f=LLtestv,lower=1e-6,upper=.999)
## integrate(f=LLtestve,lower=1e-6,upper=.999)
## integrate(f=function(x)exp(-x^2/2),lower=-Inf,upper=Inf) #sqrt(2*pi)=2.58

## LLtest1 <- function(x){
##     ans <- 2/3
##     if(abs(x-.5)<.25)
##         ans <- 2*ans
##     ans
## }

## LLtest <- function(x){
##     sum(sapply(x,function(x)log(LLtest1(x))))
## }

## LLtest(c(.5,.5))
         

## test <- runif(1e4)
## for(i in 1:1e4)test[i] <- LLtest1(test[i])
## sum(test)/1e4

## testch <- enest(Nlive=1e3,niter=3e3,DIM=20,
##                      toteval = 1.5e4,rejtol=1e-3,
##                      LL=LLtest)


## str(testch)
## plot(-na.omit(-testch$LLall),type='s')#,log='y')

## plot(testch$acc,log='y')

## tcho <- postprocess(testch)
## str(tcho)
## tcho$Z/(2*pi)^(20/2)
## pks <- sample(nrow(tcho$X),1e3,replace=TRUE,prob=tcho$w/sum(tcho$w))

## pairs((tcho$X[pks,]),cex=.1)
## plot(logit(tcho$X[pks,1]),logit(tcho$X[pks,2]))
## var(logit(tcho$X[pks,1]))

## tmp <- logit(tcho$X[pks,])
## colMeans(tmp)
## pairs(tmp[,1:10],cex=.1)
## sum(tcho$w)^2/sum(tcho$w^2)
