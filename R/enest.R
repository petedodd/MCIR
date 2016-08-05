## todo:

inelipse <- function(x,U,m){
    if(!is.null(dim(x))){
        m <- matrix(m,nrow=nrow(x),ncol=ncol(x),byrow=TRUE)
        x <- x-m                            #change for matrix
        x <- x %*% solve(U)
        ans <- rowSums(x^2)<1
    } else {
        x <- x-m                            #change for matrix
        x <- x %*% solve(U)
        ans <- sum(x^2)<1    
    }
    ans    
}

rsph <- function(n,DIM){
    ans <- matrix(rnorm(n=DIM*n),nrow=n,ncol=DIM)
    ans <- ans / sqrt(rowSums(ans^2))
    ans * runif(n)^(1/DIM)
}

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
    INSL <- INSX <- INSV <- INSM <- INSI <- list()
    cat('initial live points calculated...\n')
    cat('...beginning iterations...\n')
    for(i in 1:niter){
        if(!i%%(ceiling(niter/10)))
            cat('iteration ',i,' / ',niter,'...\n')
        k <- 0
        accept <- FALSE
        LIVE <- liveall[(i):(Nlive+i-1),]
        mn <- colMeans(LIVE)
        V <- chol(cov(LIVE)) * expand
        INSM[[i]] <- mn
        INSV[[i]] <- V
        INSX[[i]] <- INSL[[i]] <- list()
        while(k<maxit & !accept){
            ## xnw <- mn + rnorm(DIM) %*% V
            xnw <- mn + rsph(1,DIM) %*% V 
            if(any(xnw>1)|any(xnw<0)){
                bad <- xnw>1 | xnw<0
                xnw[bad] <- runif(sum(bad))
            } #stop('out of U!')
            LLnew <- LL(xnw)
            k <- k + 1
            kmax <- max(kmax,k)
            neval <- neval + 1
            INSL[[i]][[k]] <- LLnew
            INSX[[i]][[k]] <- xnw
            INSI[[neval]] <- i
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
         INSM=INSM,INSV=INSV,INSX=INSX,INSL=INSL,INSI=INSI)
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
        vZ=Z^2*sum((wz-1/length(wz)^2)^2) #maybe?? todo
    } else {                            #INS version
        MNZ <- inpt$INSM
        IZ <- do.call('rbind',inpt$INSI)
        DIM <- ncol(inpt$X)
        XZ <-  matrix(unlist(inpt$INSX),nrow=length(IZ),ncol=DIM,byrow=TRUE) 
        LZ <- unlist(inpt$INSL)
        VZ <- inpt$INSV
        ni <- unlist(lapply(inpt$INSL,length))
        print(sum(ni))
        print(nrow(XZ))
        Vi <- unlist(lapply(inpt$INSV,function(x) prod(diag(x)))) * (pi^(DIM/2)) / gamma(DIM/2+1) #vols
        plot(Vi,type='s',ylim=c(0,max(Vi)));abline(h=0,lty=2)
        gtheti <- rep(0,nrow(XZ))
        for(i in 1:max(IZ))
            gtheti <- gtheti + ni[i]*inelipse(XZ,VZ[[i]],MNZ[[i]]) / Vi[i]
        LLmax <- max(LZ)
        wz <- exp(LZ-LLmax) / gtheti
        Z <- exp(LLmax) * sum(wz)
        wz <- wz / sum(wz)
        P <- XZ
        vZ <- sum((exp(LZ) / gtheti-Z / sum(ni)^2)^2)
    }
    list(Z=Z,w=wz,X=P,LLmax=LLmax,vZ=vZ)
}


## prod(eigen(SS)$values)

## ## prod(diag(V)^2)
## DIM <- 2

## VZ <- testch$INSV

## Vi <- unlist(lapply(testch$INSV,function(x) prod(diag(x)^2))) * (pi^(DIM/2)) / gamma(DIM/2+1)

## plot(Vi,type='s')

## names(testch)
## testch$neval                            #15115
## ## MNZ <- do.call('rbind',testch$INSM)     #ellipse means, i, 1995
## MNZ <- testch$INSM     #ellipse means, i, 1995
## length(testch$INSV)     #ellipse V, i, 1995
## length(testch$INSX)                     #1995
## length(testch$INSI)                     #15115
## IZ <- do.call('rbind',testch$INSI)
## testch$neval-testch$Nlive == length(testch$INSI) #?

## XZ <-  matrix(unlist(testch$INSX),nrow=length(IZ),ncol=2,byrow=TRUE) #todo:change DIM
## LZ <- unlist(testch$INSL)
## ni <- unlist(lapply(testch$INSL,length))

## gtheta <- function(theta){
##     ans <- 0
##     for(i in 1:max(IZ)){
##         bit <- ni[i] * inelipse(theta,VZ[[i]],MNZ[[i]]) / Vi[i]
##         if(bit>0)cat('i=',i,';',bit,'\n')
##         ans <- ans + bit
##     }
##     ans/max(IZ)
## }


## gtheta <- function(theta){
##     ans <- unlist(mapply(function(X,Y) inelipse(theta,X,Y), VZ, MNZ))
##     ans <- ans * ni / Vi
##     sum(ans)/max(IZ)
## }


## ni[1] * inelipse(XZ[1,],VZ[[1]],MNZ[[1]]) / Vi[1]

## gtheta(XZ[4,])

## GZ <- apply(XZ[1:1e3,],1,gtheta)                #oof! slow...
## system.time({GZ <- apply(XZ[1:1e1,],1,gtheta)}) #0.5 s for 10
## ## looking like 10 mins for full 14e3... :-(


## gtheti <- ni[1]*inelipse(XZ,VZ[[1]],MNZ[[1]]) / Vi[1]

## ## good timing example
## system.time({                           #500x faster
##     gtheti <- rep(0,nrow(XZ))
##     for(i in 1:max(IZ)){
##         bit <- ni[i] * inelipse(XZ,VZ[[i]],MNZ[[i]]) / Vi[i]
##         gtheti <- gtheti + bit
##         ## if(i==620)cat('i=',i,';',bit[4],'\n')
##         ## if(i==620){
##         ##     print(XZ[4,])
##         ##     print(inelipse(XZ[4,],VZ[[i]],MNZ[[i]]))
##         ##     print(ni[i])
##         ##     print(Vi[i])
##         ## }
##     }
## })
## gtheti <- gtheti / max(IZ)

## 42% are 0 -- going to be a problem when dividing!! TODO
## GZ very different but also have zeros...

## summary(gtheti)
## summary(GZ)
## head(GZ)
## head(gtheti)


## difference is late i for 4th coord post i=618

## kk <- 620
## ni[kk] * inelipse(XZ[4,],VZ[[kk]],MNZ[[kk]]) / Vi[kk]

## inelipse(XZ[4,],VZ[[kk]],MNZ[[kk]])

## elz <- inelipse(XZ,VZ[[kk]],MNZ[[kk]]) / Vi[kk]
## elz[4]

## ## save(testch,file='../testch.Rdata')

## ## 1% different!?!
## ## BUG FOUND

## inelipse <- function(x,U,m){
##     if(!is.null(dim(x))){
##         m <- matrix(m,nrow=nrow(x),ncol=ncol(x),byrow=TRUE)
##         x <- x-m                            #change for matrix
##         print(x[4,])
##         x <- x %*% solve(U)
##         print(x[4,])
##         print(rowSums(x^2)[4])
##         ans <- rowSums(x^2)<1
##     } else {
##         x <- x-m                            #change for matrix
##         print(x)
##         x <- x %*% solve(U)
##         print(x)
##         print(sum(x^2))
##         ans <- sum(x^2)<1    
##     }
##     ans    
## }


## summary(gtheti)
## hist(gtheti)

## gtheti[4]
## gtheta(XZ[4,])

## ni[1] * inelipse(XZ[4,],VZ[[1]],MNZ[[1]]) / (Vi[1]*max(IZ)) +
##     ni[2] * inelipse(XZ[4,],VZ[[2]],MNZ[[2]]) / (Vi[2]*max(IZ)) 


## plot(rsph(1e5,2),cex=.1)

## SS <- matrix(c(1,-.4,-.4,2),ncol=2)
## test <- mvtnorm::rmvnorm(1e3,mean=c(1,1),sigma=SS)
## plot(test)
## ## SS <- diag(2)
## V <- chol(SS)
## test2 <- rep(1,2) + rnorm(2) %*% V
## test2 <- rep(1,2) + matrix(rnorm(2*1e3),ncol=2) %*% V
## test2 <- rep(1,2) + rsph(1e3,2) %*% V
## points(test2,col=2,pch=4)

## ## area
## RJ <- matrix(runif(2e5),ncol=2,nrow=1e5)
## RJ <- RJ-.35
## RJ <- RJ*6

## plot(RJ,col=4,cex=.1)
## points(test2,col=2,pch=4)
## points(RJ[ie,],col=3,pch=4)

## ie <- inelipse(RJ,V,c(1,1))
## mie <- mean(ie)
## 36*mie
## prod(diag(V)) * (pi^(DIM/2)) / gamma(DIM/2+1)


## out <- inelipse(test,m=c(1,1),U=V)
## mean(out)
## plot(test)
## points(test[out,],col=2)

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

## ## integrate(f=LLtestv,lower=1e-6,upper=.999)
## integrate(f=LLtestve,lower=1e-6,upper=.999)
## integrate(f=function(x)exp(-x^2/2),lower=-Inf,upper=Inf) #sqrt(2*pi)=2.58

## ## LLtest1 <- function(x){
## ##     ans <- 2/3
## ##     if(abs(x-.5)<.25)
## ##         ans <- 2*ans
## ##     ans
## ## }

## ## LLtest <- function(x){
## ##     sum(sapply(x,function(x)log(LLtest1(x))))
## ## }

## ## LLtest(c(.5,.5))
         

## ## test <- runif(1e4)
## ## for(i in 1:1e4)test[i] <- LLtest1(test[i])
## ## sum(test)/1e4

## testch <- enest(Nlive=5e2,niter=3e3,DIM=2, expand=2,
##                      toteval = 5e3,rejtol=1e-3,
##                      LL=LLtest)


## ## str(testch)
## plot(-na.omit(-testch$LLall),type='s')#,log='y')

## plot(testch$acc,log='y')


## tcho <- postprocess(testch,INS=1)
## str(tcho)
## tcho$Z/(2*pi)^(2/2)                     #WRONG!!
## pks <- sample(nrow(tcho$X),1e3,replace=TRUE,prob=tcho$w/sum(tcho$w))

## ## pairs((tcho$X[pks,]),cex=.1)
## plot(logit(tcho$X[pks,1]),logit(tcho$X[pks,2]))
## ## var(logit(tcho$X[pks,1]))

## pairs(logit(tcho$X[pks,]),cex=2e2*tcho$w)

## ## tmp <- logit(tcho$X[pks,])
## ## colMeans(tmp)
## ## pairs(tmp[,1:10],cex=.1)
## ## sum(tcho$w)^2/sum(tcho$w^2)
