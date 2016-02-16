## an affine-invariant MCMC sampler, like emcee etc

## to sample K/sqrt(z) in [1/a,a], 1/K=2*(sqrt(a)-1/sqrt(a))
##  F(z) = 2*K*sqrt(z)
random_g <- function(n,a){
  y <- runif(n)
  ans <- (  y*sqrt(a) + (1-y)/sqrt(a) )^2
  ans[ans<1/a] <- 0
  ans[ans>a] <- 0
  return(ans)
}


## testm <- matrix(runif(10),ncol=2)
## lv <- getlist(testm)
## lv
## unlist(lapply( getlist(testm), logprob))
## mclapply(getlist(testm), sum, mc.cores=2)

## rep(1:nrow(testm), each = ncol(testm)) ##NB
## rep(1:nrow(testm),ncol(testm))

## splitting matrix into list of rows for parallelisation
getlist <- function(x)split((x), rep(1:nrow(x), ncol(x)))




## z <- random_g(1e4,2)
## hist(z,freq=FALSE)
## xz <- seq(from=.5,to=2,len=100)
## yz <- 1/(2*(sqrt(2)-1/sqrt(2)) * sqrt(xz))
## lines(xz,yz,col=2)


update_step <- function(X1,X2,L1,L2,lnprob,a,CPUs=1){
  ## X1 rows are walkers and cols are parms
  ## L1 is a vector of the lnprob(X1)s
  n <- dim(X1)[1]                       #no walkers
  p <- dim(X1)[2]                       #no dims
  X1old <- X1
  X2old <- X2
  ## update 1st ensemble
  sel <- sample(n,replace=TRUE) #random selection
  z <- random_g(n,a)
  Y <- X2old[sel,] + z * (X1old-X2old[sel,])
  if(CPUs==1)
    LY <- apply(Y,1,lnprob)               #calculate the LLs
  else{
    LY <- unlist(parLapply(NULL,getlist(Y), lnprob))
    ## LY <- unlist(mclapply(getlist(Y), lnprob, mc.cores=CPUs)) #NOT ON WINDOWS!
  }
  lq <- (p-1)*log(z) + LY - L1          #prob ratio
  lq[is.nan(lq)] <- -Inf; lq[is.na(lq)] <- -Inf                 #safety
  lr <- log(runif(n))
  news <- lr<lq
  X1[news,] <- Y[news,]
  L1[news] <- LY[news]
  acc <- sum(news)                      #acceptance
  ## update 2nd ensemble
  sel <- sample(n,replace=TRUE) #random selection
  z <- random_g(n,a)
  Y <- X1old[sel,] + z * (X2old-X1old[sel,])
  if(CPUs==1)
    LY <- apply(Y,1,lnprob)               #calculate the LLs
  else{
    ##clusterExport(NULL,varlist=c('Y'))
    LY <- unlist(parLapply(NULL,getlist(Y), lnprob)) 
    ## LY <- unlist(mclapply(getlist(Y), lnprob, mc.cores=CPUs)) #NOT ON WINDOWS!
  }
  lq <- (p-1)*log(z) + LY - L2          #prob ratio
  lq[is.nan(lq)] <- -Inf; lq[is.na(lq)] <- -Inf                 #safety
  lr <- log(runif(n))
  news <- lr<lq
  X2[news,] <- Y[news,]
  L2[news] <- LY[news]
  acc <- acc + sum(news)                      #acceptance
  return(list(X1=X1,X2=X2,L1=L1,L2=L2,acc=acc))
}


getAIESamples <- function(halfN,nosteps,lnprob,x0,a=2,eps=1e-2,
                       CPUs=1,dynlib=NULL,echo=FALSE,lgf=NULL){
    ## eps is the pterturbation around initial state
    if(!is.null(lgf)){
        if(file.exists(lgf)) file.remove(lgf)
        file.create(lgf)
    }
    ## parallelisation prep
    if(CPUs>1){
        if(CPUs>detectCores()) {
            CPUs <- detectCores()
            print(paste0('CPUs>cores! Using ',CPUs))
        }
        
        print("preparing parallel environment...")
        cl <- makeCluster(CPUs)
        on.exit(stopCluster(cl))
        varlist <- unique(c(ls(), ls(envir=.GlobalEnv),ls(envir=parent.env(environment()))))
        clusterExport(cl, varlist=varlist, envir=environment())
        if(!is.null(dynlib))
            dummy <- parLapply(cl, as.list(1:CPUs), function(x) dyn.load(dynlib))
        ## clusterExport(cl,varlist=varlist)
        clusterSetRNGStream(cl)
        setDefaultCluster(cl)
        ## wdir <- getwd()
        ## clusterExport(cl, varlist=c("packages", "dyn.libs", "wdir"),envir=environment())
    }
  
    ## INITIALIZE
    np <- length(x0)
    ## initial states
    X1 <- matrix( (1+eps*runif(halfN*np))*x0, nrow=halfN, ncol=np, byrow=TRUE )
    X2 <- matrix( (1+eps*runif(halfN*np))*x0, nrow=halfN, ncol=np, byrow=TRUE )
    ## compute LLS
    L1 <- apply(X1,1,lnprob)
    L2 <- apply(X2,1,lnprob)

    ## DO STEPS
    U <- list(X1=X1,X2=X2,L1=L1,L2=L2)
    R1 <- R2 <- list()
    LL1 <- LL2 <- list()
    acc <- 0
    accv <- rep(NA,nosteps)
    if(!is.null(lgf)){
        cat(Sys.time(),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        cat(summary(U$L1),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        cat(summary(U$L2),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        if( any(U$L1==-Inf) | any(U$L2==-Inf) ){ #a problem with the initial set up
            for(i in which(U$L1==-Inf)){
                cat(U$X1[i,],file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            }
            for(i in which(U$L2==-Inf)){
                cat(U$X2[i,],file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            }
        }
    }
    if( any(U$L1==-Inf) | any(U$L2==-Inf) ){ stop('Bailing! Initial walker positions have infinite likelihood!') }
    for(i in 1:nosteps){
        if(echo & !i%%50 ) print(paste0('---',i,'---'))
        U <- update_step(X1=U$X1,X2=U$X2,L1=U$L1,L2=U$L2,lnprob=lnprob,a=a,CPUs=CPUs)
        R1[[i]] <- U$X1                     #records
        R2[[i]] <- U$X2
        LL1[[i]] <- U$L1                     #records
        LL2[[i]] <- U$L2
        acc <- acc + U$acc                  #acceptance
        accv[i] <- U$acc/(2*halfN)          #acceptance by step
        if(!is.null(lgf)){
            cat(paste0('---------',i,'---------'),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(summary(unlist(LL1)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(summary(unlist(LL2)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(paste0('==acc== ',accv[i]),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        }
    }
    acc <- acc/(2*halfN*nosteps)
    return(list(chains=list(R1=R2,R2=R2),acc=acc,accv=accv,a=a,
                llike=list(LL1=LL1,LL2=LL2),lnprob=lnprob))
}

## carry on this for another notsteps
continueSamples <- function(nosteps,Y,CPUs=1,dynlib=NULL,echo=FALSE,lgf=NULL){
    if(!is.null(lgf)){
        if(file.exists(lgf)) file.remove(lgf)
        file.create(lgf)
    }
    
    ## parallelisation prep
    if(CPUs>1){
        if(CPUs>detectCores()) {
            CPUs <- detectCores()
            print(paste0('CPUs>cores! Using ',CPUs))
        }
        print("preparing parallel environment...")
        cl <- makeCluster(CPUs)
        on.exit(stopCluster(cl))
        varlist <- unique(c(ls(), ls(envir=.GlobalEnv),ls(envir=parent.env(environment()))))
        clusterExport(cl, varlist=varlist, envir=environment())
        if(!is.null(dynlib))
            dummy <- parLapply(cl, as.list(1:CPUs), function(x) dyn.load(dynlib))
        ## clusterExport(cl,varlist=varlist)
        clusterSetRNGStream(cl)
        setDefaultCluster(cl)
    }
    n <- length(Y$chains$R1)              #previous nosteps
    w <- length(Y$llike$LL1[[n]])         #halfn
    a <- Y$a
    accv <- Y$accv
    accv <- c(accv,rep(NA,nosteps))
    lnprob <- Y$lnprob
    ## DO STEPS
    U <- list(X1 = Y$chains$R1[[n]],X2= Y$chains$R2[[n]],
              L1 = Y$llike$LL1[[n]], L2 = Y$llike$LL2[[n]])
    R1 <- Y$chains$R1
    R2 <- Y$chains$R2
    LL1 <- Y$llike$LL1
    LL2 <- Y$llike$LL2
    acc <- Y$acc * 2 * n * w
    if(!is.null(lgf)){
        cat(Sys.time(),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        cat(summary(unlist(LL1)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        cat(summary(unlist(LL2)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
    }
    for(i in 1:nosteps){
        if(echo & !i%%50 ) print(paste0('---',i,'---'))
        U <- update_step(X1=U$X1,X2=U$X2,L1=U$L1,L2=U$L2,lnprob=lnprob,a=a,CPUs=CPUs)
        R1[[n+i]] <- U$X1                     #records
        R2[[n+i]] <- U$X2
        LL1[[n+i]] <- U$L1                     #records
        LL2[[n+i]] <- U$L2
        acc <- acc + U$acc                  #acceptance
        accv[n+i] <- U$acc/(2*w)          #acceptance by step
        if(!is.null(lgf)){
            cat(paste0('---------',i,'---------'),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(summary(unlist(LL1)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(summary(unlist(LL2)),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
            cat(paste0('==acc== ',accv[n+i]),file=lgf,append = TRUE);cat('\n',file=lgf,append = TRUE)
        }
    }
    acc <- acc/( 2 * w * (nosteps + n) )
    return(list(chains=list(R1=R2,R2=R2),acc=acc,accv=accv,a=a,
                llike=list(LL1=LL1,LL2=LL2),lnprob=lnprob))  
}

## todo: 
##  - DIC?


processAIESamples <- function(Y,burnin=1,filename=NULL){
    chains <- Y$chains
    n <- length(chains$R1)                #number of steps
    p <- dim(chains$R1[[1]])[2]           #number of parameters
    w <- dim(chains$R1[[1]])[1]           #half n: particles
    print(mean(Y$accv[burnin:n]))         #acceptance for this
    if(!is.null(filename))cat(mean(Y$accv[burnin:n]),file=filename) #write this out to file, optionally
    dat <- matrix(nrow=2*w*(n-burnin),ncol=p)
    dat[1:((n-burnin)*w),] <- do.call(rbind,chains$R1[(burnin+1):n])
    dat[1:((n-burnin)*w) + dim(dat)[1]/2,] <- do.call(rbind,chains$R2[(burnin+1):n])
    print(dim(dat))                  
    return(dat)                           #this will need rejigging for autocorrelation?
}


## processSamples <- function(Y,burnin=1,thresh=-Inf){
##   chains <- Y$chains
##   n <- length(chains$R1)
##   p <- dim(chains$R1[[1]])[2]
##   w <- dim(chains$R1[[1]])[1]
##   dat <- matrix(nrow=2*(n-burnin)*w,ncol=p)
##   dat[1:((n-burnin)*w),] <- do.call(rbind,chains$R1[(burnin+1):n])
##   dat[1:((n-burnin)*w) + dim(dat)[1]/2,] <- do.call(rbind,chains$R2[(burnin+1):n])
##   if(thresh > -Inf){
##     ## thresholding
##     x <- seq(from=0,to=(n-burnin-1)*w,by=w)
##     w1 <- which(Y$llike$LL1[[n]] > thresh) #OK walkers
##     w1 <- matrix(w1,ncol=length(w1),nrow=length(x),byrow=TRUE)
##     w1 <- c(w1+x)
##     w2 <- which(Y$llike$LL2[[n]] > thresh)
##     w2 <- matrix(w2,ncol=length(w2),nrow=length(x),byrow=TRUE)
##     w2 <- c(w2+x) + dim(dat)[1]/2
##     print(dim(dat))                     #THIS IS NOT CORRECT - overlap at end!!
##     print(summary(w1))
##     print(summary(w2))
##     dat <- dat[c(w1,w2),]
##   }
##   return(dat)                           #this will need rejigging for autocorrelation?
## }

plotAIEchains <- function(Y,labels=NULL,file='',ncol=1,main='',lwdadj=1){
  chains <- Y$chains 
  ## max/min
  n <- length(chains$R1)
  p <- dim(chains$R1[[1]])[2]
  w <- dim(chains$R1[[1]])[1]
  if(is.null(labels))labels <- rep('',p)
  print(c(n,p,w))
  if(ncol!=1)                           #number of cols in plot
    p1 <- ceiling( p / ncol )
  else
    p1 <- p
  lt <- .1                 #todo:this and pointsize should be made adaptive
  if(file!='')png(file,width=480*ncol,height=160*p1,res=600,pointsize=3)
  oldpar <- par(mfrow=c(p1,ncol),mar=c(.3,4,.3,2),oma=c(2,1,1,1))
  ## loop over parz
  for(i in 1:p){
    M1 <- max(unlist(lapply(chains$R1,function(X)max(X[,i]))))
    M2 <- max(unlist(lapply(chains$R2,function(X)max(X[,i]))))
    m1 <- min(unlist(lapply(chains$R1,function(X)min(X[,i]))))
    m2 <- min(unlist(lapply(chains$R2,function(X)min(X[,i]))))
    M <- max(M1,M2)
    m <- min(m1,m2)
    plot(c(1,n),c(m,M),type='n',xlab='',ylab=labels[i],axes=FALSE,main=main)
    if(i<p)
      axis(1,labels=FALSE,lwd=lt)
    else
      axis(1,labels=TRUE,lwd=lt)
    axis(2,lwd=lt);box(lwd=lt);
    for(j in 1:w){                      #plot chains by walker
      lines(1:n,unlist(lapply(chains$R1,function(X)X[j,i])),type='s',
            lwd=10*lt*lwdadj/w)
      lines(1:n,unlist(lapply(chains$R2,function(X)X[j,i])),type='s',
            lwd=10*lt*lwdadj/w)
    }
  }
  if(file!='')dev.off()
  par(oldpar)
}


plotAIEposts <- function(Y,labels=NULL,file='',ncol=1,burnin=0,main=''){
  chains <- Y$chains 
  ## max/min
  n <- length(chains$R1)
  p <- dim(chains$R1[[1]])[2]
  w <- dim(chains$R1[[1]])[1]
  X <- processAIESamples(Y,burnin=burnin)
  if(is.null(labels))labels <- rep('',p)
  print(c(n,p,w))
  if(ncol!=1)                           #number of cols in plot
    p1 <- ceiling( p / ncol )
  else
    p1 <- p
  if(file!='')pdf(file)#,width=480*ncol,height=160*p1)
  oldpar <- par(mfrow=c(p1,ncol),mar=c(2,4,.3,2),oma=c(2,1,1,1))
  ## loop over parz
  for(i in 1:p){
    dd <- density(X[,i]) 
    plot(dd,xlab='',ylab='',main=main)
    if(!is.null(labels)) text(labels[i],x=0.9*min(dd$x)+.1*max(dd$x),
                              y=0.1*min(dd$y)+.9*max(dd$y))
    abline(v=median(X[,i]),col=2)
  }
  if(file!='')dev.off()
  par(oldpar)
}

## plotting the LL chain
plotAIELLchain <- function(Y,file='',main='',topy=NULL){
  chains <- Y$llike 
  ## max/min
  n <- length(chains$LL1)
  w <- length(chains$LL1[[1]])
  if(file!='')png(file,width=480,height=360)
  M1 <- unlist(chains$LL1)
  M2 <- unlist(chains$LL2)
  M1[is.infinite(abs(M1))] <- NA; M2[is.infinite(abs(M2))] <- NA
  M1 <- matrix(M1,ncol=n)
  M2 <- matrix(M2,ncol=n)
  if(is.null(topy)) topy <- min(0,max(c(M1),c(M2),na.rm=TRUE))
  plot(c(1,n),c(min(c(M1),c(M2),na.rm=TRUE),topy),type='n',xlab='',ylab='log likelihood', axes=FALSE,main=main)
  axis(1,labels=TRUE)
  axis(2);box();
  for(j in 1:w){                      #plot chains by walker
    lines(1:n,M1[j,],type='s',lwd=.1)
    lines(1:n,M2[j,],type='s',lwd=.1)
    if(all(diff(M1[j,])==0))print(paste0('1:',j))
    if(all(diff(M2[j,])==0))print(paste0('2:',j))
  }
  if(file!='')dev.off()
}


## logprob <- function(x) -x^2/2
## test <- getSamples(100,200,logprob,x0=1)
## plotchains(test$chains)
## dd <- processSamples(test$chains,burnin=50)
## hist(dd,freq=FALSE);curve(dnorm,from=-5,to=5,col='red',add=TRUE)
## acf(dd[seq(from=3,by=2*100,to=length(dd))],lag.max=100)


## ## harder test
## logprob <- function(x) -(x[1]+x[2])^2/2 -100*(x[1]-x[2])^2/2
## test <- getSamples(100,200,logprob,x0=rep(1,2),CPUs=5)


## plotchains(test$chains)
## dd <- processSamples(test$chains,burnin=50)
## plot(dd[,1],dd[,2],pch='.')
## acf(dd[seq(from=3,by=2*100,to=dim(dd)[1]),1],lag.max=100) #acf of single walker



## testfun <- function(){
##   cat('hi\n')
## }

## testfun <- function(x, CPUs=1,packages=NULL,dyn.libs=NULL){
##   ## setup
##   detectedCores <- detectCores()
##   cat("\n\nCPUs Detected:", detectedCores, "\n")
##   if(CPUs > detectedCores) {
##     cat("\nOnly", detectedCores, "will be used.\n")
##     CPUs <- detectedCores
##   }
##   cat("preparing environments for CPUs...\n")
##   cl <- makeCluster(CPUs)
##   on.exit(stopCluster(cl))
##   varlist <- unique(c(ls(), ls(envir=.GlobalEnv),ls(envir=parent.env(environment()))))
##   clusterExport(cl, varlist=varlist, envir=environment())
##   clusterSetRNGStream(cl)
##   wdir <- getwd()
##   clusterExport(cl, varlist=c("packages", "dyn.libs", "wdir"),envir=environment())
                
##   ## work
##   cat('\nwork!')
##   ans <- unlist(parLapply(cl, list(x),function(x) x^2))
##   return(sum(ans))
## }

## testfun(1:100,CPUs=3)

## a

## ## parallel code....
## ## prep...
## detectedCores <- detectCores()
## cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
##     append=TRUE)
## if(CPUs > detectedCores) {
##   cat("\nOnly", detectedCores, "will be used.\n",
##       file=LogFile, append=TRUE)
##   CPUs <- detectedCores}
## cat("\nLaplace's Demon is preparing environments for CPUs...",
##     file=LogFile, append=TRUE)
## cat("\n##################################################\n",
##     file=LogFile, append=TRUE)
## cl <- makeCluster(CPUs)
## cat("\n##################################################\n",
##     file=LogFile, append=TRUE)
## on.exit(stopCluster(cl))
## varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
##                     ls(envir=parent.env(environment()))))
##           clusterExport(cl, varlist=varlist, envir=environment())
## clusterSetRNGStream(cl)
## wd <- getwd()
## clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
##               envir=environment())

## ## use
## Mo1 <- parLapply(cl, 1:G,
##           function(x) Model(prop[x,], Data))
