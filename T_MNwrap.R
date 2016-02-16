
## ## todo: PROBLEM with gaussian prior
## ## todo: sine example
## ## todo: not sure I understand/believe posterior construction
## ## todo: actually think this is correct - just wrong likelihood
## check paper for tolerance and sample weights - not as above! see F code for tol
## profile
## points into mesh
## AFSS

library(lattice)
## setwd('/Users/pjd/Documents/Python/matlabmultinestR/src')

## testing the outputs
DEBUG <- FALSE

## ## get the source files
## srcfiles <- dir(path='../src/',pattern = '\\.R')
## srcfiles <- srcfiles[srcfiles!='testing.R']
## srcfiles <- srcfiles[srcfiles!='wrap.R']

## ## load em
## for(sf in srcfiles)
##     source(paste0('../src/',sf))

## source('~/Documents/Rwork/useful/mysplom3.R')
## source('~/Documents/Rwork/useful/mysplom2.R')



nestwrap <- function(Nlive = 500,tolerance = .5,likelihood,prior, verbose=TRUE){
    ans <- nested_sampler(data=data.frame(),Nlive = Nlive,tolerance = tolerance,likelihood = likelihood, model=function(x,parnames,parvals){},prior=prior,verbose=verbose)
    nc <- ncol(ans$nest_samples)
    wts <- ans$nest_samples[,nc]
    wts <- wts-max(wts)
    wts <- exp(wts)
    wts <- wts/sum(wts)
    list(logZ=ans$logZ,samps=ans$nest_samples[,1:(nc-1)],wts=wts)
}



## ---------- testings ----------


rosenbrock <- function(x){
    f <- (1-x[1])^2 + 100*(x[2] - x[1]^2)^2
    g <- c(100*2*(x[2] - x[1]^2)*(-2)*x[1] - 2*(1-x[1]),
     200*(x[2] - x[1]^2) )
    return(list(logp=-f,grad=-g))
}

arosen <- function(x) return(rosenbrock(x)$logp)

g2d <- function(x) -.5*t(x) %*% matrix(c(1,.2,2,.2),ncol=2) %*% x

pspec22 <- data.frame(name=c('x','y'),type='uniform',bottom=-2,top=4)

nstep <- 1

inc <- function(){
    assign("nstep", nstep+1, envir = .GlobalEnv)
}

like22 <- function(data, model, parnames, parvals){
    inc()
    arosen(parvals)
    ## g2d(parvals)
}

## run the sampler -- VERY slow
system.time({
    test <- nestwrap(likelihood = like22, prior=pspec22)
})


## psmp <- test$samps[,1:2]
## twts <- test$samps[,3]; twts <- twts-max(twts)
## twts <- exp(twts); twts <- twts/sum(twts)

corplot(test$samps,labels=c('x','y'),points=TRUE)
corplot(test$samps,labels=c('x','y'))



## mysplom2(psmp,labels=c('x','y'),weights=twts)

## efficiency about 67%!


## ## ------2D simple
## nullmodel <- function(x,parnames,parvals){
## }


## like2d <- function(data, model, parnames, parvals){
##     dy <- data[,1] - parvals[1]
##     loglike <- -0.5*sum(dy^2/parvals[2]) - parvals[2] - parvals[1]^2/1e4
##     return(exp(loglike))
## }

## gdata <- data.frame(x=rnorm(20))

## pspec2d <- data.frame(name=c('mu','sig'),type=c('gaussian','uniform'),p1=c(0,0.1),p2=c(10,5))

## pspec2d <- data.frame(name=c('mu','sig'),type=c('uniform','uniform'),p1=c(-5,0.1),p2=c(5,5))

## like2d(gdata,nullmodel,parnames=pspec2d[,1],parvals=c(0,1))

## DEBUG <- TRUE

## ## system('rm ../dbg/*pdf')

## testo <- nested_sampler(data=gdata,Nlive = 500,tolerance = .01,likelihood = like2d,
##                         maxit = 1e4,
##                         verbose=TRUE, model=nullmodel,h=1.1,prior=pspec2d)

## ## testm <- testo

## ## look at samples
## clz <- heat.colors(n=nrow(testo$nest_samples)); clz <- rev(clz)
## plot(testo$nest_samples[,1:2],col=clz,pch='+') #by wave

## qplot(x=testo$nest_samples[,1],y=testo$nest_samples[,2],color=testo$nest_samples[,3],geom='point')


## ## posterior?
## ## clz <- heat.colors(n=nrow(testo$post_samples)); clz <- rev(clz)
## ## plot(testo$post_samples[,1:2],col=clz,pch='+') #by wave

## ## psmp <- testo$post_samples
## ## twts <- exp(psmp[,4]); twts <- twts/sum(twts)

## psmp <- testo$nest_samples[,1:2]
## twts <- testo$nest_samples[,3]; twts <- twts-max(twts)
## twts <- exp(twts); twts <- twts/sum(twts)

## plot(density(psmp[,1],weights = twts)); rug(x=psmp[,1]); grid()
## plot(density(psmp[,2],weights = twts)); rug(x=psmp[,2]); grid()


## ## densityplot(psmp[,1],weights = twts)    #doesn't like 0??
## ## densityplot(psmp[,2],weights=twts)      #over by 5

