library(lattice)
## testing the outputs
DEBUG <- TRUE

## get the source files
srcfiles <- dir(path='../src/',pattern = '\\.R')
srcfiles <- srcfiles[srcfiles!='testing.R']

## load em
for(sf in srcfiles)
    source(paste0('../src/',sf))

## source('sinusoid_model.R')              #load the sine model

## make some data

## write a likelihood

## --likelihood(data, model, parnames, parvals) #??

## estimate
## NB a likelihood can take a null model as only used inside likelihood
trueparms <- data.frame(name=c('amp','f','phi'),values=c(1,.1,0.2))

testdata <- data.frame(time=1:10,y=0)

for(i in 1:nrow(testdata))
    testdata[i,2] <- sinusoid_model(testdata[i,1],trueparms$name,trueparms$values) 

eps <- .1
testdata[,2] <- testdata[,2] + eps * rnorm(nrow(testdata))
plot(testdata)


likefun <- function(data, model, parnames, parvals){
    dy <- data[,2] - model(testdata[,1],parnames,parvals)
    ## loglike <- sum(dnorm(dy,sd=eps,log=TRUE))
    loglike <- -0.5*sum(dy^2/eps)
    return(exp(loglike))
}

pspec <- data.frame(name=trueparms$name,type='uniform',bottom=0,top=5)

likefun(testdata,sinusoid_model,trueparms$name,trueparms$values)

DEBUG <- TRUE

testo <- nested_sampler(data=testdata,Nlive = 500,tolerance = .02,likelihood = likefun,model=sinusoid_model,prior=pspec)

psmp <- testo$nest_samples[,1:3]
twts <- testo$nest_samples[,4]; twts <- twts-max(twts)
twts <- exp(twts); twts <- twts/sum(twts)

## densityplot(psmp[,1:3],weights=twts)

plot(density(psmp[,1],weights = twts)); rug(x=psmp[,1]); grid(); abline(v=trueparms$values[1])
plot(density(psmp[,2],weights = twts)); rug(x=psmp[,2]); grid(); abline(v=trueparms$values[2])
plot(density(psmp[,3],weights = twts)); rug(x=psmp[,3]); grid(); abline(v=trueparms$values[3])



## -------- 2D

testo <- nested_sampler(data=testdata,Nlive = 400,tolerance = .01,likelihood = likefun,model=sinusoid_model2D,prior=pspec[1:2,])


psmp <- testo$nest_samples[,1:2]
twts <- testo$nest_samples[,3]; twts <- twts-max(twts)
twts <- exp(twts); twts <- twts/sum(twts)

plot(density(psmp[,1],weights = twts)); rug(x=psmp[,1]); grid(); abline(v=trueparms$values[1])
plot(density(psmp[,2],weights = twts)); rug(x=psmp[,2]); grid(); abline(v=trueparms$values[2])

clz <- heat.colors(n=nrow(testo$nest_samples)); clz <- rev(clz)
plot(testo$nest_samples[,1:2],col=clz,pch='+') #by wave
abline(v=trueparms$values[1])
abline(h=trueparms$values[2])

## ---------2D

source('~/Documents/Rwork/useful/mysplom3.R')

rosenbrock <- function(x){
    f <- (1-x[1])^2 + 100*(x[2] - x[1]^2)^2
    g <- c(100*2*(x[2] - x[1]^2)*(-2)*x[1] - 2*(1-x[1]),
     200*(x[2] - x[1]^2) )
    return(list(logp=-f,grad=-g))
}

arosen <- function(x) return(rosenbrock(x)$logp)

pspec22 <- data.frame(name=trueparms$name,type='uniform',bottom=-2,top=4)

pspec22 <- data.frame(name=c('x','y'),type='uniform',bottom=-2,top=4)

nstep <- 1

inc <- function(){
    assign("nstep", nstep+1, envir = .GlobalEnv)
}

like22 <- function(data, model, parnames, parvals){
    inc()
    arosen(parvals)
}

testo <- nested_sampler(data=testdata,Nlive = 400,tolerance = .01,likelihood = like22, model=nullmodel,prior=pspec22)


psmp <- testo$nest_samples[,1:2]
twts <- testo$nest_samples[,3]; twts <- twts-max(twts)
twts <- exp(twts); twts <- twts/sum(twts)


plot(density(psmp[,1],weights = twts)); rug(x=psmp[,1]); grid(); abline(v=trueparms$values[1])
plot(density(psmp[,2],weights = twts)); rug(x=psmp[,2]); grid(); abline(v=trueparms$values[2])

clz <- heat.colors(n=nrow(testo$nest_samples)); clz <- rev(clz)
plot(testo$nest_samples[,1:2],col=clz,pch='+') #by wave

qplot(x=testo$nest_samples[,1],y=testo$nest_samples[,2],geom='point',alpha=twts)

mysplom2(testo$nest_samples[,1:2],labels=c('x','y'),weights=twts,file='../graphs/Rosen.pdf')


nestwrap <- function(Nlive = 500,tolerance = .01,likelihood,prior, verbose=TRUE){
    ans <- nested_sampler(data=data.frame(),Nlive = Nlive,tolerance = tolerance,likelihood = likelihood, model=function(x,parnames,parvals){},prior=prior,verbose=verbose)
    nc <- ncol(ans$nest_samples)
    wts <- ans$nest_samples[,nc]
    wts <- wts-max(wts)
    wts <- exp(wts)
    wts <- wts/sum(wts)
    list(logZ=ans$logZ,samps=ans$nest_samples[,1:(nc-1)],wts=wts)
}

test <- nestwrap(likelihood = like22,prior=pspec22)

psmp <- test$samps[,1:2]
twts <- test$samps[,3]; twts <- twts-max(twts)
twts <- exp(twts); twts <- twts/sum(twts)

mysplom2(test$samps,labels=c('x','y'),weights=test$wts)

mysplom2(psmp,labels=c('x','y'),weights=twts)

## efficiency about 67%!
## check paper for tolerance and sample weights - not as above! see F code for tol
## profile
## points into mesh
## AFSS


## ------2D simple
nullmodel <- function(x,parnames,parvals){
}


like2d <- function(data, model, parnames, parvals){
    dy <- data[,1] - parvals[1]
    loglike <- -0.5*sum(dy^2/parvals[2]) - parvals[2] - parvals[1]^2/1e4
    return(exp(loglike))
}

gdata <- data.frame(x=rnorm(20))

pspec2d <- data.frame(name=c('mu','sig'),type=c('gaussian','uniform'),p1=c(0,0.1),p2=c(10,5))

pspec2d <- data.frame(name=c('mu','sig'),type=c('uniform','uniform'),p1=c(-5,0.1),p2=c(5,5))

like2d(gdata,nullmodel,parnames=pspec2d[,1],parvals=c(0,1))

DEBUG <- TRUE

## system('rm ../dbg/*pdf')

testo <- nested_sampler(data=gdata,Nlive = 500,tolerance = .01,likelihood = like2d,
                        maxit = 1e4,
                        verbose=TRUE, model=nullmodel,h=1.1,prior=pspec2d)

## todo: PROBLEM with gaussian prior
## todo: sine example
## todo: not sure I understand/believe posterior construction
## todo: actually think this is correct - just wrong likelihood

## testm <- testo

## look at samples
clz <- heat.colors(n=nrow(testo$nest_samples)); clz <- rev(clz)
plot(testo$nest_samples[,1:2],col=clz,pch='+') #by wave

qplot(x=testo$nest_samples[,1],y=testo$nest_samples[,2],color=testo$nest_samples[,3],geom='point')


## posterior?
## clz <- heat.colors(n=nrow(testo$post_samples)); clz <- rev(clz)
## plot(testo$post_samples[,1:2],col=clz,pch='+') #by wave

## psmp <- testo$post_samples
## twts <- exp(psmp[,4]); twts <- twts/sum(twts)

psmp <- testo$nest_samples[,1:2]
twts <- testo$nest_samples[,3]; twts <- twts-max(twts)
twts <- exp(twts); twts <- twts/sum(twts)

plot(density(psmp[,1],weights = twts)); rug(x=psmp[,1]); grid()
plot(density(psmp[,2],weights = twts)); rug(x=psmp[,2]); grid()


## densityplot(psmp[,1],weights = twts)    #doesn't like 0??
## densityplot(psmp[,2],weights=twts)      #over by 5


## Error in nest2pos(list(nest_samples), Nlive) (from nest2pos.R!257hlt#101) : object 'Ns.' not found

## why does tolerance effect when singularity occurs?
