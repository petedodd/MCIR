## library(MCMCpack);

## source('batchmeans.R')# Scripts to compute Effective Sample Size
## source('factorSS.R');
## source('stepOut.R');
## source('heuristicOpt.R');
## source('estimateFactors.R');



getSliceSamples <- function(
    LL, x0,                             #log-like and starting
    nBurnin = 10 * 1024,
    nSample = 10000,
    verbose = TRUE, useFactor = TRUE,
    useTrueFactors = FALSE              #to delete?
    ){

    K <- length(x0)                     #dimension of problem
    sliceParams <- list( nExpands = array(0,c(K,1)),
                        nShrinks       = array(0,c(K,1)),
                        intervalWidths = array(1,c(K,1)),
                        factors        = diag(K),
                        nProposals     = 0,
                        nIterPerAdapt       = 1,
                        nIterPerFactorAdapt = 2048 );

    if ( useTrueFactors ){
        stop('This useTrueFactors is not implemented!')
	sliceParams$factors = eigen( t(X) %*% X )$vectors; #won't work
    }
    if ( !useFactor || useTrueFactors ){
	nBurnin = 1024;
    }


    beta <- x0
    betaMean      <- array( 0, c(K,1));
    betaSampleCov <- array( 0, c(K,K));
    nStored <- 0;


    ##  Initial Tuning Phase
    printScreenThin = ceiling(nBurnin / 10.0)

    if ( verbose  ){
        cat('Starting burnin... \n\n');
    }
    
    for ( I in 1:nBurnin ){
	##  Output Iteration Number
	if ( verbose && (0 == (I %% printScreenThin)) ){
		cat('Burnin Iteration ', I, '\n\n');
	}

	##  Run Iteration of Slice Sampler
	if ( useFactor ){
            output <- factorSliceSampler( beta, LL, sliceParams );
	} else {
            output <- stepOutSS( beta, LL, sliceParams );
	}

## cat('Iter=',I,' \n\n');
        
	beta        <- output$beta;
	sliceParams <- output$sliceParams;

	##  Store Summary Statistics
	betaMean      <- betaMean      + beta;
	betaSampleCov <- betaSampleCov + beta %*% t(beta);
	nStored <- nStored + 1;

	##  Tune Initial Interval Widths
	if ( 0 == (sliceParams$nProposals %% sliceParams$nIterPerAdapt) ){
            sliceParams <- heuristicTuning( sliceParams , K=K);
	}

	##  Tune Sampling Factors
	if ( useFactor && !useTrueFactors &&
	     (0 == (I %% sliceParams$nIterPerFactorAdapt)) &&
	     (I < (nBurnin - sliceParams$nIterPerFactorAdapt)) ){
            sliceParams <- tuneFactors( sliceParams, betaMean, betaSampleCov, nStored, K=K );

            ## Clear Summary Statistics
            betaMean      <- array( 0, c(K,1));
            betaSampleCov <- array( 0, c(K,K));
            nStored <- 0;
	}

	if ( useFactor && !useTrueFactors &&
	     (0 == (I %% sliceParams$nIterPerFactorAdapt)) ){
		sliceParams$nIterPerAdapt <- 1;
		sliceParams$nProposals    <- 0;
	}
    }                                   #end burnin

    
    ##  Sampling Phase
    printScreenThin = nSample / 10.0;
    betaDraws <- array( 0, c(nSample,K));

    for ( I in 1:nSample ){
	##  Output Iteration Number
	if ( verbose && (0 == (I %% printScreenThin)) ){
            cat('Sampling Iteration ', I, '\n\n');
	}

	##  Run Iteration of Slice Sampler
	if ( useFactor ){
            output <- factorSliceSampler( beta, LL, sliceParams );
	}
	else{
            output <- stepOutSS( beta, LL, sliceParams );
	}
	beta        <- output$beta;
	sliceParams <- output$sliceParams;

	##  Store Draw
	betaDraws[I,] <- beta;
    }

    ##  Update Interval Widths to Test Optimality
    if ( verbose ){
	tmp <- heuristicTuning( sliceParams, K=K );
    }

    ## ## from SSDriver.R
    ## ##  Print Sampling Report
    ## betaMoments <- array(0, c(K,2));
    ## betaESS     <- array(0, c(K,1));
    ## betaESsec   <- array(0, c(K,1));
    ## for ( I in 1:K ){
    ##     temp <- bm(betaDraws[,I]);
    ##     betaMoments[I,1] <- temp$est;
    ##     betaMoments[I,2] <- temp$se;
    ##     ## Ignore Negative Autocorrelation (minimal)
    ##     betaESS[I] <- min(nSample,ess(betaDraws[,I]));
    ## }

    return(betaDraws)
    
}


## ##  --------------  tests -------------
## ##  run this on banana
## source('~/Documents/Rwork/useful/mysplom2.R')

## rosenbrock <- function(x){
##     f <- (1-x[1])^2 + 100*(x[2] - x[1]^2)^2
##     g <- c(100*2*(x[2] - x[1]^2)*(-2)*x[1] - 2*(1-x[1]),
##      200*(x[2] - x[1]^2) )
##     return(list(logp=-f,grad=-g))
## }

## ## rosen
## arosen <- function(x) return(rosenbrock(x)$logp)

## runss <- getSliceSamples(LL=arosen,x0=c(0,0),useFactor=TRUE)
## mysplom(runss)

## ## gaussian
## gL <- function(x) -0.5*sum(x^2)

## runss <- getSliceSamples(LL=gL,x0=c(0,0),useFactor=TRUE)
## mysplom(runss)

## ## medium...big...

## N <- 5

## testf2 <- function(x){
##     f <- 0
##     df <- rep(0,2*N)
##     for(i in 1:N) f <- f - 100*(x[1+2*(i-1)]-x[2+2*(i-1)])^2 - (x[1+2*(i-1)]+x[2+2*(i-1)])^2
##     for(i in 1:N) df[c(1+2*(i-1),2+2*(i-1))] <- c(-200*(x[1+2*(i-1)]-x[2+2*(i-1)]) -2*(x[1+2*(i-1)]+x[2+2*(i-1)]), +200*(x[1+2*(i-1)]-x[2+2*(i-1)]) -2*(x[1+2*(i-1)]+x[2+2*(i-1)]))
##     return(list(logp=f,grad=df))
## }

## testf <- function(x) testf2(x)$logp

## runss <- getSliceSamples(LL=testf,x0=rep(0,2*N),useFactor=TRUE)

## colnames(runss) <- letters[1:(2*N)]

## chains <- mcmc(runss[,1:10])

## xyplot(chains,ask=FALSE,layout=c(2,5))

## mysplom(runss[,1:10])


