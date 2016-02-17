##' Function to run automated-factor slice sampling
##'
##' todo: more about provenance and optional arguments...
##' @title getSliceSamples
##' @param LL the loglikelihood function
##' @param x0 a vector specifying the start for the chain
##' @return a matrix with columns as parameters and rows as post-burnin samples
##' @author Pete Dodd
##' @examples
##' rosen <- function(x) -(1-x[1])^2 - 100*(x[2] - x[1]^2)^2
##' runss <- getSliceSamples(LL=rosen,x0=c(0,0),useFactor=TRUE)
##' corplot(runss)
getSliceSamples <- function(
    LL, x0,                             #log-like and starting
    nBurnin = 10 * 1024,
    nSample = 10000,
    verbose = FALSE, useFactor = TRUE,
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

    pb <- txtProgressBar(min=1,max=(nBurnin+nSample),char='.',style=3) #progress bar
    
    for ( I in 1:nBurnin ){
        setTxtProgressBar(pb, I)        #progress
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
        setTxtProgressBar(pb, I+nBurnin)        #progress
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
    close(pb)                           #close progress bar
    return(betaDraws)
    
}
