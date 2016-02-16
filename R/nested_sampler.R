## function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
##           Nlive, tolerance, likelihood, model, prior, extraparams, ...
##           varargin)

## % function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
## %           Nlive, Nmcmc, tolerance, likelihood, model, prior, extraparams)
## %
## % This function performs nested sampling of the likelihood function from
## % the given prior (given a set of data, a model, and a set of extra model
## % parameters).
## %
## % By default the algorithm will draw new samples from a set of bounding
## % ellipsoids constructed using the MultiNest algorithm for partitioning
## % live points. However, if the optional 'Nmcmc' argument is set and
## % Nmcmc > 0, new samples will be drawn from a proposal using an MCMC. This
## % method is based on that of Veitch & Vecchio. For both methods the
## % sampling will stop once the tolerance critereon has been reached.
## %
## % The likelihood should be the function handle of a likelihood function to
## % use. This should return the log likelihood of the model parameters given
## % the data.
## %
## % The model should be the function handle of the model function to be
## % passed to the likelihood function.
## %
## % The prior should be a cell array with each cell containing five values:
## %   parameter name (string)
## %   prior type (string) e.g. 'uniform', 'gaussian' of 'jeffreys'
## %   minimum value (for uniform prior), or mean value (for Gaussian prior)
## %   maximum value (for uniform prior), or width (for Gaussian prior)
## %   parameter behaviour (string):
## %       'reflect' - if the parameters reflect off the boundaries
## %       'cyclic'  - if the parameter space is cyclic
## %       'fixed'   - if the parameters have fixed boundaries
## %       ''        - for gaussian priors
## %   e.g., prior = {'h0', 'uniform', 0, 1, 'reflect';
## %                  'r', 'gaussian', 0, 5, '';
## %                  'phi', 'uniform', 0, 2*pi, 'cyclic'};
## %
## % extraparams is a cell array of fixed extra parameters (in addition
## % to those specified by prior) used by the model
## % e.g.  extraparams = {'phi', 2;
## %                      'x', 4};
## %
## % Optional arguments:
## %  Set these via e.g. 'Nmcmc', 100
## %   Nmcmc - if this is set then MultiNest will not be used as the sampling
## %           algorithm. Instead an MCMC chain with this number of iterations
## %           will be used to draw the number nested sample point.
## %   Nsloppy - if this is set then during the MCMC the likelihood will only
## %             be evaluted once every Nsloppy points rather than at every
## %             iteration of the chain.
## %   covfrac - the relative fraction of the iterations for which the MCMC
## %             proposal distribution will be based on a Students-t
## %             distribution defined by the covariance of the current live
## %             points.
## %   diffevfrac - the relative fraction of the iterations that will use
## %                differential evolution to draw the new sample.
## %   stretchfrac - the relative fraction of the iterations that will use the
## %                 affine invariant ensemble stretch method for drawing a
## %                 new sample
## %   walkfrac - the relative fraction of the iterations that will use the
## %              affine invariant ensemble walk method for drawing a new
## %              sample
## %   propscale - the scaling factor for the covariance matrix used by the
## %               'covfrac' Students-t distribution proposal. This defaults
## %               to 0.1.
## %
## % E.g. if covfrac = 10 then diffevfrac = 5 the Students-t proposal will be
## % used 2/3s of the time and differential evolution 1/3. The default is to
## % use the affine invariant samplers with the stretch move 75% of the time
## % and the walk move 25% of the time.
## %
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nested_sampler <- function(data, Nlive, tolerance, likelihood, model, prior,h=1.1,
                           extraparams=data.frame(), varargin=NULL, verbose=TRUE,
                           maxit = Inf,
                           Nmcmc = 0, Nsloppy = 0, covfrac = 0, diffevfrac = 0,
                           walkfrac = 25, stretchfrac = 75,
                           propscale = 0.1 #optional MCMC parms
                           ){
    ## % get optional input arguments
    if( ! Nmcmc > 0 ){
        cat('Using MultiNest algorithm\n')
    } 

    
    ## % get the number of parameters from the prior array
    D <- nrow(prior)

    ## % get all parameter names
    parnames <- prior[,1]

    if(nrow(extraparams)>0){
        extraparnames <- extraparams[,1]
        extraparvals <- extraparams[,2]
        parnames <- c(parnames, extraparnames)
    } else {
        extraparvals <- c()
    }

    ## % draw the set of initial live points from the prior
    livepoints <- matrix(0,nrow=Nlive,ncol=D)

    ## setting up prior
    for(i in 1:D){
        priortype <- prior[i,2]
        p3 <- prior[i,3]
        p4 <- prior[i,4]

        ## % currently only handles uniform or Gaussian priors
        if (priortype == 'uniform')
            livepoints[,i] <- p3 + (p4-p3)*runif(Nlive)
        if (priortype == 'gaussian')
            livepoints[,i] <- p3 + p4*rnorm(Nlive)
        if (priortype == 'jeffreys') ## % uniform in log space
            livepoints[,i] <- 10^(log10(p3) + (log10(p4)-log10(p3))*runif(Nlive))
    }

    ## % calculate the log likelihood of all the live points
    logL <- rep(0,Nlive)

    for(i in 1:Nlive){                  #todo: consider vectorizing
        parvals <- c(livepoints[i,],extraparvals)
        logL[i] <- likelihood(data, model, parnames, parvals) #
    }

    if(any(is.na(logL))) stop('Likelihood contains NAs! Bailing...')
    if(any(is.nan(logL))) stop('Likelihood contains NaNs! Bailing...')

    ## % now scale the parameters, so that uniform parameters range from 0->1,
    ## % and Gaussian parameters have a mean of zero and unit standard deviation
    livepoints <- scale_parameters(prior,livepoints) #pjd vectorized
    
    ## % initial tolerance
    tol <- Inf

    ## % initial width of prior volume (from X_0=1 to X_1=exp(-1/N))
    logw <- log(1 - exp(-1/Nlive))

    ## % initial log evidence (Z=0)
    logZ <- -Inf

    ## % initial information
    H <- 0

    ## % initialize array of samples for posterior
    ## nest_samples <- matrix(0,nrow=1,ncol=D+1)
    nest_samples <- list()              #pjd a list

    ## %%%%%%%%%%%%%%%%
    ## % some initial values if MultiNest sampling is used
    ## h <- 1.1 ## % h values from bottom of p. 1605 of Feroz and Hobson
    FS <- h ## % start FS at h, so ellipsoidal partitioning is done first time
    K <- 1 ## % start with one cluster of live points

    ## % get maximum likelihood
    logLmax <- max(logL)

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## % initialize iteration counter
    j <- 1
    ## %figure;

    ## % MAIN LOOP
    while( (tol > tolerance || j <= Nlive) & j < maxit ){ #pjd: maximum iterations
        
        ## % expected value of true remaining prior volume X
        VS <- exp(-j/Nlive)

        ## % find minimum of likelihoods
        logLmin <- min(logL)
        idx <- which.min(logL)

        ## % set the sample to the minimum value
        ## nest_samples[j,] <- c(livepoints[idx,],logLmin)
        nest_samples[[j]] <- c(livepoints[idx,],logLmin) #pjd: a list

        ## % get the log weight (Wt = L*w)
        logWt <- logLmin + logw

        ## % save old evidence and information
        logZold <- logZ
        Hold <- H

        ## % update evidence, information, and width
        logZ <- logplus(logZ, logWt)    
        H <- exp(logWt - logZ)*logLmin + exp(logZold - logZ)*(Hold + logZold) - logZ
        
        ## %logw = logw - logt(Nlive)
        logw <- logw - 1/Nlive

        ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (Nmcmc > 0){
            ## % do MCMC nested sampling

            ## % get the Cholesky decomposed covariance of the live points
            ## % (do every 100th iteration - CAN CHANGE THIS IF REQUIRED)
            if( !(j-1)%%100 ){
                ## % NOTE that for numbers of parameters >~10 covariances are often
                ## % not positive definite and cholcov will have "problems".
                ## %cholmat = cholcov(propscale*cov(livepoints));

                ## % use modified Cholesky decomposition, which works even for
                ## % matrices that are not quite positive definite
                ## % from http://infohost.nmt.edu/~borchers/ldlt.html
                ## % (via http://stats.stackexchange.com/questions/6364
                ## % /making-square-root-of-covariance-matrix-positive-definite-matlab
                cv <- cov(livepoints)
                cholmat <- chol(propscale*cv) #todo: needs checking
            
            } 
            
            ## % draw a new sample using mcmc algorithm
            tmp <-  draw_mcmc(livepoints, cholmat,
                              logLmin, prior, data, likelihood, model, Nmcmc, Nsloppy, 
                              covfrac, diffevfrac, walkfrac, stretchfrac, parnames, 
                              extraparvals) #todo: fn
            livepoints[idx,] <- tmp$smpl
            logL[idx] <- tmp$logL
        } else {

            ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ## % do MultiNest nested sampling

            ## % separate out ellipsoids
            if( FS >= h ){
                ## % NOTE: THIS CODE IS GUARANTEED TO RUN THE 1ST TIME THROUGH
                ## % calculate optimal ellipsoids
                tmp <- optimal_ellipsoids(livepoints, VS) #todo: protection against bollixed livepoints
                Bs <- tmp$Bs
                mus <- tmp$mus
                VEs <- tmp$VEs
                ns <- tmp$ns
                K <- length(VEs) ## % number of ellipsoids (subclusters)    
            }

            if( FS < h || !length(VEs) ){ 
                if( !length(VEs) ) {
                    VEs <- VEstmp ## % revert to previous VEs value    
                }
                ## % simply rescale the bounding ellipsoids
                for( k in 1:K){
                    scalefac <- max(1, (exp(-(j+1)/Nlive)*ns[k]/Nlive)/VEs[k]) #todo: check the max statement
                    ## % scale bounding matrix and volume
                    if( scalefac != 1 ){
                        Bs[((k-1)*D+1):(k*D),] <- Bs[((k-1)*D+1):(k*D),]*scalefac^(2/D)
                        VEs[k] <- scalefac*VEs[k]
                    } 
                }
            } 

            VEstmp <- VEs

            if( DEBUG && D==2 ){
                ## % plot 2-dimensionsal live points and bounding ellipses
                plot_2d_livepoints_with_ellipses(livepoints, Bs, mus, j)     
            } 

            ## % calculate ratio of volumes (FS>=1) and cumulative fractional volume
            Vtot <- sum(VEs)
            FS <- Vtot/VS
            fracvol <- cumsum(VEs)/Vtot
            ## % draw a new sample using multinest algorithm
            tmp <- draw_multinest(fracvol,
                                  Bs, mus, logLmin, prior, data, likelihood, model,
                                  parnames, extraparvals) #
            livepoints[idx,] <- tmp$smpl
            logL[idx] <- tmp$logL


        }
    
        ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ## % update maximum likelihood if appropriate
        if( logL[idx] > logLmax )
            logLmax <- logL[idx]
    

        ## % work out tolerance for stopping criterion
        tol <- logplus(logZ, logLmax - (j/Nlive)) - logZ 

        ## % display progress (optional)
        if( verbose || tol < tolerance ) 
            cat(sprintf('log(Z): %.5e, tol = %.5e, K = %d, iter = %d\n',logZ, tol, K, j))

        ## % update counter
        j <- j+1
        
    }                                   #end main loop
    
    ## % sort the remaining points (in order of likelihood) and add them on to
    ## % the evidence
    isort <- order(logL)
    logL_sorted <- logL[isort]
    livepoints_sorted <- livepoints[isort,]

    for( i in 1:Nlive){
        logZ <- logplus(logZ, logL_sorted[i] + logw) 
    }

    ## % append the additional livepoints to the nested samples
    nest_samples <- do.call('rbind',nest_samples) #pjd: convert from list
    nest_samples <- rbind(nest_samples, cbind(livepoints_sorted, logL_sorted))

    ## % rescale the samples back to their true ranges
    end <- ncol(nest_samples)
    nest_samples[,1:(end-1)] <- rescale_parameters(prior,nest_samples[,1:(end-1)]) #pjd vectorized
    
    ## % convert nested samples into posterior samples - nest2pos assumes that the
    ## % final column in the sample chain is the log likelihood
    post_samples <- list()## nest2pos(list(nest_samples), Nlive) #pjd: now takes list

    
    return(list(logZ=logZ, nest_samples=nest_samples, post_samples=post_samples))## [logZ, nest_samples, post_samples] 
}






