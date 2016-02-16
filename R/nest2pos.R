nest2pos <- function(chains, Nlives, Lcol=-1){
    ## pjd: takes list of chains
    ## % Convert a set of nested sample chains (each containing the same
    ## % number of parameters), computed with given numbers of live points
    ## % into a posterior array.
    ## %
    ## % The inputs are:
    ## %	chains: A cell array (or if just one set of nested samples is being
    ## %       used this can be just an array) containing a set of nested sample
    ## %		parameter chains, including a column with the log likelihood
    ## %       values. The nested samples from each chains will be combined when
    ## %       forming the posteriors. Each nested sample chain must have the same
    ## %       number of parameters.
    ## %	Nlives: An array, or single value, containing the number of live
    ## %       points, used for each parameter chain. If this is a single value
    ## %       then it will be assumed that all chains used the same number of
    ## %		live points.
    ## %	Lcol:   The column within the chains that contains the log likelihood
    ## %       values. If not set then this defaults to using the final column.
    ## %
    ## % This is based on the draw_posterior_many python function in
    ## % nest2pos.py at
    ## % https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/tree/lalapps/src/inspiral/posterior/nest2pos.py

    ## % check for Lcol value
    if(Lcol==-1){
        logLcol <- ncol(chains[[1]]) ## % use final column        
    } else {
        logLcol <- Lcol        
    }

    nchains <- length(chains) ## % get number of chains
    logevs <- matrix(0,nrow=nchains, ncol=1) ## % log evidence values for each chain
    ## logwts <- matrix(nchains, 1) ## % log weights for each chain
    logwts <- list() ## % log weights for each chain, pjd: list

    ## % check whether only one live point value has been given
    if( length(Nlives) == 1 ) {
        Nlives = Nlives*matrix(1,nrow=nchains, ncol=1)        
    }

    
    ## % compute weights and evidence for each chain
    for(i in 1:nchains){
        ## % get log likelihood values
        if( logLcol > 0 ) {
            loglikes <- chains[[i]][,logLcol]
        } else{
            loglikes <- chains[[i]][,ncol(chains[[i]])] #never runs
        }
    
        Nlive <- Nlives[i]

        end <- length(loglikes)         #pjd
        start_data <- loglikes[1:(end-Nlive)] #todo: these seem to overlap, check matlab syntax and logic
        end_data <- loglikes[(end-Nlive+1):end] #chng
        
        log_wts <- matrix(0,nrow=length(loglikes), ncol=1)
        log_vol_factor <- log1p(-1/Nlive)
        log_dvol <- -1/Nlive
        log_vol <- 0
        log_ev <- -Inf
        for(j in 1:length(start_data)){
            log_this_vol <- log_vol + log_dvol
            log_wts[j] <- start_data[j] + log_this_vol
            log_ev <- logplus(log_ev, log_wts[j])
            log_vol <- log_vol + log_vol_factor    
        } 

        avg_log_like_end <- -Inf
        for(j in 1:length(end_data)){
            avg_log_like_end <- logplus(avg_log_like_end, end_data[j])      
        } 

        avg_log_like_end <- avg_log_like_end - log(Nlive)
        end <- length(log_wts)          #pjd
        log_wts[(end-Nlive+1):end] <- end_data + log_vol #chng
        log_ev <- logplus(log_ev, avg_log_like_end + log_vol)
        log_wts <- log_wts - log_ev

        logevs[i] <- log_ev
        logwts[[i]] <- log_wts
    }

    
    ## % get total evidence
    log_total_evidence <- logevs[1]
    for( i in 2:nchains){
        log_total_evidence <- logplus(log_total_evidence, logevs[i])
    }

    log_max_evidence <- max(logevs)

    ## % get relative weights of each set of nested samples
    Ns <- matrix(0,ncol=nchains, nrow=1)
    for(i in 1:nchains){
        frac <- exp(logevs[i]-log_max_evidence)
        Ns[i] <- frac/length(chains[[i]])    
    } 
    Ntot <- max(Ns)
    fracs <- Ns/Ntot

    
    ## % get posterior samples from each chain
    posts <- list()                     #pjd orig: cell(nchains, 1);
    for(i in 1:nchains) {
        Ncol <- ncol(chains[[i]]) ## % number of columns
        data <- chains[[i]]
        logwt <- logwts[[i]]            #pjd: list
        maxWt <- max(logwt)
        normalised_wts <- logwt - maxWt
        rands <- log(runif(length(logwt)))
        idxs <- normalised_wts > rands
        pos <- data[idxs,]
        wts <- normalised_wts[idxs]

        ## % append normalised weights to posterior
        ## pos[,Ncol+1] <- wts
        pos <- cbind(pos,wts)           #pjd

        ## % now weight samples according to the relative weights of each set of nested samples
        ## if(nchains > 1){
        nrands <- runif(nrow(pos))      #chng
        posts[[i]] <- pos[nrands < fracs[i],]
        ## } else {
        ## posts{i} = pos;    
        ## }
        
    }
    
    ## % combine posteriors
    post <- do.call('rbind',posts)

    ## % sort in ascending order based on the weights
    ## %post = sortrows(post, Ncol+1);

    return(post)
}



