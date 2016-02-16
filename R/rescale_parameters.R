rescale_parameters <- function(prior, params){
    ## function scaled = rescale_parameters(prior, params)
    ## % scaled = rescale_parameters(prior, params)
    ## %
    ## % This function will do the reverse of scale_parameters.

    lp <- ncol(params)                  #pjd: use df for params
    scaled <- matrix(NA,nrow=nrow(params),ncol=ncol(params)) #pjd - vectorising

    for(i in 1:lp){
        priortype <- prior[i,2]
        p3 <- prior[i,3]
        p4 <- prior[i,4]

        ## % currently only handles uniform or Gaussian priors
        if( priortype == 'uniform' )
            scaled[,i] <- params[,i]*(p4 - p3) + p3
        if( priortype == 'gaussian' )
            scaled[,i] <- params[,i]*p4 + p3
        if( priortype == 'jeffreys' )
            scaled[,i] <- 10^(params[,i]*(log10(p4) - log10(p3)) + log10(p3))

    }
    return(scaled)
}
