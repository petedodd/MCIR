## function [sample, logL] = draw_multinest(fracvol, Bs, mus, ...
##     logLmin, prior, data, likelihood, model, parnames, extraparvals)

## % function [sample, logL] = draw_multinest(fracvol, Bs, mus, ...
## %     logLmin, prior, data, likelihood, model, parnames, extraparvals)
## %
## % This function draws a multi-dimensional sample from the prior volume
## % for use in the nested sampling algorithm. The new point will have a
## % likelihood greater than the value logLmin. The new point will be found by
## % drawing a random multi-dimensional sample from within the set of optimal
## % ellipsoids constructed using the MultiNest algorithm.  The bounding
## % ellipsoids are defined by their bounding matrices Bs and centroids mus.
## % extraparvals is a vector of additional parameters needed by the model.
## %
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


draw_multinest <- function(fracvol, Bs, mus, logLmin, prior, data, likelihood, model, parnames, extraparvals){

    ## % extra number of ellipsoids, number of dimensions
    K <- nrow(mus)                          #
    ndims <- ncol(mus)
    
    while( TRUE ){

        ## % find the ellipsoid from which to draw a new point
        rval <- runif(1)

        for(k in 1:K){
            if(!rval < fracvol[k])          #
                break  
        }

        k0 <- k

        ## % extract bounding matrix and centroid for that ellipsoid
        B <- Bs[((k0-1)*ndims+1):(k0*ndims),]
        mu <- mus[k0,]

        ## % draw points from that ellipsoid until logL >= logLmin
        logL <- -Inf
        while( logL < logLmin ){
            in_range <- 1 ## % default value

            ## % draw one point from the ellipsoid
            pnt <- draw_from_ellipsoid(B, mu, 1) #
            ## if(any(pnt<0)){print(pnt);print(mu)}
            ## pnt[pnt<0] <- -pnt[pnt<0]   #reflect
            ## pnt[pnt>1] <- 2-pnt[pnt>1]  #reflect
            ## pnt[pnt<0] <- 1+pnt[pnt<0]   #wrap
            ## pnt[pnt>1] <- pnt[pnt>1]-1  #wrap
            
            ## % make sure that the point lies in unit hypercube
            if( any(pnt>1) | any(pnt<0) ){
                in_range <- 0
                if(DEBUG) cat('new point not in range!!!!\n')
            }
                    
            if( in_range ){
                ## % assign as candidate replacement live point
                smpl <- pnt                #nb: changed var spelling

                ## % rescale point back to full range
                rescaledpnt <- rescale_parameters(prior, pnt) #

                ## % get new likelihood
                logL <- likelihood(data, model, parnames, c(rescaledpnt,extraparvals))
            } 
        }                                   #end while
    

        ## % check how many ellipsoids this point lies in
        inN <- in_ellipsoids(pnt, Bs, mus)  #

        ## % only accept sample with 1/inN probability
        if( runif(1) < 1/inN) break

    }
    
    return(list(smpl=smpl,logL=logL)) ## function [sample, logL] =    
}
