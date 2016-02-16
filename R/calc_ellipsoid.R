## library(Matrix)

calc_ellipsoid <- function(u, VS){
    ## function [B, mu, VE, flag] =
    ## %
    ## % calculate properties of ellipsoid given a set of points u
    ## %
    ## % Inputs:
    ## %    u:  Nxndims array where N is the number point and ndims is the
    ## %        number of dimensions
    ## %    VS: minimum volume that the bounding ellipsoid should have
    ## %
    ## % Outputs:
    ## %    B:    bounding matrix for ellipsoid including scale factor
    ## %          for mininimum volume
    ## %    mu:   centroid
    ## %    VE:   volume of ellipsoid
    ## %    flag: = 1 if number of points too small or bounding matrix
    ## %          has bad condition number; otherwise = 0
    ## %
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## % default values
    B <- c()
    mu <- c()
    VE <- c()
    flag <- 0

    ## % extract number of points and number of dimensions
    if(is.null(dim(u))){                #pjd: cope with vector input
        N <- ndims <- 1
        goahead <- FALSE
        flag <- 1                       #don't do below
    } else {
        N <- nrow(u)
        ndims <- ncol(u)
        goahead <- TRUE
    }

    ## % check that total number of points is large enough
    if( goahead & N < ndims+1){ #pjd no 
        if(DEBUG)
            cat('number of samples too small to calculate bounding matrix for ellipsoid\n')
        flag <- 1
    } 

    if( goahead ){        #pjd: not vector
        ## % constant factor for volume of ellipsoid
        const <- pi^(ndims/2)/gamma(ndims/2 + 1)

        ## % calculate covariance matrix and centroid
        C <- cov(u)
        mu <- colMeans(u)

        ## % check condition number of C (eps = 2.2204e-16)
        eps <- 1e-10
        if( any(is.na(C)) || rcond(C) < eps || is.nan(rcond(C)) ){ #pjd mod 1st term
            if (DEBUG)
                cat('bad condition number!\n')
            flag <- 1
        } 
    }
    
    ## % find scale factor for bounding ellipsoid E
    ## fB <- 0
    ## print(dim(C));print(dim(u[i,]-mu));print(dim(u[i,]))
    ## for( i in 1:N){
    ##     f <-  (u[i,]-mu)  %*% solve(C)  %*% t(u[i,]-mu) #todo: check
    ##     if( f > fB )
    ##         fB <- f
    ## }

    if( !flag ){
        ## pjd R version
        u <- u - matrix(mu,nrow = N,ncol = ndims,byrow = TRUE)
        fz <- u %*% solve(C)                #100 x 5
        fz <- u * fz
        fB <- max( rowSums(fz) )                   #100
    
        ## % calculate volume of bounding ellipsoid E
        VE <- const*sqrt(det( fB * C ))
        
        ## % expand volume of bounding ellipsoid to VS if necessary
        fV <- 1
        if(VE < VS){
            fV <- (VS/VE)^(2/ndims)
            VE <- VS
        }
    
        ## % scale C to get bounding matrix B
        B <- fV * fB * C
    } else {
        mu <- rep(0,ndims)
        B <- diag(mu)
        VE <- 0
    }
    
    return(list(B=B, mu=mu, VE=VE, flag=flag))## function [B, mu, VE, flag] =
}


## uz <- matrix(runif(1e2*5),nrow=1e2,ncol=5)
## testo <- calc_ellipsoid( uz, 1 )
