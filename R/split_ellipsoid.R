split_ellipsoid <- function(u, VS){

    ## % function [u1, u2, VE1, VE2, nosplit] = split_ellipsiod(u, VS)
    ## %
    ## % This function takes in a set of multi-dimensional data points u and the
    ## % sample volume (VS) that they occupy. It uses the k-means algorthim to
    ## % split the points into two sub-clusters and uses an optimisation scheme to
    ## % re-assign points, if necessary, between the sub-clusters. This is based
    ## % on the description in Algorithm 1 of the MULTINEST paper by Feroz,
    ## % Hobson, and Bridges, MNRAS, 398, 1601-1614 (2009).
    ## %
    ## % The function returns the points in the two sub-cluster u1 and u2, and
    ## % the volumes of the ellipsoid subclusters VE1 and VE2.  The flag nosplit
    ## % is set to 1 if the splitting cannot be done; otherwise = 0.
    ## %
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    max_attempt <- 50 ## % maximum number of attempts to recluster points

    ## % default return values
    nosplit <- 0
    VE1 <- VE2 <- u1 <- u2 <- c()
    
    ## % extract number of samples and number of dimensions
    N <- nrow(u)
    D <- ncol(u)

    ## % check total number of samples
    if(N < 2*(D+1)){
        if (DEBUG)
            cat(sprintf('CANT SPLIT: total number of samples is too small!  N = %d\n', N))
        nosplit <- 1
        return(list(u1=u1, u2=u2, VE1=VE1, VE2=VE2, nosplit=nosplit)) 
    } 

    ## % use kmeans to separate the data points into two sub-clusters
    tmp <- kmeans(u,2)
    idx <- tmp$cluster
    mu <- tmp$centres
    u1 <- u[idx==1,]
    u2 <- u[idx==2,]
    
    n1 <- nrow(u1) ## % number of samples in S1
    n2 <- nrow(u2) ## % number of samples in S2

    ## % check number of points in subclusters
    if( is.null(n1) || is.null(n2) || n1 < D+1 || n2 < D+1 ){ #pjd nulls for nrow=1
        if (DEBUG)
            cat(sprintf('CANT SPLIT: number of samples in subclusters is too small! n1 = %d, n2 = %d\n', n1, n2))
        nosplit <- 1
        return(list(u1=u1, u2=u2, VE1=VE1, VE2=VE2, nosplit=nosplit))     
    } 


    ## % preallocate temp arrays
    temp_u1 <- list() ## cell(max_attempt,1);
    temp_u2 <- list() ## cell(max_attempt,1);
    temp_VE1 <- temp_VE2<- rep(0,max_attempt)
    FS <- Inf
    FSidx <- 1

    
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%
    numreassigned <- 0
    counter <- 1
    while( TRUE ){

        ## % calculate minimum volume of ellipsoids
        VS1 <- VS*n1/N
        VS2 <- VS*n2/N

        ## % calculate properties of bounding ellipsoids for the two subclusters
        tmp <- calc_ellipsoid(u1, VS1)
        B1 <- tmp$B
        mu1 <- tmp$mu
        VE1 <- tmp$VE
        flag1 <- tmp$flag
        tmp <- calc_ellipsoid(u2, VS2)
        B2 <- tmp$B
        mu2 <- tmp$mu
        VE2 <- tmp$VE
        flag2 <- tmp$flag

        ## % check flags
        if(flag1 || flag2){
            if( DEBUG )
                cat(sprintf('CANT SPLIT!!\n'))
            nosplit <- 1
            return(list(u1=u1, u2=u2, VE1=VE1, VE2=VE2, nosplit=nosplit))     
        } 

        
        ## % construct temporary arrays and cell arrays containing results for
        ## % each pass through the loop
        temp_u1[[counter]] <- u1
        temp_u2[[counter]] <- u2
        temp_VE1[counter] <- VE1
        temp_VE2[counter] <- VE2

        
        if((VE1+VE2)/VS < FS){
            FS <- (VE1+VE2)/VS
            FSidx <- counter
        } 

        ## % DEBUG print statement
        if( DEBUG ){
            cat(sprintf('SPLIT ELLIPSOID: counter = %d, FS = %f, numreassigned = %d\n',counter, (VE1+VE2)/VS, numreassigned))
        }

        ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ## % check if points need to be reassigned to the other subcluster
        reassign <- 0
        m1 <- 0
        m2 <- 0
        u1new <- u2new <- matrix(0,nrow=N,ncol=D)

        
        ## % for all points get the Mahalanobis distance between each point and
        ## % the centroid of each ellipse and assign accordingly
        numreassigned <- 0
        u1m <- u - matrix(mu1,nrow=nrow(u),ncol=ncol(u),byrow=TRUE) #pjd vn
        u2m <- u - matrix(mu2,nrow=nrow(u),ncol=ncol(u),byrow=TRUE) #pjd vn
        u1mb <- u1m %*% solve(B1)                                   #pjd vn
        u2mb <- u2m %*% solve(B2)                                   #pjd vn
        for(i in 1:N){
            ## % get d = (u-mu)^T * B^-1 * (u-mu)
            du1 <- sum( u1m[i,] * u1mb[i,] )
            du2 <- sum( u2m[i,] * u2mb[i,] )
        
            ## % calculate hk = VEk * duk / VSk;
            h1 <- VE1 * du1 / VS1
            h2 <- VE2 * du2 / VS2
            if( h1 < h2 ){
                m1 <- m1 + 1
                u1new[m1,] <- u[i,]

                ## % check if point has been reassigned or not
                if( idx[i]!=1){
                    reassign <- 1
                    idx[i] <- 1
                    numreassigned <- numreassigned + 1
                }                 
            } else {
                m2 <- m2 + 1
                u2new[m2,] <- u[i,]

                ## % check if point has been reassigned or not
                if( idx[i] != 2 ) {
                    reassign <- 1
                    idx[i] <- 2
                    numreassigned <- numreassigned + 1
                }
            }        
        }                               #end for loop

        
        rm(u1, u2, mu1, mu2) 
        n1 <- m1
        n2 <- m2

        u1 <- u1new[1:n1,]
        u2 <- u2new[1:n2,]

        rm(u1new, u2new)

        ## % update counter
        counter <- counter + 1

        if(reassign && counter <= max_attempt){
            
        } else {
            ## % DEBUG print statement
            if(DEBUG){
                cat(sprintf('SPLIT ELLIPSOID: counter = %d, FS = %f, numreassigned = %d\n', counter, (VE1+VE2)/VS, numreassigned))
                if(counter > max_attempt) 
                    cat(sprintf('SPLIT ELLIPSOID: exceeded maximum attempts; take min F(S).\n'))
            } 
            break
        }

        
    }                                   #end while loop

    
    u1 <- temp_u1[[FSidx]]
    u2 <- temp_u2[[FSidx]]
    VE1 <- temp_VE1[FSidx]
    VE2 <- temp_VE2[FSidx]

    if(DEBUG)
        cat(sprintf('SPLIT ELLIPSOID: min F(S) = %f\n', FS))

    return(list(u1=u1, u2=u2, VE1=VE1, VE2=VE2, nosplit=nosplit))  #function [u1, u2, VE1, VE2, nosplit] = 
}


## testo2 <- split_ellipsoid( uz, 1 )
