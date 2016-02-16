optimal_ellipsoids <- function(u, VS){
    ## % function [Bs, mus, VEs, ns] = optimal_ellipsoids(u, VS)
    ## %
    ## % This function attempts to optimally partition the multi-dimensional
    ## % samples u (uniformly distributed within the sample volume VS), into
    ## % a set of subclusters enclosed by bounding ellipsoids.  The algorithm
    ## % is based on Algorithm 1 of the MULTINEST paper by Feroz, Hobson,
    ## % and Bridges, MNRAS, 398, 1601-1614 (2009).
    ## %
    ## % Output:
    ## %   Bs:  an array of bounding matrices for the ellipsoids enclosing
    ## %        the subclusters, scaled to have at least the minimum volume
    ## %        required by the subclusters. ( (K x ndims) x ndims )
    ## %   mus: an array of centroids for the bounding ellipsoids (K x ndims)
    ## %   VEs: an array of volumes for the bounding ellipsoids   (K x 1)
    ## %   ns:  an array containing the number of points for each subcluster (K x 1)
    ## %
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N <- nrow(u) ## % number of samples in multi-dimensional space
    ndims <- ncol(u) ## % number of dimensions

    ## % calculate bounding matrix, etc. for bounding ellipsoid associated
    ## % with the original set of points u
    tmp <- calc_ellipsoid(u, VS)        #
    B <- tmp$B
    mu <- tmp$mu
    VE <- tmp$VE
    tmp <- tmp$flag

    ## % attempt to split u into two subclusters
    tmp <- split_ellipsoid(u, VS)       #
    u1 <- tmp$u1
    u2 <- tmp$u2
    VE1 <- tmp$VE1
    VE2 <- tmp$VE2
    nosplit <- tmp$nosplit
    n1 <- nrow(u1)
    n2 <- nrow(u2)

    if( nosplit || n1<ndims+1 || n2<ndims+1 ){
        ## % couldn't split the cluster
        Bs <- B
        mus <- matrix(mu,nrow=1,ncol=length(mu)) #pjd slight mod here
        VEs <- VE
        ns <- N    
    } else {
        ## % check if we should keep the partitioning of S
        if (VE1 + VE2 < VE || VE > 2*VS) {
            if (DEBUG){
                cat(sprintf('PARTITION ACCEPTED: N=%d splits to n1=%d, n2=%d\n', N, n1, n2))
            } 
 
            VS1 <- n1 * VS / N
            VS2 <- n2 * VS / N

            tmp <- optimal_ellipsoids(u1, VS1) #todo: fn
            B1 <- tmp$Bs
            mu1 <- tmp$mus
            VE1 <- tmp$VEs
            n1 <- tmp$ns

            tmp <- optimal_ellipsoids(u2, VS2) #todo: fn
            B2 <- tmp$Bs
            mu2 <- tmp$mus
            VE2 <- tmp$VEs
            n2 <- tmp$ns

            Bs <- rbind(B1 , B2)
            mus <- rbind(mu1 , mu2)
            VEs <- rbind(VE1 , VE2)
            ns <- rbind(n1 , n2)
        } else {
            if(DEBUG) 
                cat(sprintf('PARTITION REJECTED: N=%d doesnt split into n1=%d and n2=%d\n', N, n1, n2))
            Bs <- B
            mus <- matrix(mu,nrow=1,ncol=length(mu)) #pjd slight mod
            VEs <- VE
            ns <- N   
        }
    }

    return(list(Bs=Bs, mus=mus, VEs=VEs, ns=ns))## function [Bs, mus, VEs, ns] = 
}

               
## testo3 <- optimal_ellipsoids(uz,1)
