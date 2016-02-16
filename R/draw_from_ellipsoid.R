draw_from_ellipsoid <- function(B, mu, N ){
    ## % function pnts = draw_from_ellipsoid(B, mu, N )
    ## %
    ## % This function draws points uniformly from an ndims-dimensional ellipsoid
    ## % with edges and orientation defined by the the bounding matrix B and
    ## % centroid mu.  The output is a Nxndims dimensional array pnts.
    ## %
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## % get number of dimensions from the bounding matrix B
    ndims <- nrow(B)

    ## % calculate eigenvalues and vectors of the bounding matrix
    tmp <- eigen(B)
    V <- tmp$vectors
    E <- tmp$values
    D <- sqrt((E))                      #no diag needed, vector

    ## % check size of mu and transpose if necessary
    mnr <- nrow(mu)
    if(!is.null(mnr))
        if( mnr >1 )
            mu <- t(mu)

    ## % generate radii of hyperspheres
    rs <- runif(N)

    ## % generate points
    pt <- matrix(rnorm(N*ndims),nrow=N,ncol=ndims)

    ## % get scalings for each point onto the surface of a unit hypersphere
    fac <- rowSums(pt^2)

    ## % calculate scaling for each point to be within the unit hypersphere
    ## % with radii rs
    fac <- (rs^(1/ndims)) / sqrt(fac)
    pnts <- matrix(0,nrow=N,ncol=ndims)

    ## % scale points to the ellipsoid using the eigenvalues and rotate with
    ## % the eigenvectors and add centroid
    ## % scale points to a uniform distribution within unit hypersphere
    pnts <- fac * pt                #row-wise loop
    ## % scale and rotate to ellipsoid
    ## for( i in 1:N){
    ##     pnts[i,] <- matrix((pnts[i,] * D),nrow=1,ncol=ncol(V)) %*% t(V) + mu
    ## }

    ## M=(ev) %*% diag(mv) %*% solve(ev )
    ## todo: fix this!
    pnts <- sweep(pnts,2,D,`*`)         #multiply cols by D, pjd
    pnts <- pnts %*% t(V)
    pnts <- pnts + matrix(mu,nrow = N, ncol=ndims,byrow=TRUE)
    
    return(Re(pnts))
}

## test <- draw_from_ellipsoid(B=matrix(c(-1,-1,1,1),ncol=2),c(0,0),1e1)

## plot(c(-1,1),c(-1,1),col=2)
## abline(v=1,col=2);abline(v=-1,col=2);
## abline(h=1,col=2);abline(h=-1,col=2);
## points(test)
