## function plot_2d_livepoints_with_ellipses(livepoints, Bs, mus)
## %
## % plot 2-d livepoints in unit square and bounding ellipses
## %
## % NOTE:  mainly used for debug purposes
## %
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## close all

plot_2d_livepoints_with_ellipses <- function(livepoints, Bs, mus, j ){
    ## % extract number of ellipsoids and dimension
    K <- nrow(mus)
    D <- ncol(mus)

    ## % plot live points and bounding ellipse for ndims=2 only
    if(D==2 ){
        pdf(paste0('../dbg/',j,'.pdf')) 
        plot(livepoints[,1],livepoints[,2],pch='+',xlim=c(-.5,1.5),ylim=c(-.5,1.5))
        for(k in 1:K){
            tmp <- ellipsoid_2d(Bs[((k-1)*D+1):(k*D),], mus[k,])
            u <- tmp$x; v <- tmp$y; a <- tmp$a; b <- tmp$b; vol <- tmp$vol
            lines(u,v,col=2)
        }
        
        dev.off()
    }

    
}
