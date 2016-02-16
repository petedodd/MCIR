## function logz  = logplus(logx, logy)
## %
## % logz = logplus(logx, logy)
## %
## % Given logx and logy, this function returns logz=log(x+y).
## % It avoids problems of dynamic range when the
## % exponentiated values of x or y are very large or small.
## %
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logplus <- function(logx, logy){
    if( is.na(logx) | is.na(logy) )
        return(-Inf)        
        
    if( abs(logx)==Inf & abs(logy)==Inf )
        return(-Inf)

    if( logx > logy ){
        logz <- logx+log(1+exp(logy-logx))        
    } else {
        logz <- logy+log(1+exp(logx-logy))        
    }
    return(logz)
}
