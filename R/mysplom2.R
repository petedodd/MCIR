##' A hexbinned pairs and density plot to examine correlations. 
##'
##' The diagonal panels are histograms of each parameter representing estimates of marginal density. The lower half shows hexbinned plots of correlations between parameters (estimates of pairwise joint densities) together with a loess smoother to indicate trend. If \code{points=TRUE} the pair plots are not hexbinned and points are plotted with alpha proportional to weight.
##' @title corplot
##' @param X a matrix whose columns refer to parameters and whose rows to samples
##' @param main the plot title (if desired)
##' @param labels a character vector of length \code{ncol(X)} to label parameters 
##' @param file a filename to save pdf output plot to (else plotted to default to device)
##' @return None
##' @author Pete Dodd
##' @examples
##' corplot(matrix(rnorm(3e4),ncol=3),labels=c('x','y','z'),main='3D isotropic Gaussian')
corplot <- function(X,main='',labels=NULL,file='',weights=1,points=FALSE){
  if(!is.null(labels)){
    colnames(X) <- labels
    f <- 1/(length(labels))^.2
  } else {f <- 1}
  if(file!='') pdf(file)
  if(!points){                          #hexbin version
      myp <- lattice::splom(X,main=main,xlab='',
                        panel=hexbin::panel.hexbinplot,
                        colramp=hexbin::BTC,
                        diag.panel = function(x, ...){
                            yrng <- lattice::current.panel.limits()$ylim
                            h <- hist(x, plot = FALSE,breaks=30)
                            breaks <- h$breaks
                            nB <- length(breaks)
                            y <- h$density
                            x <- h$mids
                            y <- yrng[1] + 0.95 * diff(yrng) * y / max(y)
                            lattice::panel.lines(x,y,type='s',col='blue')
                            lattice::diag.panel.splom(x, ...)
                        },
                        lower.panel = function(x, y, ...){
                            hexbin::panel.hexbinplot(x, y, ...)
                            lattice::panel.loess(x, y, ..., col = 'red')
                        },
                        upper.panel = function(x, y, ...){
                            pl <- lattice::current.panel.limits()
                            lattice::panel.text(x=mean(pl$xlim),y=mean(pl$ylim),sprintf('%1.2f',cor(x,y,use="complete")),col = 2,cex=0.75)
                        },
                        pscale=5,varname.cex=1, varname.col='red',axis.text.cex=.5*f
                        )
      
  } else{                               #nonhexbin version
      myp <- lattice::splom(X,main=main,xlab='',
                        colramp=hexbin::BTC,
                        diag.panel = function(x, ...){
                            yrng <- lattice::current.panel.limits()$ylim
                            h <- hist(x, plot = FALSE,breaks=30)
                            breaks <- h$breaks
                            nB <- length(breaks)
                            y <- h$density
                            x <- h$mids
                            y <- yrng[1] + 0.95 * diff(yrng) * y / max(y)
                            lattice::panel.lines(x,y,type='s',col='blue')
                            lattice::diag.panel.splom(x, ...)
                        },
                            lower.panel = function(x, y, ...){
                                clz <- heat.colors(1e3)
                                az <- as.numeric(cut(weights,breaks=length(clz)))
                                ptclz <- clz[az]
                                lattice::panel.xyplot(x,y,pch='+',col=1,alpha=az/length(clz))
                        },
                        upper.panel = function(x, y, ...){
                            pl <- lattice::current.panel.limits()
                            lattice::panel.text(x=mean(pl$xlim),y=mean(pl$ylim),sprintf('%1.2f',cor(x,y,use="complete")),col = 2,cex=0.75)
                        },
                        pscale=5,varname.cex=1, varname.col='red',axis.text.cex=.5*f
                        )
  }
  
  print(myp)
  if(file!='')dev.off()
}
