tuneFactors <- function( sliceParams, sampleMean, sampleCov, nStored, K )
{
	normalizedSampleCov  <- (1.0 / (nStored-1)) *
	(sampleCov - sampleMean %*% t(sampleMean) * (1.0 / nStored));

	# Update Factors using Eigenvectors of Sample Covariance Estimate
	sliceParams$factors <- eigen( normalizedSampleCov )$vectors;

	#  Reset Interval Width Counters
	sliceParams$nExpands <- array(0,c(K,1));
	sliceParams$nShrinks <- array(0,c(K,1));

	sliceParams$nIterPerAdapt <- 1;
	sliceParams$nProposals    <- 0;

	return ( sliceParams );
}
