heuristicTuning <- function ( sliceParams, targetRatio = 0.5, K )
{
	for ( I in 1:length(sliceParams$nExpands) )
	{
		denom <- sliceParams$nExpands[I] + sliceParams$nShrinks[I];

		if ( denom > 0 )
		{
			ratio <- sliceParams$nExpands[I] / denom;

			if ( 0.0 == ratio )
				ratio <- 1.0 / denom;

			multiplier <- (ratio / targetRatio);

			# Modify Initial Interval Width
			sliceParams$intervalWidth[I] <-
				sliceParams$intervalWidth[I] * multiplier;
		}
	}

	#  Reset Interval Width Counters
	sliceParams$nExpands   <- array(0,c(K,1));
	sliceParams$nShrinks   <- array(0,c(K,1));
	sliceParams$nProposals <- 0;

	# Double Number of Iterations to Next Adaptation
	sliceParams$nIterPerAdapt <- sliceParams$nIterPerAdapt * 2;

	return ( sliceParams );
}
