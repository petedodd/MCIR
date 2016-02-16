stepOutSS <- function( currValue, LogLikelihood, sliceParams )
{
	proposedParam <- currValue;

	for ( I in 1:length(currValue) )
	{
		# Current Initial Interval Width
		width <- sliceParams$intervalWidth[I];

		#####  Sample Random Height under Density  #####
		height <- LogLikelihood( proposedParam ) - rexp(1);

		#####  Construct Slice Approximation  #####
		lowerBound <- proposedParam[I] - width * runif(1);
		upperBound <- lowerBound       + width;

		#####  Step Out Procedure  #####
		proposedParam[I] <- lowerBound;
		while ( height < LogLikelihood( proposedParam ) )
		{
			sliceParams$nExpands[I] <- sliceParams$nExpands[I] + 1;
			lowerBound       <- lowerBound - width;
			proposedParam[I] <- lowerBound;
		}

		proposedParam[I] <- upperBound;
		while ( height < LogLikelihood( proposedParam ) )
		{
			sliceParams$nExpands[I] <- sliceParams$nExpands[I] + 1;
			upperBound       <- upperBound + width;
			proposedParam[I] <- upperBound;
		}

		#####  Sample from Approximate Interval  #####
		proposal <- proposedParam;

		while ( TRUE )
		{
			#  Sample Uniformly from [ lowerBound, upperBound ]
			proposal[I] <- lowerBound +
				runif(1) * ( upperBound - lowerBound );

			#  Test to determine if proposal falls within target slice
			if ( height < LogLikelihood( proposal ) )
			{
				proposedParam <- proposal;
				break;
			}

			#####  Shrink Interval if Proposal Failed  #####
			sliceParams$nShrinks[I] <- sliceParams$nShrinks[I] + 1;

			if ( proposal[I] < currValue[I] )
			{
				lowerBound <- proposal[I];
			}
			else
			{
				upperBound <- proposal[I];
			}
		}
	}

	sliceParams$nProposals <- sliceParams$nProposals + 1;
	return ( list( beta = proposedParam, sliceParams = sliceParams ) );
}