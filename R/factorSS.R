factorSliceSampler <- function( currValue, LogLikelihood, sliceParams )
{
	for ( I in 1:length(currValue) )
	{
		width  <- sliceParams$intervalWidth[I];
		factor <- sliceParams$factors[,I];

		#####  Sample Random Height under Density  #####
		height <- LogLikelihood( currValue ) - rexp(1);
 ## cat('SS =',I,' height =',height,' \n\n');
                
		#####  Construct Slice Approximation  #####
		lowerBnd <- -1.0 * width * runif(1);
		upperBnd <- lowerBnd + width;
 ## cat('lb =',lowerBnd,' ub =',upperBnd,' \n\n');
		#####  Step Out Procedure  #####
		while ( height < LogLikelihood( currValue + lowerBnd * factor ) )
		{
			sliceParams$nExpands[I] <- sliceParams$nExpands[I] + 1;
			lowerBnd <- lowerBnd - width;
		}
## cat('SS1=',I,' \n\n');
                
		while ( height < LogLikelihood( currValue + upperBnd * factor ) )
		{
			sliceParams$nExpands[I] <- sliceParams$nExpands[I] + 1;
			upperBnd <- upperBnd + width;
		}
## cat('SS2=',I,' \n\n');
                
		#####  Sample from Approximate Interval  #####
		while ( TRUE )
		{
			#  Sample Uniformly from [ lowerBnd, upperBnd ]
			proposal <- lowerBnd + runif(1) * ( upperBnd - lowerBnd );

			#  Test to determine if proposal falls within target slice
			if ( height < LogLikelihood( currValue + proposal * factor ) )
			{
				currValue <- currValue + proposal * factor;
				break;
			}

			#####  Shrink Interval if Proposal Failed  #####
			sliceParams$nShrinks[I] <- sliceParams$nShrinks[I] + 1;

			if ( proposal < 0 ) { lowerBnd <- proposal; }
			if ( proposal > 0 ) { upperBnd <- proposal; }
		}
	}

	sliceParams$nProposals <- sliceParams$nProposals + 1;
	return ( list( beta = currValue, sliceParams = sliceParams ) );
}
