// Description: Estimation of multivariate transfer entropy.

#ifndef TIM_TRANSFER_ENTROPY_H
#define TIM_TRANSFER_ENTROPY_H

#include "tim/core/signal.h"

#include <vector>

namespace Tim
{

	//! Computes the multivariate transfer entropy.
	/*!
	Preconditions:
	flowFrom >= 0
	flowFrom < signalSet.height()
	flowTo >= 0
	flowTo < signalSet.height()
	fromDimension > 0
	timeWindowRadius >= 0
	kNearest > 0
	Signal_Iterator dereferences to SignalPtr.

	embeddedSet:
	An array containing signals for which the multivariate
	transfer entropy (MTE) is computed.
	Each row of the array represents a signal from the 
	same source, and each column represents a trial. I.e. the
	same test has been repeated and measured for a number of
	signal sources and the signals have been measured in parallel
	for each trial. The signals must already have been delay-embedded.

	xIndex, yIndex:
	The row indices of the signals x and y in 'embeddedSet' 
	between which the directed information flow from x to y is of interest. 
	The effect of other signals will be removed from the end result 
	so that, for example, if the connectivity between signals x, y, 
	and z is x->y->z, then x->z won't report information flow.

	xFutureBegin, xFutureEnd:
	The future of signal x _before_ the delay-embedding, for each trial.
	Let x : R -> R^n. The future of x is then defined by:
	xFuture[t] = last n components of xEmbedded[t + 1].
	If you used the 'delayEmbed' function to do the delay-embedding,
	you can use the 'delayEmbedFuture' function to compute the future.

	timeWindowRadius:
	The radius of the time-window in samples to use for MTE estimation
	at each time instant. Smaller 'timeWindowRadius's give better temporal resolution,
	but greater errors. In the other direction, if the underlying
	connectivity pattern stays fixed, then the error can be made 
	smaller by using a larger 'timeWindowRadius'. In the general case, the MTE 
	estimate is averaged over the time window and thus its correspondence 
	to the actual temporal MTE is dependent on how rapidly the underlying
	connectivity changes.

	kNearest:
	The k:th neighbor to use in the estimation.

	estimateSet (output):
	For each time instant, an MTE estimate.
	
	The embedding and shifting can result in signals that have slightly 
	different number of samples in them. This is handled by using
	the minimum sample count and ignoring the rest of the samples.

	If the 'signalSet' contains only two signals, this function
	reduces to the bivariate transfer entropy.
	*/
	template <typename Signal_Iterator>
	void transferEntropy(
		const Array<SignalPtr, 2>& embeddedSet,
		integer xIndex,
		integer yIndex,
		const Signal_Iterator& xFutureBegin,
		const Signal_Iterator& xFutureEnd,
		integer timeWindowRadius,
		integer kNearest,
		std::vector<real>& estimateSet);

}

#include "tim/core/transfer_entropy.hpp"

#endif
