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
	sigma >= 0
	kNearest > 0	

	Usage:

	Assume that you want to analyze the information
	flow from signal B to signal A while
	excluding the effect of signals {C_i}:

	A : Z -> R^a
	B : Z -> R^b
	C_i : Z -> R^c_i

	One way to do this is via Multivariate
	Transfer Entropy (MTE). To use this function, you are expected 
	to delay-embed all of the signals A, B and {C_i} properly 
	before giving them as input to this function (see
	'tim/core/embed.h' for this). This results in the signals:

	A' : Z -> R^a'
	B' : Z -> R^b'
	C'_i : Z -> R^c'_i

	In addition to these signals, you are required to pass
	the 'future' of the signal A. Assume you want to do the delay-
	embedding of A with embedding factor k. Then the future
	of A is the last component of the signal A' if you
	delay-embed it with embedding factor (k + 1).
	You can obtain the future conveniently as follows:
	
	SignalPtr aEmbeddedBig = embed(aSignal, k + 1, shift, step);
	SignalPtr aEmbedded = slice(aEmbeddedBig, 0, k);
	SignalPtr aFuture = slice(aEmbeddedBig, k, k + 1);

	The embedding can result in signals that have slightly 
	different number of samples in them. This is handled by using
	the minimum sample count and ignoring the rest of the samples.

	If cEmbeddedSet is empty, this function effectively 
	computes bivariate transfer entropy.
	*/
	TIMCORE void transferEntropy(
		const std::vector<SignalPtr>& aEnsemble,
		const std::vector<SignalPtr>& aFutureEnsemble,
		const std::vector<SignalPtr>& bEnsemble,
		const Array<2, SignalPtr>& cEnsembleSet,
		integer sigma,
		integer kNearest,
		std::vector<real>& estimateSet);

}

#endif
