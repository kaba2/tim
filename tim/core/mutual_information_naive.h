// Description: Mutual information estimation using naive algorithms
// Detail: Includes computation from differential entropies and from a binning.

#ifndef TIM_MUTUAL_INFORMATION_NAIVE_H
#define TIM_MUTUAL_INFORMATION_NAIVE_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

#include <pastel/sys/smallset.h>

namespace Tim
{

	//! Computes pairwise 1d mutual information by binning.
	/*!
	Preconditions:
	bins > 0

	signal:
	Signal that contains n

	bins:
	Number of bins to use to estimate a 1d probability
	distribution function (pdf). For the joint pdf, a
	2d array of extents bins x bins is used.

	result (output):
	The element (i, j) contains the mutual information
	between the i:th and j:th 1d marginal signals of
	the 'signal'.

	The approximation of probability distribution functions
	using binning does not generalize practically to higher
	dimensions than 2, because of the exponential explosion
	of the number of needed bins. Therefore, using this
	technique, the mutual information can only be computed
	between two 1d signals. Given a multi-dimensional 
	signal, this function computes pairwise mutual
	information between the 1d marginal signals of the
	given signal.

	This technique is not recommended because it tends
	to give large errors in estimation. The intent of
	this function is to demonstrate the non-applicability
	of the technique.
	*/
	TIM Array<real> mutualInformationFromBinning(
		const SignalPtr& signal,
		integer bins);

	//! Computes mutual information from entropies.
	/*!
	Preconditions:
	kNearest > 0

	signalSet:
	A set of signals between which the mutual information
	is computed.

	kNearest:
	The k:th neighbor to use in the estimation of the
	differential entropies.

	normBijection:
	The norm to use to do the estimations.
	See 'pastel/math/normbijections.h'.

	This technique is not recommended because it tends
	to give large errors in estimation. The intent of
	this function is to demonstrate the non-applicability
	of the technique.
	*/
	template <typename NormBijection>
	real mutualInformationFromEntropy(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		const NormBijection& normBijection);

}

#include "tim/core/mutual_information_naive.hpp"

#endif
