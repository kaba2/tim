// Description: Generic entropy estimation
// Detail: Encapsulates properties common to k-nn entropy estimators.

#ifndef TIM_GENERIC_ENTROPY_H
#define TIM_GENERIC_ENTROPY_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	// Temporal generic entropy
	// -----------------------------

	//! Computes temporal generic entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	maxRelativeError >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	entropyAlgorithm:
	Encapsulates the specifics of the used
	entropy algorithm.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	result:
	A real output iterator, denoting the start
	of the region where the sequence of temporal
	differential entropies are to be stored.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate generic entropy.

	Returns:
	True, if all estimates were succesfully estimated.
	False, if at least one estimate could not be 
	estimated and was given NaN.
	One can later apply the theory of irregular 
	sampling to reconstruct these missing values
	assuming continuity.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename EntropyAlgorithm,
		typename Real_OutputIterator>
	bool temporalGenericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes temporal generic entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalGenericEntropy(
		forwardRange(constantIterator(signal)),
		entropyAlgorithm,
		timeWindowRadius, result, maxRelativeError, 
		kNearest);

	See the documentation for that function.
	*/

	template <
		typename EntropyAlgorithm, 
		typename Real_OutputIterator>
	bool temporalGenericEntropy(
		const SignalPtr& signal,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError = 0,
		integer kNearest = 1);

	// Generic entropy
	// --------------------

	//! Computes generic entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	maxRelativeError >= 0
	signalSet contains SignalPtr's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	entropyAlgorithm:
	Encapsulates the specifics of the used
	entropy algorithm.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate generic entropy.

	Returns:
	A generic entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/

	template <
		typename SignalPtr_Iterator,
		typename EntropyAlgorithm>
	real genericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes generic entropy of a signal.
	/*!
	This is a convenience function that calls:

	genericEntropy(
		forwardRange(constantIterator(signal)),
		entropyAlgorithm,
		maxRelativeError,
		kNearest);

	See the documentation for that function.
	*/

	template <typename EntropyAlgorithm>
	real genericEntropy(
		const SignalPtr& signal,
		const EntropyAlgorithm& entropyAlgorithm,
		real maxRelativeError = 0,
		integer kNearest = 1);

}

#include "tim/core/generic_entropy.hpp"

#endif
