// Description: Tsallis entropy estimation
// Detail: Leonenko-Pronzato-Savani nearest-neighbor estimator

#ifndef TIM_TSALLIS_ENTROPY_LPS_H
#define TIM_TSALLIS_ENTROPY_LPS_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	// Temporal Tsallis entropy
	// ------------------------

	//! Computes temporal Tsallis entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	maxRelativeError >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	result:
	A real output iterator, denoting the start
	of the region where the sequence of temporal
	Tsallis entropies are to be stored.

	q:
	The exponent in the definition of Tsallis entropy.
	In case q == 1, the result of 
	temporalDifferentialEntropyKl() is returned instead.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate Tsallis entropy.

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
		typename Real_OutputIterator>
	bool temporalTsallisEntropyLps(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real q = 2,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes temporal Tsallis entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalTsallisEntropyLps(
		forwardRange(constantIterator(signal)),
		timeWindowRadius, q, result, maxRelativeError, 
		kNearest);

	See the documentation for that function.
	*/

	template <typename Real_OutputIterator>
	bool temporalTsallisEntropyLps(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real q = 2,
		real maxRelativeError = 0,
		integer kNearest = 1);

	// Tsallis entropy
	// ---------------

	//! Computes Tsallis entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	maxRelativeError >= 0
	signalSet contains SignalPtr's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	q:
	The exponent in the definition of Tsallis entropy.
	In case q == 1, the result of 
	differentialEntropyKl() is returned instead.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate Tsallis entropy.

	Returns:
	A Tsallis entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/

	template <typename SignalPtr_Iterator>
	real tsallisEntropyLps(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		real q = 2,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes Tsallis entropy of a signal.
	/*!
	This is a convenience function that calls:

	tsallisEntropyLps(
		signal, q, 
		maxRelativeError,
		kNearest);

	See the documentation for that function.
	*/
	TIM real tsallisEntropyLps(
		const SignalPtr& signal,
		real q = 2,
		real maxRelativeError = 0,
		integer kNearest = 1);

}

#include "tim/core/tsallis_entropy_lps.hpp"

#endif
