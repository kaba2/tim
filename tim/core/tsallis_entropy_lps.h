// Description: Tsallis entropy estimation
// Detail: Leonenko-Pronzato-Savani nearest-neighbor estimator

#ifndef TIM_TSALLIS_ENTROPY_LPS_H
#define TIM_TSALLIS_ENTROPY_LPS_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>

namespace Tim
{

	//! Returns the actual k:th neighbor to use given a suggestion.
	/*!
	Preconditions:
	q > 0
	kNearestSuggestion >= 0

	Returns:
	If kNearestSuggestion == 0, then k = 2 * ceil(q).
	If 0 < kNearestSuggestion < q - 1, then k = ceil(q - 1).
	If kNearestSuggestion == q - 1, then k = kNearestSuggestion + 1.
	Otherwise, k = kNearestSuggestion.

	The k in the Leonenko-Pronzato-Savani estimator can't be set 
	freely because the algorithm is not defined when k <= q - 1.
	This functions helps to decide a proper k.
	*/
	TIM real tsallisDecideK(real q, integer kNearestSuggestion);

	// Temporal Tsallis entropy
	// ------------------------

	//! Computes temporal Tsallis entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearestSuggestion >= 0
	q > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	q:
	The exponent in the definition of Tsallis entropy.
	If q == 1, the result of temporalDifferentialEntropyKl() 
	is returned instead.
	If q < 1, the results have huge errors: you
	should not use this estimator for those values.

	kNearestSuggestion:
	A suggestion for the k:th nearest neighbor that should be
	used for estimation. The k can't be set	freely because the 
	estimation algorithm is only defined for k > q - 1. 
	Value zero means an accurate (q-dependent) default is used.
	The actual k that is used is given by tsallisDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.
	*/

	template <
		typename SignalPtr_Range, 
		typename Real_Filter_Iterator>
	Signal temporalTsallisEntropyLps(
		const SignalPtr_Range& signalSet,
		integer timeWindowRadius,
		real q,
		integer kNearestSuggestion,
		const boost::iterator_range<Real_Filter_Iterator>& filter);

	//! Computes temporal Tsallis entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalTsallisEntropyLps(
		signalSet, timeWindowRadius,
		q, kNearestSuggestion,
		constantRange((real)1, 1));

	See the documentation for that function.
	*/

	template <typename SignalPtr_Range>
	Signal temporalTsallisEntropyLps(
		const SignalPtr_Range& signalSet,
		integer timeWindowRadius,
		real q = 2,
		integer kNearestSuggestion = 0);

	// Tsallis entropy
	// ---------------

	//! Computes Tsallis entropy of a signal.
	/*!
	Preconditions:
	signalSet contains Signal's.
	q > 0
	kNearestSuggestion >= 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	q:
	The exponent in the definition of Tsallis entropy.
	If q == 1, the result of differentialEntropyKl() 
	is returned instead.
	If q < 1, the results have huge errors: you
	should not use this estimator for those values.

	kNearestSuggestion:
	A suggestion for the k:th nearest neighbor that should be
	used for estimation. The k can't be set	freely because the 
	estimation algorithm is only defined for k > q - 1. 
	Value zero means an accurate (q-dependent) default is used.
	The actual k that is used is given by tsallisDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.

	Returns:
	A Tsallis entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/

	template <typename SignalPtr_Range>
	real tsallisEntropyLps(
		const SignalPtr_Range& signalSet,
		real q = 2,
		integer kNearestSuggestion = 0);

}

#include "tim/core/tsallis_entropy_lps.hpp"

#endif
