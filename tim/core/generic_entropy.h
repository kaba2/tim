// Description: Generic entropy estimation
// Detail: Encapsulates properties common to k-nn entropy estimators.

#ifndef TIM_GENERIC_ENTROPY_H
#define TIM_GENERIC_ENTROPY_H

#include "tim/core/signal.h"

#include <pastel/sys/iteratorrange.h>

namespace Tim
{

	//! Computes generic entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	signalSet contains SignalPtr's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	entropyAlgorithm:
	Encapsulates the specifics of the used
	entropy algorithm.

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
		integer kNearest = 1);

}

#include "tim/core/generic_entropy.hpp"

#endif
