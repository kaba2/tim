// Description: Temporal generic entropy estimation
// Detail: Encapsulates properties common to k-nn entropy estimators.

#ifndef TIM_GENERIC_ENTROPY_T_H
#define TIM_GENERIC_ENTROPY_T_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes temporal generic entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
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

	kNearest:
	The k:th nearest neighbor that is used to
	estimate generic entropy.

	filter:
	An array of coefficients by which to weight the results
	in the time-window. The center of the array corresponds 
	to the current time instant. The width of the array can 
	be arbitrary but must be odd. The coefficients must sum
	to a non-zero value.

	Returns:
	The number of time instants that had an
	undefined estimate. If not all estimates
	were undefined, they were reconstructed from 
	the defined estimates using interpolation.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename EntropyAlgorithm,
		typename Real_OutputIterator,
		typename Real_Filter_Iterator>
	integer temporalGenericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter);

	//! Computes temporal generic entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalGenericEntropy(
		signalSet,
		entropyAlgorithm,
		timeWindowRadius,
		result,
		kNearest,
		constantRange((real)1, 1));

	See the documentation for that function.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename EntropyAlgorithm,
		typename Real_OutputIterator>
	integer temporalGenericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest = 1);

}

#include "tim/core/generic_entropy_t.hpp"

#endif
