// Description: Differential entropy estimation
// Detail: Kozachenko-Leonenko nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_H
#define TIM_DIFFERENTIAL_ENTROPY_KL_H

#include "tim/core/signal.h"

#include <pastel/sys/iterator_range.h>

namespace Tim
{

	// Temporal differential entropy
	// -----------------------------

	//! Computes temporal differential entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	A measure of distance to use. See 
	'pastel/math/normbijection.txt'	for documentation.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection,
		typename Real_Filter_Iterator>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection,
		const ForwardIterator_Range<Real_Filter_Iterator>& filter);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropyKl(
		signalSet, timeWindowRadius,
		kNearest, Default_NormBijection());

	See the documentation for that function.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropyKl(
		signalSet, timeWindowRadius,
		kNearest, Default_NormBijection());

	See the documentation for that function.
	*/

	template <typename SignalPtr_Iterator>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest = 1);

	// Differential entropy
	// --------------------

	//! Computes differential entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	signalSet contains SignalPtr's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	A measure of distance to use. See 
	'pastel/math/normbijection.txt'	for documentation.

	Returns:
	A differential entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropyKl(
		signalSet, kNearest, 
		Default_NormBijection());

	See the documentation for that function.
	*/

	template <typename SignalPtr_Iterator>
	real differentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer kNearest = 1);

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
