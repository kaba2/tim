// Description: Differential entropy estimation
// Detail: Kozachenko-Leonenko nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_H
#define TIM_DIFFERENTIAL_ENTROPY_KL_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

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

	result:
	A real output iterator, denoting the start
	of the region where the sequence of temporal
	differential entropies are to be stored.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	A measure of distance to use. See 
	'pastel/math/normbijection.txt'	for documentation.

	Returns:
	The number of time instants that had an
	undefined estimate. If not all estimates
	were undefined, they were reconstructed from 
	the defined estimates using interpolation.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator,
		typename NormBijection>
	integer temporalDifferentialEntropyKl(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropyKl(
		signalSet, timeWindowRadius, result, 
		kNearest, Default_NormBijection());

	See the documentation for that function.
	*/
	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator>
	integer temporalDifferentialEntropyKl(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest = 1);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropyKl(
		signalSet, timeWindowRadius, result, 
		kNearest, Default_NormBijection());

	See the documentation for that function.
	*/

	template <
		typename Real_OutputIterator,
		typename NormBijection>
	integer temporalDifferentialEntropyKl(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropyKl(
		constantRange(signal),
		timeWindowRadius, result, 
		kNearest, Default_NormBijection());

	See the documentation for that function.
	*/

	template <typename Real_OutputIterator>
	integer temporalDifferentialEntropyKl(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
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
		const ForwardRange<SignalPtr_Iterator>& signalSet,
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
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer kNearest = 1);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropyKl(
		constantRange(signal),
		kNearest,
		normBijection);

	See the documentation for that function.
	*/

	template <typename NormBijection>
	real differentialEntropyKl(
		const SignalPtr& signal,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropyKl(
		signal, kNearest, 
		Default_NormBijection());

	See the documentation for that function.
	*/
	TIM real differentialEntropyKl(
		const SignalPtr& signal,
		integer kNearest = 1);

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
