// Description: Differential entropy estimation via Kozachenko-Leonenko 

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
	maxRelativeError >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	A measure of distance to use. See 
	'pastel/math/normbijection.txt'	for documentation.

	Returns:
	An array of temporal differential entropy
	estimates.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator,
		typename NormBijection>
	void temporalDifferentialEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropy(
		signalSet, timeWindowRadius, result, 
		maxRelativeError, kNearest, Euclidean_NormBijection());

	See the documentation for that function.
	*/
	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator>
	void temporalDifferentialEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropy(
		signalSet, timeWindowRadius, result, 
		maxRelativeError, kNearest, Euclidean_NormBijection());

	See the documentation for that function.
	*/

	template <
		typename Real_OutputIterator,
		typename NormBijection>
	void temporalDifferentialEntropy(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes temporal differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalDifferentialEntropy(
		forwardRange(constantIterator(signal)),
		timeWindowRadius, result, maxRelativeError, 
		kNearest, Euclidean_NormBijection());

	See the documentation for that function.
	*/

	template <typename Real_OutputIterator>
	void temporalDifferentialEntropy(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError = 0,
		integer kNearest = 1);

	// Differential entropy
	// --------------------

	//! Computes differential entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	maxRelativeError >= 0
	signalSet contains SignalPtr's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	maxRelativeError:
	The maximum relative error allowed for
	distance in nearest neighbor searching.
	Zero gives exact matches. Higher values can
	result in improved performance.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	A measure of distance to use. See 
	'pastel/math/normbijection.txt'	for documentation.

	Returns:
	A differential entropy estimate.
	*/

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropy(
		signalSet, maxRelativeError, kNearest, 
		Euclidean_NormBijection());

	See the documentation for that function.
	*/

	template <typename SignalPtr_Iterator>
	real differentialEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropy(
		forwardRange(constantIterator(signal)),
		maxRelativeError,
		kNearest,
		normBijection);

	See the documentation for that function.
	*/

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	This is a convenience function that calls:

	differentialEntropy(
		signal, maxRelativeError,
		kNearest, Euclidean_NormBijection());

	See the documentation for that function.
	*/
	TIM real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError = 0,
		integer kNearest = 1);

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
