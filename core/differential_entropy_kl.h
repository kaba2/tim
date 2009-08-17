// Description: Kozachenko-Leonenko differential entropy estimation

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
	kNearest > 0
	maxRelativeError >= 0

	Let X be a real random variable with
	probability density function p(x).
	The differential entropy H of X is defined by:
	H(X) = -int(p(x) ln(p(x)), x = -Infinity..Infinity)
	where 'int' denotes integration.
	Note that differential entropy is _not_
	a generalization of (Shannon's discrete) 
	entropy and thus does not measure information content.
	*/

	template <
		typename Signal_Iterator, 
		typename Real_OutputIterator,
		typename NormBijection>
	void temporalDifferentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	template <
		typename Signal_Iterator, 
		typename Real_OutputIterator>
	void temporalDifferentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes temporal differential entropy of a signal.

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
	differentialEntropy(signal, timeWindowRadius, result, 
		maxRelativeError, kNearest, Euclidean_NormBijection());
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
	*/

	template <
		typename Signal_Iterator, 
		typename NormBijection>
	real differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	differentialEntropy(
		signalSet, maxRelativeError, kNearest, 
		Euclidean_NormBijection());
	*/

	template <typename Signal_Iterator>
	real differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		real maxRelativeError = 0,
		integer kNearest = 1);

	//! Computes differential entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	maxRelativeError >= 0

	differentialEntropy(
		forwardRange(constantIterator(signal)),
		maxRelativeError,
		kNearest,
		normBijection);
	*/

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection);

	//! Computes differential entropy of a signal.
	/*!
	differentialEntropy(
		signal, maxRelativeError,
		kNearest, Euclidean_NormBijection());
	*/
	TIMCORE real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError = 0,
		integer kNearest = 1);

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
