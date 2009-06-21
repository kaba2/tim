#ifndef TIM_DIFFERENTIAL_ENTROPY_H
#define TIM_DIFFERENTIAL_ENTROPY_H

#include "tim/core/signal.h"

namespace Tim
{

	//! Computes the differential entropy of a signal.
	/*!
	Let X be a real random variable with
	probability density function p(x).
	The differential entropy H of X is defined by:
	H(X) = -int(p(x) ln(p(x)), x = -Infinity..Infinity)
	where 'int' denotes integration.
	Note that differential entropy is _not_
	a generalization of (Shannon's discrete) 
	entropy and thus does not measure information content.
	*/

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	//! Returns the differential entropy of a gaussian random variable.
	/*!
	Preconditions:
	dimension > 0
	variance >= 0
	*/
	TIMCORE real gaussianDifferentialEntropy(
		integer dimension, real varianceDeterminant);

	//! Returns the differential entropy of a uniform random variable.
	/*!
	Preconditions:
	supportVolume >= 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})
	*/
	TIMCORE real uniformDifferentialEntropy(
		real supportVolume);

	TIMCORE real generalizedGaussianDifferentialEntropy(
		integer dimension, real shape, real scale);

}

#include "tim/core/differential_entropy.hpp"

#endif
