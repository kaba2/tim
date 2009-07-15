// Description: Kozachenko-Leonenko differential entropy estimation

#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_H
#define TIM_DIFFERENTIAL_ENTROPY_KL_H

#include "tim/core/signal.h"

namespace Tim
{

	//! Computes the differential entropy of a signal.
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

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
