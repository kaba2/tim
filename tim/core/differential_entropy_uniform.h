// Description: Differential entropy of a uniform distribution

#ifndef TIM_DIFFERENTIAL_ENTROPY_UNIFORM_H
#define TIM_DIFFERENTIAL_ENTROPY_UNIFORM_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real/real_concept.h>
#include <cmath>

namespace Tim
{

	//! Differential entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume > 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})
	*/
	template <typename Real>
	Real differentialEntropyUniform(
		const NoDeduction<Real>& supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		// Let X be a random variable in R^n 
		// with a probability density function p.
		// Let S subset R^n be the support set of p.
		// The measure of S is m(S). Then because
		// of uniformity p(x) = 1 / m(S).
		// Now
		// H(X) = -int_S p(x) ln(p(x)) dx
		// = -int_S (1 / m(s)) ln(1 / m(S)) dx
		// = -int_S (1 / m(s)) (ln(1) - ln(m(S))) dx
		// = -int_S (1 / m(s)) (-ln(m(S))) dx
		// = (ln(m(S)) / m(S)) int_S 1 dx
		// = (ln(m(S)) / m(S)) m(S)
		// = ln(m(S))

		return std::log(supportVolume);
	}

}

#endif
