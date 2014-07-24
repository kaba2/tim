#ifndef TIM_DIFFERENTIAL_ENTROPY_UNIFORM_HPP
#define TIM_DIFFERENTIAL_ENTROPY_UNIFORM_HPP

#include "tim/core/differential_entropy_uniform.h"

#include <cmath>

namespace Tim
{

	template <typename Real>
	Real uniformDifferentialEntropy(
		const PASTEL_NO_DEDUCTION(Real)& supportVolume)
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

