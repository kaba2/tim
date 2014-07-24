#ifndef TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_HPP

#include "tim/core/differential_entropy_generalized_normal.h"

#include <cmath>

namespace Tim
{

	template <typename Real>
	Real differentialEntropyGeneralizedNormal(
		integer dimension, 
		const PASTEL_NO_DEDUCTION(Real)& shape, 
		const PASTEL_NO_DEDUCTION(Real)& scale)
	{
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(shape, >, 0);
		PENSURE_OP(scale, >, 0);

		// Let
		// a = scale
		// b = shape
		//
		// Then
		//
		// differential entropy 
		// = (1 / b) - log(b / (2a gamma(1 / b)))
		// = (1 / b) - (log(b / (2a)) - log(gamma(1 / b)))
		// = (1 / b) - log(b / (2a)) + log(gamma(1 / b))
		// = (1 / b) + log((2a) / b) + log(gamma(1 / b))
		// = (1 / b) + log((2a) / b) + log(gamma(1 / b))
		//
		// The point of this derivation is to evaluate lnGamma 
		// instead of gamma for better numerical behaviour.

		Real invShape = inverse(shape);

		return dimension * (
			inverse(shape) + 
			std::log(2 * scale * invShape) +
			lnGamma<real>(invShape));
	}

}

#endif
