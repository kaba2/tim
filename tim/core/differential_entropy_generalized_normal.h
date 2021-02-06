// Description: Differential entropy of a generalized normal distribution

#ifndef TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_H
#define TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Differential entropy of a generalized normal distribution.
	/*!
	Preconditions:
	dimension > 0
	shape > 0
	scale > 0

	dimension:
	The dimension of the distribution.

	shape:
	The shape parameter of the distribution.

	scale:
	The scale parameter of the distribution.
	*/
	template <typename Real>
	Real differentialEntropyGeneralizedNormal(
		integer dimension, 
		const NoDeduction<Real>& shape, 
		const NoDeduction<Real>& scale)
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
			lnGamma<dreal>(invShape));
	}

}

#endif
