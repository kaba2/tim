// Description: Differential entropy of a generalized normal distribution

#ifndef TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_H
#define TIM_DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Differential entropy of a generalized gaussian distribution.
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
	Real generalizedGaussianDifferentialEntropy(
		integer dimension, 
		const PASTEL_NO_DEDUCTION(Real)& shape, 
		const PASTEL_NO_DEDUCTION(Real)& scale);

}

#include "tim/core/differential_entropy_generalized_normal.hpp"

#endif
