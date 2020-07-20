// Description: Differential entropy of a normal distribution

#ifndef TIM_DIFFERENTIAL_ENTROPY_NORMAL_H
#define TIM_DIFFERENTIAL_ENTROPY_NORMAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real/real_concept.h>
#include <cmath>

namespace Tim
{

	//! Differential entropy of a normal distribution.
	/*!
	Preconditions:
	dimension > 0
	covarianceDeterminant >= 0

	dimension:
	The dimension of the distribution.

	covarianceDeterminant:
	The determinant of the covariance matrix of the distribution.
	*/
	template <typename Real>
	Real differentialEntropyNormal(
		integer dimension, 
		const NoDeduction<Real>& covarianceDeterminant)
	{
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(covarianceDeterminant, >=, 0);

		static const Real ConstantFactor = std::log(
			2 * constantPi<Real>()) + 1;

		return 0.5 * (
			std::log(covarianceDeterminant) + 
			dimension * ConstantFactor);
	}

}

#endif
