#ifndef TIM_DIFFERENTIAL_ENTROPY_NORMAL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_NORMAL_HPP

#include "tim/core/differential_entropy_normal.h"

#include <cmath>

namespace Tim
{

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
