#ifndef TIM_RENYI_ENTROPY_NORMAL_HPP
#define TIM_RENYI_ENTROPY_NORMAL_HPP

#include "tim/core/renyi_entropy_normal.h"
#include "tim/core/differential_entropy_normal.h"

#include <cmath>

namespace Tim
{

	template <typename Real>
	Real gaussianRenyiEntropy(
		const PASTEL_NO_DEDUCTION(Real)& q,
		integer dimension,
		const PASTEL_NO_DEDUCTION(Real)& covarianceDeterminant)
	{
		PENSURE_OP(q, >, 0);
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(covarianceDeterminant, >=, 0);

		if (q == 1)
		{
			return normalDifferentialEntropy<Real>(
				dimension,
				covarianceDeterminant);
		}

		// H_q^*(X) = (n / 2)[log(2 pi) - (log(q) / (1 - q))] + 
		//            (1 / 2)log(det(C))

		return ((Real)dimension / 2) * 
			(std::log(2 * constantPi<Real>()) - std::log(q) / (1 - q)) +
			std::log(covarianceDeterminant) / 2;
	}

}

#endif
