#ifndef TIM_TSALLIS_ENTROPY_NORMAL_HPP
#define TIM_TSALLIS_ENTROPY_NORMAL_HPP

#include "tim/core/tsallis_entropy_normal.h"
#include "tim/core/differential_entropy_normal.h"

namespace Tim
{

	template <typename Real>
	Real normalTsallisEntropy(
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

		// A = 1 / ((2 pi)^(n / 2) det(C)^(1 / 2))
		
		Real A = inverse(
			std::pow(2 * constantPi<Real>(), (Real)dimension / 2) *
			std::sqrt(covarianceDeterminant));

		// I = A^(q - 1) / q^(n / 2)

		Real I = std::pow(A, q - 1) / 
			std::pow(q, (Real)dimension / 2);

		// H_q(X) = (1 - I) / (q - 1)

		return (1 - I) / (q - 1);
	}

}

#endif
