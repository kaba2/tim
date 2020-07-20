// Description: Tsallis entropy of a normal distribution
// Documentation: tsallis_entropy_analytic.txt

#ifndef TIM_TSALLIS_ENTROPY_NORMAL_H
#define TIM_TSALLIS_ENTROPY_NORMAL_H

#include "tim/core/mytypes.h"
#include "tim/core/differential_entropy_normal.h"

#include <pastel/sys/real/real_concept.h>

namespace Tim
{

	//! Tsallis entropy of a normal distribution.
	/*!
	Preconditions:
	dimension > 0
	covarianceDeterminant >= 0

	q:
	The power in the definition of Renyi entropy.
	If q = 1, differentialEntropyNormal() is returned
	instead.

	dimension:
	The dimension of the distribution.

	covarianceDeterminant:
	Determinant of the covariance matrix of the distribution.
	*/
	template <typename Real>
	Real normalTsallisEntropy(
		const NoDeduction<Real>& q,
		integer dimension,
		const NoDeduction<Real>& covarianceDeterminant)
	{
		PENSURE_OP(q, >, 0);
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(covarianceDeterminant, >=, 0);

		if (q == 1)
		{
			return differentialEntropyNormal<Real>(
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

#include "tim/core/tsallis_entropy_normal.hpp"

#endif
