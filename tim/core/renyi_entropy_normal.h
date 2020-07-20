// Description: Analytical solutions for Renyi entropies
// Documentation: renyi_entropy_analytic.txt

#ifndef TIM_RENYI_ENTROPY_NORMAL_H
#define TIM_RENYI_ENTROPY_NORMAL_H

#include "tim/core/mytypes.h"
#include "tim/core/differential_entropy_normal.h"

#include <pastel/sys/real/real_concept.h>

#include <cmath>

namespace Tim
{

	//! Renyi entropy of a gaussian distribution.
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
	Real normalRenyiEntropy(
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

		// H_q^*(X) = (n / 2)[log(2 pi) - (log(q) / (1 - q))] + 
		//            (1 / 2)log(det(C))

		return ((Real)dimension / 2) * 
			(std::log(2 * constantPi<Real>()) - std::log(q) / (1 - q)) +
			std::log(covarianceDeterminant) / 2;
	}

}

#endif
