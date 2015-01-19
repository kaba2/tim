// Description: Analytical solutions for Renyi entropies
// Documentation: renyi_entropy_analytic.txt

#ifndef TIM_RENYI_ENTROPY_NORMAL_H
#define TIM_RENYI_ENTROPY_NORMAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real/real_concept.h>

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
		const NoDeduction<Real>& covarianceDeterminant);

}

#include "tim/core/renyi_entropy_normal.hpp"

#endif
