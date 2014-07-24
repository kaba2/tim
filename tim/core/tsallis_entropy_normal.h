// Description: Tsallis entropy of a normal distribution

#ifndef TIM_TSALLIS_ENTROPY_NORMAL_H
#define TIM_TSALLIS_ENTROPY_NORMAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real_concept.h>

namespace Tim
{

	//! Tsallis entropy of a normal distribution.
	/*!
	Preconditions:
	dimension > 0
	covarianceDeterminant >= 0

	q:
	The power in the definition of Renyi entropy.
	If q = 1, normalDifferentialEntropy() is returned
	instead.

	dimension:
	The dimension of the distribution.

	covarianceDeterminant:
	Determinant of the covariance matrix of the distribution.
	*/
	template <typename Real>
	Real normalTsallisEntropy(
		const PASTEL_NO_DEDUCTION(Real)& q,
		integer dimension,
		const PASTEL_NO_DEDUCTION(Real)& covarianceDeterminant);

}

#include "tim/core/tsallis_entropy_normal.hpp"

#endif