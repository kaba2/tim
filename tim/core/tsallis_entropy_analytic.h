// Description: Analytical solutions for Tsallis entropies

#ifndef TIM_TSALLIS_ENTROPY_ANALYTIC_H
#define TIM_TSALLIS_ENTROPY_ANALYTIC_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Tsallis entropy of a gaussian distribution.
	/*!
	Preconditions:
	dimension > 0
	covarianceDeterminant >= 0

	dimension:
	The dimension of the distribution.

	q:
	The power in the definition of Renyi entropy.
	If q = 1, gaussianDifferentialEntropy() is returned
	instead.

	covarianceDeterminant:
	Determinant of the covariance matrix of the distribution.

	Returns:
	The Tsallis entropy of the gaussian distribution.
	*/
	TIM real gaussianTsallisEntropy(
		integer dimension,
		real q,
		real covarianceDeterminant);

	//! Tsallis entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume >= 0

	q:
	The power in the definition of Tsallis entropy.
	If q == 1, uniformDifferentialEntropy() is
	returned instead.

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})

	Returns:
	The Tsallis entropy of the uniform distribution.
	*/
	TIM real uniformTsallisEntropy(
		real q,
		real supportVolume);

}

#endif
