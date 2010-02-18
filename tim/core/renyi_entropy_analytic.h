// Description: Analytical solutions for Renyi entropies

#ifndef TIM_RENYI_ENTROPY_ANALYTIC_H
#define TIM_RENYI_ENTROPY_ANALYTIC_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Renyi entropy of a gaussian distribution.
	/*!
	Preconditions:
	dimension > 0
	covarianceDeterminant >= 0

	q:
	The power in the definition of Renyi entropy.
	If q = 1, gaussianDifferentialEntropy() is returned
	instead.

	dimension:
	The dimension of the distribution.

	covarianceDeterminant:
	Determinant of the covariance matrix of the distribution.

	Returns:
	The Renyi entropy of the gaussian distribution.
	*/
	TIM real gaussianRenyiEntropy(
		real q,
		integer dimension,
		real covarianceDeterminant);

	//! Renyi entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume > 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})

	Returns:
	The Renyi entropy of the uniform distribution.

	Note:
	The Renyi entropy of a uniform distribution
	is independent of q.
	*/
	TIM real uniformRenyiEntropy(
		real supportVolume);

}

#endif
