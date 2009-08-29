// Description: Analytical solutions for differential entropies

#ifndef TIM_DIFFERENTIAL_ENTROPY_ANALYTIC_H
#define TIM_DIFFERENTIAL_ENTROPY_ANALYTIC_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Differential entropy of a gaussian distribution.
	/*!
	Preconditions:
	dimension > 0
	covariance >= 0
	*/
	TIMCORE real gaussianDifferentialEntropy(
		integer dimension, real covarianceDeterminant);

	//! Differential entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume >= 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})
	*/
	TIMCORE real uniformDifferentialEntropy(
		real supportVolume);

	//! Differential entropy of a generalized gaussian distribution.
	/*!
	Preconditions:
	dimension > 0
	shape > 0
	scale > 0
	*/
	TIMCORE real generalizedGaussianDifferentialEntropy(
		integer dimension, real shape, real scale);

}

#endif
