// Description: Tools to generate various signals.

#ifndef TIM_SIGNAL_GENERATE_H
#define TIM_SIGNAL_GENERATE_H

#include "tim/core/signal.h"

#include <pastel/math/cholesky_decomposition.h>

namespace Tim
{

	//! Generates uniform random variables in [0, 1]^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateUniform(
		integer samples,
		integer dimension);

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0
	*/

	TIMCORE SignalPtr generateGaussian(
		integer samples,
		integer dimension);

	//! Generates correlated gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0

	The correlated gaussian random variable is given by
	multiplying a standard gaussian random variable
	with the lower triangular part of the cholesky decomposition 
	of the correlation matrix.

	If the given correlation matrix turns out not to
	be numerically positive definite then
	the function call is equivalent to calling
	generateGaussian() (resulting in the 
	correlation matrix being identity).
	*/

	TIMCORE SignalPtr generateCorrelatedGaussian(
		integer samples,
		integer dimension,
		const CholeskyDecompositionD& covarianceCholesky);

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale);

}

#endif
