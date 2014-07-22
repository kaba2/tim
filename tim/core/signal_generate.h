// Description: Generation of common signals

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

	TIM Signal generateUniform(
		integer samples,
		integer dimension);

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0
	*/

	TIM Signal generateGaussian(
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

	TIM Signal generateCorrelatedGaussian(
		integer samples,
		integer dimension,
		const CholeskyDecomposition<real>& covarianceCholesky);

	TIM Signal generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale);

	//! Generates a signal with time-varying coupling.
	/*!
	Preconditions:
	samples >= 0
	yzShift >= 0
	zyShift >= 0

	The signals are divide into three time regions.
	In the first and the third time regions, there is
	no coupling between x, y, and z. However, in
	the second time region x drives y, and y drives z.
	The amplitudes of these drives are given by a sine
	wave for the x->y, and by a cosine wave for the
	y->z. Thus, those estimators which are sensitive
	to temporal changes in coupling (e.g. partial 
	transfer entropy) should give similar coupling curves.	
	*/
	TIM void generateTimeVaryingCoupling(
		integer samples,
		integer yxShift,
		integer zyShift,
		Signal& xSignal,
		Signal& ySignal,
		Signal& zSignal);

}

#endif
