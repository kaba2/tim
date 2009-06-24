#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <iostream>

namespace Tim
{

	TIMCORE std::ostream& operator<<(std::ostream& stream, const Signal& signal);

	//! Generates uniform random variables in [0, 1]^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateUniform(
		integer size,
		integer dimension);

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateGaussian(
		integer size,
		integer dimension);

	//! Generates correlated gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0

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
		integer size,
		integer dimension,
		const CholeskyDecomposition<Dynamic, real>& correlationCholesky);

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer size,
		integer dimension,
		real shape,
		real scale);

	/*!
	Preconditions:
	!signal.empty()
	*/

	TIMCORE void splitDimensions(
		const SignalPtr& signal,
		std::vector<SignalPtr>& signalSet);

	TIMCORE SignalPtr mergeDimensions(
		const std::vector<SignalPtr>& signalList);

}

#endif
